#!/usr/bin/env python3
"""
mod_ncdata.py -- modify a CAM initial-condition file (the user_nl_cam `ncdata`).

Author: Eric T. Wolf (IDL originals: changevert_cesm.pro, changepress_cesm.pro,
        mod_cesm_files.pro). Python port of the level change by Russell Deitrick
        (mod_cam.i.file.py, which this replaces). Merged, made importable, and
        the pressure rescale finished 2026-09-13 (Claude, with E.T. Wolf).

Operations, applied in this order:

  1. Change the number of vertical levels.  The new hybrid-sigma grid is cut
     from the bottom of a 66-level WACCM template
     (cesm1.2.1/initial_files/other/oxygen_CE.cam2.avg.nc); every lev-dependent
     field is interpolated column by column in pressure from the old grid onto
     the new one (scipy interp1d, linear, extrapolating at the ends).  Skipped
     when the requested nlev equals the template's.

     KNOWN DIVERGENCE FROM THE IDL -- winds.  changevert_cesm.pro does NOT
     interpolate US/VS: it zeroes them outright (the interpol calls are
     present but commented out, at lines 183-184 for VS and 207-208 for US).
     This port interpolates them like any other lev-dependent field, a
     behavior inherited from mod_cam.i.file.py, which diverged from the IDL
     the same way and was never checked against it.  Audited and confirmed
     2026-09-18; E.T. Wolf's decision that day was to KEEP interpolating for
     now and revisit later, so a level-changed IC made here carries a
     spun-up wind field where an IDL-made one starts from rest.  That is an
     O(1) difference in the wind field, not a rounding effect.  Everything
     else in the level change -- the cut-from-bottom index, linear-in-
     pressure interpolation, coordinate arrays, unbounded extrapolation and
     the PS/P0 ordering -- was checked against the IDL and agrees.

  2. Change the surface pressure.  This reproduces changepress_cesm.pro exactly:
     the fields are LEFT ALONE in hybrid-sigma space and only PS and P0 are
     rewritten, either
        flat   : PS(:,:) = P0 = pstd            (the ExoCAM convention -- every
                 ic_*bar_* file in initial_files/ was made this way; T is
                 byte-identical to its template)
        scaled : PS(:,:) = PS_template * pstd/P0_template, P0 = pstd
     In ExoCAM the composition (and so the pressure) is set in exoplanet_mod.F90;
     ncdata must match it to within a few percent or the run crashes, and the
     model uses the exoplanet_mod value thereafter.

  3. Dry: zero Q, CLDLIQ, CLDICE.

Library use (what exocam-casemgr's `build.py prep` calls):

    from mod_ncdata import modify_ncdata
    modify_ncdata(template, out, nlev=51, pstd_bar=2.0, dry=True)

CLI (the old mod_cam.i.file.py flags, plus --exocam-root and --ps-mode):

    python mod_ncdata.py out.nc -ic control_L40.cam.i.0009-01-01-00000.nc \\
           -n 51 -ps 2.0 --dry -w

The template may be an absolute path or a bare filename searched under the
usual ExoCAM initial_files/<config> directories (and CESM inputdata if
--ccsm-inputdata is given).  The ExoCAM root defaults to the checkout this
file lives in (tools/py_progs/../..), overridable with --exocam-root or $EXOCAM.
"""
import argparse
import os
import pathlib
import shutil
import sys

import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exocampy_tools as exotools  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_THIS_DIR = pathlib.Path(__file__).resolve().parent


def default_exocam_root():
    """$EXOCAM if set, else the ExoCAM checkout containing this file."""
    env = os.environ.get('EXOCAM')
    if env:
        return pathlib.Path(env)
    return _THIS_DIR.parent.parent


# Relative to the ExoCAM root: the 66-level WACCM grid the level change cuts
# its new hybrid coordinate from.
LEVEL_GRID_FILE = 'cesm1.2.1/initial_files/other/oxygen_CE.cam2.avg.nc'

IC_SUBDIRS = ['cam_aqua_fv', 'cam_aqua_se', 'cam_land_fv', 'cam_mixed_fv',
              'mars/atm']

# Fields interpolated in the level change (everything else is copied). US is
# on the staggered latitude grid.
#
# US/VS are here deliberately: the IDL (changevert_cesm.pro) zeroes both
# instead of interpolating. Keeping them interpolated is an explicit owner
# decision (2026-09-18), revisitable -- see the docstring. Removing them from
# this tuple is all that is needed to match the IDL.
LEV_FIELDS = ('T', 'Q', 'CLDLIQ', 'CLDICE', 'US', 'VS')


def find_template(name, exocam_root=None, ccsm_inputdata=None):
    """Resolve an IC template: absolute/relative path as given, else a bare
    filename searched under the ExoCAM initial_files tree."""
    p = pathlib.Path(name)
    if p.is_file():
        return p
    root = pathlib.Path(exocam_root) if exocam_root else default_exocam_root()
    candidates = [root / 'cesm1.2.1' / 'initial_files' / sub / name
                  for sub in IC_SUBDIRS]
    if ccsm_inputdata:
        ci = pathlib.Path(ccsm_inputdata)
        candidates += [ci / 'atm/cam/inic/fv' / name, ci / 'atm/waccm/ic' / name]
    for c in candidates:
        if c.is_file():
            return c
    raise FileNotFoundError(
        f"IC template '{name}' not found; searched: "
        + ', '.join(str(c.parent) for c in candidates))


# ---------------------------------------------------------------------------
# Level change
# ---------------------------------------------------------------------------

def _read_grid_template(path):
    with nc.Dataset(path, 'r') as d:
        g = {k: np.array(d[k][:]) for k in ('lev', 'ilev', 'hyai', 'hybi', 'hyam', 'hybm')}
    return g


def _interp_columns(P_old, F_old, P_new):
    """Interpolate F_old(lev, lat, lon) from P_old(lev, lat, lon) onto
    P_new(lev', lat, lon), column by column, linear with extrapolation --
    the original mod_cam.i.file.py loop, kept explicit so the
    numerics are exactly the original's)."""
    nlev_new, nlat, nlon = P_new.shape
    out = np.zeros((nlev_new, nlat, nlon))
    for y in range(nlat):
        for x in range(nlon):
            out[:, y, x] = interp1d(P_old[:, y, x], F_old[:, y, x],
                                    fill_value='extrapolate')(P_new[:, y, x])
    return out


def change_levels(template, out, nlev_out, exocam_root=None, grid_file=None,
                  verbose=True):
    """Write `out` = `template` re-gridded to nlev_out levels.

    Every variable, dimension and attribute of the template is copied; the
    lev/ilev dimensions and coordinates come from the bottom nlev_out levels of
    the WACCM grid file; LEV_FIELDS are interpolated in pressure.
    """
    root = pathlib.Path(exocam_root) if exocam_root else default_exocam_root()
    grid_path = pathlib.Path(grid_file) if grid_file else root / LEVEL_GRID_FILE
    if not grid_path.is_file():
        raise FileNotFoundError(f"level grid template not found: {grid_path}")
    g = _read_grid_template(grid_path)
    nlev_g = len(g['lev'])
    nilev_out = nlev_out + 1
    if nlev_out > nlev_g:
        raise ValueError(f"nlev={nlev_out} exceeds the {nlev_g}-level grid template")

    with nc.Dataset(template, 'r') as src:
        lat = np.array(src['lat'][:]); lon = np.array(src['lon'][:])
        nlat, nlon = len(lat), len(lon)
        P0 = float(np.array(src['P0'][:]).squeeze())
        PS = np.array(src['PS'][:]).squeeze()
        old = {k: np.array(src[k][:]) for k in ('hyai', 'hybi', 'hyam', 'hybm')}
        nlev_old = len(src['lev'])
        if verbose:
            print(f"change_levels: {template}: L{nlev_old} -> L{nlev_out} "
                  f"(grid from {grid_path.name})")

        # Pressure of the old grid and of the FULL new grid, both under the
        # template's own PS -- the level change is done at the template's
        # pressure; any pressure rescale happens afterwards in sigma space.
        lev_P_old, _ = exotools.hybrid2pressure(
            nlon, nlat, nlev_old, PS, P0, old['hyam'], old['hybm'], old['hyai'], old['hybi'])
        lev_P_new, _ = exotools.hybrid2pressure(
            nlon, nlat, nlev_g, PS, P0, g['hyam'], g['hybm'], g['hyai'], g['hybi'])
        n = nlev_g - nlev_out          # first index of the kept (bottom) levels
        lev_P_new = lev_P_new[n:, :, :]

        # New coordinate arrays: the bottom nlev_out levels of the grid template
        coord = {
            'lev':  g['lev'][nlev_g - nlev_out:],
            'ilev': g['ilev'][nlev_g + 1 - nilev_out:],
            'hyam': g['hyam'][nlev_g - nlev_out:],
            'hybm': g['hybm'][nlev_g - nlev_out:],
            'hyai': g['hyai'][nlev_g + 1 - nilev_out:],
            'hybi': g['hybi'][nlev_g + 1 - nilev_out:],
        }

        interp = {}
        for name in LEV_FIELDS:
            if name not in src.variables:
                continue
            F = np.array(src[name][:])[0]                 # drop time
            if name == 'US':                              # staggered lat
                nslat = F.shape[1]
                interp[name] = _interp_columns(
                    lev_P_old[:, :nslat, :], F, lev_P_new[:, :nslat, :])[None]
            else:
                interp[name] = _interp_columns(lev_P_old, F, lev_P_new)[None]

        with nc.Dataset(out, 'w', format=src.data_model) as dst:
            for dname, dim in src.dimensions.items():
                if dname == 'lev':
                    size = nlev_out
                elif dname == 'ilev':
                    size = nilev_out
                else:
                    size = None if dim.isunlimited() else len(dim)
                dst.createDimension(dname, size)
            for vname, var in src.variables.items():
                fill = getattr(var, '_FillValue', None)
                v = dst.createVariable(vname, var.dtype, var.dimensions,
                                       fill_value=fill)
                v.setncatts({k: var.getncattr(k) for k in var.ncattrs()
                             if k != '_FillValue'})
                if vname in coord:
                    v[:] = coord[vname]
                elif vname in interp:
                    v[:] = interp[vname]
                else:
                    v[:] = var[:]
            dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
    return out


# ---------------------------------------------------------------------------
# Pressure change + dry (in place)
# ---------------------------------------------------------------------------

def apply_pressure_and_dry(path, pstd_bar=None, dry=False, ps_mode='flat',
                           verbose=True):
    """Rewrite PS/P0 (changepress_cesm.pro semantics) and/or zero the water
    fields, in place. Fields other than PS/P0/Q/CLDLIQ/CLDICE are untouched."""
    if ps_mode not in ('flat', 'scaled'):
        raise ValueError("ps_mode must be 'flat' or 'scaled'")
    with nc.Dataset(path, 'r+') as d:
        if pstd_bar is not None:
            new_ps = float(pstd_bar) * 1.0e5
            P0_in = float(np.array(d['P0'][:]).squeeze())
            PS_in = np.array(d['PS'][:])
            if ps_mode == 'flat':
                PS_out = np.full_like(PS_in, new_ps)
                P0_out = new_ps
            else:
                scalefac = new_ps / P0_in
                PS_out = PS_in * scalefac
                P0_out = P0_in * scalefac
            if verbose:
                print(f"change_pressure ({ps_mode}): P0 {P0_in:.1f} -> {P0_out:.1f} Pa, "
                      f"PS mean {PS_in.mean():.1f} -> {PS_out.mean():.1f} Pa")
            d['P0'][:] = P0_out
            d['PS'][:] = PS_out
        if dry:
            for name in ('Q', 'CLDLIQ', 'CLDICE'):
                if name in d.variables:
                    d[name][:] = 0.0
            if verbose:
                print("dry: Q, CLDLIQ, CLDICE set to zero")
    return path


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def modify_ncdata(template, out, nlev=None, pstd_bar=None, dry=False,
                  ps_mode='flat', exocam_root=None, grid_file=None,
                  overwrite=False, verbose=True):
    """Produce `out` from `template` with the requested level count, surface
    pressure and dryness. Returns a provenance dict.

    Writes to a temporary sibling and renames at the end, so a killed run never
    leaves a partial file under the final name.
    """
    template = find_template(template, exocam_root)
    out = pathlib.Path(out)
    if out.exists() and not overwrite:
        raise FileExistsError(f"{out} exists (overwrite=True to replace)")
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + '.partial')

    with nc.Dataset(template, 'r') as d:
        nlev_in = len(d['lev'])
    level_changed = nlev is not None and int(nlev) != nlev_in
    if level_changed:
        change_levels(template, tmp, int(nlev), exocam_root=exocam_root,
                      grid_file=grid_file, verbose=verbose)
    else:
        shutil.copyfile(template, tmp)
        if verbose:
            print(f"copy: {template} (L{nlev_in}, levels unchanged)")
    apply_pressure_and_dry(tmp, pstd_bar=pstd_bar, dry=dry, ps_mode=ps_mode,
                           verbose=verbose)
    os.replace(tmp, out)
    if verbose:
        print(f"wrote {out}")
    return {
        'template': str(template), 'out': str(out),
        'nlev_in': nlev_in, 'nlev': int(nlev) if nlev is not None else nlev_in,
        'level_changed': level_changed,
        'pstd_bar': pstd_bar, 'ps_mode': ps_mode if pstd_bar is not None else None,
        'dry': bool(dry),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fname_out', help='Name of new IC file')
    p.add_argument('-ic', '--input_IC_file', required=True,
                   help='template IC file (path, or bare name under initial_files/)')
    p.add_argument('-n', '--num_lev', type=int, default=None,
                   help='Number of levels in new IC file (default: keep the template\'s)')
    p.add_argument('-ps', '--surface_pressure', type=float, default=None,
                   help='Surface pressure in bars (default: keep the template\'s)')
    p.add_argument('--ps-mode', choices=('flat', 'scaled'), default='flat',
                   help='flat: PS=P0=pstd everywhere (ExoCAM convention); '
                        'scaled: scale the template PS field')
    p.add_argument('-d', '--dry', action='store_true',
                   help='zero out water and clouds from new ic')
    p.add_argument('-w', '--overwrite', action='store_true',
                   help='force overwrite of output files')
    p.add_argument('--exocam-root', default=None,
                   help='ExoCAM checkout (default: $EXOCAM or this file\'s checkout)')
    p.add_argument('--grid-file', default=None,
                   help=f'66-level grid template (default: <exocam-root>/{LEVEL_GRID_FILE})')
    args = p.parse_args(argv)

    print("\n-----------------------------------------------------------------------------")
    print("Entering mod_ncdata.py ...")
    try:
        modify_ncdata(args.input_IC_file, args.fname_out, nlev=args.num_lev,
                      pstd_bar=args.surface_pressure, dry=args.dry,
                      ps_mode=args.ps_mode, exocam_root=args.exocam_root,
                      grid_file=args.grid_file, overwrite=args.overwrite)
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print("Exiting ...")
    print("-----------------------------------------------------------------------------\n")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
