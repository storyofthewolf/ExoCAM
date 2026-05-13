# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What ExoCAM is (and isn't)

ExoCAM is **a patch/overlay for NCAR CESM 1.2.1**, not a standalone model. The repo contains SourceMods (Fortran files that override CESM source), namelist templates, machine config files, initial-condition file references, and post-processing tools. It does not build on its own — files here are copied into a checked-out CESM 1.2.1 tree, and runtime radiative transfer is supplied by a sister repo, [ExoRT](https://github.com/storyofthewolf/ExoRT).

There is no build system, test suite, lint config, or CI in this repo. "Running the code" means: install CESM 1.2.1 separately, install ExoRT separately, copy files from here into both, then use CESM's `create_newcase` / `cesm_setup` / `$CASE.build` / `$CASE.run` workflow on a supported HPC machine. Full procedure is in `cesm1.2.1/instructions/general_instructions.txt`.

## Repository layout

- `cesm1.2.1/configs/<config>/` — one directory per supported model configuration. Each contains a `SourceMods/` tree (copied into `$CASEROOT/SourceMods/`) and a `namelist_files/` directory (copied into `$CASEROOT/`). Configs: `cam_aqua_fv`, `cam_aqua_se`, `cam_land_fv`, `cam_mixed_fv`, `circumbinary`, `carma` (add-on), plus `experimental/` and `extras/` (gitignored, not for general use).
- `cesm1.2.1/ccsm_utils_files/` — machine/compset/grid XML and batch templates that get copied into `cesm1_2_1/scripts/ccsm_utils/Machines/` and `Case.template/`. Supported machines: `hyak`, `discover`, `summit`, `computecanada`.
- `cesm1.2.1/initial_files/` — paths to NetCDF initial condition / topography / ozone / gravity-wave files referenced from `user_nl_*` namelists. Per-config subdirectories mirror the configs above.
- `cesm1.2.1/instructions/` — authoritative how-to docs. Read these before changing build/run procedure: `general_instructions.txt`, `supported_planet_mods.txt`, `adding_oxygen_ozone.txt`, `computecanada_instructions.txt`.
- `tools/py_progs/`, `tools/idl_progs/` — post-processing utilities; users add these directories to `PYTHONPATH` / `IDL_PATH`. `tools/spectral_albedos/` holds spectral surface albedo files for `broadband_albedo_calculator.py`.
- `tools/externals` — pointers to companion GitHub repos (`trend`, `analysis`, docker image, inverse-CO2 tools).

## The central file: `exoplanet_mod.F90`

`cesm1.2.1/configs/<config>/SourceMods/src.share/exoplanet_mod.F90` is the **single source of truth for planet/atmosphere/run parameters**. Almost all "what kind of planet/star/atmosphere" choices live here as `parameter` constants and are baked in at compile time — they override CESM namelist GHG settings when `do_exo_atmconst = .true.` (the default). Editing this file requires rebuilding the case.

Key parameter groups (see file comments for full docs):
- Run options: `do_exo_synchronous`, `do_exo_rt`, `do_exo_atmconst`, `do_carma_exort`, `do_exo_gw`, `exo_convect_plim`
- Radiation: `exo_rad_step`, `do_exo_rt_clearsky`, `do_exo_rt_spectral`, `do_exo_rt_optimize_bands`, `Tmax`
- Planet: `exo_planet_radius`, `exo_surface_gravity`, `exo_ndays`, `exo_porb`, `exo_sday`, `exo_eccen`, `exo_obliq`, `exo_mvelp`
- Stellar: `exo_scon`, `exo_solar_file` (must match RT spectral resolution)
- Atmosphere: `exo_n2bar`, `exo_co2bar`, `exo_ch4bar`, `exo_h2bar`, `exo_o2bar`, `exo_c2h6bar`, `exo_nh3bar`, `exo_cobar`. Derived `vmr`/`mmr`/`exo_mwdair`/`exo_cpdair` are computed from these — **do not edit the derived block**. Note: NH₃ and CO default to `0.0` and contribute to `exo_cpdair` via `cpnh3 = 2.175e3` and `cpco = 1.040e3` J/kg/K.

Critical consistency rules enforced only by convention (the model will not check):
- `exo_pstd` (sum of partial pressures × 10⁵) must match the total pressure of the `ncdata` initial condition file in `user_nl_cam`.
- `exo_solar_file` must match the spectral binning of the linked ExoRT package (e.g. `n68equiv` solar file with `src.cam.n68equiv`).
- For synchronous rotators set `exo_sday = 86400 * exo_ndays`; for asynchronous set `exo_sday = 86400 * exo_ndays / (1 + exo_ndays/exo_porb)`. Verify by inspecting output `FDS` (incident stellar flux) for the expected diurnal cycle.

## Configuration vs. ExoRT pairing

`do_exo_rt = .true.` requires a matching ExoRT package, linked via the CESM build:

```
xmlchange CAM_CONFIG_OPTS="-nlev 40 -phys cam4 -usr_src $EXORT_PATH/ExoRT/3dmodels/src.cam.n<spectra>"
```

Common spectra: `n68equiv` (default, recommended since 2020), `n84equiv` (extends shortward for F-dwarfs), `n28archean` (high-CO₂ Archean), `n42h2o` (H₂O/N₂/H₂ only). Match `exo_solar_file`, the ozone file, and (for `circumbinary`) the in-tree RT to the chosen package. Per-config quickstart in each `cesm1.2.1/configs/<config>/README`.

For "MG" (CAM5 Morrison-Gettelman) clouds instead of default "RK" (CAM4 Rasch-Kristjansson), append `-chem none -microphys mg1` to `CAM_CONFIG_OPTS` **and** uncomment the MG block (`eddy_scheme`, `macrop_scheme`, `shallow_scheme`, `uwshcu_rpen`) in `user_nl_cam`. The two go together.

## SourceMods layering

Each config's `SourceMods/` mirrors CESM component subdirs and only the listed components are present:
- `src.cam` — atmosphere (always present)
- `src.share` — shared, **always contains `exoplanet_mod.F90`**, plus `shr_orb_mod.F90`, `shr_const_mod.F90`
- `src.cice`, `src.docn` — sea ice + slab ocean (aqua/mixed configs)
- `src.clm` — land model (land/mixed configs)
- `src.drv` — coupler (`ccsm_comp_mod.F90`, `seq_flux_mct.F90`)

When propagating a change to `exoplanet_mod.F90` or other shared SourceMods, **apply it to every config that contains the file** — there is no shared copy. Recent commit `019c2d3` ("propagate changes from cam_mixed_fv to cam_aqua_fv and cam_land_fv") and the 2026-04-29 NH₃/CO additions show the expected pattern. Configs are kept in sync by hand.

Adding a new gas absorber to the ExoRT 3-D interface requires parallel edits in this repo: `mwXXX` in `src.cam/physconst.F90`, and `exo_xxxbar`/`exo_xxxmmr`/`cpXXX`/`exo_cpdair` in `src.share/exoplanet_mod.F90`, applied to every config. See `ExoRT/CLAUDE.md` §"Connecting a New Gas to the 3-D Interface" for the full checklist.

## Namelist file conventions

Each `cesm1.2.1/configs/<config>/namelist_files/` holds `user_nl_cam`, `user_nl_cpl`, plus `user_nl_cice` / `user_nl_docn` / `user_docn.streams.txt.som` for ocean configs. Paths inside (`ncdata`, `bnd_topo`, `gw_drag_file`, `prescribed_ozone_datapath`) are written for Discover (`/discover/nobackup/etwolf/...`) — **users must edit these to point at their local ExoCAM checkout** before running. Same applies to `exo_solar_file` in `exoplanet_mod.F90` and `exort_rootdir` in ExoRT's `sys_rootdir.F90`.

## Adding O₂/O₃

See `cesm1.2.1/instructions/adding_oxygen_ozone.txt`. Two-step: set `exo_o2bar` in `exoplanet_mod.F90`; for O₃ set `prescribed_ozone_datapath`/`prescribed_ozone_file`/`prescribed_ozone_cycle_yr` in `user_nl_cam` (must point to a NetCDF with O₃ as `(lat,lon,lev,time)` and the cycle year must exist in the file or the run fails immediately). With non-zero ozone, raise `exo_convect_plim` from the ExoCAM default 5 Pa toward CESM's original 4000 Pa or the convection scheme will crash.

## Working in this repo

- This is research code maintained primarily by one author. Configurations drift; not all combinations are tested. Cross-config compatibility is not guaranteed.
- Prefer editing existing config files over adding new configs. If a change touches several configs, apply it everywhere consistently.
- Do not commit machine-specific paths back to upstream namelists — but the file convention is to keep Discover paths there as templates. Don't "fix" them away.
- `experimental/` and `extras/` are gitignored intentionally; don't restore them.
- Companion repo ExoRT must be checked out separately for any actual model run.

## Active work: CO2 condensation debugging (branch `co2_claude_diag`)

As of 2026-05-13, active development is on branch `co2_claude_diag`. This targets the `experimental/co2condense` configuration only — a Mars-targeted experimental build. Do not apply changes from this branch to other configs without explicit review.

**Full session notes:** `DEVELOPER_NOTES.md` (project root)  
**Full bug inventory:** `Opus4.7_Bug_Report.md` (project root)

### Key design constraints for `experimental/co2condense`

- The model currently runs **pure CO2** (`exo_co2vmr = 1.0`, compile-time constant). Local composition feedback is intentionally disabled. Do not add local vmr tracking until explicitly requested.
- `SHR_CONST_CPDAIR` is overridden to the CO2 value in `shr_const_mod.F90`. Uses of `SHR_CONST_CPDAIR` in `exo_condense_mod.F90` are correct — do not change them to Earth's 1004.64 J/kg/K.
- `cam_out%co2srf_snow` carries kg/m² (not kg/m²/s) despite being registered as `a2x_fluxes`. This is intentional and relies on `cpl_dt == atm_dt`. Both ends are commented. Do not change `cpl_dt` without revisiting this field.

### What was fixed (commits `eea21c9`, `05695cf`)

1. **P1 — Atmospheric mass budget:** `state%pdel` and `state%ps` now updated in `exo_condense_co2` after condensation. Bottom-up layer-peel in dp_coupling replaces layer-averaging. Non-`local_dp_map` branch now has identical Mars adjustment.
2. **P2 — Forget-1998 pot/heat correction:** `fxcld_in` (incoming flux per layer) now returned from `exo_cloud_sediment_tend` and used in place of `sed_tend` for pot/heat corrections. Eliminates positive feedback over thin-layer high-terrain cells.
3. **P3 — k=1 top-layer loop:** Typo `do i=i,ncol` fixed to `do i=1,ncol`. Missing `cld_tmp`/`tmid_tmp` updates and diagnostics added for k=1.
4. **P4 — CLM double-counted latent energy:** Full-sublimation branch in `exo_clm_condense.F90` no longer both cools the soil and sends energy to atmosphere from the same latent release.
5. **Budget diagnostics:** `exo_condense_diag_calc` now outputs `CO2_COL_GAS`, `CO2_COL_ICE`, `CO2_COL_SRF`, `CO2_COL_TOT` (all kg/m²). `physpkg.F90` call updated to pass `cam_in`.

### Next steps (in priority order)

1. **Run the model and inspect `CO2_COL_TOT` globally.** Flat global mean = P1 working. Negative drift = likely double-subtraction of snowfall path (see A4 in bug report).
2. **If budget closes, check Olympus Mons column** for `CO2_COL_GAS` decreasing when frost forms. If still runaway, address B3: recompute `tcond` each substep using evolving pressure rather than fixed `state%pmid`.
3. **Update derived pressure fields** (`pmid`, `pint`, `lnpmid`) after the P1 `pdel` change so downstream parameterizations see the post-condensation pressure (B4 in bug report).
4. **Investigate double-subtraction (A4):** If `CO2_COL_TOT` drifts, the snowfall path may subtract mass from `pdel` twice — once in physics via `ptend%q` and once in dp_coupling via `co2_mass_change`. Fix by either excluding sedimented mass from the physics `pdel` update, or by zeroing the snowfall component of `co2_mass_change` before dp_coupling acts on it.
