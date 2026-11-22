#!/usr/bin/env python3
"""gen_exoplanet_mod.py -- write a case's exoplanet_mod.F90 from user_nl_exo.

ExoCAM 2.0 (EXOCAM2_PLAN.md, Phase 2a). The planet is defined by namelist-style
settings instead of hand-edited Fortran:

    exoplanet_mod.F90.template      every parameter, documented, aqua defaults
  + configs/<config>/user_nl_exo    the config's departures from those defaults
  + $CASEROOT/user_nl_exo           the case's own settings (optional)
  = SourceMods/src.share/exoplanet_mod.F90

Later layers override earlier ones. The values are still compile-time Fortran
`parameter`s, so changing one still means rebuilding the case.

user_nl_exo syntax: one `name = value` per line; `!` starts a comment; values
are Fortran literals (.true./.false., 1.0, 4.e3, 1.7e-6_r8, 'a/string').
Only parameters declared in the template's editable section may be set.

N2 is the fill gas and is always computed here, never by hand: unless the
settings give exo_n2bar explicitly, exo_n2bar = target - (sum of the other
gases), where target is 1 bar, or the sum itself if the other gases already
exceed 1 bar. This is the rule exocam-casemgr has always used, and the value is
written in the same format, so a case generated here compiles to the same
numbers as one casemgr rendered. The derived block at the end of the template
(mixing ratios, molecular weight, heat capacity) stays as Fortran expressions
for the compiler to evaluate.

Usage:
  gen_exoplanet_mod.py CONFIG [--user-nl FILE ...] [--out FILE]
  gen_exoplanet_mod.py CONFIG --show          # print the resolved settings
Without --out the file is written to stdout.
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CESM_DIR = os.path.dirname(HERE)
TEMPLATE = os.path.join(HERE, 'exoplanet_mod.F90.template')

# The editable section ends at this marker; nothing below it may be set.
END_EDITABLE = 'No modifications below this point'

# Radiatively active gases; N2 fills the remainder (exocam-casemgr GAS_BAR_PARAMS).
GAS_BARS = ('exo_co2bar', 'exo_ch4bar', 'exo_c2h6bar', 'exo_nh3bar',
            'exo_cobar', 'exo_h2bar', 'exo_o2bar')

# Computed in the template from other parameters; never set directly.
COMPUTED = ('exo_pstd',)

# `<type>, public, parameter :: name = value   ! comment`
RE_PARAM = re.compile(
    r'^(\s*(real|integer|logical|character)\b[^:]*\bparameter\b[^:]*::\s*)'
    r'(\w+)(\s*=\s*)(.*?)(\s*!.*)?$', re.IGNORECASE)


class ExoError(Exception):
    pass


def read_template(path=TEMPLATE):
    """Return (lines, decls). decls maps lower-case name -> (line index, type)
    for each active parameter line in the editable section."""
    with open(path) as fh:
        lines = fh.readlines()
    decls = {}
    for i, line in enumerate(lines):
        if END_EDITABLE in line:
            break
        if line.lstrip().startswith('!'):
            continue
        m = RE_PARAM.match(line.rstrip('\n'))
        if m:
            decls[m.group(3).lower()] = (i, m.group(2).lower())
    return lines, decls


def parse_user_nl(path):
    """Parse a user_nl_exo file into an ordered {name: literal} dict."""
    out = {}
    with open(path) as fh:
        for n, raw in enumerate(fh, 1):
            line = _strip_comment(raw).strip()
            if not line:
                continue
            if '=' not in line:
                raise ExoError(f"{path}:{n}: expected 'name = value', got: {raw.strip()}")
            name, value = (s.strip() for s in line.split('=', 1))
            if not re.fullmatch(r'[A-Za-z]\w*', name) or not value:
                raise ExoError(f"{path}:{n}: expected 'name = value', got: {raw.strip()}")
            out[name.lower()] = value
    return out


def _strip_comment(line):
    """Drop a trailing ! comment, leaving ! inside quotes alone."""
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
        elif ch == '!':
            return line[:i]
    return line


def _real(literal, where):
    s = re.sub(r'_(r8|R8|r4|R4|\d+)$', '', literal.strip()).lower().replace('d', 'e')
    try:
        return float(s)
    except ValueError:
        raise ExoError(f"{where}: not a real number: {literal}")


def format_real(value):
    """exocam-casemgr's _fortran_value for reals, so both render identically."""
    f = float(value)
    if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e8):
        s = f"{f:.10e}"
    else:
        s = f"{f:.12g}"
        if '.' not in s and 'e' not in s:
            s += '.0'
    return f"{s}_r8"


def format_value(name, ftype, literal):
    where = f"user_nl_exo {name}"
    lit = literal.strip()
    if ftype == 'logical':
        v = lit.lower().strip('.')
        if v in ('true', 't'):
            return '.true.'
        if v in ('false', 'f'):
            return '.false.'
        raise ExoError(f"{where}: not a logical (.true./.false.): {literal}")
    if ftype == 'integer':
        if not re.fullmatch(r'[+-]?\d+', lit):
            raise ExoError(f"{where}: not an integer: {literal}")
        return str(int(lit))
    if ftype == 'character':
        if len(lit) >= 2 and lit[0] == lit[-1] and lit[0] in ("'", '"'):
            lit = lit[1:-1]
        if "'" in lit:
            raise ExoError(f"{where}: a single quote inside the string is not supported")
        return f"'{lit}'"
    return format_real(_real(lit, where))


def resolve(config, user_nl_files=()):
    """Merge the config's user_nl_exo and the given files, later winning."""
    layers = []
    cfg_nl = os.path.join(CESM_DIR, 'configs', config, 'user_nl_exo')
    if not os.path.isdir(os.path.join(CESM_DIR, 'configs', config)):
        raise ExoError(f"no config '{config}' under {os.path.join(CESM_DIR, 'configs')}")
    if os.path.isfile(cfg_nl):
        layers.append(cfg_nl)
    layers.extend(user_nl_files)
    settings, origin = {}, {}
    for path in layers:
        for k, v in parse_user_nl(path).items():
            settings[k] = v
            origin[k] = path
    return settings, origin


def render(settings, template=TEMPLATE):
    lines, decls = read_template(template)
    for name in settings:
        if name in COMPUTED:
            raise ExoError(f"{name} is computed from the gas inventories; set the gases instead")
        if name not in decls:
            raise ExoError(f"unknown parameter '{name}' (not in the editable section of "
                           f"{os.path.basename(template)})")

    values = dict(settings)
    if 'exo_n2bar' not in values:
        gases = 0.0
        for g in GAS_BARS:
            lit = values.get(g)
            if lit is None:
                i, _ = decls[g]
                lit = RE_PARAM.match(lines[i].rstrip('\n')).group(5)
            gases += _real(lit, f"template {g}")
        target = 1.0 if gases <= 1.0 else gases
        values['exo_n2bar'] = repr(target - gases)

    out = list(lines)
    for name, lit in values.items():
        i, ftype = decls[name]
        m = RE_PARAM.match(lines[i].rstrip('\n'))
        prefix, decl_name, eq, comment = m.group(1), m.group(3), m.group(4), m.group(6) or ''
        out[i] = f"{prefix}{decl_name}{eq}{format_value(name, ftype, lit)}{comment}\n"
    return ''.join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], allow_abbrev=False)
    ap.add_argument('config', help='config under cesm1.2.1/configs/ (e.g. cam_aqua_fv)')
    ap.add_argument('--user-nl', action='append', default=[], metavar='FILE',
                    help="a user_nl_exo file; repeatable, later files win")
    ap.add_argument('--out', help='write here instead of stdout')
    ap.add_argument('--show', action='store_true', help='print the resolved settings and exit')
    args = ap.parse_args(argv)
    try:
        for f in args.user_nl:
            if not os.path.isfile(f):
                raise ExoError(f"user_nl_exo file not found: {f}")
        settings, origin = resolve(args.config, args.user_nl)
        if args.show:
            for k, v in settings.items():
                print(f"{k} = {v}    ! {origin[k]}")
            return 0
        text = render(settings)
    except ExoError as e:
        print(f"gen_exoplanet_mod: ERROR: {e}", file=sys.stderr)
        return 1
    if args.out:
        with open(args.out, 'w') as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
