# ExoCAM 2.0 — plan

**Status:** agreed 2026-11-21. Phase 0 baselines recorded for aqua and mixed; Phase 1 steps 1–4 done and BFB for aqua and mixed (2026-11-22); steps 5–6 (one `exoplanet_mod.F90`, shared namelists) added 2026-11-23. Land baseline (Ben1 hybrid) recorded and BFB (2026-11-23). Phase 2 deferred past 2.0 (2026-11-23): `exoplanet_mod.F90` stays compile-time; 2a reverted. This is the working plan; edit it as decisions land.

## Decisions

- **ExoCAM 2.0 is an infrastructure release.** The model keeps its current
  shape (an overlay on CESM 1.2.1), but:
  - one SourceMods set replaces the per-config copies;
  - ~~`user_nl_exo` replaces hand-edited `exoplanet_mod.F90`~~ (deferred past 2.0,
    2026-11-23: see Phase 2);
  - ExoRT 2.0 (`src.cam.exort`) is hooked up;
  - regression baselines guard every step.
- **2.0 ships as one presentable package** (author, 2026-11-23): ExoCAM
  `refactor`, ExoRT `refactor` and casemgr `v2` are released together, and a
  new user must be able to install and run them from the docs, as with the
  original ExoCAM. Working code alone is not enough. Phase 4 holds the items.
- **Released together with ExoRT 2.0.** Until that day, ExoCAM 1.x with
  ExoRT 1.x (`v1.1.0`, `n68equiv`) stays the public pair, and ExoRT `main`
  stays v1.
- **The proposal physics lands afterwards, as 2.x releases** on the single
  tree: the colleague's convection PR, working CO₂ condensation, and a
  cleaner CARMA. Then the test/benchmark suite grows, then an agent harness.
- **Out of 2.0:**
  - `circumbinary`: frozen at tag `circumbinary_jgrplanets_frozen`;
  - `cam_aqua_se`: 2.0 is FV only;
  - the unified-repo redesign: later, a 3.0 question.

## Where things stand (2026-11-21)

- **ExoCAM `main`** is `v1.0.0` plus 14 commits: the AM conservation fixes
  (PR #14), the dycore robustness patches and the `mod_ncdata` tooling.
  `tau4` changed, so it is not bit-for-bit with 1.0 even with the AM flags
  off.
- **The three FV configs share most of their SourceMods.**

  | Files (59 paths in total) | Status |
  |---|---|
  | 51 | identical wherever present (37 in all three, 14 in two) |
  | `exoplanet_mod.F90` | one version per config: the planet definition |
  | `convect_shallow.F90` | 3 versions |
  | `zm_conv_intr.F90` | 2 versions (mixed differs) |
  | `docn_comp_mod.F90` | aqua vs mixed |
  | `clm_varcon`, `clm_driver`, `Biogeophysics1Mod`, `BareGroundFluxesMod` | land vs mixed |

  The two convection files are the ones the convection PR touches.
- **`exoplanet_mod.F90` values are compile-time Fortran `parameter`s.** They
  feed seven `shr_const_mod` parameters: `SCON`, `MWDAIR`, `CPDAIR`, `PSTD`,
  `SDAY`, `REARTH` and `G`. Through those they also feed the derived
  `OMEGA`, `RDAIR`, `ZVIR` and `RHODAIR`. Twelve SourceMods files import
  `exoplanet_mod` directly.
- **ExoRT 2.0 has run in ExoCAM** (sandbox `sbx_exort_v2_n1`, 2026-11-21), at
  ~1.9× the original throughput. What a v2 case needs is in ExoRT
  `REFACTOR_PLAN.md` §0.
- **A 2024 precursor exists:** `/gpfsm/dnb33/etwolf/models/ExoCAM2` on
  Discover (not a git repo), a hand-rolled first try at the same idea. It has a
  single `SourceMods/` tree; `config_exocam.py <case> --newcase -config aqua
  -rad n68equiv`, which reads machine paths from a 5-line `sys.in` and runs
  create_newcase/clone, the XML changes and the copies (the seed of casemgr);
  `config_ccsm_utils.py`, a stub to install the machine files into a CESM tree;
  and `diff_exocam.py` (now casemgr `diff.py`). Two of its notes are still open
  here: shared namelists as well as shared SourceMods, and copying the ExoRT
  bundle into each case's SourceMods so ExoRT drift cannot break old cases
  (`-usr_src` takes a comma-separated list, so subdirectories can stay apart).
- **nsplit:** the configs ship `nsplit = 32`. Measured on pt04 (126 ranks):
  - 41.7 sim yr/day at 32, 49.0 at 16, 53.7 at 8; 4 crashed.
  - It is a per-campaign science choice, so the default stays at 32.
  - Detail: `~/Research/models/DYCORE_NEXT_PHASE.md`.

---

## Phase 0 — Baselines, before anything moves

Implement the Tier 0–2 parts of `REGRESSION_HARNESS_DESIGN.md` (as a
casemgr `regress.py`) far enough to:
1. **Build one master case per FV config** (aqua, land, mixed) from current
   `main` with `n68equiv`, and record the ExoCAM, ExoRT and CESM SHAs,
   compiler, machine and the rendered `exoplanet_mod.F90`.
2. **Run each for a short fixed period and archive the restart set.** About
   5 days is enough for bit-for-bit checks; a 30-day pt04 run takes ~5 min on
   the debug queue. Also record the timing file.
3. **The comparison:** an `hpc-nc diff` of the restarts (`cam.r`, plus
   `cice.r`/`docn.r`/`clm.r` as present). It reports which fields differ and
   by how much.

Every later phase must reproduce these baselines bit for bit, except Phase 3,
which changes the answer on purpose.

**Status (2026-11-21):** matrix drafted as casemgr
`exp_matrices/reg0_baselines.yaml` (cases `reg0_{aqua,land,mixed}_fv`;
nlev 40/40/51; 126 tasks; 5 days, restart at day 5). The baselines use
existing casemgr verbs (`build.py make --send-it`); `regress.py` wraps the
same steps later. Two things found while drafting:
- casemgr does not inherit template gases, so the matrix states the shipped
  composition (400 ppm CO₂, 1.7 ppm CH₄; mixed also O₂ 0.2095 bar and
  `exo_convect_plim` 4000 Pa).
- casemgr renders `exo_n2bar` as a rounded literal (`0.9995983_r8`), not the
  template's expression (`1.0 - exo_co2bar - ...`), which can differ in the
  last bit. The baseline is therefore defined by the **rendered** file, and
  the Phase 2a generator must render `exo_n2bar` the same way casemgr does
  (or both switch to the expression before the baselines are recorded).

**Done 2026-11-21:** `main` tagged `v1.1.0` at `1f81c30` (paired with ExoRT
`v1.1.0`), so public users get the AM and robustness fixes now rather than at
2.0. The baselines are recorded against this tag.

## Phase 1 — One SourceMods set (bit-for-bit per config)

**Target layout (proposal):**
```
cesm1.2.1/
  SourceMods/              # the single set: src.cam src.share src.drv src.cice src.docn src.clm
  configs/aqua/            # namelist_files/ + README only
  configs/land/
  configs/mixed/
  addons/carma/            # or folded in behind a switch (decide)
```

**Steps:**
1. **Move the 51 identical files into the shared tree.** One commit; every
   config's file list must be unchanged.
   **Done 2026-11-22:** shared tree at `cesm1.2.1/SourceMods/` (config
   directory names kept, so casemgr's `config_type` values still apply);
   `configs/<cfg>/COMPONENTS` plus `cesm1.2.1/install_sourcemods.sh`
   assemble a case. All three assemble byte-identical to `v1.1.0`, so this
   step is BFB by construction.
2. **Reconcile the 7 differing files one at a time.** Diff each and classify:
   - **drift** (one config missed an update): take the right version;
   - **a genuine config difference** (land vs ocean surface): a runtime
     branch on something the model already knows (the compset's active
     components, a surface-type flag), never a CPP switch.

   Each reconciliation must pass all three baselines.

   **Diff results (2026-11-21, `d2c8720`):**

   | File | Configs | Real difference | Class |
   |---|---|---|---|
   | `zm_conv_intr.F90` | aqua=land / mixed | trailing whitespace only | cosmetic |
   | `convect_shallow.F90` | aqua / land=mixed\* | one space before `)` (\*land and mixed also differ only in trailing whitespace) | cosmetic |
   | `clm_driver.F90` | land / mixed | spacing inside a comment | cosmetic |
   | `docn_comp_mod.F90` | aqua / mixed | aqua has `+ shr_const_stebol*t_int**4`; mixed lacks `t_int` | drift: propagate to mixed (adds exactly `+0.0` at `t_int = 0`; BFB under `-fp-model source`, confirm against the baseline) |
   | `clm_varcon.F90` | land / mixed | `watmin` = 1e-15 (land) vs 1e-10 (mixed) | science choice: per-config parameter |
   | `Biogeophysics1Mod.F90` | land / mixed | ground emissivity 1.0 (land) vs stock 0.96/0.97 (mixed) | science choice: per-config parameter |
   | `BareGroundFluxesMod.F90` | land / mixed | `q_ref2m` = grid-box `forc_q` (land) vs stock formula (mixed) | science choice: per-config parameter |

   `exoplanet_mod.F90`: the three parameter sets differ only in `t_int`
   (aqua only, drift), `exo_o2bar` (mixed 0.2095) and `exo_convect_plim`
   (mixed 4000 Pa); every other value is identical.

   **Decided 2026-11-21, the three CLM files:** keep both behaviours as
   parameters in `exoplanet_mod.F90`, set per config (for example
   `exo_clm_watmin`, `exo_ground_emis`, `do_exo_q2m_gridbox`). This is not a
   CPP switch, and it keeps each config BFB. They become `user_nl_exo` entries
   when the runtime namelist is built (Phase 2, deferred).
3. **Check that SourceMods for an inactive component is harmless.** Does
   `src.clm` copied into an aqua case compile into anything? Check it against
   the compsets; if not, have casemgr copy per component.
   **Answer (2026-11-22): harmless.** Each component's build reads only
   `SourceMods/src.<its model>`: a stub reads `src.slnd`/`src.sice`/`src.socn`
   (`models/*/s*/bld/*.buildexe.csh`), so `src.clm` in an aqua case (SLND) or
   `src.cice`/`src.docn` in a land case (SICE, SOCN) is never compiled. The
   `COMPONENTS` lists are kept anyway, so a case's SourceMods file list stays
   exactly what 1.x installed (handy for casemgr `diff.py`); they could be
   dropped later.
4. **Update casemgr and `general_instructions.txt`** for the new layout.
   **Done 2026-11-22:** `general_instructions.txt`, the config READMEs and
   `cesm1.2.1/README` (`98b666f`); casemgr `c378430` runs
   `install_sourcemods.sh` when the config has `COMPONENTS`, otherwise the old
   `cp -r`, so `main` is unaffected.
5. **One shared `exoplanet_mod.F90`** (added 2026-11-23, author). The three
   copies differ only in values:

   | Parameter | aqua | land | mixed |
   |---|---|---|---|
   | `exo_o2bar` | 0.0 | 0.0 | 0.2095 |
   | `exo_convect_plim` | 5 Pa | 5 Pa | 4000 Pa |
   | `exo_lnd_watmin` | 1e-10 | 1e-15 | 1e-10 |
   | `exo_lnd_emis_soil` / `_snow` | 0.96 / 0.97 | 1.0 / 1.0 | 0.96 / 0.97 |
   | `do_exo_lnd_q2m_gridbox` | .false. | .true. | .false. |

   The land values are read only by CLM, so aqua never sees them. Move one
   file to `SourceMods/src.share/`, with the aqua values as defaults. Each
   config's README lists its overrides for a manual setup. casemgr applies them
   from the registry (`cesm_config.<type>`), under the matrix values, which
   win. Whether mixed should simply take land's CLM values is a science choice
   for the author; it is not part of this step, since it changes mixed's answer.
   **BFB:** all three reg0 baselines, built through casemgr.
6. **Shared namelists** (added 2026-11-23, author; the 2024 precursor's other
   half). Today each config ships its own `user_nl_*` files. They differ for
   real reasons and from drift:
   - **real, per surface:** `ncdata`, `bnd_topo`, `finidat`, the sea-ice
     grid/IC files, and the ozone files (mixed has Earth ozone) in `user_nl_cam`;
     `orb_iyear` in `user_nl_cpl` (mixed 2000, the others -1);
   - **drift, to classify one by one as in step 2:** `eddy_scheme`'s position,
     `cldfrc_rhminl = 0.88` in mixed only, small `fincl1` differences, and
     commented-out spectral-output blocks for the old n68/n28 band counts (update
     to n84 or drop);
   - **a defect:** land's `user_nl_cam` has `prescribed_ozone_datapath =
     /discover/...'` with no opening quote. casemgr overwrites the line, so its
     cases are unaffected, but a manual land setup cannot parse it.

   Target: `cesm1.2.1/namelist_files/` holds the common files; each config
   keeps a short overlay of its surface inputs only. The installer combines
   them and refuses a key set in both layers, as `install_sourcemods.sh` does for
   files (CESM 1.2.1's build-namelist may reject a duplicate key anyway). casemgr
   uses the installer in place of its `cp namelist_files/*`.
   **BFB:** each config's generated `atm_in`, `drv_in`, `docn_in`, `ice_in`,
   `lnd_in` must match its reg0 case's (ignoring comments, order and the
   case-name strings), and the reg0 restarts must match bit for bit.

   With steps 5 and 6 done, a config is its `COMPONENTS` line, its namelist
   overlay and its README.

**Phase 1 status (2026-11-22):** steps 1–4 committed on `refactor`
(`98b666f` shared tree, `b71d48e` whitespace merges, `88be16f` docn/`t_int`,
`bcfc3ce` CLM parameters). Every FV config's own SourceMods is now just
`exoplanet_mod.F90`, and the three differ only in values.

**BFB result (2026-11-22, casemgr `exp_matrices/sbx_reg0r_phase1.yaml`,
built from a Discover clone of `refactor` at `d08510d`):** aqua and mixed
restarts (`cam.r`, `cam.rs`, `cice.r`, `docn.r`, `cpl.r`, and `clm2.r` for
mixed) are **identical in every data variable** to `reg0`; only the case-name
attributes differ. Land cannot be checked: `reg0_land_fv` crashes on stock
`v1.1.0` (below), and the refactor build crashes identically, with the same
`nstep, te` values to all digits for steps 0–3, the same 97 fill messages,
and the same line.

**Open: `cam_land_fv` crashes on stock v1.1.0.** SIGSEGV on 2 of 126 ranks at
the first radiation step (nstep 3), ExoRT `aerad_driver` at the call to
`rad_precalc` (`exo_radiation_mod.F90:652`), bad address
`0xfffffffc06e05698`, a hugely negative offset, as from a NaN converted to an
array index. Logs in `exocam-baselines/reg0/reg0_land_fv_FAILED/`. Needs a
debug build (`DEBUG=TRUE`) of the land case before land has a baseline.
**Land baseline (2026-11-23):** `reg0_land_hyb`, a 5-day hybrid from the spun-up
`ben1_L51_p6.1` restart (year 156) with Ben1's planet and the shipped land SourceMods;
archived in `exocam-baselines/reg0/reg0_land_hyb/`. The refactor build
(`sbx_reg0r_land_hyb`, `400aa22`: Phase 1 + 2a) is **identical in every data variable**
(`cam.r`, `cam.rs`, `clm2.r`, `cpl.r`). Phases 1 and 2a are now BFB for all three FV configs.
casemgr `9bbdbee` drops `finidat` from `user_nl_clm` for land/mixed hybrids (CLM takes the
refcase restart), as the exovolc_* hybrids have it.

**Fresh-start crash, cause found 2026-11-23:** a dry IC with Q = 1e-12 instead of exactly 0
starts and runs 5 days clean (exocam-operations `dev/dry-start`, trial 9). Fix not chosen.

**Shelved 2026-11-22 by the author.** Findings so far: the shipped
`landplanet_dry.clm2.r.nc` is byte-identical to the TRAPPIST-1e Ben1 restart
`ben1_L51_p6.1.clm2.r.0003-01-01` (commit `29b77bf`); a start from the L51 dry
IC with soil colour 5 crashes the same way at nstep 6. A `DEBUG=TRUE` trial
(`sbx_landdbg_debug`, job 59376732) was submitted and left unread.

`exoplanet_mod.F90` stays per config in this phase; Phase 2 replaces it.

## Phase 2 — the planet definition: deferred out of 2.0

**Decision (author, 2026-11-23): `exoplanet_mod.F90` stays as it is, compile-time
`parameter`s, one file per FV config, for ExoCAM 2.0.** Revisit after 2.0.

**Why.** Two designs were tried or drafted on 2026-11-22 and rejected:
- *2a, `user_nl_exo` generating Fortran* (implemented in `960a18a`, bitwise-verified,
  reverted). A "namelist" that is really compile-time input is unlike every other
  `user_nl_*` file, needs a regenerate step that is easy to forget, and has to be
  kept in lock-step with casemgr's own rendering.
- *2b, a two-tier split* (run options at runtime, planet constants compiled). One
  file holding two kinds of setting with different rules is bad organisation, and it
  optimises the wrong thing.

casemgr is changed in tandem with whatever replaces `exoplanet_mod.F90`, so the
design is settled first and casemgr rebuilt once.

**Target design for the revisit: a true runtime `user_nl_exo`.** One independent
namelist that behaves like the others: read at model start, no rebuild to change a
value, every setting runtime (planet constants included). That means turning
`exoplanet_mod`'s parameters, and everything computed from them at compile time,
into module variables set once at initialisation.

**Audit (2026-11-23): what would have to change.** Every compile-time context
(`parameter` declaration or declaration initialiser, which Fortran requires to be
constant) that reaches a settable `exoplanet_mod` value, followed across modules
through `use`, in the 599 Fortran files a `cam_mixed_fv` case compiles (atm, lnd,
ice, docn, driver, csm_share; source lists from `reg0_mixed_fv`'s build; the aqua
and land configs compile subsets):

| Where | File | Compile-time uses |
|---|---|---|
| ExoCAM SourceMods (already owned) | `exoplanet_mod.F90` | `exo_n2bar`, `exo_pstd`, the vmr/mmr/`exo_mwdair`/`exo_cpdair` block (21) |
| | `shr_const_mod.F90` | `SCON MWDAIR CPDAIR PSTD SDAY REARTH G` + derived `OMEGA DLAPSE RDAIR ZVIR RHODAIR CPVIR` (13) |
| | `physconst.F90` | `scon cpair pstd` parameters (3); `gravit sday mwdry rearth rga ra omega rair epsilo zvir cpvir rhodair cappa cpd_on_cpv` initialisers (14) |
| | `hycoef.F90` | `ps0 psr` = `SHR_CONST_PSTD` (2) |
| | `seq_flux_mct.F90` | `albdif albdir` = `exo_albdif/dir`, twice (4) |
| | `shr_orb_mod.F90` | `dayspy ve` = `exo_porb exo_ve` (2) |
| | `ice_constants.F90` | `cp_air gravit omega radius zvir` (5) |
| | `clm_varcon.F90` | `grav rair roverg zvir cpair re` initialisers (6) |
| **Stock CESM, new SourceMods copies** | `csm_share/shr/shr_flux_mod.F90` | `g cpdair cpvir zvir` (4) |
| | `atm/.../chemistry/mozart/mo_neu_wetdep.F90` | `rearth gravit` (2) |
| | `atm/.../physics/cam/clubb_intr.F90` | `grav` (1) |
| | `atm/.../physics/cam/co2_cycle.F90` | `c_cp = (/cpair,…/)` (1) |
| | `ice/cice/src/source/ice_atmo.F90` | `cpvir` (1) |
| | `ice/cice/src/source/ice_mechred.F90` | `cp` (uses `gravit`) (1) |
| ExoRT n68equiv | — | none: ExoRT reads `exoplanet_mod` values at run time only |

So the compile-time surface is small: 8 files ExoCAM already overrides, plus 6 stock
files to copy into the shared tree (several are compiled but not used by these
configs, e.g. `clubb_intr`, `mo_neu_wetdep`, and still have to build). CAM itself
already holds `gravit`, `rearth`, `omega`, `rair`, `cappa` … as initialised
variables, not parameters.

**Not covered by the audit, to settle in the design:**
- **Initialisation order.** Components copy constants into their own saved
  variables during their init (FV dycore, CLM, CICE, the orbit). `user_nl_exo` has
  to be read and every derived value set before any component initialises: in the
  driver's init (`ccsm_comp_mod.F90`, already overridden), on the root task and
  broadcast. One executable (`cesm.exe`), so module variables are shared.
- **Write protection.** Declare the runtime values `protected` (Fortran 2003) so
  only the reader can set them.
- **Other constant contexts** the scan did not parse: `data` statements, `case`
  selectors, CPP macros. A compile of the converted tree is the real test.
- **BFB.** Values identical → results should be, but a compile-time constant
  folded by the compiler and the same value computed at run time can differ in the
  last bit (e.g. `x/g` versus `x*(1/g)`). Verify against `reg0` on all three
  configs.
- **The rundir copy and the log**, as drafted before: how `user_nl_exo` reaches
  the run directory (casemgr or the case scripts), print every resolved value in
  `atm.log`, and the checks the model has never made (`exo_pstd` vs the `ncdata`
  surface pressure; the solar file grid vs the ExoRT bundle; `exo_sday` vs the
  synchronous-rotation formula).
- **Template defects** found on the way, to fix in that work: the shipped
  `exo_n2bar` expression omits H₂, NH₃ and CO, and `exo_pstd` omits NH₃ and CO
  (wrong whenever either is non-zero).

## Phase 3 — ExoRT 2.0 hook-up (changes the answer on purpose)

1. **ExoCAM side.** `do_exo_clouds`/`do_exo_haze` go into `exoplanet_mod.F90`
   as compile-time parameters, and the ExoRT data root stays in ExoRT's compiled
   `sys_rootdir.F90`, as in `sbx_exort_v2_n1` (both move to the runtime
   `user_nl_exo` when Phase 2 is revisited). The `exo_solar_file` defaults move to
   the n84 grid (`data/stellar/*_n84.nc`).
   **Done 2026-11-23:** the gates (`do_exo_clouds = .true.`,
   `do_exo_haze = do_carma_exort`) and the n84 default
   `data/stellar/G2V_SUN_n84.nc` in all three FV configs. Checked locally with
   ExoRT's CAM compile gate built against ExoCAM's *real* `exoplanet_mod.F90` (each
   config) and shared `shr_const_mod.F90` instead of its stubs: the mainline
   (no-macro) bundle compiles for aqua, land and mixed, and the pre-change file fails
   on exactly the two missing gates. BFB for n68equiv cases (the flags are unused
   there; casemgr sets `exo_solar_file` per case).

   **Data root, decided 2026-11-23 (author):** ExoRT owns it, in its
   `src.cam.exort/sys_rootdir.F90`, as the 1.x instructions already say. ExoCAM
   does not carry it (an `exo_exort_root` in `exoplanet_mod.F90` was proposed and
   declined: ExoRT must stay a model that runs on its own). `exo_solar_file` stays a
   full path in `exoplanet_mod.F90`. casemgr will therefore need the ExoRT checkout's
   `sys_rootdir.F90` to name that checkout, not a per-case edit. Note: ExoRT
   `REFACTOR_PLAN.md` §3.11 still says CAM should take the root from `user_nl_exo`;
   that line conflicts with this decision.
   **Docs, 2026-11-23:** `general_instructions.txt` STEPs 4–5 and the three FV
   config READMEs now name `src.cam.exort`, the `sys_rootdir.F90` edit, the n84
   stellar file and `NTHRDS_ATM = 1`, so the manual path matches the configs.
2. **casemgr.**
   - an `exort` package with the n84 stem;
   - solar files from `data/stellar/`;
   - an ExoRT `refactor` checkout alongside `main`'s.

   **Done 2026-11-23 on casemgr branch `v2`** (`139e1f1`; `main` stays the
   production 1.x tool until 2.0, author's call). Registry `package_roots.exort`
   names `models/ExoCAM-refactor` and `models/ExoRT-refactor` (fresh Discover
   clones of `refactor`, made the same day); an `exort` case resolves to them and
   counts as production, a 1.x case resolves as before (its newcase script is
   byte-identical to `main`'s). An `exort` newcase gets the n84 sun from
   `data/stellar/`, `NTHRDS_ATM = 1`, and a guard that refuses an ExoRT checkout
   whose `sys_rootdir.F90` names any data root but itself; `generate` also refuses
   a 1.x `exoplanet_mod.F90` (no `do_exo_clouds`/`do_exo_haze`). `a7dfd5f`
   reverts casemgr's dead Phase 2a path (`deb5be0`). The `n84equiv` and `exort`
   n84 grids are identical, so `_n84` star products serve both.
   The Discover clone's `sys_rootdir.F90` names itself (hand edit, uncommitted,
   2026-11-23). **Still open:** no Discover checkout of casemgr `v2` to run it
   from yet.
3. **Validation (ExoRT `REFACTOR_PLAN.md` §1.1):**
   - ~~the clear-sky n68 vs exort comparison~~: dropped 2026-11-23, the RT was
     validated on ExoRT `refactor` (its `CHANGES.md`, `REGRESSION_STATUS.md`);
   - the optional physics (`-DEXORT_CO2CLD`, `-DEXORT_CARMA`,
     `-DEXORT_SRF_EMISS`) in a real case;
   - the OpenMP checks.
   - Optional: a 3-D climate comparison of one published configuration for
     the release notes.
4. **Re-record the masters on ExoRT 2.0.** The answer changes by design, so
   these become the 2.x baselines.
5. **The ExoRT build switches stay as they are for 2.0.** They are replaced
   by runtime flags during the 2.x physics work (below), not before.

## Phase 4 — Release (same day for both repos)

1. **Docs.** Rewrite `README.md`, `instructions/general_instructions.txt`,
   the config READMEs, `CLAUDE.md` and `AGENTS.md` for the single tree. Add a compatibility table to both READMEs (ExoCAM 1.x ↔
   ExoRT 1.x; 2.0 ↔ 2.0).
2. **ExoRT.** Branch `v1` from `v1.1.0` first, so ExoCAM 1.x users keep a
   working target. On `refactor` the frozen `src.cam.n68equiv` reads
   `data/kdist/n68h2o/...`, which no longer exists. Then merge `refactor` →
   `main` and tag `v2.0.0`.
3. **ExoCAM.** Tag `v2.0.0`.
4. **casemgr.** `v2` and `main` stay separate until the release (author,
   2026-11-23: `main` carries ongoing 1.x work). At the release: merge `v2` →
   `main`, delete `package_roots.exort`, tag `v2.0.0` with the other two. 1.x
   fixes made on `main` meanwhile are merged into `v2` as they land.
5. **The package: what a new user gets** (decision above). Known gaps:
   - **One front door.** A getting-started page that goes from three clones
     plus CESM 1.2.1 to a built, running case, with the compatibility table.
     Today the procedure is split across `general_instructions.txt`, three
     config READMEs, ExoRT's README and casemgr's README.
   - **casemgr's registry is the author's.** `config_registry.yaml` carries
     Discover paths and the author's IC tables in the repo. Ship a template
     (`config_registry.example.yaml`) and keep the live one out of git, as
     `sys.in` did in the precursor.
   - **Machine files.** Installing `ccsm_utils_files` into the CESM tree is
     still a manual copy (the precursor's `config_ccsm_utils.py` was a stub).
   - **Install-time edits.** ExoRT's `sys_rootdir.F90` and the Discover paths
     in the namelist templates: state them in the front door, in one place.
   - **A worked example** that a new user runs first: one aquaplanet matrix,
     shipped with casemgr (the `sbx_v2smoke` matrix is the model for it).

---

## After 2.0: the proposal work, as 2.x releases

| Release | Content | Notes |
|---|---|---|
| 2.1 | Convection PR (colleague) | Rebased onto the single tree. Ask the colleague to wait for Phase 1 rather than land it across the current versions of `convect_shallow`/`zm_conv_intr`. |
| 2.2 | CO₂ condensation fully working | ExoRT `-DEXORT_CO2CLD` becomes a runtime flag in `user_nl_exo` (`do_exo_condense_co2`). |
| 2.3 | Cleaner CARMA | ExoRT `-DEXORT_CARMA` becomes its own build directive. `SRF_EMISS` becomes standard in `camsrfexch`. |
| ongoing | Regression tiers 2–3 and benchmarks | A throughput benchmark per config (sim yr/day, timing file), using the November method. |
| later | Agent harness | Needs the stable layout and a pass/fail suite underneath it; casemgr is most of it already. |

## Open decisions

- ~~Tag ExoCAM `v1.1.0` now, or go straight to 2.0?~~ Tagged 2026-11-21.
- CARMA: an add-on directory, or part of the single tree behind a switch?
- The runtime `user_nl_exo` (Phase 2, deferred past 2.0): design, then casemgr in tandem.
- When the convection PR lands, relative to Phase 1.
- **The front door for a new user:** casemgr's matrix workflow (generate →
  make → submit), or a thin single-case command on top of it, like the
  precursor's `config_exocam.py`?
- **Backwards compatibility after 2.0** (author unsure, 2026-11-23). casemgr's
  `package_roots` works in either direction: after the release it can send
  1.x packages to a v1 checkout. Open: how long 1.x cases stay supported,
  whether ExoCAM keeps a `v1` branch as ExoRT will, and whether a case should
  snapshot its ExoRT bundle (precursor note above).
