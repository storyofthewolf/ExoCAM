# ExoCAM 2.0 — plan

**Status:** agreed 2026-11-21. Phase 0 baselines recorded for aqua and mixed; Phase 1 done and BFB for aqua and mixed (2026-11-22). Land baseline (Ben1 hybrid) recorded and BFB (2026-11-23). Phase 2 deferred past 2.0 (2026-11-23): `exoplanet_mod.F90` stays compile-time; 2a reverted. This is the working plan; edit it as decisions land.

## Decisions

- **ExoCAM 2.0 is an infrastructure release.** The model keeps its current
  shape (an overlay on CESM 1.2.1), but:
  - one SourceMods set replaces the per-config copies;
  - ~~`user_nl_exo` replaces hand-edited `exoplanet_mod.F90`~~ (deferred past 2.0,
    2026-11-23: see Phase 2);
  - ExoRT 2.0 (`src.cam.exort`) is hooked up;
  - regression baselines guard every step.
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
3. **Validation (ExoRT `REFACTOR_PLAN.md` §1.1):**
   - the clear-sky n68 vs exort comparison;
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
