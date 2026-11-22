# ExoCAM 2.0 — plan

**Status:** agreed 2026-11-21; Phase 0 drafted, Phase 1 diffs classified (2026-11-21). This is the working plan; edit it as decisions land.

## Decisions

- **ExoCAM 2.0 is an infrastructure release.** The model keeps its current
  shape (an overlay on CESM 1.2.1), but:
  - one SourceMods set replaces the per-config copies;
  - `user_nl_exo` replaces hand-edited `exoplanet_mod.F90`;
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
   CPP switch, it keeps each config BFB, and the parameters become
   `user_nl_exo` entries in Phase 2.
3. **Check that SourceMods for an inactive component is harmless.** Does
   `src.clm` copied into an aqua case compile into anything? Check it against
   the compsets; if not, have casemgr copy per component.
4. **Update casemgr and `general_instructions.txt`** for the new layout.

`exoplanet_mod.F90` stays per config in this phase; Phase 2 replaces it.

## Phase 2 — `user_nl_exo` replaces `exoplanet_mod.F90`

**The constraint.** The planet constants are compile-time `parameter`s in
`shr_const_mod`. Stock CESM code uses `SHR_CONST_*` (and CAM's `physconst`
values) inside `parameter` declarations and variable initialisers, and Fortran
requires both to be compile-time constants. Making a planet constant truly
runtime therefore means overriding every stock file that does this, across CAM,
CLM, CICE, docn and the driver. That set needs an audit before committing to
it. Two steps:

**2a. `user_nl_exo` as the interface, Fortran generated from it.**
- A build step (in ExoCAM, called by casemgr) reads `user_nl_exo` and writes
  the case's `exoplanet_mod.F90`.
- Users never edit Fortran, and changing a value means a rebuild, as today.
- The derived block (`exo_n2bar`, `exo_pstd`, `mw`/`cp`) is computed by the
  generator, never by hand.
- Bit-for-bit against the baselines when `user_nl_exo` reproduces each
  config's current values.

**2b. Move the runtime-safe settings to a true runtime read.** These are the
ones that never reach `shr_const_mod`:
- run options: `do_exo_rt`, `exo_rad_step`, `do_exo_rt_clearsky`,
  `do_exo_rt_spectral`, the AM flags, `exo_convect_plim`;
- `exo_solar_file`, `t_int`, `Tmax`;
- the ExoRT data root and `do_exo_clouds`/`do_exo_haze` (for Phase 3);
- the orbital parameters, if `shr_orb_mod` allows it (check).

The planet constants (gravity, radius, rotation, pressure, composition) stay
generated at compile time, unless the audit finds that the set of files to
override is small.

**Mechanics to settle:**
- **The reader.** ExoCAM reads `user_nl_exo` itself, in its own module with
  its own namelist group, much like ExoRT's `user_nl_exort`. CESM's
  build-namelist cannot register new variables without editing the stock
  `namelist_definition.xml`, which is off limits.
- **The rundir copy.** CESM copies only the `user_nl_*` files it knows, so
  decide how `user_nl_exo` reaches the rundir: casemgr, or the case's build or
  run scripts.
- **The log.** Print every resolved value in `atm.log`. Optionally write them
  as global attributes in `cam.h0`, so the output records the planet.
- **Runtime checks for today's conventions** (the model checks none of these
  now):
  - `exo_pstd` matches the `ncdata` surface pressure;
  - the solar file grid matches the ExoRT bundle;
  - `exo_sday` matches the synchronous-rotation formula when
    `do_exo_synchronous` is set.

## Phase 3 — ExoRT 2.0 hook-up (changes the answer on purpose)

1. **ExoCAM side.** `do_exo_clouds`/`do_exo_haze` and the ExoRT data root
   go into `user_nl_exo`; this replaces ExoRT's compiled `sys_rootdir.F90`,
   ExoRT §3.11. The `exo_solar_file` defaults move to the n84 grid
   (`data/stellar/*_n84.nc`).
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
   the config READMEs, `CLAUDE.md` and `AGENTS.md` for the single tree and
   `user_nl_exo`. Add a compatibility table to both READMEs (ExoCAM 1.x ↔
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
- How `user_nl_exo` reaches the rundir.
- Which constants become runtime in 2b, after the stock-file audit.
- When the convection PR lands, relative to Phase 1.
