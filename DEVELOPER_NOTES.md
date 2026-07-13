# Developer Notes — CO2 Condensation Debugging Session
**Date:** 2026-05-13  
**Branch:** `co2_claude_diag` (branched from `main`)  
**Reviewer/Developer:** Eric Wolf + Claude Opus 4.7  

---

## Session Goals

1. Diagnose why CO2 mass balance is not correctly maintained during condensation and re-evaporation in the `experimental/co2condense` configuration.
2. Diagnose why physically implausible CO2 condensation occurs over high-altitude terrain (Olympus Mons) when using Mars topography.

---

## Files Reviewed

All review was confined to `cesm1.2.1/configs/experimental/co2condense/SourceMods/`:

- `src.cam/exo_condense_mod.F90` — atmospheric CO2 condensation, sedimentation, tendencies
- `src.cam/dp_coupling.F90` — physics-to-dynamics coupler, Mars `delpxy` adjustment
- `src.clm/exo_clm_condense.F90` — CLM surface CO2 frost accumulation/sublimation
- `src.share/shr_condense_mod.F90` — CO2 frost-point curve
- `src.share/exoplanet_mod.F90` — compile-time planet/atmosphere parameters
- `src.cam/camsrfexch.F90`, `src.cam/atm_comp_mct.F90`, `src.clm/lnd_comp_mct.F90` — coupler field definitions and packing
- `src.cam/physpkg.F90` — call site for condense routines
- `src.drv/seq_flds_mod.F90` — MCT field registration

---

## Key Design Context Established During Review

- The model is currently configured as **pure CO2** (`exo_co2vmr = 1.0`, compile-time constant in `exoplanet_mod.F90`). Local composition feedback is intentionally disabled for this development phase — this simplifies the physics but means CO2 partial pressure always equals total pressure.
- `SHR_CONST_CPDAIR` is overridden to the CO2 value in `shr_const_mod.F90`. Do not "fix" uses of `SHR_CONST_CPDAIR` in `exo_condense_mod.F90` — it is already correct for CO2.
- `cam_out%co2srf_snow` is in kg/m² (accumulated over `dtime`), not kg/m²/s, despite being registered as `a2x_fluxes`. This is intentional: `cpl_dt == atm_dt` so no coupler time-averaging occurs. Both ends have been annotated with comments. Do not change `cpl_dt` without revisiting this field.

---

## Changes Made

### Commit 1: `eea21c9` — "fix CO2 condensation mass budget and numerical issues (P1-P4)"

#### P1 — Close atmospheric mass budget (`exo_condense_mod.F90`, `dp_coupling.F90`)

**Problem:** `exo_condense_co2` computed condensation tendencies but never subtracted the condensed mass from `state%pdel`. The FV dynamics received an unchanged column mass every step regardless of how much CO2 had condensed. No pressure gradient developed in response to condensation, so there was no negative feedback to limit runaway frost accumulation anywhere, and global atmospheric collapse was impossible.

**Fix in `exo_condense_mod.F90`:** After the substep loop, compute `ptend%q` as before, then immediately update:
```fortran
state%pdel(i,k)  = state%pdel(i,k) - ptend%q(i,k,ixcldice_co2) * dtime * state%pdel(i,k)
state%rpdel(i,k) = 1.0_r8 / state%pdel(i,k)
```
Then recompute `state%ps = sum(state%pdel(i,1:pver))`. This must happen before `ptend` is applied by `physics_update` in `physpkg.F90` (it currently does — the update is inside `exo_condense_co2`, called before `physics_update`).

**Fix in `dp_coupling.F90` (local_dp_map branch, lines 772–803 original):** The Mars block that adjusts `delpxy` for CLM's `co2_mass_change` had a depletion fallback that averaged layer pressures across `nlay` layers (`delpxy = delptemp/nlay`). This destroyed the hydrostatic pressure profile. Replaced with a bottom-up peel: subtract from layer `km`; if that layer is exhausted, continue to `km-1`, etc. Sign convention is handled correctly for both condensation (positive `co2_mass_change`, subtract from `delpxy`) and sublimation (negative `co2_mass_change`, add back).

**Fix in `dp_coupling.F90` (non-local_dp_map branch):** The Mars adjustment was entirely absent from the block-transpose MPI path (lines 819–899 original). Added equivalent bottom-up peel logic after the unpack loop, iterating over chunks with `get_lon_all_p`/`get_lat_all_p` indexing. This makes mass conservation decomposition-independent.

**Known open issue — potential double-subtraction (see A4 in bug report):** Ice that sediments to the surface goes through `ptend%q` (subtracted from `pdel` in physics) → `co2_snowfall` → `cam_out%co2srf_snow` → CLM → `cam_in%co2_mass_change` → dp_coupling (subtracted from `delpxy` again). This may double-count the snowfall path. Use `CO2_COL_TOT` to diagnose; if the budget drifts negative, this is the cause.

#### P2 — Fix Forget-1998 pot/heat correction sign (`exo_condense_mod.F90`)

**Problem:** In the `do_forget1998` condensation block (lines 522–527 original), `pot_tend` and `heat_tend` were multiplied by `sed_tend(i,k)`, which is the **net** flux divergence `(fxcld(k) - fxcld(k+1)) / (dt*pdel)`. For layers that are net losers (e.g., the bottom layer over high terrain that is dumping its entire ice content to the ground), `sed_tend < 0`. This made `pot_tend` and `heat_tend` negative, and since `cld_tend = cond_tend - pot_tend - heat_tend`, subtracting negative corrections increased `cld_tend`. This was a positive feedback driving more condensation in the exact layers responsible for Olympus-type runaway.

**Fix:** `exo_cloud_sediment_tend` now returns a new output argument `fxcld_in(pcols,pver)` — the incoming downward flux through the top interface of each layer:
```fortran
fxcld_in(i,k) = fxcld(i,k) / (dtime * ext_pdel(i,k))   ! [kg/kg/s]
```
The pot/heat corrections in the Forget-1998 block now use `fxcld_in(i,k)` instead of `sed_tend(i,k)`. The call to `exo_cloud_sediment_tend` was updated to pass this new argument.

#### P3 — Fix top-layer (k=1) loop typo (`exo_condense_mod.F90`, line 507 original)

**Problem:** `do i=i,ncol` — the lower bound `i` was the loop variable itself, which held its value from the previous loop iteration (typically `ncol` at the end of the `k=2,pver` loop). So the k=1 loop executed for only one column, not all `ncol` columns. Additionally, `cld_tmp(i,1)` and `tmid_tmp(i,1)` were never updated even when condensation occurred at k=1, and no diagnostics were accumulated for k=1.

**Fix:** Corrected to `do i=1,ncol`. Added `cld_tmp`, `tmid_tmp` updates and `sed_tend_diag`, `cond_tend_diag`, `temp_tend_diag` accumulation for k=1.

#### P4 — Fix double-counted latent energy on full sublimation (`exo_clm_condense.F90`)

**Problem:** In the branch where `co2dp < 0` (all frost sublimed in one step), the latent energy was both sent to the atmosphere (`eflx_lh_co2 = L/dtime * co2dp_ip`) and used to cool the soil (`temp_change = L / (rho*cp*thick) * co2dp_ip; t_grnd -= temp_change`). The same energy appeared twice in the surface energy balance. Additionally, `thick_layer1` included the frost depth but the heat capacity was pure soil, making the implied thermal mass inconsistently larger than reality.

**Fix:** Removed `temp_change` soil cooling from the full-sublimation branch. All latent energy goes to the atmosphere. The soil cools on the next step via conduction once `t_grnd` is free to evolve normally.

---

### Commit 2: `05695cf` — "add CO2 column mass budget diagnostics"

Added four new CAM history fields to `exo_condense_diag_calc` to diagnose the mass budget:

| Field | Units | Description |
|-------|-------|-------------|
| `CO2_COL_GAS` | kg/m² | Bulk atmospheric gas column: `sum(pdel/g)` |
| `CO2_COL_ICE` | kg/m² | Airborne CLDICE_CO2: `sum(q_ice * pdel/g)` |
| `CO2_COL_SRF` | kg/m² | Surface frost from CLM: `cam_in%co2dp` |
| `CO2_COL_TOT` | kg/m² | Sum of all three — globally conserved quantity |

`exo_condense_diag_calc` now takes `cam_in` as a third argument (type `cam_in_t`). The call site in `physpkg.F90` was updated accordingly.

**How to use these fields:**
- Global mean `CO2_COL_TOT` should be flat after P1. Slope reveals leak rate and sign.
- `CO2_COL_GAS` at Olympus Mons should now decrease when frost forms. Runaway means it stays at initial value despite frost accumulation.
- `CO2_COL_SRF` at polar caps: should show seasonal cycle without secular drift.
- `CO2_COL_ICE`: should be small and transient everywhere. Large sustained values indicate inefficient sedimentation.

---

---

# Session 2 — 2026-07-12 (Claude Fable 5)

## A4 static trace: sign was OPPOSITE the original hypothesis

The feared double-subtraction of the snowfall path does not exist. In
`exo_clm_condense.F90`, snowfall is added to the frost reservoir
(`co2dp(c) = co2dp_i + forc_co2srf_snow(c)`, before the energy-balance logic)
but `co2_mass_change` is set **only** from the energy-balance term
(`dtime/L * co2ice_lh * wtcol`). dp_coupling therefore never subtracts snowfall
from `delpxy`.

The real bug was on the physics side: the P1 update used the **net** tracer
change `ptend%q*dtime`, which includes sedimentation transport. Layers that
lost ice by sedimentation (Δq < 0 from transport) had their `pdel` *refunded*
gas that never sublimed. Column-integrated, physics removed only
`(column condensation − snowfall)` while the snowfall mass simultaneously
appeared in CLM's frost — **mass creation at the global snowfall rate;
`CO2_COL_TOT` would have drifted positive.**

## Changes made (all in `src.cam/exo_condense_mod.F90`)

### A4 — pdel update now uses phase conversion only

New accumulator `cond_accum(i,k)` [kg/kg] sums `cld_tend*subtime` (gas↔ice
conversion, including Forget pot/heat sublimation of falling ice) across
substeps, in all three branches (Forget k=1, Forget k≥2, LMDG). It excludes
`sed_tend` transport. When the negative-cloud clamp fires, the clamped deficit
is backed out of `cond_accum` (that ice never existed, so the gas phase is not
credited for subliming it). Final update: `pdel = pdel*(1 − cond_accum)`.

Budget closure per column per step:
- physics: removes net atmospheric gas→ice conversion, in the layers where it occurred
- dp_coupling: removes direct surface exchange only (CLM energy balance; snowfall excluded)
- snowfall: pure ice transport atm→surface, already paid for in gas at the condensing layer

### B3 — tcond/tnuc recomputed every substep from evolving pressure

`pdel_tmp`/`pmid_tmp` are rebuilt at the end of each substep from `cond_accum`
(running interface-pressure integration starting at `state%pint(i,1)`);
`tcond`/`tnuc` are recomputed from `pmid_tmp` at the top of each substep. This
restores the negative feedback: as a column thins, its frost point falls and
condensation self-limits. Guard: `cond_accum` clamped at 0.99 with a warning
(prevents non-positive pressures if a layer would totally deplete).
Sedimentation calls deliberately stay on the fixed original `state%` pressure
basis, consistent with the fixed-basis mixing-ratio bookkeeping.

### B4 — derived pressure fields rebuilt after the pdel update

`pint`, `lnpint`, `pmid`, `lnpmid` rebuilt with the same formulas
`d_p_coupling` uses (`pint(k+1)=pint(k)+pdel(k)`, `pmid=0.5*(pint(k)+pint(k+1))`),
`pint(1)` unchanged, and `ps = pint(pverp)` — which also fixes the original
P1 `ps = sum(pdel)` omitting ptop.

Two facts discovered while verifying downstream consistency:
- `exo_condense_tend` is called from **tphysbc** (not tphysac as Session 1
  notes said), immediately before radiation — so B4 directly changes what
  radiation sees the same step.
- `physics_update` (physpkg.F90:2156) calls `geopotential_dse` with
  `lnpint/lnpmid/pint/pmid/pdel/rpdel`, so `t`, `zi`, `zm` are automatically
  reconciled with the new pressures right after our routine returns.
  `cam_export` runs later in tphysbc, so `forc_pbot` (used for CLM's surface
  frost point) is also post-condensation.

## Open Items for Next Session

### High priority

1. **Run model, check `CO2_COL_TOT`** — global mean should now be flat. If a
   drift remains, suspects are (in order): coupler lag between `co2srf_snow`
   delivery and `co2_mass_change` return (one-step offset is fine, secular
   drift is not); CLM energy-balance mass terms; the `cond_accum` 0.99 clamp
   firing (grep logs for `cond_accum > 0.99 clamped`).
2. **Olympus Mons column** — `CO2_COL_GAS` should decrease as frost forms and
   condensation should self-limit (B3). Compare `CLDICE_CO2_COND_TEND` before/after.

### Low priority

3. **Fix C2 (pft loop)** — `co2_mass_change` overwritten not summed in pft loop; only relevant for multi-pft columns (vegetated land surfaces). Not an issue for bare-ground Mars runs.
4. **Fix C3 (`t_grnd = tcond` clamp)** — unconditional clamping of surface temperature to frost point when frost is present removes thermal inertia. Worth revisiting for thin-frost sublimation at sunrise.
5. **`check_energy_chng` noise** — column mass changes are not representable in CAM's fixed-mass energy check; expect logged imbalances after `exo_condense_tend`. Log-only, but worth a comment or a bypass if it becomes distracting.

---

## Architecture Notes for Future Reference

### Mass flow diagram for CO2 in this model

```
Atmospheric bulk gas (pdel)
    |  condensation (exo_condense_co2)
    v
CLDICE_CO2 tracer (ptend%q)       <-- pdel reduced by P1 fix
    |  sedimentation (exo_cloud_sediment_tend)
    v
co2_snowfall [kg/m2]
    |  cam_out%co2srf_snow → coupler → forc_co2srf_snow
    v
CLM surface frost (co2dp)         <-- co2_mass_change sent back to atm
    |  cam_in%co2_mass_change → dp_coupling → delpxy adjustment
    v
Atmospheric bulk gas (delpxy → next dynamics step)

Direct surface condensation (no CLDICE_CO2 path):
Atmospheric bulk gas → [CLM energy balance] → co2dp → co2_mass_change → delpxy
```

### Call sequence for CO2 condensation (within one `tphysbc` call — Session 2 correction: this is tphysbc, not tphysac; radiation follows in the same call)

```
tphysbc
  └─ exo_condense_tend (physpkg.F90:2154)
       └─ exo_condense_co2
            ├─ substep loop (20x)
            │    ├─ [B3] recompute tcond/tnuc from evolving pmid_tmp
            │    ├─ calc_co2_reff
            │    ├─ exo_cloud_sediment_vel
            │    ├─ exo_cloud_sediment_tend  [returns fxcld_in]
            │    ├─ Forget-1998 condensation [uses fxcld_in for pot/heat;
            │    │                            accumulates cond_accum (A4)]
            │    └─ [B3] rebuild pdel_tmp/pmid_tmp from cond_accum
            ├─ set ptend%q, ptend%s
            ├─ [P1/A4] state%pdel *= (1 - cond_accum); rpdel
            └─ [B4] rebuild pint, lnpint, pmid, lnpmid; ps = pint(pverp)
  └─ physics_update (applies ptend; geopotential_dse reconciles t, zi, zm
                     with the updated pressure fields)
  └─ exo_condense_diag_calc (physpkg.F90:2181)
       └─ [outputs CO2_COL_GAS/ICE/SRF/TOT]
  └─ radiation (sees post-condensation pressures)
  └─ cam_export (surface forcing, incl. forc_pbot -> CLM frost point)

p_d_coupling (dp_coupling.F90, called from stepon.F90)
  └─ copy phys_state%pdel → delpxy  [includes P1 update]
  └─ [Mars block] adjust delpxy for cam_in%co2_mass_change
  └─ p_d_adjust → recompute pe, ps, pk from delpxy
```
