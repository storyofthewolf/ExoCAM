# CO2 Condensation Bug Report
**Reviewer:** Claude Opus 4.7  
**Date:** 2026-05-13  
**Scope:** `cesm1.2.1/configs/experimental/co2condense/SourceMods/` only  
**Branch:** `co2_claude_diag`

---

## Context and Symptoms

The CO2 condensation configuration (`experimental/co2condense`) is a Mars-targeted experimental build of ExoCAM. Two primary failure modes were reported:

1. **CO2 mass is not correctly removed from the atmosphere** during condensation and re-evaporation — total atmospheric mass does not decrease as expected during atmospheric collapse.
2. **Physically implausible amounts of CO2 condensation** occur over high-altitude terrain (Olympus Mons analog) when Mars topography is used.

The model assumes a **pure CO2 atmosphere** (`exo_co2vmr ≡ 1.0`, compile-time constant) to avoid local composition changes during condensation. This simplifies some issues but means the local CO2 partial pressure equals total pressure at all times.

---

## Bug Inventory

### A. Atmospheric Mass Budget (Primary — causes both symptoms)

#### A1 — Atmospheric condensation never reduces `pdel` or `ps`  
**File:** `src.cam/exo_condense_mod.F90`  
**Severity:** Critical  
**Status:** Fixed in commit `eea21c9`

`exo_condense_co2` reports its result only through `ptend%q` (CLDICE_CO2 tracer tendency) and `ptend%s` (temperature tendency). The mass that changes phase from bulk gas to ice is never subtracted from `state%pdel`. The FV dynamics therefore sees no column mass change, no pressure gradient develops, and there is no negative feedback to limit condensation at any grid cell. This is the root cause of both reported symptoms.

**Fix applied:** After the substep loop and before setting tendencies, subtract the net condensed mass from each layer:
```fortran
state%pdel(i,k)  = state%pdel(i,k) - ptend%q(i,k,ixcldice_co2) * dtime * state%pdel(i,k)
state%rpdel(i,k) = 1.0_r8 / state%pdel(i,k)
```
Then update `state%ps` as `sum(state%pdel(i,1:pver))`.

#### A2 — `dp_coupling.F90` depletion branch averages layer pressures  
**File:** `src.cam/dp_coupling.F90`, lines 780–798 (original)  
**Severity:** High  
**Status:** Fixed in commit `eea21c9`

The Mars-tagged block in `p_d_coupling` handles the surface-direct condensation path by adjusting `delpxy` using `cam_in%co2_mass_change`. When a single layer's `delpxy` was insufficient, the original code spread the remaining mass deficit uniformly across `nlay` layers (`delpxy = delptemp/nlay`). This destroys the hydrostatic pressure profile, producing spurious vertical motion on the next dynamics step.

**Fix applied:** Replace layer-averaging with a bottom-up peel: subtract from layer `km`, if exhausted move to `km-1`, etc.

#### A3 — Mars CO2 mass adjustment absent from non-`local_dp_map` branch  
**File:** `src.cam/dp_coupling.F90`  
**Severity:** High  
**Status:** Fixed in commit `eea21c9`

The Mars `delpxy` adjustment existed only in the `local_dp_map` branch (single-PE or shared-memory runs). The block-transpose MPI decomposition branch (lines 819–899) had no equivalent code. Mass conservation therefore depended on PE count, which is a latent but serious portability bug.

**Fix applied:** Added equivalent bottom-up peel logic after the block-transpose unpack, iterating over `cam_in(lchnk)%co2_mass_change` with chunk/lon/lat indexing.

#### A4 — Potential double-subtraction for snowfall path (unresolved)  
**File:** `src.cam/exo_condense_mod.F90` + `src.cam/dp_coupling.F90`  
**Severity:** Medium — may or may not manifest depending on test case  
**Status:** Open — needs budget diagnostic to confirm

`exo_condense_co2` now subtracts airborne ice mass from `pdel` via `ptend%q`. Ice that subsequently sediments to the surface accumulates in `co2_snowfall` → `cam_out%co2srf_snow` → CLM's `co2dp` → `cam_in%co2_mass_change`, which dp_coupling then also subtracts from `delpxy`. If the sedimented ice was already removed from `pdel` in the physics step, this is a double subtraction. Magnitude depends on the ratio of direct-surface condensation to snowfall. The `CO2_COL_TOT` diagnostic will reveal any drift.

---

### B. Olympus Mons / High-Terrain Runaway

#### B1 — Top-layer (k=1) loop uses uninitialized loop bound  
**File:** `src.cam/exo_condense_mod.F90`, line 507 (original)  
**Severity:** Medium  
**Status:** Fixed in commit `eea21c9`

```fortran
do i=i,ncol   ! was "i=i" — uninitialized lower bound
```
After execution of the `k=2,pver` loop, `i` equals `ncol`, so the k=1 loop ran for only one column (or zero). Also, `cld_tmp` and `tmid_tmp` were not updated for k=1 even when condensation occurred, and no diagnostics were accumulated for k=1.

**Fix applied:** Corrected to `do i=1,ncol`; added `cld_tmp`, `tmid_tmp`, and diagnostic accumulation for k=1.

#### B2 — `pot_tend`/`heat_tend` use net `sed_tend` instead of incoming flux  
**File:** `src.cam/exo_condense_mod.F90`, lines 523, 527 (original)  
**Severity:** High — direct cause of local runaway condensation  
**Status:** Fixed in commit `eea21c9`

The Forget (1998) potential-energy and sensible-heat corrections are meant to account for ice falling *into* a layer. They were multiplied by `sed_tend(i,k)`, which is the **net** flux divergence `(fxcld(k) - fxcld(k+1)) / (dt * pdel)`. For layers that are net losers (e.g., the bottom layer over high terrain that dumps ice onto the ground), `sed_tend` is negative, making `pot_tend` and `heat_tend` negative. Subtracting negative corrections increases `cld_tend`, creating a positive feedback that continuously drives more condensation in exactly the layers most prone to Olympus-type runaway.

**Fix applied:** `exo_cloud_sediment_tend` now returns a new output `fxcld_in(i,k)` — the incoming flux through the top interface of layer k, `fxcld(i,k) / (dt * pdel)`. The pot/heat corrections use this instead of `sed_tend`.

#### B3 — `tcond` computed once before substep loop using initial `pmid`  
**File:** `src.cam/exo_condense_mod.F90`, lines 450–456 (original)  
**Severity:** Medium  
**Status:** Open

`tcond(i,k)` is computed before the substep loop using `state%pmid(i,k)` — the pressure at the start of the physics timestep. Now that A1 is fixed and `pdel` decreases as condensation proceeds, the actual mid-layer pressure during substeps is lower than `state%pmid`. The frost-point temperature therefore does not drop as the column thins, removing a physical negative feedback. Over Olympus Mons, where thin layers are cold and the local pressure is already low, this makes the model continue condensing at the same `tcond` even as atmospheric mass is being depleted.

**Recommended fix:** Recompute `tcond` each substep using the current `tmid_tmp` and a pressure estimated from the evolving `pdel`, or at minimum recompute once after the `pdel` update.

#### B4 — `pmid`, `pint` not updated after `pdel` change  
**File:** `src.cam/exo_condense_mod.F90`  
**Severity:** Low-medium  
**Status:** Open

After the P1 fix updates `state%pdel` and `state%ps`, the derived fields `state%pmid`, `state%pint`, `state%lnpmid`, etc. are stale for the remainder of that physics column step. Parameterizations called after `exo_condense_tend` (radiation, boundary layer) see the pre-condensation pressure. This is a one-timestep lag that accumulates. Consider calling a pressure-recompute utility after the `pdel` update, or moving the pdel adjustment to after all other parameterizations.

---

### C. CLM Surface Condensation (`src.clm/exo_clm_condense.F90`)

#### C1 — Double-counted latent energy on full sublimation  
**File:** `src.clm/exo_clm_condense.F90`, lines 309–320 (original)  
**Severity:** Medium  
**Status:** Fixed in commit `eea21c9`

When all available frost sublimes in one timestep (the `co2dp < 0` branch), the latent energy was sent **both** to the atmosphere via `eflx_lh_co2` *and* used to cool the soil via `temp_change = L / (rho_soil * cp_soil * thick_layer1) * co2dp_ip`. The same energy was counted twice, overcooling the soil and over-warming the atmosphere.

Additionally, `thick_layer1` included the CO2 frost depth, but the heat capacity used was pure soil (`SHR_CONST_SOILBD * SHR_CONST_SOILCP`), making the implied thermal mass too large.

**Fix applied:** In the full-sublimation branch, removed the `temp_change` soil cooling. All latent energy goes to the atmosphere via `eflx_lh_co2`. The soil cools on the next timestep via conduction once `t_grnd` is free to evolve.

#### C2 — `co2_mass_change` overwritten rather than summed in pft loop  
**File:** `src.clm/exo_clm_condense.F90`  
**Severity:** Low (only affects multi-pft columns)  
**Status:** Open

`co2_mass_change(c)` is a column-level variable reset and set inside a pft loop. For columns with more than one pft in the condense filter, each pft iteration overwrites the column value rather than accumulating weighted contributions. Only the last pft's contribution reaches the atmosphere. For bare-ground Mars configurations this is typically one pft per column and does not bite.

#### C3 — `t_grnd` unconditionally set to `tcond` when frost is present  
**File:** `src.clm/exo_clm_condense.F90`  
**Severity:** Low  
**Status:** Open

Whenever frost exists or condensation occurs, `t_grnd` is clamped to `tcond` regardless of frost layer thickness. For very thin frost the underlying soil is warmer and will conduct heat into the frost; the artificial clamp prevents this. Combined with A1 (now fixed), this contributed to runaway frost at high terrain by continuously presenting the surface at the condensation temperature with no self-limiting thermal inertia.

---

### D. Additional Notes (Non-bugs)

- **`SHR_CONST_CPDAIR` usage throughout `exo_condense_mod.F90`** — initially flagged as a bug (Earth's cp vs. CO2's cp). Confirmed non-issue: `SHR_CONST_CPDAIR` is overridden in `src.share/shr_const_mod.F90` to the CO2 value for this build.
- **`co2srf_snow` units (kg/m² vs. kg/m²/s)** — `cam_out%co2srf_snow` is registered as a flux (`a2x_fluxes`) but populated with kg/m² accumulated over `dtime`. CLM consumes it as kg/m². This is intentional: `cpl_dt == atm_dt` means no coupler time-averaging occurs. A comment has been added at both ends to flag this reliance. Do not change `cpl_dt` without revisiting this.
- **`nlay` declared as `integer` initialized with `1.` (real literal)** — cosmetic Fortran issue in dp_coupling, not a runtime bug.

---

## Files Changed in This Session

| File | Changes |
|------|---------|
| `src.cam/exo_condense_mod.F90` | P1: pdel/ps update; P2: fxcld_in for pot/heat; P3: k=1 typo; budget diagnostics |
| `src.cam/dp_coupling.F90` | P1: depletion branch rewrite; non-local_dp_map mirror |
| `src.clm/exo_clm_condense.F90` | P4: double-counted latent energy |
| `src.cam/physpkg.F90` | Updated call to `exo_condense_diag_calc` |

## Open Items (Next Session)

1. Run model and inspect `CO2_COL_TOT` globally — confirm mass budget closes
2. Investigate potential double-subtraction (A4) if budget drifts
3. Recompute `tcond` each substep using evolving pressure (B3)
4. Update `pmid`/`pint` after `pdel` change (B4)
5. Fix C2 (pft loop) if running with vegetated land surfaces
