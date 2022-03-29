module BareGroundFluxesMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: BareGroundFluxesMod
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BareGroundFluxes
!
! !INTERFACE:
  subroutine BareGroundFluxes(lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_varpar         , only : nlevgrnd
    use clm_varcon         , only : cpair, vkc, grav, denice, denh2o, istsoil, zvir
    use clm_varcon         , only : istcrop
    use clm_varctl         , only : use_c13
    use shr_const_mod      , only : SHR_CONST_RGAS, SHR_CONST_CPDAIR, SHR_CONST_KVISC
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! This routine originally had a long list of parameters, and also a reference to
! the entire clm derived type.  For consistency, only the derived type reference
! is passed (now pointing to the current column and pft), and the other original
! parameters are initialized locally. Using t_grnd instead of tg (tg eliminated
! as redundant).
! 1/23/02, PET: Added pft reference as parameter. All outputs will be written
! to the pft data structures, and averaged to the column level outside of
! this routine.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! pft's column index
    integer , pointer :: pgridcell(:)      ! pft's gridcell index
    integer , pointer :: plandunit(:)      ! pft's landunit index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(r8), pointer :: t_grnd(:)         ! ground surface temperature [K]
    real(r8), pointer :: thm(:)            ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(r8), pointer :: qg(:)             ! specific humidity at ground surface [kg/kg]
    real(r8), pointer :: thv(:)            ! virtual potential temperature (kelvin)
    real(r8), pointer :: dqgdT(:)          ! temperature derivative of "qg"
    real(r8), pointer :: htvp(:)           ! latent heat of evaporation (/sublimation) [J/kg]
    real(r8), pointer :: beta(:)           ! coefficient of conective velocity [-]
    real(r8), pointer :: zii(:)            ! convective boundary height [m]
    real(r8), pointer :: forc_u(:)         ! atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)         ! atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_th(:)        ! atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q(:)         ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rho(:)       ! density (kg/m**3)
    real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft level [m]
    real(r8), pointer :: psnsun(:)         ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha(:)         ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: z0mg_col(:)       ! roughness length, momentum [m]
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: dz(:,:)           ! layer depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: soilbeta(:)       ! soil wetness relative to field capacity
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: z0hg_col(:)       ! roughness length, sensible heat [m]
    real(r8), pointer :: z0qg_col(:)       ! roughness length, latent heat [m]
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: dlrad(:)         ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad(:)         ! upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: cgrnds(:)        ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)        ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: cgrnd(:)         ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: taux(:)          ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy(:)          ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: t_ref2m(:)       ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref1m(:)       ! 2 m height surface air temperature (Kelvin) ! mars
    real(r8), pointer :: q_ref2m(:)       ! 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: t_ref2m_r(:)     ! Rural 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: rh_ref2m_r(:)    ! Rural 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m(:)      ! 2 m height surface relative humidity (%)
    real(r8), pointer :: t_veg(:)         ! vegetation temperature (Kelvin)
    real(r8), pointer :: btran(:)         ! transpiration wetness factor (0 to 1)
    real(r8), pointer :: rssun(:)         ! sunlit stomatal resistance (s/m)
    real(r8), pointer :: rssha(:)         ! shaded stomatal resistance (s/m)
    real(r8), pointer :: ram1(:)          ! aerodynamical resistance (s/m)
    real(r8), pointer :: fpsn(:)          ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: rootr(:,:)       ! effective fraction of roots in each soil layer
    real(r8), pointer :: rresis(:,:)      ! root resistance by layer (0-1)  (nlevgrnd)	
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l            ! indices
    integer  :: filterp(ubp-lbp+1)     ! pft filter for vegetated pfts
    integer  :: fn                     ! number of values in local pft filter
    integer  :: fp                     ! lake filter pft index
    integer  :: iter                   ! iteration index
    real(r8) :: zldis(lbp:ubp)         ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(lbp:ubp)        ! displacement height [m]
    real(r8) :: zeta                   ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                     ! convective velocity [m/s]
    real(r8) :: dth(lbp:ubp)           ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                   ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(lbp:ubp)           ! diff of humidity between ref. height and surface
    real(r8) :: obu(lbp:ubp)           ! Monin-Obukhov length (m)
    real(r8) :: ur(lbp:ubp)            ! wind speed at reference height [m/s]
    real(r8) :: um(lbp:ubp)            ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(lbp:ubp)         ! relation for potential temperature profile
    real(r8) :: temp12m(lbp:ubp)       ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp11m(lbp:ubp)       ! relation for potential temperature profile applied at 1-m ! mars
    real(r8) :: temp2(lbp:ubp)         ! relation for specific humidity profile
    real(r8) :: temp22m(lbp:ubp)       ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(lbp:ubp)         ! friction velocity [m/s]
    real(r8) :: tstar                  ! temperature scaling parameter
    real(r8) :: qstar                  ! moisture scaling parameter
    real(r8) :: thvstar                ! virtual potential temperature scaling parameter
    real(r8) :: cf                     ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                    ! aerodynamical resistance [s/m]
    real(r8) :: rah                    ! thermal resistance [s/m]
    real(r8) :: raw                    ! moisture resistance [s/m]
    real(r8) :: raih                   ! temporary variable [kg/m2/s]
    real(r8) :: raiw                   ! temporary variable [kg/m2/s]
    real(r8) :: fm(lbp:ubp)            ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_pft(lbp:ubp)
    real(r8) :: z0hg_pft(lbp:ubp)
    real(r8) :: z0qg_pft(lbp:ubp)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(r8) :: www                    ! surface soil wetness [-]
! + VLH 07/01/13  ! mars
    integer :: dust
    real(r8) :: eg                     ! water vapor saturation pressure at T [Pa] !mars
    real(r8) :: e                      ! water vapor pressure [Pa]
    real(r8) :: qsatg                  ! saturated humidty [kg/kg]
    real(r8) :: degdT                  ! d(eg)/dT 
    real(r8) :: qsatgdT                ! d(qsatg)/dT 
    real(r8) :: d                      ! diffusion coefficient of H2O in CO2 [cm2/s] 
    real(r8) :: rhos                   ! saturation vapor density [g/cm3]
    real(r8) :: rho                    ! vapor density [g/cm3]
    real(r8) :: rho_ratio              ! ratio in delta density over density [unitless]
    real(r8) :: vis                    ! kinematic visocity [cm2/s]  
    real(r8) :: qflx_evap_conv         ! free convective mass loss [kg/m2/s=mm/s] 
! - VLH ! mars

!------------------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_u     => clm_a2l%forc_u
    forc_v     => clm_a2l%forc_v
!    forc_th    => clm_a2l%forc_th  ! mars
!    forc_t     => clm_a2l%forc_t   ! mars
!    forc_pbot  => clm_a2l%forc_pbot ! mars
!    forc_rho   => clm_a2l%forc_rho ! mars
!    forc_q     => clm_a2l%forc_q   ! mars

    ! Assign local pointers to derived type members (landunit-level)

    ltype      => lun%itype

    ! Assign local pointers to derived type members (column-level)

    forc_th    => ces%forc_th  ! mars
    forc_t     => ces%forc_t   ! mars
    forc_pbot  => cps%forc_pbot ! mars
    forc_rho   => cps%forc_rho   ! mars
    forc_q     => cws%forc_q   ! mars
    pcolumn    => pft%column
    pgridcell  => pft%gridcell
    frac_veg_nosno => pps%frac_veg_nosno
    dlrad  => pef%dlrad
    ulrad  => pef%ulrad
    t_grnd => ces%t_grnd
    qg     => cws%qg
    z0mg_col => cps%z0mg
    z0hg_col => cps%z0hg
    z0qg_col => cps%z0qg
    thv    => ces%thv
    beta   => cps%beta
    zii    => cps%zii
    ram1   => pps%ram1
    cgrnds => pef%cgrnds
    cgrndl => pef%cgrndl
    cgrnd  => pef%cgrnd
    dqgdT  => cws%dqgdT
    htvp   => cps%htvp
    watsat         => cps%watsat
    h2osoi_ice     => cws%h2osoi_ice
    dz             => cps%dz
    h2osoi_liq     => cws%h2osoi_liq
    frac_sno       => cps%frac_sno
    soilbeta       => cws%soilbeta

    ! Assign local pointers to derived type members (pft-level)

    taux => pmf%taux
    tauy => pmf%tauy
    eflx_sh_grnd => pef%eflx_sh_grnd
    eflx_sh_tot => pef%eflx_sh_tot
    qflx_evap_soi => pwf%qflx_evap_soi
    qflx_evap_tot => pwf%qflx_evap_tot
    t_ref2m => pes%t_ref2m
    t_ref1m => pes%t_ref1m ! mars
    q_ref2m => pes%q_ref2m
    t_ref2m_r => pes%t_ref2m_r
    rh_ref2m_r => pes%rh_ref2m_r
    plandunit => pft%landunit
    rh_ref2m => pes%rh_ref2m
    t_veg => pes%t_veg
    thm => pes%thm
    btran => pps%btran
    rssun => pps%rssun
    rssha => pps%rssha
    rootr => pps%rootr
    rresis => pps%rresis
    psnsun => pcf%psnsun
    psnsha => pcf%psnsha
    fpsn => pcf%fpsn
    forc_hgt_u_pft => pps%forc_hgt_u_pft

    ! Filter pfts where frac_veg_nosno is zero

    fn = 0
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (frac_veg_nosno(p) == 0) then
          fn = fn + 1
          filterp(fn) = p
       end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       ! Initialization variables

       displa(p) = 0._r8
       dlrad(p)  = 0._r8
       ulrad(p)  = 0._r8

       ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p) = thm(p)-t_grnd(c)
       dqh(p) = forc_q(c) - qg(c)
       dthv = dth(p)*(1._r8+zvir*forc_q(c))+zvir*forc_th(c)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p)

       ! Copy column roughness to local pft-level arrays

       z0mg_pft(p) = z0mg_col(c)
       z0hg_pft(p) = z0hg_col(c)
       z0qg_pft(p) = z0qg_col(c)

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_pft(p), um(p), obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1, niters

       call FrictionVelocity(lbp, ubp, fn, filterp, &
                             displa, z0mg_pft, z0hg_pft, z0qg_pft, &
                             obu, iter, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, temp11m, fm)

       do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
          g = pgridcell(p)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)
          z0hg_pft(p) = z0mg_pft(p)/exp(0.13_r8 * (ustar(p)*z0mg_pft(p)/SHR_CONST_KVISC)**0.45_r8)
          z0qg_pft(p) = z0hg_pft(p)
          thvstar = tstar*(1._r8+zvir*forc_q(c)) + zvir*forc_th(c)*qstar
          zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0._r8) then                   !stable
             zeta = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                                      !unstable
             zeta = max(-100._r8,min(zeta,-0.01_r8))
             wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p) + wc*wc)
          end if
          obu(p) = zldis(p)/zeta
       end do

    end do ! end stability iteration

     do j = 1, nlevgrnd
       do f = 1, fn
          p = filterp(f)
          rootr(p,j) = 0._r8
          rresis(p,j) = 0._r8
        end do
     end do

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Determine aerodynamic resistances

       ram     = 1._r8/(ustar(p)*ustar(p)/um(p))
       rah     = 1._r8/(temp1(p)*ustar(p))
       raw     = 1._r8/(temp2(p)*ustar(p))
       raih    = forc_rho(c)*cpair/rah

       ! Soil evaporation resistance
       www     = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
       www     = min(max(www,0.0_r8),1._r8)

       !changed by K.Sakaguchi. Soilbeta is used for evaporation
       if (dqh(p) .gt. 0._r8) then   !dew  (beta is not applied, just like rsoil used to be)
          raiw    = forc_rho(c)/(raw)
       else
       ! Lee and Pielke 1992 beta is applied
          raiw    = soilbeta(c)*forc_rho(c)/(raw)
       end if

       ram1(p) = ram  !pass value to global variable

       ! Output to pft-level data structures
       ! Derivative of fluxes with respect to ground temperature

       cgrnds(p) = raih
       cgrndl(p) = raiw*dqgdT(c)
       cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       taux(p)          = -forc_rho(c)*forc_u(g)/ram
       tauy(p)          = -forc_rho(c)*forc_v(g)/ram
       eflx_sh_grnd(p)  = -raih*dth(p)
       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       qflx_evap_soi(p) = -raiw*dqh(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)

       
       ! mars  
       ! +VLH: add free convective mass loss [eqn 1 Ingersoll 1970], cgs units!
       ! E = 0.17*delta_nu*denisty*D*(delta_rho/rho * g/v^2)^(1/3)
       ! delta_nu = difference b/w water vapor concentration by mass of gas at
       ! evaporatiing surface and gas away from the surface (at surface - 
       ! saturated) = delta_nu = rhos/rho 
       ! d = diffusion coefficient of water vapor in carbon dioxide [cm2/s]
       ! rho, rhos = vapor and saturation vapor densities [g/cm3]
       ! vis = kinematic visocity [cm2/s] = dyn vis/air density, dyn visc calculated 
       ! from sutherlands equation for co2 [uPa s = 1e-5 dyn/cm2 s] (To=293.15, vo=14.8, Co=240) 
       ! rho_ratio = difference between the density of ambient gas and gas at 
       ! the surface, divided by density of gas at the surface [unitless] = 
       ! (mc-mw)es/ (mcPo-(mc-mw)es) where mc = 44, mw=18 
       ! need to convert cgs units to model SI [g/cm2/s] * 10 = [kg/m2/s=mm/s]
       ! mars wolf.. there are a bunch of imbeed hard numbers here
       if (frac_sno(c) > 0.) then
          call QSat(forc_t(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
          ! dynamic viscosity (sutherlands formula) / air denisty                                                                     
           vis = 1.e-5 * 14.8_r8 * ((293.15_r8 + 240._r8) / (forc_t(c) + 240._r8))*(forc_t(c)/293.15_r8)**(3./2.)*1./(1.e-3*forc_rho(c))
           rho = forc_q(c) * forc_rho(c) * (1.e-3)   ! water vapor air density [g/cm3]
           rhos = qsatg * forc_rho(c) * (1.e-3)      ! saturation density [g/cm3]
           d = 32._r8*t_grnd(c)**(3./2.) / (10._r8*forc_pbot(c))   ! [cm2/s]
           rho_ratio = (26._r8 * eg) / (44._r8  * forc_pbot(c) - 26._r8 * eg)
           qflx_evap_conv = frac_sno(c) * 1.7_r8 * rhos/rho * rho * d * (rho_ratio * (100._r8*grav)/vis**2.)**(1./3.)
           qflx_evap_soi(p) = qflx_evap_soi(p) + qflx_evap_conv
        endif

        qflx_evap_tot(p) = qflx_evap_soi(p)

        if (qflx_evap_soi(p) > .5) then
                write(6,*)'qflx high in bgf_mod:qflx=',qflx_evap_soi(p), &
                        'raiw=',raiw,'dqh=',dqh(p),'forc_q=',forc_q(c),'qg=',qg(c)
        endif
        if (eflx_sh_grnd(p) > 3000.) then
                write(6,*)'high sh_grnd in bgf_mod:sh=',eflx_sh_grnd(p), &
                        'raih=',raih,'dth=',dth(p),'thm=',thm(p),'tgrnd=',t_grnd(c)
        endif
        if (eflx_sh_grnd(p) < -100.) then
                write(6,*)'low sh_grnd in bgf_mod:sh=',eflx_sh_grnd(p), &
                        'raih=',raih,'dth=',dth(p),'thm=',thm(p),'tgrnd=',t_grnd(c)
        endif

! \mars

       ! 1 m height air temperature

        t_ref1m(p) = thm(p) + temp1(p)*dth(p)*(1./temp11m(p) - 1./temp1(p))  ! mars

       ! 2 m height air temperature

       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity

       !q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))
       q_ref2m(p) = forc_q(c) !set 2-m specific humidity equal to grid box specific humidity.

       ! 2 m height relative humidity
                                                                                
       call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         rh_ref2m_r(p) = rh_ref2m(p)
         t_ref2m_r(p) = t_ref2m(p)
       end if

       ! Variables needed by history tape

       t_veg(p) = forc_t(c)
       btran(p) = 0._r8
       cf = forc_pbot(c)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
       rssun(p) = 1._r8/1.e15_r8 * cf
       rssha(p) = 1._r8/1.e15_r8 * cf

       ! Add the following to avoid NaN

       psnsun(p) = 0._r8
       psnsha(p) = 0._r8
       fpsn(p) = 0._r8
       pps%lncsun(p) = 0._r8
       pps%lncsha(p) = 0._r8
       pps%vcmxsun(p) = 0._r8
       pps%vcmxsha(p) = 0._r8
       ! adding code for isotopes, 8/17/05, PET
       pps%cisun(p) = 0._r8
       pps%cisha(p) = 0._r8
       if (use_c13) then
          pps%alphapsnsun(p) = 0._r8
          pps%alphapsnsha(p) = 0._r8
          pepv%rc13_canair(p) = 0._r8
          pepv%rc13_psnsun(p) = 0._r8
          pepv%rc13_psnsha(p) = 0._r8
          pc13f%psnsun(p) = 0._r8
          pc13f%psnsha(p) = 0._r8
          pc13f%fpsn(p) = 0._r8
       end if

    end do

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
