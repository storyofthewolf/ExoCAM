
module exo_radiation_mod

!---------------------------------------------------------------------
! Purpose:
!
! Contains radiative transfer algorithm
!
! Adapts correlated K radiative transfer code for use with CAM/WACCM.
! Thie code uses the two stream radiative transfer method described in 
! Toon et al (1989).  Quadrature is used for shortwave, hemispheric mean
! is used for longwave.  Gas phase optical depths are calculate using a 
! correlated K-distribution method (Mlawer, 1997), with overlapping bands 
! treated via an amount weighted scheme (Shi et al, 2009). 
!
! Cloud optics treated using mie scattering for both liquid and ice clouds.
! Cloud overlap is treated using Monte Carlo Independent Column Approximation.
!
! Water vapor and CO2 continuum from MTCKD
!
!
! Revision history
! September 2010, E. T. Wolf, R. Urata CAM3
! March     2014, E. T. Wolf --- decoupled solar and IR streams
!                            --- merged with CESM1.2.1
! February  2015, E.T. Wolf  --- new build
!---------------------------------------------------------------------
!  
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_G, SHR_CONST_PSTD, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, SHR_CONST_CDAY, &
                              SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
			      SHR_CONST_LOSCHMIDT
  use physconst,        only: scon,mwn2, mwco2, mwch4, mwh2o, mwo2, mwh2, mwdry, cpair, cappa
  use ppgrid            ! pver, pverp is here
  use pmgrid            ! ?masterproc is here?
  use spmd_utils,       only: masterproc
  use rad_interp_mod    
  use radgrid
  use kabs
  use time_manager,     only: get_nstep
  use calc_opd_mod
  use exo_init_ref
 
  implicit none
  private
  save

!------------------------------------------------------------------------
!
! Public interfaces
!  
  public :: init_planck
  public :: aerad_driver

!------------------------------------------------------------------------
!
! private data
!
  
  ! Default values for namelist variables

  integer :: openstatus

  !
  ! Planck function interpolation table temperatures, etc. [K]:
  !
  integer, parameter  :: tpft_inc = 10                      ! increment [# pts/K]
  integer, parameter  :: tpft_beg = 1                       ! start temp (1 K)
  integer, parameter  :: tpft_end = 1000*tpft_inc           ! end temp (1000 K)
  integer, parameter  :: tpft_nt = tpft_end                 ! table temp dimension
  real(r8), parameter :: tpft_invfinc = 1.d0/dble(tpft_inc) ! factor dble(1/tpft_inc)
  real(r8), parameter :: tpft_dinc = dble(tpft_inc)         ! factor dble(tpft_inc)
  real(r8), parameter :: tpft_finc = real(tpft_inc)         ! factor real(tpft_inc)
  real(r8), dimension(tpft_nt,ntot_gpt) :: ptemp_itp        ! table [tpft_nt,ntot_gpt]
 

!============================================================================
contains
!============================================================================

!============================================================================
!
! Public subroutines 
!
!============================================================================

   subroutine init_planck
!------------------------------------------------------------------------
!
! Purpose: Initial reference value arrays
!
!------------------------------------------------------------------------

  implicit none

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: iq
    integer :: iw
    integer :: ip
    integer :: ig
    integer :: iv
    integer :: it

!    real(r8) :: wn1 
!    real(r8) :: wn2 
!    real(r8) :: tempsol 
!    real(r8) :: tempsol_iw
    real(r8) :: w1 
    real(r8) :: w2 
    real(r8) :: tempk
 !   logical :: found
 
!------------------------------------------------------------------------
!
! Start Code
! 

    !
    ! Calculate Planck function quantities as a function of wavelength point and
    !  temperature: 
    !

    if (masterproc) then
      write(*,*)   "INIT_PLANCK: CREATING PLANCK FUNCTION TABLE"
    endif
    ip = lw_ipbeg-1

    do iw=1,ntot_wavlnrng   ! Loop over all wavelength intervals

      ! Set these so that difference below is
      !   f(longer_wavlen)-f(shorter_wavlen),
      !   w1 = longer_wavlen, w2 = shorter_wavlen:
      !   h*c*nu/kc 
      w1 = dble(1.439*wavenum_edge(iw))      ! 1.439 ~ ((h*c)/k) ,
      w2 = dble(1.439*wavenum_edge(iw+1))    ! convert nu cm-1 to lambda (um

      do ig=1,ngauss_pts(iw)
        ip = ip+1

        do it=tpft_beg,tpft_end
          tempk = tpft_invfinc*dble(it)   ! Temperature [K]

          !ptemp_itp, spectrally integrated radiance (isotropic) in each longwave guass interval
          !value weighted according to gauss weights, (B outside mapping)
          ptemp_itp(it,ip) = (PLANCKf(w1,tempk)-PLANCKf(w2,tempk))* &
                              g_weight(ip)*SHR_CONST_STEBOL/SHR_CONST_PI
!write(*,*) "PLANCK", (PLANCKf(w1,tempk)-PLANCKf(w2,tempk)), g_weight(ip)
        enddo
        ptemp_itp(1:tpft_beg,ip) = ptemp_itp(tpft_beg,ip)
      enddo
    enddo

    return

  end subroutine init_planck


!============================================================================

  subroutine aerad_driver(ext_H2O, ext_CO2, ext_CH4, ext_H2, ext_N2, &
      ext_cicewp, ext_cliqwp, ext_cfrc, ext_rei, ext_rel, &
      ext_sfcT, ext_sfcP, ext_pmid, ext_pdel, ext_pdeldry, ext_tmid, &
      ext_pint, ext_pintdry, ext_cosZ, ext_msdist, ext_asdir,  & 
      ext_aldir, ext_asdif, ext_aldif,  &
      ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang,  &
      ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon,  &
      ext_TCx_obstruct, ext_TCz_obstruct, ext_zint, &
      sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, & 
      lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
      vis_dir, vis_dif, nir_dir, nir_dif, dzc )


!------------------------------------------------------------------------
!
! Purpose: Driver for correlated K radiative transfer code.
!          Recieves column data from radiation_tend
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    ! {intent IN}:
   
    integer, intent(in) :: ext_tslas_tog
    integer, intent(in) :: ext_tshadow_tog
    integer, intent(in) :: ext_nazm_tshadow

    real(r8), intent(in) :: ext_msdist
    real(r8), intent(in) :: ext_asdir          ! direct albedo (0.2-0.7 um) (from cam_in%asdir)
    real(r8), intent(in) :: ext_aldir          ! direct albedo (0.7-4.0 um) (from cam_in%aldir)
    real(r8), intent(in) :: ext_asdif          ! diffuse albedo (0.2-0.7 um) (from cam_in%asdif)
    real(r8), intent(in) :: ext_aldif          ! diffuse albedo (0.7-4.0 um) (from cam_in%aldif)
    real(r8), intent(in) :: ext_sfcT           ! surface temperature radiative  (from srfflx_state2d%ts)
    real(r8), intent(in) :: ext_sfcP           ! surface pressre (from state%ps) 
    real(r8), intent(in) :: ext_cosZ           ! cosine of the zenith angle 
    real(r8), intent(in) :: ext_rtgt           ! scaling used by models grid coordinate?
    real(r8), intent(in) :: ext_solar_azm_ang  ! solar azimuthal angle [rad]
    real(r8), intent(in) :: ext_tazm_ang       ! topographic slope and aspect angles from our model
    real(r8), intent(in) :: ext_tslope_ang     
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_cosz_horizon
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
    real(r8), intent(in), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct
    real(r8), intent(in), dimension(pverp) :: ext_zint     ! geopotential height at interfaces [m]
    real(r8), intent(in), dimension(pver) :: ext_pmid      ! pressure at midpoints [Pa]
    real(r8), intent(in), dimension(pver) :: ext_pdel      ! layer thickness [Pa]
    real(r8), intent(in), dimension(pver) :: ext_pdeldry   ! layer thickness dry atmosphere [Pa]
    real(r8), intent(in), dimension(pver) :: ext_tmid      ! temperature at midpoints [K]
    real(r8), intent(in), dimension(pverp) :: ext_pint     ! wet pressure at interfaces
    real(r8), intent(in), dimension(pverp) :: ext_pintdry     ! dry pressure at interfaces

    real(r8), intent(in), dimension(pver) :: ext_H2O       ! specific humidy (from state%q < q) at midlayer [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CO2       ! CO2 mass mixing ratio from state%q < co2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_CH4       ! CH4 mass mixing ratio from state%q < ch4mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_H2        ! H2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_N2        ! N2 mass mixing ratio from state%q < h2mmr)   [kg/kg]
    real(r8), intent(in), dimension(pver) :: ext_cicewp    ! in cloud ice water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver) :: ext_cliqwp    ! in cloud liquid water path at layer midpoints [g/m2]
    real(r8), intent(in), dimension(pver) :: ext_cFRC      ! cloud fraction]
    real(r8), intent(in), dimension(pver) :: ext_rei       ! ice cloud particle effective drop size ice [microns]
    real(r8), intent(in), dimension(pver) :: ext_rel       ! liquid cloud drop effective drop size liquid [micron   

    real(r8), intent(out), dimension(pver) ::  sw_dTdt     
    real(r8), intent(out), dimension(pver) ::  lw_dTdt     

    real(r8), intent(out), dimension(pverp) ::  sw_upflux   
    real(r8), intent(out), dimension(pverp) ::  sw_dnflux   
    real(r8), intent(out), dimension(pverp) ::  lw_upflux   
    real(r8), intent(out), dimension(pverp) ::  lw_dnflux

    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_dnflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_dnflux_spec

    real(r8), intent(out) ::  vis_dir
    real(r8), intent(out) ::  vis_dif
    real(r8), intent(out) ::  nir_dir
    real(r8), intent(out) ::  nir_dif

    real(r8), intent(out), dimension(pver) ::  dzc          ! [kg m-2], column amount of mass 

!------------------------------------------------------------------------
!
! Local Variables       
!
     real(r8), dimension(pverp) :: coldens      ! [molec m-2] wet columen amount per layer
     real(r8), dimension(pverp) :: coldens_dry  ! [molec m-2] dry columen amount per layer
     real(r8), dimension(pverp) :: qH2O         ! [kg/kg] H2O  mass mixing ratio mid layers 
     real(r8), dimension(pverp) :: qCO2         ! [kg/kg] CO2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qCH4         ! [kg/kg] CH4 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qO2          ! [kg/kg] O2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qO3          ! [kg/kg] O3 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qH2          ! [kg/kg] H2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: qN2          ! [kg/kg] H2 mass mixing ratio at mid layers   
     real(r8), dimension(pverp) :: cICE         ! [g/m2] in cloud ice water path at mid layers
     real(r8), dimension(pverp) :: cLIQ         ! [g/m2] in cloud liquid water path at mid layers
     real(r8), dimension(pverp) :: cfrc         ! cloud fraction at mid layers
     real(r8), dimension(pverp) :: REI          ! [microns] ice cloud particle effective radii at mid layers
     real(r8), dimension(pverp) :: REL          ! [microns] liquid cloud drop effective radii at mid layers
     real(r8), dimension(pverp) :: zlayer       ! [m] thickness of each vertical layer

     integer  :: swcut
     real(r8) :: tmp 
     real(r8) :: sinz
     real(r8) :: cosai 
     real(r8) :: atmp 
     real(r8) :: aint 
     real(r8) :: cosz_h 
     real(r8) :: x_obst
     real(r8) :: z_obst 
     real(r8) :: sw_cutoff
     real(r8) :: t1
     real(r8) :: t2
     real(r8), external :: tag_CPUtime

     integer :: iv
     integer :: ig
     integer :: ik
     integer :: k 
     integer :: i
     integer :: j
     integer :: iw
     integer :: ip
     integer :: ia0
     integer :: ia1
     integer :: ip_ibeg, ip_iend

     logical :: found
     logical :: horizon_extension 

     ! local variables used for computation (see Toon, 1989) annotate?
     real(r8), dimension(ntot_gpt,pverp) :: CK1sol, CK1ir
     real(r8), dimension(ntot_gpt,pverp) :: CK2sol, CK2ir       
     real(r8), dimension(ntot_gpt,pverp) :: CPBsol, CPBir
     real(r8), dimension(ntot_gpt,pverp) :: CMBsol, CMBir       
     real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y3   
     real(r8), dimension(ntot_gpt,pverp) :: EM1sol, EM1ir
     real(r8), dimension(ntot_gpt,pverp) :: EM2sol, EM2ir       
     real(r8), dimension(ntot_gpt,pverp) :: EL1sol, EL1ir
     real(r8), dimension(ntot_gpt,pverp) :: EL2sol, EL2ir
     real(r8), dimension(ntot_gpt,2*pverp) :: AFsol, AFir
     real(r8), dimension(ntot_gpt,2*pverp) :: BFsol, BFir
     real(r8), dimension(ntot_gpt,2*pverp) :: EFsol, EFir
     real(r8), dimension(ntot_gpt,pverp) :: AKsol, AKir
     real(r8), dimension(ntot_gpt,pverp) :: GAMIsol, GAMIir
     real(r8), dimension(ntot_gpt,pverp) :: EE1sol, EE1ir       
     real(r8), dimension(ntot_gpt,pverp) :: B1sol, B1ir  ! gamma 1  \
     real(r8), dimension(ntot_gpt,pverp) :: B2sol, B2ir  ! gamma 2   two stream parameters
     real(r8), dimension(ntot_gpt,pverp) :: B3sol, B3ir  ! gamma 3  /
     real(r8), dimension(ntot_gpt,pverp) :: DIRECTsol, DIRECTir       
     real(r8), dimension(ntot_gpt,pverp) :: DIRECTU   
     real(r8), dimension(ntot_gpt,pverp) :: DIREC     
     real(r8), dimension(ntot_gpt,pverp) :: TAUL       ! optical depth of each layer
     real(r8), dimension(ntot_gpt,pverp) :: OPD        ! cumulative optical depth (top down)
     real(r8), dimension(ntot_gpt,pverp) ::  tau_gas    ! gas optical depth array
     real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_ray    ! rayleigh optical depth
     real(r8), dimension(ntot_gpt,pverp) :: SOL       
     real(r8), dimension(ntot_gpt,pverp) :: W0         ! single scattering albedo        
     real(r8), dimension(ntot_gpt,pverp) :: G0         ! asymmetry parameter        

     ! Cloud optical properties
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  singscat_cld_mcica
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  asym_cld_mcica       
     real(r8), dimension(ncld_grp,ntot_gpt,pverp) ::  tau_cld_mcica
  
     ! stochastic bulk cloud properties (MCICA)
     real(r8), dimension(ntot_gpt,pverp) :: cFRC_mcica         
     real(r8), dimension(ntot_gpt,pverp) :: cICE_mcica
     real(r8), dimension(ntot_gpt,pverp) :: cLIQ_mcica

     real(r8), dimension(ntot_gpt,pverp) :: PTEMP     ! Planck function evaluated at each level  
     real(r8), dimension(ntot_gpt) :: PTEMPG          ! Planck function evaluated at ground
     real(r8), dimension(ntot_gpt,pverp) :: SLOPE     

     logical  :: part_in_tshadow
     real(r8), dimension(ntot_gpt) :: EMIS       ! Surface emissivity
     real(r8), dimension(ntot_gpt) :: RSFXdir    ! Surface reflectivity, direct radiation, gauss point grid
     real(r8), dimension(ntot_gpt) :: RSFXdif    ! Surface reflectivity, diffuse radiation, gauss point grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dir   ! Surface albedo, direct radiation, wavenumber grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_albedo_dif   ! Surface albedo, diffuse radiation, wavenumber grid
     real(r8), dimension(ntot_wavlnrng) :: sfc_emiss
     real(r8) :: sflux_frac
     real(r8) :: sfc_tempk  
     real(r8) :: sfc_press
     real(r8) :: cos_mu           
     logical  :: sw_on     ! switch for togopography and day/night
     logical  :: beamSolar ! switch for beam

     real(r8), dimension(pverp) ::  tint        ! [K] temperatures at level interfaces 
     real(r8), dimension(pverp) ::  tmid        ! [K] temperatures at level at mid layers + top (isothermal) 
     real(r8), dimension(pverp) ::  pmid        ! [Pa] pressure at level at mid layers + top (isothermal) 

     real(r8) :: dy     
!------------------------------------------------------------------------
!
! Start Code
!

    ! initialize internal RT variables
    CK1sol(:,:) = 0.0    ;    CK1ir(:,:) = 0.0
    CK2sol(:,:) = 0.0    ;    CK2ir(:,:) = 0.0
    CPBsol(:,:) = 0.0    ;    CPBir(:,:) = 0.0
    CMBsol(:,:) = 0.0    ;    CMBir(:,:) = 0.0
    Y3(:,:,:) = 0.0
    EM1sol(:,:) = 0.0    ;    EM1ir(:,:) = 0.0
    EM2sol(:,:) = 0.0    ;    EM2ir(:,:) = 0.0
    EL1sol(:,:) = 0.0    ;    EL1ir(:,:) = 0.0
    EL2sol(:,:) = 0.0    ;    EL2ir(:,:) = 0.0
    AFsol(:,:) = 0.0     ;    AFir(:,:) = 0.0
    BFsol(:,:) = 0.0     ;    BFir(:,:) = 0.0
    EFsol(:,:) = 0.0     ;    EFir(:,:) = 0.0
    AKsol(:,:) = 0.0     ;    AKir(:,:) = 0.0
    GAMIsol(:,:) = 0.0   ;    GAMIir(:,:) = 0.0
    EE1sol(:,:) = 0.0    ;    EE1ir(:,:) = 0.0
    B1sol(:,:) = 0.0     ;    B1ir(:,:) = 0.0
    B2sol(:,:) = 0.0     ;    B2ir(:,:) = 0.0
    B3sol(:,:) = 0.0     ;    B3ir(:,:) = 0.0
    DIRECTsol(:,:) = 0.0 ;    DIRECTir(:,:) = 0.0
    DIRECTU(:,:) = 0.0
    DIREC(:,:) = 0.0
    TAUL(:,:) = 0.0
    OPD(:,:) = 0.0
    tau_gas(:,:) = 0.0
    tau_cld_mcica(:,:,:) = 0.0
    singscat_cld_mcica(:,:,:) = 0.0
    asym_cld_mcica(:,:,:) = 0.0

    ! Fraction of the interplanetary solar flux at top of atmosphere:
    sflux_frac = dble(1./ext_msdist)    ! [1/AU^2]
   
    ! Set amount in layer above top defined model boundary
    qH2O(1) = ext_H2O(1)       ! H2O vapor mass concentration (specific humdity) [kg/kg]
    qCO2(1) = ext_CO2(1)       ! CO2 mass mixing ratio [kg/kg]
    qCH4(1) = ext_CH4(1)       ! CH4 mass mixing ratio [kg/kg]
    qH2(1) = ext_H2(1)         ! H2 mass mixing ratio [kg/kg]
    qN2(1) = ext_N2(1)         ! N2 mass mixing ratio [kg/kg]
    cICE(1) = ext_cicewp(1)    ! in cloud ice water path [g/m2]
    cLIQ(1) = ext_cliqwp(1)    ! in cloud liquid water path [g/m2]
    cFRC(1) = ext_cfrc(1)      ! cloud fraction
    REI(1) = ext_rei(1)        ! ice cloud particle effective radii [microns]
    REL(1) = ext_rel(1)        ! liquid cloud dropeffective radii [microns]
    tmid(1) = ext_tmid(1)      ! temperatures [K]
    pmid(1) = ext_pint(1)      ! pressure [Pa]
   
    ! Set amount in midlayers elsewhere 
    do k=2, pverp
      qH2O(k) = ext_H2O(k-1)
      qCO2(k) = ext_CO2(k-1) 
      qCH4(k) = ext_CH4(k-1) 
      qH2(k) = ext_H2(k-1)
      qN2(k) = ext_N2(k-1)
      cICE(k) = ext_cicewp(k-1) 
      cLIQ(k) = ext_cliqwp(k-1) 
      cFRC(k) = ext_cfrc(k-1) 
      REI(k) = ext_rei(k-1)
      REL(k) = ext_rel(k-1)
      tmid(k) = ext_tmid(k-1)
      pmid(k) = ext_pmid(k-1)
    enddo    

    ! Set ground (surface) values:
    sfc_tempk = ext_sfcT   ! [K]
    sfc_press = ext_sfcP   ! [pa]

    ! Set interface temperatures
    tint(pverp) = sfc_tempk
    tint(1) = ext_tmid(1)

    do k = 2, pver 
      dy = (log10(ext_pint(k)) - log10(ext_pmid(k))) / (log10(ext_pmid(k-1)) - log10(ext_pmid(k)))
      tint(k) = ext_tmid(k) - dy * (ext_tmid(k) - ext_tmid(k-1))
    enddo    

    ! Define molecular column density at each layer [molec m-2]  
  
    ! Set column density in layer above top boundary
    ! Note: ext_pintdry(1) is wet, ext_pint(>1) is dry
    coldens_dry(1) = (ext_pint(1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)*(1.0-qh2o(1))
    coldens(1) = (ext_pint(1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)
    ! Set column density for other mid layers
    do k=2, pverp   
      coldens(k) = (ext_pdel(k-1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)    !defined from wet air mass
!      coldens_dry(k) = (ext_pdeldry(k-1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)  !defined from dry air mass only
 coldens_dry(k) = (ext_pdel(k-1)*SHR_CONST_AVOGAD)/(mwdry*SHR_CONST_G)*(1.0-qh2o(k))  !kludge
    enddo    

    ! Define mass column density in each layer [kg m-2]
    dzc(:) = ext_pdel(:)/SHR_CONST_G  

    ! Define height of each layer [m]
    zlayer(1) = 0.0   !thickness of layer with lower boundary at model top is zero 
    do k=2, pverp
      zlayer(k-1) = (ext_zint(k-1) - ext_zint(k))    
    enddo

    !NOTES: surface albedo set from cam_in variable
    !Current implementation uses a gray albedo for shortwave
    !albedo and a gray albedo for longwave. The demarcation 
    !between to the two regimes is  ~5 um radiation.   

    ! Set surface direct albedo, diffuse albedo, and emissivity
    do iw=1,ntot_wavlnrng    ! Loop over relevant wavelength intervals 
      if (wavenum_edge(iw) .le. 2000) then  ! "thermal"
        sfc_albedo_dir(iw) = ext_aldir
        sfc_albedo_dif(iw) = ext_aldif
        sfc_emiss(iw) = 1.0
      endif
      if (wavenum_edge(iw) .gt. 2000 .and. wavenum_edge(iw) .le. 13000) then   ! "infrared"
        sfc_albedo_dir(iw) = ext_aldir       
        sfc_albedo_dif(iw) = ext_aldif
        sfc_emiss(iw) = 1.0 - sfc_albedo_dir(iw)
      endif
      if (wavenum_edge(iw) .ge. 13000) then     ! "solar" 
        sfc_albedo_dir(iw) = ext_asdir       
        sfc_albedo_dif(iw) = ext_asdif       
        sfc_emiss(iw) = 1.0 - sfc_albedo_dir(iw)
      endif
    enddo
    
    !write(*,*) "sfc_albedo", sfc_albedo, "sfc_emiss", sfc_emiss

    !*** Set cosine of incident solar angle (zenith angle); take into account
    !     (crudely) any horizon extension:
    cos_mu = dble(ext_cosZ)     ! [none]

    if(ext_tshadow_tog == 1) then   ! Need to check for horizon extension
    
      ! Calculate cosz_horizon array location/interpolative indices (using
      !  solar azimuth values):
    
      atmp = ext_solar_azm_ang/(2.0*SHR_CONST_PI/real(ext_nazm_tshadow)) 
      ia0 = int(atmp)
      aint = atmp-real(ia0)   ! Location crucial to avoid ia0 complications
    
      if(ia0 == 0) then
        ia0 = ext_nazm_tshadow
        ia1 = 1
      elseif(ia0 == ext_nazm_tshadow) then
        ia1 = 1
      else
        ia1 = ia0+1
      endif

      ! Calc horizon cosz, etc. for current solar azimuth via interpolation:
      z_obst = (1.-aint)*ext_TCz_obstruct(ia0)+ aint*ext_TCz_obstruct(ia1)      
      if(z_obst < 0.0) then   ! Horizon extension
        horizon_extension = .TRUE.
        cosz_h = (1.-aint)*ext_cosz_horizon(ia0)+aint*ext_cosz_horizon(ia1)
        if(real(cos_mu) > cosz_h .and. cos_mu < 4.d-2) then
          ! Sun is above extended horizon, but still low in sky;
          !  cosZ=0.04 =~ 87.71deg
          cos_mu = min(cos_mu-dble(cosz_h),4.d-2)
        endif
      else
        horizon_extension = .FALSE.
      endif
    endif

    ! Only do shortwave calculation if sun is sufficiently above horizon;
    !  set probability point interval indices for use in radiative transfer
    !  calcs:

    if(cos_mu > 1.d-6) then
      sw_on = .TRUE.       ! Sun sufficiently high, do  shortwave calculation

      ! Topographic slope/aspect correction:

      if(ext_tslas_tog == 1) then
        sinz = sin(acos(cos_mu))
        ! If cosai < 0., point not sunlit
        cosai = (cos_mu*cos(ext_tslope_ang))+  &
                 (sinz*sin(ext_tslope_ang)*  &
                 cos(ext_solar_azm_ang-ext_tazm_ang))
      endif

      ! Topographic shadowing:

      part_in_tshadow  = .FALSE.   ! Initialize   ! False, no part in shadow

      if(ext_tshadow_tog == 1 .and. .NOT.horizon_extension) then           ! possible shadowing
        cosz_h = (1.-aint)*ext_cosz_horizon(ia0)+aint*ext_cosz_horizon(ia1)
        x_obst = (1.-aint)*ext_TCx_obstruct(ia0)+aint*ext_TCx_obstruct(ia1)
        sw_cutoff = z_obst-(x_obst*tan(asin(real(cos_mu))))
      
        if(sw_cutoff > 0. .and. sw_cutoff < ext_zint(pverp)*ext_rtgt) then
          part_in_tshadow = .TRUE.
          swcut = 0
          atmp = 0.
          shadow_top: do ik=1,pver        ! Find top of shadow
            aint = ext_zint(ik)*ext_rtgt
            if(sw_cutoff < aint) then
              if(ik /= 1) then
                if(abs(aint-sw_cutoff) < abs(atmp-sw_cutoff)) then
                  swcut = pver-ik+1
                else
                  swcut = pver-ik+2
                endif
              else
                swcut = pver-ik+1
              endif
                exit shadow_top
            endif
            atmp = aint
          enddo shadow_top   ! Column sunlit above layer 'swcut'
       
          if(swcut == 0) then
            write(*,*) 'ERROR: unable to find shadow top.'
            stop '*** ERROR: aerad_driver01 ***'
          endif
      
        endif
      endif

    else     ! (cos_mu < 1.d-6) 
      sw_on = .FALSE.       ! Sun below horizon, do only longwave
    endif


    call calc_gasopd(tmid, pmid/100.0, ext_pdel/100.0, coldens, coldens_dry, qH2O, qCO2, qCH4, qO2, qO3, qH2, qN2, &
                     zlayer*100.0, tau_gas, tau_ray)

    !call calc_aeropd( )

    call calc_cldopd(ext_pint, cICE, cLIQ, REI, REL, cFRC, tau_cld_mcica, singscat_cld_mcica, & 
                     asym_cld_mcica, cFRC_mcica, cICE_mcica, cICE_mcica ) 

    call rad_precalc(pmid/100.0, tmid, tint, swcut, tau_gas, tau_ray, &
                     tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                     part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, sfc_emiss, & 
                     sflux_frac, sfc_tempk, cos_mu, sw_on, &                 
                     Y3, TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif)

    beamSolar = .true.  ! do solar calculation, all wavenlengths, two-stream quadrature
    ip_ibeg = sw_ipbeg
    ip_iend = sw_ipend
    call two_stream(TAUL, W0, G0, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                    EM1sol, EM2sol, EL1sol, EL2sol, &
                    AFsol, BFsol, EFsol, AKsol, &
                    GAMIsol, B1sol, B2sol, EE1sol)

    call add_txrad(EM1sol, EM2sol, EL1sol, EL2sol, &
                   AFsol, BFsol, EFsol, B1sol, B2sol, AKsol, &
                   TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                   ip_ibeg, ip_iend, beamSolar, CK1sol, CK2sol, CPBsol, CMBsol, B3sol, DIRECTsol)    

    beamSolar = .false.  ! do thermal calculation, all wavelengths, two-stream hemispheric mean
    ip_ibeg = lw_ipbeg
    ip_iend = lw_ipend
    call two_stream(TAUL, W0, G0, RSFXdir,RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                    EM1ir, EM2ir, EL1ir, EL2ir, &
                    AFir, BFir, EFir, AKir, &
                    GAMIir, B1ir, B2ir, EE1ir)

    call add_txrad(EM1ir, EM2ir, EL1ir, EL2ir, &
                   AFir, BFir, EFir, B1ir, B2ir, AKir, &
                   TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                   ip_ibeg, ip_iend, beamSolar, CK1ir, CK2ir, CPBir, CMBir, B3ir, DIRECTir)    

    call refine_lwflux(CK1ir, CK2ir, Y3, AKir, GAMIir, B3ir, EE1ir, &
                       TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir, RSFXdif, cos_mu, DIRECTU, DIREC)

    ! Calculate final fluxes / heating rates
    call rad_postcalc(CK1sol, CK2sol, CPBsol, CMBsol, &
                      EM1sol, EM2sol, EL1sol, EL2sol, &
                      DIRECTsol, DIRECTU, DIREC, dzc, swcut, part_in_tshadow, sw_on, &
                      sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                      lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                      vis_dir, vis_dif, nir_dir, nir_dif) 

    return

  end subroutine aerad_driver


!============================================================================

  subroutine rad_precalc(pmid, tmid, tint, swcut, tau_gas, tau_ray, &
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, &
                         part_in_tshadow, sfc_albedo_dir, sfc_albedo_dif, & 
                         sfc_emiss, sflux_frac, sfc_tempk, cos_mu, sw_on, &
                         Y3, TAUL, OPD, PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif)

!------------------------------------------------------------------------
!
! Purpose: Calculates quantities needed before the main radiative transfer 
!          calculation can be performed
!                                        
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in), dimension(pverp) :: pmid       ! [mb] pressures at mid layers)
    real(r8), intent(in), dimension(pverp) :: tmid       ! [K] temperatures at mid layers + top (isothermal)
    real(r8), intent(in), dimension(pverp) :: tint       ! [K] temperatures at level interfaces 
    integer, intent(in) :: swcut
    real(r8), intent(in), dimension(ntot_gpt,pverp) ::  tau_gas 
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) ::  tau_ray   
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  tau_cld_mcica
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  singscat_cld_mcica
    real(r8), intent(in), dimension(ncld_grp,ntot_gpt,pverp) ::  asym_cld_mcica
    logical, intent(in) :: part_in_tshadow 
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dir
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_albedo_dif
    real(r8), intent(in), dimension(ntot_wavlnrng) :: sfc_emiss
    real(r8), intent(in) :: sflux_frac
    real(r8), intent(in) :: sfc_tempk
    real(r8), intent(in) :: cos_mu           
    logical, intent(in) :: sw_on

    real(r8), intent(out), dimension(ntot_gpt,ngangles,pverp) :: Y3   
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: TAUL      
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: OPD       
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: PTEMP     
    real(r8), intent(out), dimension(ntot_gpt) :: PTEMPG
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: SLOPE      
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: SOL       
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: W0
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: G0        
    real(r8), intent(out), dimension(ntot_gpt) :: EMIS
    real(r8), intent(out), dimension(ntot_gpt) :: RSFXdir
    real(r8), intent(out), dimension(ntot_gpt) :: RSFXdif

!------------------------------------------------------------------------
!
! Local Variables
! 
    real(r8) :: den
    real(r8) :: x
    real(r8) :: f0_ig
    real(r8) :: w0_ig
    real(r8) :: g0_ig
    real(r8) :: taul_ig
    real(r8) :: tmp

    integer :: it
    integer :: ig
    integer :: k
    integer :: ip
    integer :: iw
    integer :: ib
    integer :: k_1
    integer :: ia

!   The vertical structure used by this radiation code is as follows:
!
!   Rad Index = Layer Boundaries
!   Atmos Index = Layer Centers
!
!   Layer                       Rad Index   Atmos Index
!
!    top -----------------------  j = 1
!         - - - - - - - - - - -   j = 2       k = 1
!        -----------------------
!         - - - - - - - - - - -               k = 2
!        -----------------------
!         - - - - - - - - - - -
!        -----------------------
!         - - - - - - - - - - -               k = NZ-1
!        -----------------------
!         - - - - - - - - - - -   j = 2*NZ-1  k = NZ
! bottom -----------------------  j = 2*NZ
!============================================================
 
!------------------------------------------------------------------------
!
! Start Code
!     

    ! Map/average important parameters to gauss points (it):
    it = 0
    do iw=1,ntot_wavlnrng    ! Loop over wavenumber bands
      do ig=1,ngauss_pts(iw)
        it = it+1
        EMIS(it) = sfc_emiss(iw)
        RSFXdir(it) = sfc_albedo_dir(iw)
        RSFXdif(it) = sfc_albedo_dif(iw)        
      enddo
    enddo    

    if(sw_on) then
      if(part_in_tshadow) then   ! Eliminate direct "shortwave" flux within shadow

        it = sw_ipbeg-1
        do iw=sw_iwbeg,sw_iwend  
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:swcut-1) = gw_solflux(it)*sflux_frac
            SOL(it,swcut:pverp) = 0.d0
          enddo
        enddo        

      else    ! No topographic shadowing in current column
       
        it = sw_ipbeg-1
        do iw=sw_iwbeg,sw_iwend
          do ig=1,ngauss_pts(iw)
            it = it+1
            SOL(it,1:pverp) = gw_solflux(it)*sflux_frac
          enddo
        enddo

      endif
    endif
    !if (sw_on) write(6,*) 'sw_on is true in precalc, solin1 is',solin1,'solinl1 is',solinl1
    
  
    OPD(1:ntot_gpt,:) = 0.d0     ! Initialize necessary array portions in
    TAUL(1:ntot_gpt,:) = 0.d0    !  anticipation of summing below

    do k=1,pverp     ! Loop over all layer BOUNDARIES  
   
      k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
      it = 0

      do iw=1,ntot_wavlnrng      ! Only necessary wavelength bands
        do ig=1,ngauss_pts(iw)
          it = it+1
 
          ! Combine the aerosol and gas optical parameters to get the total
          !  "effective" layer parameters:

          taul_ig = tau_gas(it,k) + tau_ray(iw,k)
          w0_ig = tau_ray(iw,k)
          g0_ig = 0.                             

          ! Add cloud optical depths
          do ip =1,2
            taul_ig = taul_ig + tau_cld_mcica(ip,it,k)      
            w0_ig = w0_ig + singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
            g0_ig = g0_ig + asym_cld_mcica(ip,it,k) * singscat_cld_mcica(ip,it,k) * tau_cld_mcica(ip,it,k)
          enddo
         
          ! Add aerosol optical depths here
          ! place holder

          if(taul_ig < SMALLd) then   ! Clip if optical depth too small
            taul_ig = SMALLd
          endif

          w0_ig = w0_ig/taul_ig        ! "total" weighted average singscat_albd
          w0_ig = min(1.d0-SMALLd,w0_ig)  !

          if(w0_ig > SMALLd) then
            g0_ig = g0_ig/(w0_ig*taul_ig)  ! "total" weighted average asym_fact
          else
            !?any point to this?    wsum = SMALLe   !? w0_ig = epsilon
            g0_ig = 0.d0           ! "total" weighted average asym_fact
          endif

          ! Apply delta-Eddington scaling to truncate the forward scattering
          !  lobe of the scattering phase function (effectively leave the
          !  energy within the truncated portion in the direct beam), as
          !  oringinally described in Joseph et al. (1976):
 
          f0_ig = g0_ig*g0_ig
          den = 1.d0-w0_ig*f0_ig
          TAUL(it,k) = taul_ig*den              ! "delta-scaled" total opt_depth
          W0(it,k) = (1.d0-f0_ig)*w0_ig/den     ! "delta-scaled" total singscat_albd
          G0(it,k) = g0_ig/(1.d0+g0_ig)         ! "delta-scaled" total asym_fact
                                                ! Actually '(g0_ig-f0_ig)/(1.d0-f0_ig)'

          ! Cumulative "Delta-scaled" optical depth, sums through layer k
          OPD(it,k) = OPD(it,k_1)+TAUL(it,k)    
                                                 
          !write (15, fmt = '(1X, I3, I3, ES15.7, ES15.7, ES15.7, ES15.7)') & 
          !  it,k,OPD(it,k),TAUL(it,k),w0(it,k),g0(it,k) 

        enddo
      enddo

      do ia=1,ngangles
        do it=1,ntot_gpt
          X = TAUL(it,k)/g_angle(ia)
          if(X <= 400.d0) then
            Y3(it,ia,k) = min(1.d0,exp(-X))
          else
            Y3(it,ia,k) = 0.d0    ! would be effectively zero anyway
          endif
        enddo 
      enddo
    enddo

    ! Calculate Planck function for current temperatures profiles
    call dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

    return

  end subroutine rad_precalc

!============================================================================

  subroutine dplanck(pmid, tmid, tint, TAUL, sfc_tempk, PTEMP, PTEMPG, SLOPE)

!------------------------------------------------------------------------
!
! Purpose: Calculate Planck function and its derivative at the ground and
!          at all altitudes.
!                                                  
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
!          
    real(r8), intent(in), dimension(pverp) :: pmid       ! [mb] pressures at mid layers
    real(r8), intent(in), dimension(pverp) :: tmid       ! [K] temperatures at mid layers (isothermal top)
    real(r8), intent(in), dimension(pverp) :: tint       ! [K] temperatures at level interfaces 
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL   ! Layer optical depth      
    real(r8), intent(in) :: sfc_tempk  
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: PTEMP   ! Planck function at each layer (isothermal top)
    real(r8), intent(out), dimension(ntot_gpt) :: PTEMPG        ! Planck function evaluated at ground
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: SLOPE     

!------------------------------------------------------------------------
!
! Local Variables
!          
    real(r8), dimension(ntot_gpt,pverp) :: ptemp_midlayer   
    real(r8) :: wl1
    real(r8) :: wl2
    real(r8) :: ft
    integer :: ip,iw,ig
    integer :: k
    integer :: k_1
    integer :: it

!------------------------------------------------------------------------
!
! Start Code
!          

    ! Calculate stuff based on the wavelength-dependent Planck function at the
    !  ground surface and for all vertical levels:

    do ip=lw_ipbeg,lw_ipend   ! Loop over all longwave gauss points

      it = int(sfc_tempk*tpft_finc)      
      ft = dble(sfc_tempk*tpft_finc-real(it))
      
      ! Interpolate between table values:
      PTEMPG(ip) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft

    enddo

    !---- Diagnostic output -----
    !---- comment out for simulations
    !if (masterproc) then
    !  write(*,*) "surface planck function [W m-2] in each spectral interval"
    !  ip=0
    !  do iw=1,ntot_wavlnrng
    !    do ig=1,ngauss_pts(iw)  
    !      ip=ip+1
    !    enddo
    !    write(*,*) iw,PTEMPG(ip)/g_weight(ip)
    !  write(*,*) iw,PTEMPG(50)/g_weight(50) 
    !  write(*,*) "TOTAL PLANCK FUNCTION:", SUM(PTEMPG)*SHR_CONST_PI
    !  write(*,*) "SURFACE PLANCK FUNCTION:",SHR_CONST_STEBOL*sfc_tempk**4
    !  write(*,*) "SURFACE TEMPERATURE:",sfc_tempk
    !endif
    !-----------

    ! Set Planck function at interfaces
    do k=1,pverp      ! Loop over all layers
      do ip=lw_ipbeg,lw_ipend   ! Loop over all longwave gauss points

        it = int(tint(k)*tpft_finc)
        ft = dble(tint(k)*tpft_finc-real(it))

        ! If temperature exceeds table range, for to be maximum
        ! Should only effect upper most atmosphere
        if (tint(k)*tpft_finc > tpft_end) then
           it = tpft_end-1
        endif
    
        if (tint(k)*tpft_finc < tpft_beg) then
           write(*,*) "negative temperature ", tint(k),k, "tint(k), k"
           it = tpft_beg-1
        endif

        ! Interpolate between table values:        
        PTEMP(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft   

        ! Determine midlayer planck function
        ! Slope derived from midlayer planck
        it = int(tmid(k)*tpft_finc)
        ft = dble(tmid(k)*tpft_finc-real(it))
        ptemp_midlayer(ip,k) = ptemp_itp(it,ip)*(1.d0-ft)+ptemp_itp(it+1,ip)*ft

        k_1 = max(1,k-1)   ! index for level above (except for k=1, of course)
        if(TAUL(ip,k) > SMALLd) then
          SLOPE(ip,k) = 2.0 * (ptemp_midlayer(ip,k)-ptemp_midlayer(ip,k_1)) / (TAUL(ip,k)+TAUL(ip,k_1))  ! DIML
        else
          SLOPE(ip,k) = 0.d0          
        endif

      enddo  

    enddo
      
    return

  end subroutine dplanck

!============================================================================

  function PLANCKf(e,t1)

!------------------------------------------------------------------------
!
! Purpose: Computes integral of the Planck function between zero and a
!          given wavelength          
!                                                  
!------------------------------------------------------------------------
!     ******************************************************
!     *  Purpose             :  Calculate Planck Function  *
!     *  Subroutines Called  :  None                       *
!     *  Input               :  WAVE, TEMP                 *
!     *  Output              :  PLANK                      *
!     * ****************************************************
!
!  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE PLANCK FUNCTION BETWEEN
!  ZERO AND THE SPECIFIED VALUE OF LAMBDA.  THUS (USING XL AS LAMBDA)
!  WE WANT TO INTEGRATE
!  R = INTEGRAL(XL=0 TO XL=XLSPEC) ( C1*XL**-5* / (EXP(C2/XL*T)-1) )*
!  SUBSTITUTING U=C2/(XL*T), THE INTEGRAL BECOMES
!  R = A CONSTANT TIMES INTEGRAL (USPEC TO INFINITY) OF
!            ( U**3 / (EXP(U) - 1) )*DU
!  THE APPROXIMATIONS SHOWN HERE ARE ON PAGE 998 OF ABRAMOWITZ AND
!  UNDER THE HEADING OF DEBYE FUNCTIONS.  C2 IS THE PRODUCT OF PLANCK'S
!  CONSTANT AND THE SPEED OF LIGHT DIVIDED BY BOLTZMANN'S CONSTANT.
!  C2 = 14390 WHEN LAMBDA IS IN MICRONS.
!  THE FACTOR 0.15399 IS THE RECIPROCAL OF SIX TIMES
!  THE SUM OF (1/N**2) FOR ALL N FROM ONE TO INFINITY.  IT IS CHOSEN TO
!  NORMALIZE THE INTEGRAL TO A MAXIMUM VALUE OF UNITY.
!  RADIATION IN REAL UNITS IS OBTAINED BY MULTIPLYING THE INTEGRAL BY
!  THE STEFAN-BOLTZMANN CONSTANT TIMES T**4.

    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!

    real(r8), intent(in)  :: e
    real(r8), intent(in)  :: t1
    

!------------------------------------------------------------------------
!
! Local Variables
!

    real(r8) :: PLANCKf
    real(r8), dimension(5) :: am
    real(r8) :: d
    real(r8) :: v1
    real(r8) :: a
    integer :: m
 
!------------------------------------------------------------------------
!
! Start Code
!
  
    d = 0.d0
    v1 = e/t1

    if(v1 <= 1.d0) then
      d = 1.0-0.15399*V1**3 *  &
      (1./3.-v1/8.+v1**2/60.-v1**4/5040.+v1**6/272160.-V1**8/13305600.)
    endif

    if(v1 > 1.d0 .and. v1 <= 50.d0) then
      do m=1,5
        a = dble(m)*v1
        am(m) = 0.15399*exp(-a)/m**4*(((a+3.)*a+6.)*a+6.)
      enddo

      d = am(1)+am(2)+am(3)+am(4)+am(5)
    endif

    PLANCKf = d*t1**4

    return

  end function PLANCKf


!============================================================================

  subroutine two_stream(TAUL, W0, G0, RSFXdir, RSFXdif, beamSolar, ip_ibeg, ip_iend, &
                        EM1, EM2, EL1, EL2, AF, BF, EF, AK, GAMI, B1, B2, EE1)

!------------------------------------------------------------------------
!
! Purpose: Defines matrix properties and sets up coefficients that do 
!          not depend on solar zenith angle or atmospheric temperature
!                                                  
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments

     real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
     real(r8), intent(in), dimension(ntot_gpt,pverp) :: W0         
     real(r8), intent(in), dimension(ntot_gpt,pverp) :: G0         
     real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
     real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
     logical, intent(in) :: beamSolar
     integer, intent(in) :: ip_ibeg
     integer, intent(in) :: ip_iend
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EM1
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EM2
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EL1
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EL2
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: AF        
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: BF        
     real(r8), intent(out), dimension(ntot_gpt,2*pverp) :: EF        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: AK        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: GAMI      
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: B1        
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: B2  
     real(r8), intent(out), dimension(ntot_gpt,pverp) :: EE1       

!------------------------------------------------------------------------
!
! Local Variables
!                                                  
   
    real(r8) :: x1
    integer :: ip
    integer :: k
    integer :: kd
    real(r8) u1i, u1i_2

!------------------------------------------------------------------------
!
! Start Code
!                                                  

    ! HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
    !  OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
    !  NEEDED FOR MATRIX:

    if (beamSolar) then 
      u1i_2 = U1I2sol
    else
      u1i_2 = U1I2ir
    endif

    do k=1,pverp
      do ip=ip_ibeg,ip_iend

        ! THESE ARE FOR QUADRATURE AND HEMISPHERIC MEAN
        B1(ip,k) = u1i_2*(2.d0-W0(ip,k)*(1.d0+G0(ip,k)))
        B2(ip,k) = u1i_2*W0(ip,k)*(1.d0-G0(ip,k))
        AK(ip,k) = SQRT(ABS(B1(ip,k)**2-B2(ip,k)**2))
        GAMI(ip,k) = B2(ip,k)/(B1(ip,k)+AK(ip,k))
        x1 = AK(ip,k)*TAUL(ip,k)         !
        if(x1 <= 400.d0) then            !
          EE1(ip,k) = exp(-x1)           !        
        else                             ! TIM add mod
          EE1(ip,k) = 0.d0               !
        endif                            !
        EL1(ip,k) = 1.d0+GAMI(ip,k)*EE1(ip,k)
        EM1(ip,k) = 1.d0-GAMI(ip,k)*EE1(ip,k)
        EL2(ip,k) = GAMI(ip,k)+EE1(ip,k)
        EM2(ip,k) = GAMI(ip,k)-EE1(ip,k)

      enddo
    
    enddo

    ! WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
    !  L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
    !  AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
    !  FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
    !  ANGLE OR TEMPERATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD:
    
    k = 0
    do kd=2,2*pverp-1,2
      k = k+1
      do ip=ip_ibeg,ip_iend
        !           HERE ARE THE EVEN MATRIX ELEMENTS
        AF(ip,kd) = EM1(ip,k+1)*EL1(ip,k)-EM2(ip,k+1)*EL2(ip,k)
        BF(ip,kd) = EM1(ip,k+1)*EM1(ip,k)-EM2(ip,k+1)*EM2(ip,k)
        EF(ip,kd) = EL1(ip,k+1)*EM2(ip,k+1)-EL2(ip,k+1)*EM1(ip,k+1)
        !           HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
        AF(ip,kd+1) = EM1(ip,k)*EL2(ip,k)-EL1(ip,k)*EM2(ip,k)
        BF(ip,kd+1) = EL1(ip,k+1)*EL1(ip,k)-EL2(ip,k+1)*EL2(ip,k)
        EF(ip,kd+1) = EL2(ip,k)*EM2(ip,k+1)-EL1(ip,k)*EM1(ip,k+1)
      enddo
    enddo

    ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
    !  NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY:

    do ip=ip_ibeg,ip_iend
      AF(ip,1) = 0.d0
      BF(ip,1) = EL1(ip,1)
      EF(ip,1) = -EM1(ip,1)
      AF(ip,2*pverp) = EL1(ip,pverp)-RSFXdif(ip)*EL2(ip,pverp)
      BF(ip,2*pverp) = EM1(ip,pverp)-RSFXdif(ip)*EM2(ip,pverp)
      EF(ip,2*pverp) = 0.d0
    enddo
  
    return    

  end subroutine two_stream


!============================================================================

  subroutine add_txrad (EM1, EM2, EL1, EL2, AF, BF, EF, B1, B2, AK, TAUL, OPD, & 
                        PTEMP, PTEMPG, SLOPE, SOL, W0, G0, EMIS, RSFXdir, RSFXdif, cos_mu, sw_on, &
                        ip_ibeg, ip_iend, beamSolar, CK1, CK2, CPB, CMB, B3, DIRECT )

!------------------------------------------------------------------------
!
! Purpose: Defines source terms, forms matrix for multiple layers and solves
!          the tri-diagonal equations to obtain mean intensity and net flux 
!                                                            
!------------------------------------------------------------------------

    implicit none

!------------------------------------------------------------------------
!
! Arguments
!          
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM1       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM2       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL1       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL2       
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: AF        
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: BF
  real(r8), intent(in), dimension(ntot_gpt,2*pverp) :: EF
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: B1        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: B2        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: AK        
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: OPD       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: PTEMP
  real(r8), intent(in), dimension(ntot_gpt) :: PTEMPG 
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: SLOPE     
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: SOL       
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: W0  
  real(r8), intent(in), dimension(ntot_gpt,pverp) :: G0 
  real(r8), intent(in), dimension(ntot_gpt) :: EMIS
  real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
  real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
  real(r8), intent(in) :: cos_mu           
  logical, intent(in) :: sw_on
  logical, intent(in) :: beamSolar
  integer, intent(in) :: ip_ibeg
  integer, intent(in) :: ip_iend
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CK1       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CK2       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CPB       
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: CMB
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: B3        
  real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIRECT       

!------------------------------------------------------------------------
!
! Local Variables
!          
    
    real(r8), dimension(ntot_gpt,pverp) :: CPP
    real(r8), dimension(ntot_gpt,pverp) :: CM
    real(r8), dimension(ntot_gpt,pverp) :: EE3
    real(r8), dimension(ntot_gpt,pverp) :: EL3
    real(r8), dimension(ntot_gpt,2*pverp) :: AS
    real(r8), dimension(ntot_gpt,2*pverp) :: DF
    real(r8), dimension(ntot_gpt,2*pverp) :: DS
    real(r8), dimension(ntot_gpt,2*pverp) :: XK

    real(r8), dimension(ntot_gpt) :: sfcs
    real(r8) :: DUo
    real(r8) :: B4
    real(r8) :: X2
    real(r8) :: X3
    real(r8) :: C1
    real(r8) :: C2
    real(r8) :: CP1
    real(r8) :: CM1
    real(r8) :: X
    real(r8) :: x_lim
    integer :: KINDEX
    integer :: k
    integer :: ip
    integer ::k_1
    integer :: kd

!------------------------------------------------------------------------
!
! Start Code
!        
    ! THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
    !  USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
    !  ATMOSPHERE.
    !

    if (beamSolar) then
    ! ******************************
    ! *   CALCULATIONS FOR SOLAR   *
    ! ******************************
    ! use delta eddington scaling (Joseph et al. 1976)

      if(.not. sw_on) return  ! if sun is below horizon, skip calculation

      DUo = 1.d0/cos_mu

      do k=1,pverp

        do ip=sw_ipbeg,sw_ipend  
    
          B3(ip,k) = 0.5d0*(1.d0-sqrt3*G0(ip,k)*cos_mu)     ! TIM mod
          B4 = 1.d0-B3(ip,k)
          X2 = TAUL(ip,k)*DUo
          if(X2 > 400.d0) then   !
            X2 = 400.d0          ! TIM add
          endif                  ! 
          X3 = OPD(ip,k)*DUo  
          if(X3 > 400.d0) then   !
            X3 = 400.d0          ! TIM add
          endif                  !         
          EL3(ip,k) = exp(-X3)*SOL(ip,k)        ! Beers law, layer bottom, tau=opd
          EE3(ip,k) = exp(-(X3-X2))*SOL(ip,k)   ! Beers law, layer top, tau=opd-taul
          DIRECT(ip,k) = cos_mu*EL3(ip,k)   
          C1 = B1(ip,k)-DUo
          if(abs(C1) < SMALLd) then
           C1 = sign(SMALLd,C1)
          endif
          C2 = AK(ip,k)*AK(ip,k)-DUo*DUo
          if(abs(C2) <= SMALLd) then
            C2 = SMALLd
          endif
          CP1 = W0(ip,k)*(B3(ip,k)*C1+B4*B2(ip,k))/C2   
          CM1 = (CP1*B2(ip,k)+W0(ip,k)*B4)/C1
          CPB(ip,k) = CP1*EL3(ip,k)     ! C+ lower boundary (tau = opd)
          CPP(ip,k) = CP1*EE3(ip,k)     ! C+ upper boundary (tau = opd-taul)   
          CMB(ip,k) = CM1*EL3(ip,k)     ! C- lower boundary (tau = opd)
          CM(ip,k) = CM1*EE3(ip,k)      ! C- upper boundary (tau = opd-taul) 

        enddo
      enddo

      ! CALCULATE SFCS, SHORTWAVE SOURCE AT THE BOTTOM (REFLECTION):
    
      do ip=sw_ipbeg,sw_ipend
        sfcs(ip) = DIRECT(ip,pverp)*RSFXdir(ip)
      enddo

    else  ! beamSolar

    ! ******************************
    ! * CALCULATIONS FOR INFRARED. *
    ! ******************************
      do k=1,pverp
        kindex = max(1,k-1)
        do ip=lw_ipbeg,lw_ipend
          B3(ip,k) = 1.d0/(B1(ip,k)+B2(ip,k))
          CPP(ip,k) = (PTEMP(ip,KINDEX)+SLOPE(ip,k)*B3(ip,k))*U1Sir  ! TIM mod
          CPB(ip,k) = CPP(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir         
          CM(ip,k) = (PTEMP(ip,KINDEX)-SLOPE(ip,k)*B3(ip,k))*U1Sir
          CMB(ip,k) = CM(ip,k)+SLOPE(ip,k)*TAUL(ip,k)*U1Sir
          EL3(ip,k)    = 0.d0
          DIRECT(ip,k) = 0.d0
          EE3(ip,k)    = 0.d0
        enddo
      enddo
 
      do ip=lw_ipbeg,lw_ipend
        sfcs(ip) = EMIS(ip)*PTEMPG(ip)*SHR_CONST_PI     
      enddo
    
    endif !beamSolar

    k = 0
    do kd=2,2*pverp-1,2
      k = k+1
      do ip=ip_ibeg,ip_iend
        ! HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
        DF(ip,kd+1) = EL2(ip,k)*(CPP(ip,k+1)-CPB(ip,k))+  &   ! TIM mod
           EL1(ip,k)*(CMB(ip,k)-CM(ip,k+1))
        ! HERE ARE THE EVEN MATRIX ELEMENTS
        DF(ip,kd) = (CPP(ip,k+1)-CPB(ip,k))*EM1(ip,k+1)-  &   ! TIM mod
           (CM(ip,k+1)-CMB(ip,k))*EM2(ip,k+1)
      enddo
    enddo

    ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
    !  BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
    !  DIFFUSE RADIATION IS INCIDENT AT THE TOP:
    do ip=ip_ibeg,ip_iend
      DF(ip,1) = -CM(ip,1)
      DF(ip,2*pverp) = SFCS(ip)+RSFXdif(ip)*CMB(ip,pverp)-CPB(ip,pverp)
      DS(ip,2*pverp) = DF(ip,2*pverp)/BF(ip,2*pverp)
      AS(ip,2*pverp) = AF(ip,2*pverp)/BF(ip,2*pverp)
    enddo

    ! ********************************************
    ! *     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
    ! ********************************************
    x_lim = tiny(x)*1.d100   ! Lower limit for 'x' below        ! TIM add
    do kd=2,2*pverp
      k = 2*pverp+1-kd     ! TIM add
      do ip=ip_ibeg,ip_iend
        x = (BF(ip,k)-EF(ip,k)*AS(ip,k+1))   ! TIM mod
        if(abs(x) < x_lim) then     !Ã
          x = sign(1.d0,x)*x_lim    ! TIM add
        endif                       !
        X = 1.d0/x                  !
        AS(ip,k) = AF(ip,k)*X
        DS(ip,k) = (DF(ip,k)-EF(ip,k)*DS(ip,k+1))*X
      enddo
    enddo

    do ip=ip_ibeg,ip_iend
      XK(ip,1) = DS(ip,1)
    enddo
    
    do kd=2,2*pverp
      do ip=ip_ibeg,ip_iend
        XK(ip,kd) = DS(ip,kd)-AS(ip,kd)*XK(ip,kd-1)
      enddo
    enddo

    ! ***************************************************************
    !    CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
    ! ***************************************************************
    
    do k=1,pverp
      do ip=ip_ibeg,ip_iend
        CK1(ip,k) = XK(ip,2*k-1)
        CK2(ip,k) = XK(ip,2*k)
        !### FNET and TMI not used in this application (similar quantities are instead
        !     computed in `postcalc`)
        !FNET(ip,k)  = CK1(ip,k)  *( EL1(ip,k) -EL2(ip,k))   +  &
        !    CK2(ip,k) *( EM1(ip,k)-EM2(ip,k) ) + CPB(ip,k) -  &
        !    CMB(ip,k) - DIRECT(ip,k)
        !TMI(ip,k)   =  EL3(ip,k) + U1I(ip) *( CK1(ip,k)  *  &
        !   ( EL1(ip,k) + EL2(ip,k))   +  &
        !   CK2(ip,k) *( EM1(ip,k)+EM2(ip,k) ) +  &
        !   CPB(ip,k) + CMB(ip,k) )
      enddo
    enddo
    
    return
  end subroutine add_txrad
                  

!============================================================================

  subroutine refine_lwflux (CK1, CK2, Y3, AK, GAMI, B3, EE1, TAUL, PTEMP, PTEMPG, SLOPE, EMIS, RSFXdir,RSFXdif, cos_mu, DIRECTU, DIREC)

!------------------------------------------------------------------------
!
! Purpose: Calculate upward and downward radiances and intensities using 
!          Gauss Quadrature angles and weights          
!          original                                                  
!------------------------------------------------------------------------

     implicit none

!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK1       
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK2       
    real(r8), intent(in), dimension(ntot_gpt,ngangles,pverp) :: Y3   
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: AK        
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: GAMI      
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: B3        
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EE1           
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: TAUL      
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: PTEMP     
    real(r8), intent(in), dimension(ntot_gpt) :: PTEMPG     
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: SLOPE     
    real(r8), intent(in), dimension(ntot_gpt) :: EMIS
    real(r8), intent(in), dimension(ntot_gpt) :: RSFXdir
    real(r8), intent(in), dimension(ntot_gpt) :: RSFXdif
    real(r8), intent(in) :: cos_mu

    real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIRECTU   
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: DIREC    

!
!     *  Input               :  PTEMP, SLOPE, Y3, B3, EE1, EE2     *
!     *  Output              :  Irradiances: DIREC, DIRECTU
!
!------------------------------------------------------------------------
!
! Local Variables
!          
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: uintent  ! [ntot_gpt,ngangles,pverp]
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: dintent  ! [ntot_gpt,ngangles,pverp]

    integer :: openstatus
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y1   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y2   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y4   
    real(r8), dimension(ntot_gpt,ngangles,pverp) :: Y8   
    real(r8), dimension(ntot_gpt,pverp) :: Y5        
    real(r8), dimension(ntot_gpt,pverp) :: A1        
    real(r8), dimension(ntot_gpt,pverp) :: A2        
    real(r8), dimension(ntot_gpt,pverp) :: A3        
    real(r8), dimension(ntot_gpt,pverp) :: A4        
    real(r8), dimension(ntot_gpt,pverp) :: A7        
    real(r8) :: X4
    real(r8) :: YA
    real(r8) :: YB
    real(r8) :: CKP
    real(r8) :: CKM
    real(r8) :: cos_mu1
    integer :: KINDEX
    integer :: k
    integer :: kd
    integer :: ip
    integer :: ia
!------------------------------------------------------------------------
!
! Start Code

    do k=1,pverp
      kindex = max(1,k-1)
      do ip=lw_ipbeg,lw_ipend
        ! HERE WE DO NO SCATTERING COEFFICIENTS
        A3(ip,k) = PTEMP(ip,KINDEX)*2.0*SHR_CONST_PI
        A4(ip,k) = 2.0*SHR_CONST_PI*SLOPE(ip,k)
        A7(ip,k) = A3(ip,k)
        Y5(ip,k) = A4(ip,k)*TAUL(ip,k)
      enddo
      ! HERE WE DO SCATTERING
      do ip=lw_ipbeg,lw_ipend
        X4 = SLOPE(ip,k)*(2.0*SHR_CONST_PI*B3(ip,k)-U1Sir)
        A1(ip,k) = U1Iir-AK(ip,k)
        A2(ip,k) = GAMI(ip,k)*(AK(ip,k)+U1Iir)
        A3(ip,k) = A3(ip,k)+X4
        A7(ip,k) = A7(ip,k)-X4
      enddo
    enddo

    ! CALCULATIONS FOR ALL GAUSS POINTS:
    do k=1,pverp   
      do ia=1,ngangles
       
        ! HERE WE DO NO SCATTERING COEFFS

        do ip=lw_ipbeg,lw_ipend
          !%waste CPU%         Y1(ip,ia,k)  =  0.d0
          !%  Y2(ip,ia,k)  =  0.d0
          Y4(ip,ia,k) = A7(ip,k)-A4(ip,k)*g_angle(ia)
          Y8(ip,ia,k) = A3(ip,k)+A4(ip,k)*g_angle(ia)
          !%      enddo

          ! HERE WE DO SCATTERING

          !% do ip=lw_ipbeg,lw_ipend
          YA = A1(ip,k)*(Y3(ip,ia,k)-EE1(ip,k))/(AK(ip,k)*g_angle(ia)-1.d0)
          YB = A2(ip,k)*(1.d0-EE1(ip,k)*Y3(ip,ia,k))/  &
               (AK(ip,k)*g_angle(ia)+1.d0)
          CKP = CK1(ip,k)+CK2(ip,k)
          CKM = CK1(ip,k)-CK2(ip,k)
          Y1(ip,ia,k) = CKP*YB+CKM*YA
          Y2(ip,ia,k) = CKP*YA+CKM*YB

        enddo

      enddo
    enddo

    ! INITIALIZE IRRADIANCES TO ZERO
    do  k=1, pverp
      do ip=lw_ipbeg,lw_ipend
        DIREC(ip,k) = 0.d0
        DIRECTU(ip,k) = 0.d0
      enddo   
    enddo

    ! DIREC IS DOWNWARD IRRADIANCE. DIRECTU IS UPWARD IRRADIANCE.
    !  CALCULATE DINTENT THE DOWNWARD RADIANCE AND DIREC THE DOWNWARD IRRADIANCE
   
    ! BOUNDARY CONDITIONS: DOWNWARD IRRADIANCE, RADIANCE AT TOA (k = camtop) was 1

    cos_mu1 = max(cos_mu,0.0)

    do ia=1,ngangles
      do ip=lw_ipbeg,lw_ipend
        
        DINTENT(ip,ia,camtop) = (1.d0-Y3(ip,ia,camtop))*Y4(ip,ia,camtop)+Y1(ip,ia,camtop)

        DIREC(ip,camtop) = DIREC(ip,camtop)+DINTENT(ip,ia,camtop)*g_ang_weight(ia)
        
      enddo
    enddo


    ! DINTENT IS DOWNWARD RADIANCE * TwoPI. DIREC IS THE DOWNWARD IRRADIANCE.   
    !  CALCULATE FOR REST OF ATMOSPHERE.
    do k=camtop+1, pverp     !was k=2
      do ia=1,ngangles
        do ip=lw_ipbeg,lw_ipend
          DINTENT(ip,ia,k) = DINTENT(ip,ia,k-1)*Y3(ip,ia,k)  & 
               +Y1(ip,ia,k)+Y5(ip,k)+(1.d0-Y3(ip,ia,k))*Y4(ip,ia,k)             

          DIREC(ip,k) = DIREC(ip,k)+DINTENT(ip,ia,k)*g_ang_weight(ia)
        enddo
      enddo             
    enddo

    ! UINTENT IS THE UPWARD RADIANCE * TwoPI. DIRECTU IS THE UPWARD IRRADIANCE.
    !  ASSUME THAT THE REFLECTIVITY IS LAMBERT.
    
    ! BOUNDARY CONDITIONS: UPWARD IRRADIANCE, RADIANCE AT BOTTOM (k = pverp)

    do ia=1,ngangles
      do ip=lw_ipbeg,lw_ipend

        UINTENT(ip,ia,pverp) = EMIS(ip)*PTEMPG(ip)*2.0*SHR_CONST_PI !+RSFXdif(ip)*DIREC(ip,pverp)*2.0d0

        DIRECTU(ip,pverp) = DIRECTU(ip,pverp)+UINTENT(ip,ia,pverp)*g_ang_weight(ia)

      enddo
    enddo

    ! CALCULATE FOR THE REST OF THE ATMOSPHERE
    ! ITERATED FROM BOTTOM UP
    kd = pverp
    do k=camtop+1,pverp   ! was k=2
      kd = kd-1
      do ia=1,ngangles
        do ip=lw_ipbeg,lw_ipend

          UINTENT(ip,ia,kd) = (UINTENT(ip,ia,kd+1)-Y5(ip,kd+1))  &
                              *Y3(ip,ia,kd+1)+Y2(ip,ia,kd+1)+  &
                              (1.d0-Y3(ip,ia,kd+1))*Y8(ip,ia,kd+1)          

          DIRECTU(ip,kd) = DIRECTU(ip,kd)+UINTENT(ip,ia,kd)*g_ang_weight(ia)

        enddo
      enddo
    enddo

    return

  end subroutine refine_lwflux

!============================================================================


  subroutine rad_postcalc (CK1, CK2, &
                           CPB, CMB, &
                           EM1, EM2, EL1, EL2, &
                           DIRECT, DIRECTU, DIREC, dzc, swcut, part_in_tshadow, sw_on, &
                           sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux, &
                           lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
                           vis_dir, vis_dif, nir_dir, nir_dif )

!------------------------------------------------------------------------
!
! Purpose: Calculate total radiative fluxes; put in a form suitable for output 
!                                                                      
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
!              

    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CK2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CPB
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: CMB
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EM2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL1
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: EL2
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIRECT       
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIRECTU   
    real(r8), intent(in), dimension(ntot_gpt,pverp) :: DIREC
    real(r8), intent(in), dimension(pver) ::  dzc          ! [kg m-2], column amount of mass 
    integer, intent(in) :: swcut
    logical, intent(in) :: part_in_tshadow
    logical, intent(in) :: sw_on
 
    real(r8), intent(out), dimension(pver) ::  sw_dTdt     
    real(r8), intent(out), dimension(pver) ::  lw_dTdt     
    real(r8), intent(out), dimension(pverp) ::  sw_upflux   
    real(r8), intent(out), dimension(pverp) ::  sw_dnflux   
    real(r8), intent(out), dimension(pverp) ::  lw_upflux   
    real(r8), intent(out), dimension(pverp) ::  lw_dnflux       
    ! spectral outputs
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  sw_dnflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_upflux_spec
    real(r8), intent(out), dimension(pverp,ntot_wavlnrng) ::  lw_dnflux_spec
    real(r8), intent(out) ::  vis_dir       
    real(r8), intent(out) ::  vis_dif       
    real(r8), intent(out) ::  nir_dir       
    real(r8), intent(out) ::  nir_dif       


!------------------------------------------------------------------------
!
! Local Variables
!          
    
    integer :: k
    integer :: ip,ipc,iw, ig
    integer :: k_1
    integer :: j
    real(r8) :: lyr_mass_fact
    real(r8) :: fdiv_sw
    real(r8) :: fdiv_lw


!------------------------------------------------------------------------
!
! Start Code
!          

    sw_dTdt(:) = 0.    ! Initialize heating rate arrays
    lw_dTdt(:) = 0.    !
 
    lw_dnflux(:) = 0.   !
    lw_upflux(:) = 0.   ! Initialize entire arrays for summing below
    sw_upflux(:) = 0.   !
    sw_dnflux(:) = 0.   !

    lw_dnflux_spec(:,:) = 0.   !
    lw_upflux_spec(:,:) = 0.   ! Initialize entire arrays for summing below
    sw_upflux_spec(:,:) = 0.   !
    sw_dnflux_spec(:,:) = 0.   !

    vis_dir = 0.     !
    vis_dif = 0.     !  Initialize solar fluxes to surface
    nir_dir = 0.     !  Pass to land model in CESM
    nir_dif = 0.     !
    
    ! Finalize fluxes: 
    do k=camtop,pverp    ! Loop over all layer BOUNDARIES, was k=1
      ip=lw_ipbeg
      do iw=lw_iwbeg,lw_iwend    ! Loop over wavenumber bands
        do ig=1,ngauss_pts(iw)

          !broadband fluxes
          lw_dnflux(k) = lw_dnflux(k)+DIREC(ip,k)
          lw_upflux(k) = lw_upflux(k)+DIRECTU(ip,k)

          !spectral fluxes
          lw_dnflux_spec(k,iw) = lw_dnflux_spec(k,iw)+DIREC(ip,k)
          lw_upflux_spec(k,iw) = lw_upflux_spec(k,iw)+DIRECTU(ip,k)

          ip=ip+1
        enddo
      enddo
    enddo

    if (sw_on) then
      do k=camtop,pverp  ! Loop over all layer BOUNDARIES, was k=1
        ip=sw_ipbeg
        do iw=sw_iwbeg,sw_iwend  ! Loop over all "shortwave" wavelength intervals
          do ig=1,ngauss_pts(iw)

            ! broadband fluxes
            sw_upflux(k) = sw_upflux(k)+CK1(ip,k)*EL1(ip,k)+  &
                           CK2(ip,k)*EM1(ip,k)+CPB(ip,k)
            sw_dnflux(k) = sw_dnflux(k)+CK1(ip,k)*EL2(ip,k)+  &
                           CK2(ip,k)*EM2(ip,k)+CMB(ip,k)+DIRECT(ip,k)

            ! spectral fluxes
            sw_upflux_spec(k,iw) = sw_upflux_spec(k,iw)+CK1(ip,k)*EL1(ip,k)+  &
                                   CK2(ip,k)*EM1(ip,k)+CPB(ip,k)
            sw_dnflux_spec(k,iw) = sw_dnflux_spec(k,iw)+CK1(ip,k)*EL2(ip,k)+  &
                                   CK2(ip,k)*EM2(ip,k)+CMB(ip,k)+DIRECT(ip,k)

            ip=ip+1
          enddo
        enddo
      enddo
    endif

    if(sw_on) then      !***** BOTH "longwave" and "shortwave" *****

      ! Calculate atmospheric heating rates (dT/dt) [K/s]:
      ! Column sunlit above layer 'swcut'

      if(part_in_tshadow) then  ! NOTES: set to false always
       
        do k=camtop+1,swcut-1  ! Above shadow, was k=2
      
          lyr_mass_fact = dzc(k-1)*cpair
        
          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact      ! "shortwave" heating rate [K/s]

          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
          
        enddo

        do camtop=swcut,pverp  ! Within shadow, no shortwave calculation

          lyr_mass_fact = dzc(k-1)*cpair
      
          sw_dTdt(k-1) = 0.0   ! "shortwave" heating rate [K/s]
          
          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact    ! "longwave" heating rate [K/s]

        enddo 

      else   ! No shadow in column
       
        do k=camtop+1,pverp  ! was k=2

          lyr_mass_fact = dzc(k-1)*cpair

          fdiv_sw = (sw_upflux(k)-sw_dnflux(k))-(sw_upflux(k-1)-sw_dnflux(k-1))
          sw_dTdt(k-1) = fdiv_sw/lyr_mass_fact      ! "shortwave" heating rate [K/s]

          fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
          lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
  
        enddo

      endif
 
    else   !***** ONLY "longwave" ***** 
     
      do k=camtop+1, pverp   ! was k=2

        lyr_mass_fact = dzc(k-1)*cpair
        fdiv_lw = (lw_upflux(k)-lw_dnflux(k))-(lw_upflux(k-1)-lw_dnflux(k-1))
        lw_dTdt(k-1) = fdiv_lw/lyr_mass_fact       ! "longwave" heating rate [K/s]
      
        !write(*,*) "--------------------------------------------------------"
        !write(*,*) k,"lw_up(k)  lw_up(k-1)  lw_dn(k)  lw_dn(k-1)"
        !write(*,*) lw_upflux(k), lw_upflux(k-1), lw_dnflux(k), lw_dnflux(k-1)
        !write(*,*) fdiv_lw,lyr_mass_fact
        !write(*,*) lw_dTdt(k-1)*24.0*60.*60.0

      enddo

    endif

    !if (camtop > 1) then
    !  lw_dTdt(:camtop-1) = 0.0
    !  sw_dTdt(:camtop-1) = 0.0
    !endif

    ! Calculate surface shortwave fluxes to land model
    ipc=sw_ipbeg 
    do iw=sw_iwbeg,sw_iwend    ! Loop over relevant wavelength intervals
      if (wavenum_edge(iw) .gt. 13000) then 
        do ip=1,ngauss_pts(iw)
          vis_dir = vis_dir + DIRECT(ipc,pverp)
          vis_dif = vis_dif + CK1(ipc,pverp)*EL2(ipc,pverp)+CK2(ipc,pverp)*EM2(ipc,pverp)+CMB(ipc,pverp)
          ipc=ipc+1
        enddo
      else
        do ip=1,ngauss_pts(iw)
          nir_dir = nir_dir + DIRECT(ipc,pverp)
          nir_dif = nir_dif + CK1(ipc,pverp)*EL2(ipc,pverp)+CK2(ipc,pverp)*EM2(ipc,pverp)+CMB(ipc,pverp)
          ipc=ipc+1
        enddo
      endif
    enddo

    return

  end subroutine rad_postcalc


!============================================================================


end module exo_radiation_mod
