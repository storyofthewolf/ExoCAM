
module exo_radiation_cam_intr

!---------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the radiation code
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
			      SHR_CONST_LOSCHMIDT, &
			      SHR_CONST_MEARTH, SHR_CONST_MSUN
  use physconst,        only: scon,mwn2, mwco2, mwch4, mwh2o, mwo2, mwh2, mwdry, cpair, cappa
  use ppgrid            ! pver, pverp is here
  use pmgrid            ! ?masterproc is here?
  use spmd_utils,       only: masterproc
  use rad_interp_mod    
  use radgrid
  use kabs
  use exoplanet_mod,    only: do_exo_rt_clearsky, exo_rad_step, do_exo_rt_spectral, &
                              do_exo_circumbinary, &
                              exo_n2mmr, exo_h2mmr, exo_co2mmr, exo_ch4mmr, &
                              cb_mass, cb_teff, cb_radius, cb_semia, cb_eccen, cb_spin, &
                              cb_man, cb_transit
  use time_manager,     only: get_nstep
  use initialize_rad_mod_cam
  use exo_radiation_mod
  use abortutils,      only: endrun
 
  implicit none
  private
  save

!------------------------------------------------------------------------
!
! Public interfaces
!  

  public :: exo_radiation_init        
  public :: exo_radiation_tend  
  public :: exo_radiation_nextsw_cday

!------------------------------------------------------------------------
!
! private data
!
  
  ! Default values for namelist variables

  integer :: openstatus

  !
  ! Physics buffer indices
  !
  integer :: qrs_idx      = 0
  integer :: qrl_idx      = 0
  integer :: cld_idx      = 0
  integer :: concld_idx   = 0
  integer :: rel_idx      = 0
  integer :: rei_idx      = 0
  integer :: cicewp_idx = -1
  integer :: cliqwp_idx = -1
  integer :: cldemis_idx = -1
  integer :: cldtau_idx = -1
  integer :: nmxrgn_idx = -1
  integer :: pmxrgn_idx = -1

  !NOTES: THESE AREN'T HOOKED UP, NO AEROSOLS IN RT
  ! Aerosol optical properties
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  singscat_aer 
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  asym_aer     
  !real(r8), dimension(naerspc_rad,ntot_wavlnrng,pverp) ::  tau_aer      

!============================================================================
contains
!============================================================================

!============================================================================
!
! Public subroutines 
!
!============================================================================

!============================================================================

  subroutine exo_radiation_init()
!-----------------------------------------------------------------------
!
! Purpose: Initialize the radiation parameterization, add fields to the 
!          history buffer
!
!-----------------------------------------------------------------------

    use cam_history,     only: addfld, add_default, phys_decomp
    use physics_buffer,  only: pbuf_get_index
    use radiation_data,  only: init_rad_data
    use exo_init_ref,    only: init_ref
    use spectral_output
    use circumbinary_mod, only: circumbinary_addfld_stdoutput, circumbinary_init_orbit

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: iv
    integer :: ig
    
!------------------------------------------------------------------------
!
! Start Code
!
    ! initialize solar spectrum, clouds, and kcoeff
    call initialize_solar
    call initialize_kcoeff
    call initialize_cldopts
    call init_ref
    call init_planck

    ! set top layer of cam for computation
    camtop = 1

    ! Get physics buffer indices
    cld_idx    = pbuf_get_index('CLD')
    !concld_idx = pbuf_get_index('CONCLD')
    rel_idx    = pbuf_get_index('REL')
    rei_idx    = pbuf_get_index('REI')
    cicewp_idx = pbuf_get_index('CICEWP')
    cliqwp_idx = pbuf_get_index('CLIQWP')
    !cldemis_idx= pbuf_get_index('CLDEMIS')
    !cldtau_idx = pbuf_get_index('CLDTAU')
    nmxrgn_idx = pbuf_get_index('NMXRGN')
    pmxrgn_idx = pbuf_get_index('PMXRGN')
    qrs_idx = pbuf_get_index('QRS')
    qrl_idx = pbuf_get_index('QRL')

    ! Add Shortwave radiation fields
    call addfld ('FDSTOA   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLSD  ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLC   ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface clearsky',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLSC   ','W/m2    ',1,    'A','Solar downward visible direct  to surface clearsky',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLDC  ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface clearsky',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLSDC  ','W/m2    ',1,    'A','Solar downward visible diffuse to surface clearsky',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRS     ','K/day   ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRSC    ','K/day   ',pver, 'A','Clearsky solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
                                                                                              sampling_seq='rad_lwsw')
    call addfld ('FSDTOA  ','W/m2    ',pverp,'A','Shortwave downward flux at top of atmosphere',phys_decomp)
    call addfld ('FSN200  ','W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSN200C ','W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
                                                                                                   sampling_seq='rad_lwsw')
    call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUS     ','W/m2    ',pverp,'A','Shortwave upward flux',phys_decomp)
    call addfld ('FDS     ','W/m2    ',pverp,'A','Shortwave downward flux',phys_decomp)
    call addfld ('FUSC    ','W/m2    ',pverp,'A','Shortwave clear-sky upward flux',phys_decomp)
    call addfld ('FDSC    ','W/m2    ',pverp,'A','Shortwave clear-sky downward flux',phys_decomp)
    call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                                                                               phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                                                                                 phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
                                                                                              sampling_seq='rad_lwsw')
    call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNR    ','W/m2    ',1,    'A','Net solar flux at tropopause',phys_decomp, sampling_seq='rad_lwsw')
    call add_default ('QRS     ', 1, ' ')
    call add_default ('FSNS    ', 1, ' ')
    call add_default ('FSNT    ', 1, ' ')
    call add_default ('FSNTOA  ', 1, ' ')
    call add_default ('FSNTOAC ', 1, ' ')
    call add_default ('FSNTC   ', 1, ' ')
    call add_default ('FSNSC   ', 1, ' ')
    call add_default ('FSDSC   ', 1, ' ')
    call add_default ('FSDS    ', 1, ' ')
    call add_default ('SWCF    ', 1, ' ')

    ! aerosol forcing-only calculations
    call addfld ('FSNT_RF ','W/m^2   ',1, 'A','Total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNTC_RF','W/m^2   ',1, 'A','Clear sky total column absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNS_RF ','W/m^2   ',1, 'A','Surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('FSNSC_RF','W/m^2   ',1, 'A','Clear sky surface absorbed solar flux (radforce)' ,phys_decomp)
    call addfld ('QRS_RF  ','K/s     ',pver, 'I','Solar heating rate (radforce)' ,phys_decomp)

    ! Longwave radiation
    call addfld ('QRL     ','K/day     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRLC    ','K/day     ',pver, 'A','Longwave heating rate clearsky',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
                                                                                                   sampling_seq='rad_lwsw')
    call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200  ','W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200C ','W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUL     ','W/m2    ',pverp,'A','Longwave upward flux',phys_decomp)
    call addfld ('FDL     ','W/m2    ',pverp,'A','Longwave downward flux',phys_decomp)
    call addfld ('FULC    ','W/m2    ',pverp,'A','Longwave clear-sky upward flux',phys_decomp)
    call addfld ('FDLC    ','W/m2    ',pverp,'A','Longwave clear-sky downward flux',phys_decomp)
    call addfld ('FLDS    ','W/m2    ',1,    'A','Downwelling longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLDSC   ','W/m2    ',1,    'A','Downwelling longwave clearsky flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNR    ','W/m2    ',1,    'A','Net longwave flux at tropopause',phys_decomp, sampling_seq='rad_lwsw')

    call add_default ('QRL     ', 1, ' ')
    call add_default ('FLNS    ', 1, ' ')
    call add_default ('FLNT    ', 1, ' ')
    call add_default ('FLUT    ', 1, ' ')
    call add_default ('FLUTC   ', 1, ' ')
    call add_default ('FLNTC   ', 1, ' ')
    call add_default ('FLNSC   ', 1, ' ')
    call add_default ('LWCF    ', 1, ' ')

    ! Heating rate needed for d(theta)/dt computation
    call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)     

    call addfld_spectral_intervals
    call circumbinary_addfld_stdoutput
    call init_rad_data()
    call circumbinary_init_orbit

  end subroutine exo_radiation_init


!============================================================================

  subroutine exo_radiation_tend(state, ptend, pbuf, &
                            cam_out, cam_in, &
                            landfrac, landm, icefrac, snowh,&
                            fsns, fsnt, flns, flnt, fsds, net_flx)
!-----------------------------------------------------------------------
! 
! Purpose: Driver for correlated K radiation computation.  Uses delta eddington 
!          two stream at short wavelengths and hemispheric mean two sream 
!          at IR wavelengths.
!
! Revision history:
! 2009?     R. Urata - adapted code from MRAMS/CARMA Mars model
! September 2010: E.T.Wolf   -- ARCHEAN CAM3 RT
! June 2014 : E.T. Wolf   -- EXOPLANET CESM RT
!-----------------------------------------------------------------------
 
    use physics_buffer,    only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use phys_grid,         only: get_rlat_all_p, get_rlon_all_p
    use physics_types,     only: physics_state, physics_ptend
    use time_manager,      only: get_curr_calday, get_nstep, get_curr_calday_rotation
    use camsrfexch,        only: cam_out_t, cam_in_t
    use cam_history,       only: outfld
    use radheat,           only: radheat_tend
    use pmgrid,            only: plev, plevp
    use pspect
    use rad_constituents,  only: rad_cnst_get_gas, rad_cnst_out
    use shr_orb_mod
    use cloud_cover_diags, only: cloud_cover_diags_out
    use cam_control_mod,   only: lambm0, obliqr, eccen, mvelpp
    use radiation_data,    only: output_rad_data
    use spectral_output,   only: outfld_spectral_flux_fullsky, outfld_spectral_flux_clearsky
    use twostars_m,         only: insolation3
    use circumbinary_mod,   only: circumbinary_set_stellar, semia_adj
 
    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in), dimension(pcols) :: landfrac         ! land fraction
    real(r8), intent(in), dimension(pcols) :: landm            ! land fraction ramp
    real(r8), intent(in), dimension(pcols) :: icefrac          ! ice fraction
    real(r8), intent(in), dimension(pcols) :: snowh            ! Snow depth (liquid water equivalent)
    real(r8), intent(inout), dimension(pcols) :: fsns          ! Net solar flux at surface 
    real(r8), intent(inout), dimension(pcols) :: fsnt          ! Net column abs solar flux at model top
    real(r8), intent(inout), dimension(pcols) :: flns          ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout), dimension(pcols) :: flnt          ! Net outgoing lw flux at model top
    real(r8), intent(inout), dimension(pcols) :: fsds            ! Surface solar down flux
    real(r8), intent(inout), dimension(pcols) :: net_flx
    type(physics_state), intent(in), target :: state
    type(physics_ptend), intent(out)        :: ptend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    type(cam_out_t),     intent(inout)      :: cam_out
    type(cam_in_t),      intent(in)         :: cam_in

!------------------------------------------------------------------------
!
!  Local Variables
!

    integer, pointer, dimension(:) :: nmxrgn             ! Number of maximally overlapped regions
    real(r8), pointer, dimension(:,:) :: pmxrgn          ! Maximum values of pressure for each
!                                                        ! maximally overlapped region.
!                                                        ! 0->pmxrgn(i,1) is range of pressure for
!                                                        !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                        !    2nd region, etc
!    real(r8), dimension(pcols,pver) :: cldemis     ! Cloud longwave emissivity
    real(r8), pointer, dimension(:,:) :: cicewp      ! in-cloud cloud ice water path (from param_cldoptics_calc)
    real(r8), pointer, dimension(:,:) :: cliqwp      ! in-cloud cloud liquid water path (from param_cldoptics_calc)
    real(r8), dimension(pcols) ::  cltot           ! Diagnostic total cloud cover
    real(r8), dimension(pcols) ::  cllow           !       "     low  cloud cover
    real(r8), dimension(pcols) ::  clmed           !       "     mid  cloud cover
    real(r8), dimension(pcols) ::  clhgh           !       "     hgh  cloud cover
    real(r8), dimension(pcols,pver) :: ftem        ! Temporary workspace for outfld variables
    real(r8), dimension(pcols,pver) :: ftem2	   ! Second temp workspace

    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cfrc  ! cloud fraction
    real(r8), pointer, dimension(:,:) :: qrs     ! shortwave radiative heating rate 
    real(r8), pointer, dimension(:,:) :: qrl     ! longwave  radiative heating rate 


    integer lchnk, ncol
    real(r8) :: calday                        ! current calendar day
    real(r8), dimension(pcols) :: clat        ! current latitudes(radians)
    real(r8), dimension(pcols) :: clon        ! current longitudes(radians)
    real(r8), dimension(pcols) ::  coszrs     ! Cosine solar zenith angle
    integer  :: nstep                         ! current timestep number
    logical  :: conserve_energy = .true.      ! flag to carry (QRS,QRL)*dp across time steps

    !
    ! Local variables from radctl
    !

    integer :: i, k ,ik             
    integer :: istat
  
    real(r8) :: vis_dir
    real(r8) :: vis_dif
    real(r8) :: nir_dir
    real(r8) :: nir_dif
    real(r8), dimension(pver) :: sw_dTdt     
    real(r8), dimension(pver) :: lw_dTdt     
    real(r8), dimension(pverp) :: sw_upflux   
    real(r8), dimension(pverp) :: sw_dnflux   
    real(r8), dimension(pverp) :: lw_upflux   
    real(r8), dimension(pverp) :: lw_dnflux   
    real(r8), dimension(pverp,ntot_wavlnrng) :: sw_upflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: sw_dnflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: lw_upflux_spec
    real(r8), dimension(pverp,ntot_wavlnrng) :: lw_dnflux_spec
    real(r8), dimension(pcols) :: fsntoa        ! Net solar flux at TOA
    real(r8), dimension(pcols) :: fsntoac       ! Clear sky net solar flux at TOA
    real(r8), dimension(pcols) :: fsnirt        ! Near-IR flux absorbed at toa
    real(r8), dimension(pcols) :: fsnrtc        ! Clear sky near-IR flux absorbed at toa
    real(r8), dimension(pcols) :: fsnirtsq      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8), dimension(pcols) :: fsntc         ! Clear sky total column abs solar flux
    real(r8), dimension(pcols) :: fsnsc         ! Clear sky surface abs solar flux
    real(r8), dimension(pcols) :: fsdsc         ! Clear sky surface downwelling solar flux
    real(r8), dimension(pcols) :: flut          ! Upward flux at top of model
    real(r8), dimension(pcols) :: lwcf          ! longwave cloud forcing
    real(r8), dimension(pcols) :: swcf          ! shortwave cloud forcing
    real(r8), dimension(pcols) :: flutc         ! Upward Clear Sky flux at top of model
    real(r8), dimension(pcols) :: flntc         ! Clear sky lw flux at model top
    real(r8), dimension(pcols) :: flnsc         ! Clear sky lw flux at srf (up-down)
    real(r8), dimension(pcols) :: fln200        ! net longwave flux interpolated to 200 mb
    real(r8), dimension(pcols) :: fln200c       ! net clearsky longwave flux interpolated to 200 mb
    real(r8), dimension(pcols,pverp) :: fns     ! net shortwave flux
    real(r8), dimension(pcols,pverp) :: fcns    ! net clear-sky shortwave flux
    real(r8), dimension(pcols) :: fsn200        ! fns interpolated to 200 mb
    real(r8), dimension(pcols) :: fsn200c       ! fcns interpolated to 200 mb
    real(r8), dimension(pcols,pverp) :: fsn     ! net shortwave flux
    real(r8), dimension(pcols,pverp) :: fln     ! net longwave flux
    real(r8), dimension(pcols,pverp) :: fnl     ! net longwave flux
    real(r8), dimension(pcols,pverp) :: fcnl    ! net clear-sky longwave flux

    ! for circumbinary
    real(r8), dimension(pcols,pver) :: qrs1, qrs2, qrs_tot     ! shortwave radiative heating rate 
    real(r8), dimension(pcols,pver) :: qrl1, qrl2, qrl_tot     ! longwave radiative heating rate 
    real(r8), dimension(pcols,pver) :: fus1, fus2, fus_tot
    real(r8), dimension(pcols,pver) :: fds1, fds2, fds_tot
    real(r8), dimension(pcols) :: fdstoa1, fdstoa2, fdstoa_tot

    real(r8), dimension(pcols) :: fsns1, fsns2
    real(r8), dimension(pcols) :: fsnt1, fsnt2
    real(r8), dimension(pcols) :: fsds1, fsds2


    real(r8) :: eccf                            ! Earth/sun distance factor
    real(r8) :: delta                           ! Solar declination angle
    
    real(r8), pointer, dimension(:,:) :: h2ommr   ! h2o   mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2mmr   ! co2   mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4mmr   ! ch4   mass mixing ratio
    real(r8), dimension(pcols,pver) :: h2mmr    ! h2    mass mixing ratio
    real(r8), dimension(pcols,pver) :: n2mmr    ! n2    mass mixing ratio

    ! null cloud place holders, for clear sky calculation
    real(r8), dimension(pcols,pver) :: cicewp_zero
    real(r8), dimension(pcols,pver) :: cliqwp_zero
    real(r8), dimension(pcols,pver) :: cfrc_zero


    !the following would have dimension nazm_tshadow if shadows were taken into effect
    integer, parameter :: ext_nazm_tshadow = 1                  ! take shadows into effect
    real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon   ! cos of zenith angle of horizon 
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct
    !------------------------------------
    integer :: ext_tslas_tog
    integer :: ext_tshadow_tog
  
    real(r8) :: ext_solar_azm_ang
    real(r8) :: ext_tazm_ang
    real(r8) :: ext_tslope_ang

    real(r8) :: ext_msdist
    real(r8) :: ext_rtgt

    ! summed fluxes
    ! broadband fluxes
    real(r8), dimension(pcols,pverp) :: lwup_rad 
    real(r8), dimension(pcols,pverp) :: lwdown_rad 
    real(r8), dimension(pcols,pverp) :: swup_rad
    real(r8), dimension(pcols,pverp) :: swdown_rad
    ! spectral fluxes
    real(r8), dimension(pcols,pverp,ntot_wavlnrng) :: lwup_rad_spec
    real(r8), dimension(pcols,pverp,ntot_wavlnrng) :: lwdown_rad_spec
    real(r8), dimension(pcols,pverp,ntot_wavlnrng) :: swup_rad_spec
    real(r8), dimension(pcols,pverp,ntot_wavlnrng) :: swdown_rad_spec

    real(r8) :: frac_day
    real(r8) :: total_days
    logical :: do_exo_rad

    ! circumbinary properties
    real(r8),dimension(1:3)::mass 
    real(r8),dimension(1:2):: tf, d2, decl, phi
    real(r8):: e2
    integer :: star    
    real(r8), dimension(pver) ::  dzc   ! layer mass
    real(r8), dimension(pcols,pver) :: dzc_matrix
    real(r8) :: lyr_mass_fact, fdiv_sw
 
!------------------------------------------------------------------------
!
! Start Code
!

    lchnk = state%lchnk
    ncol = state%ncol

    calday = get_curr_calday()

    itim = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cfrc,    start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, qrs_idx, qrs)
    call pbuf_get_field(pbuf, qrl_idx, qrl)
    call pbuf_get_field(pbuf, rel_idx, rel)
    call pbuf_get_field(pbuf, rei_idx, rei)
    call pbuf_get_field(pbuf, cicewp_idx, cicewp)
    call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
    call pbuf_get_field(pbuf, pmxrgn_idx, pmxrgn)  !! Kludge for model integration
    call pbuf_get_field(pbuf, nmxrgn_idx, nmxrgn)  !! REMOVE?

    !
    ! Cosine solar zenith angle for current time step
    !
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)

    ! Wolf, length of day scaling for zenith angle calculation
    call get_curr_calday_rotation(frac_day, total_days)
    !if(masterproc) write(*,*) "Frac day :", frac_day, calday, total_days

    ! calculate Earth-Sun distance factor, scaled by eccentricity factor
    ! in the next version the orbital details will be more heavily modulated.
    ! NOTE:  SHR_CONST_MSDIST = normalized planet-star distance squared

    if (do_exo_circumbinary) then
      !this needs to be called everytime step.
      mass(:) = cb_mass(:)
      mass(3) = mass(3) *SHR_CONST_MEARTH/SHR_CONST_MSUN
      call insolation3(mass,cb_teff,cb_radius,semia_adj,cb_eccen,cb_man,cb_spin,total_days, & ! inputs
                  d2, tf, decl, phi, e2)             ! outputs

!      if(masterproc) write(*,*) "circumbinary, d2 :", d2
!      if(masterproc) write(*,*) "circumbinary, tf :", tf
!      if(masterproc) write(*,*) "circumbinary, decl :", decl
!      if(masterproc) write(*,*) "circumbinary, phi :", phi
!      if(masterproc) write(*,*) "circumbinary: ", calday, phi, phi/(2.*SHR_CONST_PI)
       phi(:) = 1.0 - phi(:)/(2.*SHR_CONST_PI) ! set hour angle
      if(masterproc) write(*,*) "circumbinary, phi :", phi
!      if(masterproc) write(*,*) "circumbinary, phi/2pi :", phi
      ext_msdist=1.0
    else
      ! do standard orbital calculation, for determining sflux_frac
      call zenith_rotation (frac_day, calday, clat, clon, coszrs, ncol)
      !call zenith (calday, clat, clon, coszrs, ncol)
      ! shr_orb_decl is also called in zenith to get delta
      ! in this formulation, I am only calling the shr_orb_decl to get eccf
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
      ext_msdist=1.0/eccf
    endif

    !! Main calculation starts here !!
    do_exo_rad = exo_radiation_do()
    !write(*,*) "do_exo_rad", do_exo_rad
    !do_exo_rad determines if RT called on a given timesteps
 
    if (do_exo_rad) then 
      !write(*,*) "inside RT if statement"
										
      !-------for now, set these topography angles to 0 (ie, no topography blocking sun)
      !-------should be moved to inside do loop if we want to use them
      !-------azimuth can be calculated from cos(azim) = (cos(hr_ang)*cos(dec)*sin(lat)-sin(dec)*cos(lat))/cos(elev)
      !-------the azimuth is the angle east of south when the hour angle, h, is negative (morning)
      !-------and the angle west of south when the hour angle, h, is positive (afternoon).
      ext_solar_azm_ang = 0.     !solar azimuthal angle? should not be zero
      ext_tazm_ang = 0.
      ext_tslope_ang = 0.
      ext_tslas_tog = 0          
      ext_tshadow_tog = 1        ! toggle shadowing 
      ext_cosz_horizon(:) = 0.
      ext_TCx_obstruct(:) = 0.
      ext_TCz_obstruct(:) = 0.
      ext_rtgt = 1.  

      sw_dTdt(:) = 0.    ! Initialize heating rate arrays
      lw_dTdt(:) = 0.    !
 
      lw_dnflux(:) = 0.   !
      lw_upflux(:) = 0.   ! Initialize entire arrays for summing below
      sw_upflux(:) = 0.   !
      sw_dnflux(:) = 0.   !

      lwup_rad_spec(:,:,:) = 0.
      lwdown_rad_spec(:,:,:) = 0.
      swup_rad_spec(:,:,:) = 0.
      swdown_rad_spec(:,:,:) = 0.

      qrs(:,:) = 0.
      qrl(:,:) = 0.
      lwup_rad(:,:) = 0.
      lwdown_rad(:,:) = 0.
      swup_rad(:,:) = 0.
      swdown_rad(:,:) = 0.

      nstep = get_nstep()

      ! Native CAM functions; returns pointer to mass mixing ratio for the gas specified        
      call rad_cnst_get_gas(0,'CO2', state, pbuf,  co2mmr)
      call rad_cnst_get_gas(0,'CH4', state, pbuf,  ch4mmr)
      call rad_cnst_get_gas(0,'H2O', state, pbuf,  h2ommr) !H2O specific humidity

      ! well mixed species from exoplanet_mod.F90
      n2mmr(:,:)  = exo_n2mmr
      h2mmr(:,:)  = exo_h2mmr

      cam_out%sols(:)  = 0.0 
      cam_out%soll(:)  = 0.0
      cam_out%solsd(:) = 0.0
      cam_out%solld(:) = 0.0


      !!!! BINARY STAR LOOP BEGIN !!!
      do star=1,2 

        if (cb_transit .eqv. .false.) then
          tf(1) = 1.0
          tf(2) = 1.0
        end if
        ! set up stellar spectra for each radiative step
        ! for each star, reset the stellar flux array (gw_solflux)
        call circumbinary_set_stellar(star, d2, tf)

        ! Compute local cosine solar zenith angle,
        do i=1,ncol
          !coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), decl(star) )
          !coszrs(i) = shr_orb_cosz( frac_day, clat(i), clon(i), decl(star) )
          coszrs(i) = shr_orb_cosz(phi(star), clat(i), clon(i), decl(star) )
        end do
!write(*,*) "Binary, coszrs: ", coszrs

        ! ------- standard radiative transfer logic ---------
        ! Do a parallel clearsky radiative calculation so we can calculate cloud forcings
        ! Setting do_exo_rt_clearsky to true, slows the code dramatically, use wisely and sparingly
        if (do_exo_rt_clearsky) then

          ! set clouds to zero everywhere
          cicewp_zero(:,:) = 0.0
          cliqwp_zero(:,:) = 0.0
          cfrc_zero(:,:) = 0.0

          do i = 1, ncol

            call aerad_driver(h2ommr(i,:), co2mmr(i,:), ch4mmr(i,:) &
                             ,h2mmr(i,:), n2mmr(i,:) &
                             ,cicewp_zero(i,:), cliqwp_zero(i,:), cfrc_zero(i,:) &
                             ,rei(i,:), rel(i,:) &
                             ,cam_in%ts(i), state%ps(i), state%pmid(i,:) &
                             ,state%pdel(i,:), state%pdeldry(i,:), state%t(i,:), state%pint(i,:), state%pintdry(i,:) &
                             ,coszrs(i), ext_msdist &
                             ,cam_in%asdir(i), cam_in%aldir(i) &
                             ,cam_in%asdif(i), cam_in%aldif(i) &
                             ,ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang  &
                             ,ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon  &
                             ,ext_TCx_obstruct, ext_TCz_obstruct, state%zi(i,:) &
                             ,sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux  &
                             ,lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec &       
                             ,vis_dir, vis_dif, nir_dir, nir_dif, dzc )
                           
            ftem(i,:) = sw_dTdt(:)       
            ftem2(i,:) = lw_dTdt(:)      

            lwup_rad(i,:) = lw_upflux(:)
            lwdown_rad(i,:) = lw_dnflux(:)
            swup_rad(i,:) = sw_upflux(:)
            swdown_rad(i,:) = sw_dnflux(:)

            if (do_exo_rt_spectral) then
              lwup_rad_spec(i,:,:)   = lw_upflux_spec(:,:)
              lwdown_rad_spec(i,:,:) = lw_dnflux_spec(:,:)
              swup_rad_spec(i,:,:)   = sw_upflux_spec(:,:)
              swdown_rad_spec(i,:,:) = sw_dnflux_spec(:,:)
            endif

            ! Fluxes sent to land model
            ! Note these values are overwritten by the full-sky rt calc
            !cam_out%sols(i) = vis_dir
            !cam_out%soll(i) = nir_dir
            !cam_out%solsd(i) = vis_dif
            !cam_out%solld(i) = nir_dif
            !cam_out%flwds(i) = lw_dnflux(pverp)
        
          enddo   ! ncol loop

          fsn(:,:) = swdown_rad(:,:) - swup_rad(:,:)
          fln(:,:) = lwup_rad(:,:) - lwdown_rad(:,:)
          fsns(:) = fsn(:,pverp)
          flns(:) = fln(:,pverp)
          fsnt(:) = fsn(:,1)
          flnt(:) = fln(:,1)
          fsds(:) = swdown_rad(:,pverp)
          qrs(:ncol,:pver) = ftem(:ncol,:pver)
          qrl(:ncol,:pver) = ftem2(:ncol,:pver)
 
          call outfld('QRSC     ',qrs*SHR_CONST_CDAY  , pcols,lchnk)    ! [K/day]
          call outfld('FSDSC    ',fsds  ,pcols,lchnk)
          call outfld('FSNTC    ',fsnt  ,pcols,lchnk)
          call outfld('FSNSC    ',fsns  ,pcols,lchnk)
          call outfld('QRLC     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
          call outfld('FLNTC    ',flnt  ,pcols,lchnk)
          call outfld('FLUTC    ',lwup_rad(:,2)  ,pcols,lchnk)
          call outfld('FLNSC    ',flns  ,pcols,lchnk)
          call outfld('FULC     ',lwup_rad, pcols, lchnk)
          call outfld('FDLC     ',lwdown_rad, pcols, lchnk)
          call outfld('FUSC     ',swup_rad, pcols, lchnk)
          call outfld('FDSC     ',swdown_rad, pcols, lchnk)
          call outfld('SOLSC    ',cam_out%sols  ,pcols,lchnk)  
          call outfld('SOLLC    ',cam_out%soll  ,pcols,lchnk)
          call outfld('SOLSDC   ',cam_out%solsd ,pcols,lchnk)
          call outfld('SOLLDC   ',cam_out%solld ,pcols,lchnk)

          if (do_exo_rt_spectral) call outfld_spectral_flux_clearsky(lchnk, lwdown_rad_spec, lwup_rad_spec, swup_rad_spec, swdown_rad_spec)
        
        endif  ! (do_exo_rt_clearsky)
!-----------------------
! Binary star updates only implemented into full-sky calculation currently.

        ! Do Column Radiative transfer calculation WITH clouds.  
        do i = 1, ncol

          call aerad_driver(h2ommr(i,:), co2mmr(i,:), ch4mmr(i,:) &
                           ,h2mmr(i,:), n2mmr(i,:) &
                           ,cicewp(i,:), cliqwp(i,:), cfrc(i,:) &
                           ,rei(i,:), rel(i,:) &
                           ,cam_in%ts(i), state%ps(i), state%pmid(i,:) &
                           ,state%pdel(i,:), state%pdeldry(i,:), state%t(i,:), state%pint(i,:), state%pintdry(i,:) &
                           ,coszrs(i), ext_msdist &
                           ,cam_in%asdir(i), cam_in%aldir(i) &
                           ,cam_in%asdif(i), cam_in%aldif(i) &
                           ,ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang  &
                           ,ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow, ext_cosz_horizon  &
                           ,ext_TCx_obstruct, ext_TCz_obstruct, state%zi(i,:) &
                           ,sw_dTdt, lw_dTdt, lw_dnflux, lw_upflux, sw_upflux, sw_dnflux  &
                           ,lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec &       
                           ,vis_dir, vis_dif, nir_dir, nir_dif, dzc ) 


          dzc_matrix(i,:) = dzc(:)

          ftem(i,:) = sw_dTdt(:)       
          ftem2(i,:) = lw_dTdt(:)      

          lwup_rad(i,:) = lw_upflux(:)
          lwdown_rad(i,:) = lw_dnflux(:)
          swup_rad(i,:) = sw_upflux(:)
          swdown_rad(i,:) = sw_dnflux(:)
 
          if (do_exo_rt_spectral) then
            lwup_rad_spec(i,:,:)   = lw_upflux_spec(:,:)
            lwdown_rad_spec(i,:,:) = lw_dnflux_spec(:,:)
            swup_rad_spec(i,:,:)   = sw_upflux_spec(:,:)
            swdown_rad_spec(i,:,:) = sw_dnflux_spec(:,:)
          endif

          ! Fluxes sent to land model
          ! Note these values overwrite those use in do_exo_rt_clearsky calc
          cam_out%sols(i) = cam_out%sols(i) + vis_dir
          cam_out%soll(i) = cam_out%soll(i) + nir_dir
          cam_out%solsd(i) = cam_out%solsd(i) + vis_dif
          cam_out%solld(i) = cam_out%solld(i) + nir_dif
          cam_out%flwds(i) = lw_dnflux(pverp)

        enddo   ! ncol loop

        fsn(:,:) = swdown_rad(:,:) - swup_rad(:,:)
        fln(:,:) = lwup_rad(:,:) - lwdown_rad(:,:)
        fsns(:) = fsn(:,pverp)
        flns(:) = fln(:,pverp)
        fsnt(:) = fsn(:,1)
        flnt(:) = fln(:,1)
        fsds(:) = swdown_rad(:,pverp)
        qrs(:ncol,:pver) = ftem(:ncol,:pver)
        qrl(:ncol,:pver) = ftem2(:ncol,:pver)
       
        if (star == 1) then 
          qrs1(:ncol,:pver) = qrs(:ncol,:pver)
          fus1(:ncol,:pverp) = swup_rad(:ncol,:pverp)
          fds1(:ncol,:pverp) = swdown_rad(:ncol,:pverp)
          fdstoa1(:ncol) = swdown_rad(:ncol,1) 
          fsns1(:) = fsns(:)  
          fsnt1(:) = fsnt(:)
          fsds1(:) = fsds(:)

          call outfld('QRS_cb1     ',qrs*SHR_CONST_CDAY  , pcols,lchnk)    ! [K/day]
          call outfld('QRL_cb1     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
          call outfld('FUL_cb1     ',lwup_rad, pcols, lchnk)
          call outfld('FDL_cb1     ',lwdown_rad, pcols, lchnk)
          call outfld('FUS_cb1     ',swup_rad, pcols, lchnk)
          call outfld('FDS_cb1     ',swdown_rad, pcols, lchnk)
          call outfld('FDSTOA_cb1     ',swdown_rad(:,1), pcols, lchnk)

!          call outfld('FSDS_cb1    ',fsds  ,pcols,lchnk)
!          call outfld('FSNT_cb1    ',fsnt  ,pcols,lchnk)
!          call outfld('FSNS_cb1    ',fsns  ,pcols,lchnk)
!          call outfld('SOLS_cb1    ',cam_out%sols  ,pcols,lchnk)
!          call outfld('SOLL_cb1    ',cam_out%soll  ,pcols,lchnk)
!          call outfld('SOLSD_cb1   ',cam_out%solsd ,pcols,lchnk)
!          call outfld('SOLLD_cb1   ',cam_out%solld ,pcols,lchnk)
        endif
        if (star == 2) then 
          qrs2(:ncol,:pver) = qrs(:ncol,:pver)
          fus2(:ncol,:pver) = swup_rad(:ncol,:pverp)
          fds2(:ncol,:pver) = swdown_rad(:ncol,:pverp)
          fsns2(:) = fsns(:)  
          fsnt2(:) = fsnt(:)
          fsds2(:) = fsds(:)
          fdstoa2(:ncol) = swdown_rad(:ncol,1) 

          call outfld('QRS_cb2     ',qrs*SHR_CONST_CDAY  , pcols,lchnk)    ! [K/day]
          call outfld('QRL_cb2     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
          call outfld('FUL_cb2     ',lwup_rad, pcols, lchnk)
          call outfld('FDL_cb2     ',lwdown_rad, pcols, lchnk)
          call outfld('FUS_cb2     ',swup_rad, pcols, lchnk)
          call outfld('FDS_cb2     ',swdown_rad, pcols, lchnk)
          call outfld('FDSTOA_cb2     ',swdown_rad(:,1), pcols, lchnk)

!          call outfld('FSDS_cb2    ',fsds  ,pcols,lchnk)
!          call outfld('FSNT_cb2    ',fsnt  ,pcols,lchnk)
!          call outfld('FSNS_cb2    ',fsns  ,pcols,lchnk)


!          call outfld('SOLS_cb2    ',cam_out%sols  ,pcols,lchnk)
!          call outfld('SOLL_cb2    ',cam_out%soll  ,pcols,lchnk)
!          call outfld('SOLSD_cb2   ',cam_out%solsd ,pcols,lchnk)
!          call outfld('SOLLD_cb2   ',cam_out%solld ,pcols,lchnk)


          qrs_tot(:ncol,:pver) = qrs1(:ncol,:pver) + qrs2(:ncol,:pver)
          fds_tot(:ncol,:pver) = fds1(:ncol,:pver) + fds2(:ncol,:pver)
          fus_tot(:ncol,:pver) = fus1(:ncol,:pver) + fus2(:ncol,:pver)
          fsns(:) = fsns1(:) + fsns2(:)
          fsnt(:) = fsnt1(:) + fsnt2(:)
          fsds(:) = fsds1(:) + fsds2(:)
          fdstoa_tot(:) = fdstoa1(:) + fdstoa2(:)

          call outfld('FDSTOA     ',fdstoa_tot, pcols, lchnk)
          call outfld('FUS     ',fus_tot, pcols, lchnk)
          call outfld('FDS     ',fds_tot, pcols, lchnk)
          call outfld('QRL     ',qrl*SHR_CONST_CDAY   ,pcols,lchnk)    ! [K/day]
          call outfld('QRS     ',qrs_tot*SHR_CONST_CDAY , pcols,lchnk)    ! [K/day]
          call outfld('FSDS    ',fsds  ,pcols,lchnk)
          call outfld('FSNT    ',fsnt  ,pcols,lchnk)
          call outfld('FSNS    ',fsns  ,pcols,lchnk)
          call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
          call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
          call outfld('SOLSD   ',cam_out%solsd   ,pcols,lchnk)
          call outfld('SOLLD   ',cam_out%solld   ,pcols,lchnk)
          call outfld('FLNT    ',flnt  ,pcols,lchnk)
          call outfld('FLUT    ',lwup_rad(:,2)  ,pcols,lchnk)
          call outfld('FLNS    ',flns  ,pcols,lchnk)
        endif
 !       if (do_exo_rt_spectral) call outfld_spectral_flux_fullsky(lchnk, lwdown_rad_spec, lwup_rad_spec, swup_rad_spec, swdown_rad_spec)
 
        ! Cloud cover diagnostics
        call cloud_cover_diags_out(lchnk, ncol, cfrc, state%pmid, nmxrgn, pmxrgn )


      enddo   ! star loop
    !!!! BINARY STAR LOOP END !!!

    ! recalculate heating and cooling rates based on summed fluxes
!    do i=1,ncol
!      do k=2,pverp      
!        lyr_mass_fact = dzc_matrix(i,k-1)*cpair
!        fdiv_sw =          (fus_tot(i,k)-fds_tot(i,k))-(fus_tot(i,k-1)-fds_tot(i,k-1))
!        qrs(i,k-1) = fdiv_sw/lyr_mass_fact      ! "shortwave" heating rate [K/s]
!      enddo
!    enddo
!    call outfld('QRS_test     ',qrs*SHR_CONST_CDAY , pcols,lchnk)    ! [K/day]

    else ! if (do_exo_rad) then 
 



      ! convert radiative heating rates to Q*dp for energy conservation
      if (conserve_energy) then
      !DIR$ CONCURRENT
        do k =1 , pver
      !DIR$ CONCURRENT
          do i = 1, ncol
            qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            qrl(i,k) = qrl(i,k)/state%pdel(i,k)
          enddo
        enddo
      endif  ! (conserve_energy)

    end if   ! (do_exo_rad)

    ! output rad inputs and resulting heating rates
    call output_rad_data(pbuf, state, cam_in, landm, coszrs(i))

    ! Compute net radiative heating tendency
    call radheat_tend(state, pbuf,  ptend, qrl*cpair, qrs_tot*cpair, fsns, &
                      fsnt, flns, flnt, cam_in%asdir, net_flx)

    ! Compute heating rate for dtheta/dt
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs_tot(i,k) + qrl(i,k))/cpair * &
                      (SHR_CONST_PSTD/state%pmid(i,k))**cappa
       enddo
    enddo
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
!DIR$ CONCURRENT
      do k =1 , pver
!DIR$ CONCURRENT
        do i = 1, ncol
          qrs(i,k) = qrs(i,k)*state%pdel(i,k)
          qrl(i,k) = qrl(i,k)*state%pdel(i,k)
        enddo
      enddo
    endif ! (conserve_energy)     

  end subroutine exo_radiation_tend

!============================================================================

  function exo_radiation_do(timestep)

!------------------------------------------------------------------------
!
! Purpose:  Returns true if the exo_rt is done this timestep
!
!------------------------------------------------------------------------

    integer, intent(in), optional :: timestep
    logical :: exo_radiation_do

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: nstep
!------------------------------------------------------------------------
!
! Start Code
!
  if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

!   write(*,*)  "exo_radiation_do ", nstep, exo_rad_step, mod(nstep-1,exo_rad_step)
!    exo_radiation_do = nstep == 0  .or.  exo_rad_step == 1                     &
!                       .or. (mod(nstep-1,exo_rad_step) == 0  .and.  nstep /= 1)

    exo_radiation_do = nstep == 0  .or.  exo_rad_step == 1                     &
                       .or. (mod(nstep,exo_rad_step) == 0  .and.  nstep /= 1)

  end function exo_radiation_do

!============================================================================

real(r8) function exo_radiation_nextsw_cday()

!-----------------------------------------------------------------------
! Purpose: Returns calendar day of next sw radiation calculation
!          This is used to ensure surface albedo claculations are sync'd
!          with exo radiation calls.
!-----------------------------------------------------------------------

   use time_manager, only: get_curr_calday, get_nstep, get_step_size


   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of
   !-----------------------------------------------------------------------

   exo_radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (exo_radiation_do(nstep)) then
         exo_radiation_nextsw_cday = get_curr_calday(offset=offset)
         dosw = .true.
      end if
   end do
   if(exo_radiation_nextsw_cday == -1._r8) then
      call endrun('error in exo_radiation_nextsw_cday')
   end if

end function exo_radiation_nextsw_cday


end module exo_radiation_cam_intr
