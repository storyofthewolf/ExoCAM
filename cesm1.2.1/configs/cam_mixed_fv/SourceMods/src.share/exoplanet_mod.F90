

module exoplanet_mod

!------------------------------------------------------------------------------
! Purpose:  Stores basic exoplanet parameters -- set before compiling
!           These options override any other NCAR switches/options, etc when
!           this package is placed in SourceMods
!------------------------------------------------------------------------------
! NOTES ::

  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  private
  save

  !! ============== RUN OPTIONS ============== !!
  logical, public, parameter :: do_exo_atmconst = .true.        !! Read gas constituents of atmosphere from this file 
                                                                !! Overrides ghg namelist options
  logical, public, parameter :: do_exo_rt = .true.              !! .true. = use correlated-k "exoplanet" RT
  	   	   	     		    		        !! .false. = use CAM4 RT, refer elsewhere for operating
                                                                !! (currently only works for true.  delete option?)
  logical, public, parameter :: do_exo_synchronous = .false.    !! Toggle to synchronous rotation mode
                                                                !! Eventually replace with better orbital computer
  logical, public, parameter :: do_carma_exort = .false.        !! Set to true only if running with the CARMA microphyiscs package
                                                                !! and linking aerosol absorption to ExoRT  
  logical, public, parameter :: do_exo_gw = .false.             !! flag to turn on gravity waves.  Note, present gw wave parameterization    
                                                                !! does not work for low pressure atmospheres.                  
  
  !! ============== RADIATION OPTIONS  ============== !!
  integer, public, parameter :: exo_rad_step = 3                !! freq. of radiation calc in time steps (positive)
                                                                !! or hours (negative).
  logical, public, parameter :: do_exo_rt_clearsky = .false.    !! Do parallel clearsky radiative calculation for exo rt
                                                                !! Slow, use sparingly
  logical, public, parameter :: do_exo_rt_spectral = .false.    !! collect and output spectrally resolved radiative fluxes
                                                                !! 
  !! Radiation  Spectral Band Optimization  !!
  !! if .false. LW and SW streams computed over fill bandpasses
  !! if .true. LW and SW bandpasses reduced to fit SED and assumed Tmax 
  logical,  public, parameter :: do_exo_rt_optimize_bands = .true.
  real(r8), public, parameter :: Tmax = 400.          !! Maximum expected temperature for thermal band optimization
  real(r8), public, parameter :: swFluxLimit = 0.999  !! Fraction of stellar flux captured in bands, rescaled
  real(r8), public, parameter :: lwFluxLimit = 0.999  !! Fraction of thermal flux captured in bands, not rescaled

  !! ==============  PLANET PARAMETERS  ============== !!
  !! It is the responsibility of the USER to set do_exo_synchronous, exo_ndays
  !! exo_sdays, exo_scon, and exo_solar_file in a self-consistent manner.  It
  !! is good practice to check that your model output incident stellar flux (FDS)
  !! looks correct for your assumed dirunal cycle.
  !!
  !! exo_sdays is the length of the sidereal period in seconds.
  !! exo_ndays is a scaler given in units of Earth days 
  !! exo_porb is the orbital period, optional input used to determine
  !! 	      the sidereal period, exo_sdays, for asynchronous rotators
  !!
  !! if (do_exo_synchronous = .true.) synchronous rotation is assumed. The
  !! zenith angle is fixed at 180 deg longitude, and the Coriolis parameter is
  !! set via, exo_sdays [sec] = 86400.0 [sec/Earth day] * exo_ndays [Earth days].
  !! There is no diurnal cycle.
  !!
  !! if (do_exo_synchronous = .false.) non-synchronous rotation occurs.  The
  !! Coriolis parameter is set via exo_sdays [sec], and the diurnal cycle
  !! (length of solar day) is then set explicitly by exo_ndays [Earth days].
  !! Note, that the length of the diurnal period does not equal to the
  !! sidereal period!  Take care in setting exo_ndays and exo_sdays!
  !!
  !! See examples below.   
  
  !! Generic
  real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8     !! radius ~ m
  real(r8), public, parameter :: exo_surface_gravity = 9.81_R8          !! gravity ~ m/s^2            
  real(r8), public, parameter :: exo_ndays           = 1.00_R8          !! scaler to number of Earth days.
  real(r8), public, parameter :: exo_porb            = 365._R8          !! orbital period, for obliquity cycle, and optionally for exo_sday
  real(r8), public, parameter :: exo_sday            = 86164.0_r8       !! sidereal period Earth
  !real(r8), public, parameter :: exo_sday            = 86400.0_r8 * exo_ndays   !! sidereal period, for synchronous rotator
  !real(r8), public, parameter :: exo_sday            = 86400.0_r8 * exo_ndays / (1._r8 + exo_ndays/exo_porb)  !! sidereal period [sec]

  !! Examples
  !! Earth
  !real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8     !! radius ~ m
  !real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8       !! gravity ~ m/s^2
  !real(r8), public, parameter :: exo_ndays           = 1.0_R8           !! scaler to number of Earth days.
  !real(r8), public, parameter :: exo_porb            = 365.0_R8
  !real(r8), public, parameter :: exo_sday            = 86164.0_R8       !! sidereal period [sec], for Earth value = 86164.0 

  !! Earth - slow asynchronous rotator
  !real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8     !! radius ~ m
  !real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8       !! gravity ~ m/s^2
  !real(r8), public, parameter :: exo_ndays           = 10.0_R8          !! scaler to number of Earth days.
  !real(r8), public, parameter :: exo_porb            = 365.0_R8         !! orbital period 
  !real(r8), public, parameter :: exo_sday            = 86400.0_r8 * exo_ndays / (1._r8 + exo_ndays/exo_porb)  !! sidereal period [sec]

  !! Trappist-1e  (Gillon et al. 2017)  --  synchronous rotator
  !real(r8), public, parameter :: exo_planet_radius    = 5.84878e6_R8    !! radius ~ m
  !real(r8), public, parameter :: exo_surface_gravity  = 7.22925_R8      !! gravity ~ m/s^2
  !real(r8), public, parameter :: exo_ndays            = 6.099615_r8               !! scaler to number of Earth days.   
  !real(r8), public, parameter :: exo_porb             = exo_ndays                   !! orbital period 
  !real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays      !! sidereal period, synchronous rotator

  !! if set user_nl_cpl::orb_iyear = -1
  real(r8), public, parameter :: exo_eccen = 0.0_r8   ! eccentricity
  real(r8), public, parameter :: exo_obliq = 0.0_r8   ! obliquity [degrees]
  real(r8), public, parameter :: exo_mvelp = 0.0_r8   ! vernal equinox

    
  !! ============== STELLAR OPTIONS ============== !!
  !! SOLAR CONSTANT
  real(r8), public, parameter :: exo_scon = 1360.0_r8         ! Solar constant (W m-2)

  ! SOLAR SPECTRAL FILE
  !! Make sure solar file matches spectral intervals for selected RT configuration !!
  character(len=256), public, parameter :: exo_solar_file = '/gpfsm/dnb53/etwolf/models/ExoRT/data/solar/G2V_SUN_n68.nc'


  !! ============== ATMOSPHERIC CONSTITUENT PARAMETERS ============== !!
  !! Activated only if (do_exo_atmconst = .true.) 
  !! Initial conditions file (ncdata) must (approximately) match the total pressure !!
  real(r8), public, parameter :: exo_co2bar  = 0.01_r8                       ! CO2 inventory (bar)
  real(r8), public, parameter :: exo_ch4bar  = 1.0e-3_r8                     ! CH4 inventory (bar)
  real(r8), public, parameter :: exo_c2h6bar = 0.0_r8                        ! C2H6 inventory (bar)
  real(r8), public, parameter :: exo_h2bar   = 0.0_r8                        ! H2 inventory (bar)
  real(r8), public, parameter :: exo_o2bar   = 0.0_r8                        ! O2 inventory (bar)
  real(r8), public, parameter :: exo_n2bar   = 1.0 - exo_co2bar - exo_ch4bar - exo_c2h6bar - exo_o2bar   ! N2 inventory (bar)
  real(r8), public, parameter :: exo_pstd    = (exo_n2bar + exo_o2bar + exo_h2bar + exo_co2bar + exo_ch4bar + exo_c2h6bar)*1.0e5  ! total pressure (Pascals)


  !! ============== OCEAN ALBEDO CONSTANTS ============== !!
  real(r8), public, parameter :: exo_albdif = 0.06 ! 60 deg reference albedo, diffuse (default = 0.06)
  real(r8), public, parameter :: exo_albdir = 0.06 ! 60 deg reference albedo, direct  (default = 0.07)


  !! ============== LAND ALBEDO CONSTANTS ============== !!
  !! Used for soil_color = 21 in the land model (if active)
  !! The input surface data set (fsurdat) contains the soil_colors
  real(r8), public, parameter :: exo_lnd_albvis_dry = 0.30    ! land albedo, visible, dry soil
  real(r8), public, parameter :: exo_lnd_albifr_dry = 0.30    ! land albedo, infrared, dry soil
  real(r8), public, parameter :: exo_lnd_albvis_sat = 0.30    ! land albedo, visible, saturated soil
  real(r8), public, parameter :: exo_lnd_albifr_sat = 0.30    ! land albedo, infrared, saturated soil


  !! ===================================================================== !!
  !! ======================= FUNDAMENTAL CONSTANTS  ====================== !!
  !! ===================================================================== !!
  !! No modifications below this point
  !! 

  ! molecular weights
  real(r8), parameter :: mwn2   = 28._r8
  real(r8), parameter :: mwh2   = 2._r8
  real(r8), parameter :: mwco2  = 44._r8
  real(r8), parameter :: mwo2   = 32._r8
  real(r8), parameter :: mwch4  = 16._r8
  real(r8), parameter :: mwc2h6 = 30._r8
  real(r8), parameter :: cpn2   = 1.039e3_r8
  real(r8), parameter :: cph2   = 14.32e3_r8
  real(r8), parameter :: cpco2  = 0.846e3_r8
  real(r8), parameter :: cpo2   = 0.918e3_r8
  real(r8), parameter :: cpch4  = 2.226e3
  real(r8), parameter :: cpc2h6 = 1.756e3

  !! DERIVED CONSTANTS -- DO NOT MODIFY
  !! automatically calculated from above inputs in bar
  ! dry volume mixing ratios  kg/kmole
  real(r8), public, parameter :: exo_n2vmr   = exo_n2bar   / (exo_pstd/1.0e5)  
  real(r8), public, parameter :: exo_o2vmr   = exo_o2bar   / (exo_pstd/1.0e5)  
  real(r8), public, parameter :: exo_h2vmr   = exo_h2bar   / (exo_pstd/1.0e5)
  real(r8), public, parameter :: exo_co2vmr  = exo_co2bar  / (exo_pstd/1.0e5)
  real(r8), public, parameter :: exo_ch4vmr  = exo_ch4bar  / (exo_pstd/1.0e5)
  real(r8), public, parameter :: exo_c2h6vmr = exo_c2h6bar / (exo_pstd/1.0e5)

  real(r8), public, parameter :: &   ! molecular weight of dry air
            exo_mwdair = exo_n2vmr*mwn2 + exo_h2vmr*mwh2 + exo_co2vmr*mwco2 + exo_ch4vmr*mwch4 + exo_c2h6vmr*mwc2h6 + exo_o2vmr*mwo2

  !! dry mass mixing ratios
  real(r8), public, parameter :: exo_n2mmr   = exo_n2vmr   * mwn2/exo_mwdair
  real(r8), public, parameter :: exo_o2mmr   = exo_o2vmr   * mwo2/exo_mwdair
  real(r8), public, parameter :: exo_h2mmr   = exo_h2vmr   * mwh2/exo_mwdair
  real(r8), public, parameter :: exo_co2mmr  = exo_co2vmr  * mwco2/exo_mwdair
  real(r8), public, parameter :: exo_ch4mmr  = exo_ch4vmr  * mwch4/exo_mwdair
  real(r8), public, parameter :: exo_c2h6mmr = exo_c2h6vmr * mwc2h6/exo_mwdair

  real(r8), public, parameter :: &   ! specific heat of dry, air J/kg/K
            exo_cpdair = exo_n2mmr*cpn2 + exo_h2mmr*cph2 + exo_co2mmr*cpco2 + exo_ch4mmr*cpch4 + exo_c2h6mmr*cpc2h6 + exo_o2mmr*cpo2

end module
