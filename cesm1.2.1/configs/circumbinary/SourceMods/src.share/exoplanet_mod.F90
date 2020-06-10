

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
  logical, public, parameter :: do_exo_rt_clearsky = .false.    !! Do parallel clearsky radiative calculation for exo rt
                                                                !! Slow, use sparingly
  logical, public, parameter :: do_exo_rt_spectral = .false.    !! collect and output spectrally resolved radiative fluxes
                                                                !! 
  logical, public, parameter :: do_exo_circumbinary = .true.    !! Stellar flux calculation for binary star systems

  integer, public, parameter :: exo_rad_step = 2                !! freq. of radiation calc in time steps (positive)
                                                                !! or hours (negative).

  
  !! ==============  planet parameters ============== !!
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
  real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_r8 !! radius ~ m
  real(r8), public, parameter :: exo_surface_gravity = 9.80616_r8          !! gravity ~ m/s^2            
  real(r8), public, parameter :: exo_ndays  = 1.0_r8                    !! scaler to number of Earth days.
  real(r8), public, parameter :: exo_porb = 353.88_r8                        !! orbital period, for obliquity cycle, and optionally for exo_sday
  !!real(r8), public, parameter :: exo_sday = 86164.0_r8                !! sidereal period Earth
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays    !! sidereal period, for synchronous rotator
  real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays / (1._r8 + exo_ndays/exo_porb)  !! sidereal period [sec]

  !! Examples
  !! Earth
  !!real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8     !! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8       !! gravity ~ m/s^2
  !!real(r8), public, parameter :: exo_ndays  = 1.0_R8                    !! scaler to number of Earth days.
  !!real(r8), public, parameter :: exo_porb = 365.0_R8
  !!real(r8), public, parameter :: exo_sday = 86164.0_R8                  !! sidereal period [sec], for Earth value = 86164.0 

  !! Earth - slow asynchronous rotator
  !!real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8      !! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8        !! gravity ~ m/s^2
  !!real(r8), public, parameter :: exo_ndays  = 10.0_R8                    !! scaler to number of Earth days.
  !!real(r8), public, parameter :: exo_porb = 365.0_R8                     !! orbital period 
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays / (1._r8 + exo_ndays/exo_porb)  !! sidereal period [sec]

  !! Trappist-1e  (Gillon et al. 2017)  --  synchronous rotator
  !!real(r8), public, parameter :: exo_planet_radius    = 5.84878e6_R8    !! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity  = 7.22925_R8      !! gravity ~ m/s^2
  !!real(r8), public, parameter :: exo_ndays  = 6.099615_r8               !! scaler to number of Earth days.   
  !!real(r8), public, parameter :: exo_porb = exo_ndays                   !! orbital period 
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays      !! sidereal period, synchronous rotator

  !! Trappist-1f  (Gillon et al. 2017)  --  synchronous rotator
  !!real(r8), public, parameter :: exo_planet_radius    = 6.65792e6_R8    !! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity  = 6.11876_R8      !! gravity ~ m/s^2
  !!real(r8), public, parameter :: exo_ndays  = 9.206690_r8               !! scaler to number of Earth days.
  !!real(r8), public, parameter :: exo_porb = exo_ndays                   !! orbital period 
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays      !! sidereal period, synchronous rotator

  !! if set user_nl_cpl::orb_iyear = -1
  real(r8), public, parameter :: exo_eccen = 0.0_r8   ! eccentricity
  real(r8), public, parameter :: exo_obliq = 23.5_r8   ! obliquity [degrees]
  real(r8), public, parameter :: exo_mvelp = 0.0_r8   ! vernal equinox

    
  !! ============== STELLAR OPTIONS ============== !!
  !! SOLAR CONSTANT
  real(r8), public, parameter :: exo_scon = 1360._r8         ! Solar constant (W m-2)

  !! SOLAR SPECTRAL FILE
  !! Make sure solar file matches spectral intervals for selected RT configuration !!
  character(len=256), public, parameter :: exo_solar_file = '/gpfsm/dnb53/etwolf/models/ExoRT/data/solar/BT_Settl_5714K_n28.nc'

  !! if (do_exo_circumbinary) 
  !! then exo_scon does not take effect
  !! make sure this makes sense with radius and surface gravity given above.
  real(r8), public, parameter :: planet_mass = 1.0_r8        ! units of M_Earth

  !! Mass (M_sun), Teff (K), Radius (R_sun), Luminosity (L_sun)
  !! self-consistency of M, T, R, L, is the responsibility of the user
  real(r8), public, parameter, dimension(4) :: stellar_params1 = (/1.0_r8, 5714._r8, 0.9989_r8, 0.958041_r8 /)   
  real(r8), public, parameter, dimension(4) :: stellar_params2 = (/0.55_r8, 3824._r8, 0.555519_r8, 0.0577802_r8 /)   
  
  real(r8), public, parameter :: cb_separation = 0.4_r8   ! binary separation in AU
  real(r8), public, parameter :: cb_semiaxis = 1.1_r8    ! planet semimajor axis in AU, no uses
  real(r8), public, parameter :: cb_eccen = 0.0_r8        ! binary eccentricity
  real(r8), public, parameter, dimension(2) :: cb_semia = (/ cb_separation, cb_semiaxis /) ! binary and planet semimajor axis.
  real(r8), public, parameter, dimension(2) :: cb_man = (/ 0.0_r8, 0.0_r8 /)    ! initial mean longitudes   
  character(len=256), public, parameter :: exo_solar_file_2 = '/gpfsm/dnb53/etwolf/models/ExoRT/data/solar/BT_Settl_3824K_n28.nc'
  logical, public, parameter :: cb_transit = .true.  

  ! Derived Arrays for Binary Star Calculations
  ! For input into twostars module
  real(r8), public, parameter :: exo_precess = 0.0_r8   !DEFINE? MOVE!
  real(r8), public, parameter, dimension(1:3) :: cb_mass   = (/ stellar_params1(1), stellar_params2(1), planet_mass  /)
  real(r8), public, parameter, dimension(1:2) :: cb_teff   = (/ stellar_params1(2), stellar_params2(2) /)
  real(r8), public, parameter, dimension(1:2) :: cb_radius = (/ stellar_params1(3), stellar_params2(3) /)
  real(r8), public, parameter, dimension(1:2) :: cb_lum    = (/ stellar_params1(4), stellar_params2(4) /)
  real(r8), public, parameter, dimension(4)   :: cb_spin   = (/ exo_obliq, exo_precess, exo_ndays, 180._r8 /) 
                                               ! obliq, precession, diurnal period, initial hour angle




  !! ============== ATMOSPHERIC CONSTITUENT PARAMETERS ============== !!
  !! Activated only if (do_exo_atmconst = .true.) 
  !! Must create matching initial conditions file (ncdata) !!
  real(r8), public, parameter :: exo_n2bar = 0.9996_r8                ! N2 inventory (bar)
  real(r8), public, parameter :: exo_h2bar = 0.0_r8                ! H2 inventory (bar)
  real(r8), public, parameter :: exo_co2bar = 0.0004_r8               ! CO2 inventory (bar)
  real(r8), public, parameter :: exo_ch4bar = 0.0_r8               ! CH4 inventory (bar)
  real(r8), public, parameter :: exo_pstd = (exo_n2bar + exo_h2bar + exo_co2bar + exo_ch4bar)*1.0e5  ! total pressure (Pascals)

  !! ============== OCEAN ALBEDO CONSTANTS ============== !!
  real(r8), public, parameter :: exo_albdif = 0.06 ! 60 deg reference albedo, diffuse (default = 0.06)
  real(r8), public, parameter :: exo_albdir = 0.06 ! 60 deg reference albedo, direct  (default = 0.07)


  !! ============== FUNDAMENTAL CONSTANTS NEEDED USED BELOW ============== !!
  !! note there is some duplication with physconst.F90, keep private
  !! No need to modifications need below this point
  real(r8), parameter :: mwn2 = 28._r8
  real(r8), parameter :: mwh2 = 2._r8
  real(r8), parameter :: mwco2 = 44._r8
  real(r8), parameter :: mwch4 = 16._r8
  real(r8), parameter :: cpn2 = 1.039e3_r8
  real(r8), parameter :: cph2 = 14.32e3_r8
  real(r8), parameter :: cpco2 = 0.846e3_r8
  real(r8), parameter :: cpch4 = 2.226e3

  !! DERIVED CONSTANTS -- DO NOT MODIFY
  !! automatically calculated from above inputs in bar
  ! dry volume mixing ratios  kg/kmole
  real(r8), public, parameter :: exo_n2vmr = exo_n2bar / (exo_pstd/1.0e5)  
  real(r8), public, parameter :: exo_h2vmr = exo_h2bar / (exo_pstd/1.0e5)
  real(r8), public, parameter :: exo_co2vmr = exo_co2bar / (exo_pstd/1.0e5)
  real(r8), public, parameter :: exo_ch4vmr = exo_ch4bar / (exo_pstd/1.0e5)

  real(r8), public, parameter :: &   ! molecular weight of dry air
            exo_mwdair = exo_n2vmr*mwn2 + exo_h2vmr*mwh2 + exo_co2vmr*mwco2 + exo_ch4vmr*mwch4

  !! dry mass mixing ratios
  real(r8), public, parameter :: exo_n2mmr = exo_n2vmr * mwn2/exo_mwdair
  real(r8), public, parameter :: exo_h2mmr = exo_h2vmr * mwh2/exo_mwdair
  real(r8), public, parameter :: exo_co2mmr = exo_co2vmr * mwco2/exo_mwdair
  real(r8), public, parameter :: exo_ch4mmr = exo_ch4vmr * mwch4/exo_mwdair

  real(r8), public, parameter :: &   ! specific heat of dry, air J/kg/K
            exo_cpdair = exo_n2mmr*cpn2 + exo_h2mmr*cph2 + exo_co2mmr*cpco2 + exo_ch4mmr*cpch4

end module
