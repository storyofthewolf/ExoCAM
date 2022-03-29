

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
  logical, public, parameter :: do_exo_rt_clearsky = .true.     !! Do parallel clearsky radiative calculation for exo rt
                                                                !! Slow, use sparingly
  logical, public, parameter :: do_exo_rt_spectral = .false.    !! collect and output spectrally resolved radiative fluxes
                                                                !! 
  integer, public, parameter :: exo_rad_step = 4 ! 4                !! freq. of radiation calc in time steps (positive)
                                                                !! or hours (negative).
  logical, public, parameter :: do_carma_exort = .false.        !! Set to true only if running with the CARMA microphyiscs package
                                                                !! and linking aerosol absorption to ExoRT  
  logical, public, parameter :: do_exo_gw = .false.             !! flag to turn on gravity waves.  Note, present gw wave parameterization
                                                                !! does not work for low pressure atmospheres.
  logical, public, parameter :: do_exo_condense_co2 = .true.    !! flag to do CO2 clouds, must add an advected constituent at build 
                                                                !! -nadv 4 (Q, CLDLIQ, CLDICE, CLDCO2)

  
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
  !!real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8 !! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8          !! gravity ~ m/s^2            
  !!real(r8), public, parameter :: exo_ndays  = 1.0_R8                    !! scaler to number of Earth days.
  !!real(r8), public, parameter :: exo_porb = 365.                        !! orbital period, for obliquity cycle, and optionally for exo_sday
  !!real(r8), public, parameter :: exo_sday = 86164.0_r8                !! sidereal period Earth
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays    !! sidereal period, for synchronous rotator
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays / (1._r8 + exo_ndays/exo_porb)  !! sidereal period [sec]

  !! Mars Orbital
  real(r8), public, parameter :: exo_planet_radius   = 3.38992e6_R8     !! radius ~ m 
  real(r8), public, parameter :: exo_surface_gravity = 3.711_R8         !! gravity ~ m/s^2
  real(r8), public, parameter :: exo_diurnal = 88800.0_R8               !! Length of diurnal period, solar-day ~ s
  real(r8), public, parameter :: exo_ndays  = exo_diurnal/86400.0_R8    !! scaler to number of Earth days
  real(r8), public, parameter :: exo_porb = 687.0_R8                    !! orbital period in Earth Days, for obliquity/eccentricity cycles
  real(r8), public, parameter :: exo_sday = 88642.663_R8                !! sidereal period 

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
  !! Modern Mars
  !real(r8), public, parameter :: exo_eccen = 0.0934_r8  ! eccentricity  
  !real(r8), public, parameter :: exo_obliq = 25.19_r8   ! obliquity [degrees]
  !real(r8), public, parameter :: exo_mvelp = 70.90_r8   ! vernal equinox
  
  !! 
  real(r8), public, parameter :: exo_eccen = 0.0_r8  ! eccentricity  
  real(r8), public, parameter :: exo_obliq = 25.0_r8   ! obliquity [degrees]
  real(r8), public, parameter :: exo_mvelp = 70.90_r8   ! vernal equinox

    
  !! ============== STELLAR OPTIONS ============== !!
  !! SOLAR CONSTANT
!!  real(r8), public, parameter :: exo_scon = 586.2_r8         ! Solar constant (W m-2)  ! Modern Mars
  real(r8), public, parameter :: exo_scon = 441.1_r8         ! Solar constant (W m-2)  ! Ancient Mars (75%)

  ! SOLAR SPECTRAL FILE
  !! Make sure solar file matches spectral intervals for selected RT configuration !!
!  character(len=256), public, parameter :: exo_solar_file = '/projects/wolfet/models/ExoRT/data/solar/G2V_SUN_n68.nc'
  character(len=256), public, parameter :: exo_solar_file = '/gpfsm/dnb53/etwolf/models/ExoRT/data/solar/G2V_SUN_n68.nc'


  !! ============== ATMOSPHERIC CONSTITUENT PARAMETERS ============== !!
  !! Activated only if (do_exo_atmconst = .true.) 
  !! Must create matching initial conditions file (ncdata) !!
  real(r8), public, parameter :: exo_n2bar = 0.0_r8                ! N2 inventory (bar)
  real(r8), public, parameter :: exo_h2bar = 0.0_r8                ! H2 inventory (bar)
  real(r8), public, parameter :: exo_co2bar = 1.0_r8               ! CO2 inventory (bar)
  !!real(r8), public, parameter :: exo_co2bar = 0.0076_r8               ! CO2 inventory (bar)
  real(r8), public, parameter :: exo_ch4bar = 0.0_r8               ! CH4 inventory (bar)
  real(r8), public, parameter :: exo_pstd = (exo_n2bar + exo_h2bar + exo_co2bar + exo_ch4bar)*1.0e5  ! total pressure (Pascals)

  !! ============== OCEAN ALBEDO CONSTANTS ============== !!
  real(r8), public, parameter :: exo_albdif = 0.06 ! 60 deg reference albedo, diffuse (default = 0.06)
  real(r8), public, parameter :: exo_albdir = 0.06 ! 60 deg reference albedo, direct  (default = 0.07)


  !! ============== Land Surface Constants ============== !!
  real(r8), public, parameter :: exo_soilw   = 0.1_r8       ! Volumetric soil water content.
  real(r8), public, parameter :: exo_soilwn  = 0.0_r8      ! Volumetric soil water content.
  real(r8), public, parameter :: exo_soilws  = 0.0_r8      ! Volumetric soil water content.
  real(r8), public, parameter :: exo_soild   = 1400.0_r8    ! avg density of Martian soil (kg/m3)
  real(r8), public, parameter :: exo_soilbd  = 1400.0_r8    ! avg bulk density of Martian soil (3060 is K.Kossacki 2000) (kg/m3)
  real(r8), public, parameter :: exo_soilcp  = 627.9_r8     ! Martian soil heat capacity ~J/kg/K (Mellon 2000 = 627.9) ! (Zent 1993 = 820)
  real(r8), public, parameter :: exo_soilpor = 0.3_r8       !(SHR_CONST_SOILD/SHR_CONST_SOILBD)           ! 1 - Martian soil porosity
  real(r8), public, parameter :: exo_soilhz  = 1.e-3_r8     ! Martian soil Hertz factor (K.Kossacki 2000)
  real(r8), public, parameter :: exo_soiltk  = 0.04_r8      ! Typical Martian soil thermal conductivity for co2 calc
  real(r8), public, parameter :: exo_thinfac = 1.0_r8       ! Factor to change the thermal inertia to match obs.

  real(r8), public, parameter :: exo_soil_roughness = 0.01_r8 ! m, constant soil roughness
  real(r8), public, parameter :: exo_soil_thermal_intertia = 250.0_r8 ! J m-2 s-1/2 K-1, soil thermal inertia

  !!================ other mars things ================== !!
  real(r8), public, parameter :: exo_tkair = 0.015_r8   ! Thermal conductivity of air at STP ~ W/m/K (CRC handbook. hbcpnetbase.com) 

  !! ============== CO2 Ice constants  ============== !!
  real(r8), public, parameter :: exo_densityco2fr  =  910.0_r8   ! Density of CO2 frost (Kg/m3) (D.Smith 2001)
!  real(r8), public, parameter :: exo_co2latsub     =  7.76e5_r8      ! latent heat of sublimation ~ J/kg   ! Richard's model
  real(r8), public, parameter :: exo_co2latsub     =  5.9e5_r8      ! latent heat of sublimation ~ J/kg    ! Used everywhere else
  real(r8), public, parameter :: exo_co2_alvis     =  0.50_r8  !Forget_2013  ! 0.60_r8    ! Visible albedo of CO2 snow (Hourdin 1993)
  real(r8), public, parameter :: exo_co2_alnir     =  0.50_r8  !Forget_2013  !0.60_r8    ! Near-IR albedo of CO2 snow
  real(r8), public, parameter :: exo_co2_als       =  0.80_r8    ! albedo of CO2 snow South (Hourdin 1993)
  real(r8), public, parameter :: exo_co2_aln       =  0.60_r8    ! albedo of CO2 snow North
  real(r8), public, parameter :: exo_co2_tk        =  1.5_r8     ! Thermal conductivity of CO2 snow (W/m/K) (J.Larsen 2000)
  real(r8), public, parameter :: exo_co2_emis      =  0.85_r8    !Forget_2013 !0.8_r8             ! Emissivity of CO2 ice (Hourdin 1993)
  real(r8), public, parameter :: exo_densityco2ice =  1620._r8    ! Density of CO2 ice (kg/m3) (this is what LMD uses for everything)
  real(r8), public, parameter :: exo_cpco2ice      =  1000._r8   ! specific heat of co2 ice [J kg-1 K-1]


  !! ============== FUNDAMENTAL CONSTANTS NEEDED USED BELOW ============== !!
  !! note there is some duplication with physconst.F90, keep private
  
  ! molecular weights
  real(r8), parameter :: mwn2 = 28._r8
  real(r8), parameter :: mwh2 = 2._r8
  real(r8), parameter :: mwco2 = 44._r8
  real(r8), parameter :: mwch4 = 16._r8

  ! specific heat of dry air
  real(r8), parameter :: cpn2 = 1.039e3_r8
  real(r8), parameter :: cph2 = 14.32e3_r8
  real(r8), parameter :: cpch4 = 2.226e3

  ! CO2 values change non-negligible with temperature
  ! may need to modify this based on atmosphere.
  real(r8), parameter :: cpco2 = 0.735e3_r8    ! 200 K    ! mars models use this
  !real(r8), parameter :: cpco2 = 0.763e3_r8    ! 225 K  
  !real(r8), parameter :: cpco2 = 0.791e3_r8    ! 250 K  
  !real(r8), parameter :: cpco2 = 0.819e3_r8    ! 275 K  
  !real(r8), parameter :: cpco2 = 0.846e3_r8    ! 300 K  


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
