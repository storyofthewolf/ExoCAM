	
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

  ! RUN OPTIONS 
  logical, public, parameter :: do_exo_atmconst = .true.        !! Read gas constituents of atmosphere from this file 
                                                                !! Overrides ghg namelist options
  logical, public, parameter :: do_exo_rt = .true.              !! .true. = use correlated-k "exoplanet" RT
  	   	   	     		    		        !! .false. = use CAM4 RT, refer elsewhere for operating
                                                                !! (currently only works for true.  delete option?)
  logical, public, parameter :: do_exo_synchronous = .false.     !! Toggle to synchronous rotation mode
                                                                !! Eventually replace with better orbital computer
  logical, public, parameter :: do_exo_rt_clearsky = .false.    !! Do parallel clearsky radiative calculation for exo rt
                                                                !! Slow, use sparingly
  logical, public, parameter :: do_exo_rt_spectral = .false.    !! collect and output spectrally resolved radiative fluxes
                                                                !! 
  integer, public, parameter :: exo_rad_step = 2                !! freq. of radiation calc in time steps (positive)
                                                                !! or hours (negative).

  !! planet parameters
  real(r8), public, parameter :: exo_planet_radius   = 1.5*6.37122e6_R8  ! radius ~ m                                                      
  real(r8), public, parameter :: exo_surface_gravity = 12.0_R8         ! gravity ~ m/s^2            
  real(r8), public, parameter :: exo_ndays  = 1.0_R8                    !! length of day scaled to number of Earth days.
  !! Earth
  !!real(r8), public, parameter :: exo_planet_radius   = 6.37122e6_R8     ! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity = 9.80616_R8       ! gravity ~ m/s^2
  !! Trappist-1e
  !!real(r8), public, parameter :: exo_planet_radius    = 5.84878e6_R8    ! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity  = 7.22925_R8      ! gravity ~ m/s^2
  !! Trappist-1f
  !!real(r8), public, parameter :: exo_planet_radius    = 6.65792e6_R8    ! radius ~ m
  !!real(r8), public, parameter :: exo_surface_gravity  = 6.11876_R8      ! gravity ~ m/s^2
  !! LHS 1140b  ;Kane DR2 values
  !!real(r8), public, parameter :: exo_planet_radius   = 1.72*6.37122e6_R8  ! radius ~ m                                                      
  !!real(r8), public, parameter :: exo_surface_gravity = 21.9734_R8         ! gravity ~ m/s^2            
  !!real(r8), public, parameter :: exo_ndays  = 24.737_r8          !! length of day scaled to number of Earth days.

  ! ORBITAL OPTIONS
  !!real(r8), public, parameter :: exo_sday = 86400.0_r8 * exo_ndays   !! sidereal period [sec], uncomment for scaling
  real(r8), public, parameter :: exo_sday = 86164.0_r8                 !! uncomment for Earth value = 86164.0
                                                                       !! becomes orbital-rotational period if (do_exo_sycnhronous)

  !! if (do_exo_synchronous) and user_nl_cpl::orb_iyear = -1
  real(r8), public, parameter :: exo_eccen = 0.0_r8   ! eccentricity
  real(r8), public, parameter :: exo_obliq = 0.0_r8   ! obliquity
  real(r8), public, parameter :: exo_mvelp = 0.0_r8   ! vernal equinox
    
  ! SOLAR CONSTANT
  real(r8), public, parameter :: exo_scon = 1360._r8         ! Solar constant (W m-2)

  ! SOLAR SPECTRAL FILE
  !! Make sure solar file matches RT configuration !!
  character(len=256), public, parameter :: exo_solar_file = '/projects/wolfet/models/ExoRT/data/solar/G2V_SUN_n28.nc '


  ! ATMOSPHERIC CONSTITUENT PARAMETERS
  ! Activated only if (do_exo_atmconst) 
  !! Must create matching initial conditions file (ncdata) !!
  real(r8), public, parameter :: exo_n2bar = 2.04700_r8                ! N2 inventory (bar)
  real(r8), public, parameter :: exo_h2bar = 0.23_r8                ! H2 inventory (bar)
  real(r8), public, parameter :: exo_co2bar = 0.023_r8               ! CO2 inventory (bar)
  real(r8), public, parameter :: exo_ch4bar = 0.0_r8               ! CH4 inventory (bar)
  real(r8), public, parameter :: exo_pstd = (exo_n2bar + exo_h2bar + exo_co2bar + exo_ch4bar)*1.0e5  ! total pressure (Pascals)

  ! OCEAN CONSTANTS
  real(r8), public, parameter :: exo_albdif = 0.06 ! 60 deg reference albedo, diffuse (default = 0.06)
  real(r8), public, parameter :: exo_albdir = 0.07 ! 60 deg reference albedo, direct  (default = 0.07)


  !! FUNDAMENTAL CONSTANTS NEEDED USED BELOW
  !! note there is some duplication with physconst.F90, keep private
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

  ! dry mass mixing ratios
  real(r8), public, parameter :: exo_n2mmr = exo_n2vmr * mwn2/exo_mwdair
  real(r8), public, parameter :: exo_h2mmr = exo_h2vmr * mwh2/exo_mwdair
  real(r8), public, parameter :: exo_co2mmr = exo_co2vmr * mwco2/exo_mwdair
  real(r8), public, parameter :: exo_ch4mmr = exo_ch4vmr * mwch4/exo_mwdair

  real(r8), public, parameter :: &   ! specific heat of dry, air J/kg/K
            exo_cpdair = exo_n2mmr*cpn2 + exo_h2mmr*cph2 + exo_co2mmr*cpco2 + exo_ch4mmr*cpch4

end module
