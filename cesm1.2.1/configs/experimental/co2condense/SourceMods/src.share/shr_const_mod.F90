!===============================================================================
! SVN $Id: shr_const_mod.F90 6749 2007-10-04 20:58:20Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n00_share3_130528/shr/shr_const_mod.F90 $
!===============================================================================
! NOTES:  Wolf, E.T.  The user should not have to touch this file. 
! shr_const_mod.F90 provides a conduit for sharing model paramter type constants throughout the code
! exoplanet_mod.F90 is where the user should collect parameters for which they want to modify for planets
! 

MODULE shr_const_mod

   use shr_kind_mod
   use exoplanet_mod

   integer(SHR_KIND_IN),parameter,private :: R8 = SHR_KIND_R8 ! rename for local readability only

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------
   public

   ! constant set by exoplanet_mod.F90 module parameters
   real(R8),parameter :: SHR_CONST_SCON    = exo_scon      ! solar constant ! W m-2
   real(R8),parameter :: SHR_CONST_MWDAIR  = exo_mwdair    ! molecular weight dry air ~ kg/kmole
   real(R8),parameter :: SHR_CONST_CPDAIR  = exo_cpdair    ! specific heat of dry air   ~ J/kg/K
   real(R8),parameter :: SHR_CONST_PSTD = exo_pstd         ! standard pressure
   real(R8),parameter :: SHR_CONST_SDAY = exo_sday         ! siderial day ~ sec
   real(R8),parameter :: SHR_CONST_REARTH = exo_planet_radius  ! planet radius
   real(R8),parameter :: SHR_CONST_G = exo_surface_gravity     ! surface gravity


   ! Orbital 
   real(R8),parameter :: SHR_CONST_SECDAYD = 88800.0_R8 ! = SHR_CONST_CDAY, exo_diurnal
   real(R8),parameter :: SHR_CONST_TPERI   = 43055876.0 ! Time (in seconds) of the orbit where periapse occurs. Measures Zero at Ls=0 
   real(R8),parameter :: SHR_CONST_LONVE   = 1.907475_SHR_KIND_R8         ! longitude of vernal equinox in radians, meas. from periapse.
   real(R8),parameter :: SHR_CONST_CDVE    = 1.0  ! 80.5 ! calendar day of vernal equinox Jan1=1                      
   real(R8),parameter :: SHR_CONST_CDAY_EARTH    = 86400.0_R8      ! sec in calendar day ~ sec, [DO NOT CHANGE] !where else in the code is SHR_CONST_CDAY used
   real(R8),parameter :: SHR_CONST_PI      = 3.14159265358979323846_R8  ! pi
   real(R8),parameter :: SHR_CONST_SECYEAR = exo_porb*SHR_CONST_CDAY_EARTH  ! seconds in one year

   real(R8),parameter :: SHR_CONST_CDAY    = exo_diurnal      ! sec in calendar day ~ sec
   real(R8),parameter :: SHR_CONST_OMEGA   = 2.0_R8*SHR_CONST_PI/SHR_CONST_SDAY ! earth rot ~ rad/sec
   real(R8),parameter :: SHR_CONST_DLAPSE = SHR_CONST_G/SHR_CONST_CPDAIR 

   ! Land constants set from exoplanet_mod.F90
   real(R8),parameter :: SHR_CONST_SOILW   = exo_soilw      ! Volumetric soil water content.
   real(R8),parameter :: SHR_CONST_SOILWN  = exo_soilwn     ! Volumetric soil water content.
   real(R8),parameter :: SHR_CONST_SOILWS  = exo_soilws     ! Volumetric soil water content.
   real(R8),parameter :: SHR_CONST_SOILD   = exo_soild      ! avg density of soil (kg/m3)
   real(R8),parameter :: SHR_CONST_SOILBD  = exo_soilbd     ! avg bulk density of soil (kg/m3)
   real(R8),parameter :: SHR_CONST_SOILCP  = exo_soilcp     ! soil heat capacity ~J/kg/K 
   real(R8),parameter :: SHR_CONST_SOILPOR = exo_soilpor    !(SHR_CONST_SOILD/SHR_CONST_SOILBD)
   real(R8),parameter :: SHR_CONST_SOILHZ  = exo_soilhz     ! soil Hertz factor (K.Kossacki 2000)
   real(R8),parameter :: SHR_CONST_SOILTK  = exo_soiltk     ! soil thermal conductivity for co2 calc
   real(R8),parameter :: SHR_CONST_THINFAC = exo_thinfac    ! Factor to change the thermal inertia to match obs.                     

   ! Constants
   real(R8),parameter :: SHR_CONST_STEBOL    = 5.67e-8_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(R8),parameter :: SHR_CONST_BOLTZ     = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
   real(R8),parameter :: SHR_CONST_AVOGAD    = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
   real(R8),parameter :: SHR_CONST_RGAS      = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ       ! Universal gas constant ~ J/K/kmole
   real(R8),parameter :: SHR_CONST_MWWV      = 18.016_R8       ! molecular weight water vapor
   real(R8),parameter :: SHR_CONST_RDAIR     = SHR_CONST_RGAS/SHR_CONST_MWDAIR        ! Dry air gas constant     ~ J/K/kg
   real(R8),parameter :: SHR_CONST_RWV       = SHR_CONST_RGAS/SHR_CONST_MWWV          ! Water vapor gas constant ~ J/K/kg
   real(R8),parameter :: SHR_CONST_ZVIR      = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_R8 ! RWV/RDAIR - 1.0
   real(R8),parameter :: SHR_CONST_LOSCHMIDT = 2.6867774e25          ! Loschmidt number ~# m^-3
   real(R8),parameter :: SHR_CONST_KARMAN    = 0.4_R8          ! Von Karman constant
   real(R8),parameter :: SHR_CONST_PDB       = 0.0112372_R8    ! ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)
   real(SHR_KIND_R8),parameter :: SHR_CONST_EPSQS  = SHR_CONST_MWWV/SHR_CONST_MWDAIR  ! molecular weight ratio of wv/dair
 
   real(R8),parameter :: SHR_CONST_TKTRIP  = 273.16_R8       ! triple point of fresh water        ~ K
   real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K 
   real(R8),parameter :: SHR_CONST_TKFRZSW = SHR_CONST_TKFRZ - 1.8_R8 ! freezing T of salt water  ~ K

   real(R8),parameter :: SHR_CONST_RHODAIR = &               ! density of dry air at STP  ~ kg/m^3
                         SHR_CONST_PSTD/(SHR_CONST_RDAIR*SHR_CONST_TKFRZ)
   real(R8),parameter :: SHR_CONST_RHOFW   = 1.000e3_R8      ! density of fresh water     ~ kg/m^3
   real(R8),parameter :: SHR_CONST_RHOSW   = 1.026e3_R8      ! density of sea water       ~ kg/m^3
   real(R8),parameter :: SHR_CONST_RHOICE  = 0.917e3_R8      ! density of ice             ~ kg/m^3
   real(R8),parameter :: SHR_CONST_CPWV    = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
   real(R8),parameter :: SHR_CONST_CPVIR   = (SHR_CONST_CPWV/SHR_CONST_CPDAIR)-1.0_R8 ! CPWV/CPDAIR - 1.0
   real(R8),parameter :: SHR_CONST_CPFW    = 4.188e3_R8      ! specific heat of fresh h2o ~ J/kg/K
   real(R8),parameter :: SHR_CONST_CPSW    = 3.996e3_R8      ! specific heat of sea h2o   ~ J/kg/K
   real(R8),parameter :: SHR_CONST_CPICE   = 2.11727e3_R8    ! specific heat of fresh ice ~ J/kg/K
   real(R8),parameter :: SHR_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
   real(R8),parameter :: SHR_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
   real(R8),parameter :: SHR_CONST_LATSUB  = &               ! latent heat of sublimation ~ J/kg
                         SHR_CONST_LATICE + SHR_CONST_LATVAP
   real(R8),parameter :: SHR_CONST_OCN_REF_SAL = 34.7_R8     ! ocn ref salinity (psu)
   real(R8),parameter :: SHR_CONST_ICE_REF_SAL =  4.0_R8     ! ice ref salinity (psu)

   real(R8),parameter :: SHR_CONST_SPVAL   = 1.0e30_R8       ! special missing value


   ! CO2 Ice Constants for exoplanet_mod.F90
   real(r8),parameter :: SHR_CONST_DENSITYCO2FR  = exo_densityco2fr    ! Density of CO2 frost (Kg/m3) (D.Smith 2001)
   real(r8),parameter :: SHR_CONST_CO2LATSUB     = exo_co2latsub       ! latent heat of sublimation ~ J/kg
   real(r8),parameter :: SHR_CONST_CO2_ALVIS     = exo_co2_alvis       ! Visible albedo of CO2 snow (Hourdin 1993)
   real(r8),parameter :: SHR_CONST_CO2_ALNIR     = exo_co2_alnir       ! Near-IR albedo of CO2 snow
   real(r8),parameter :: SHR_CONST_CO2_ALS       = exo_co2_als         ! albedo of CO2 snow South (Hourdin 1993)
   real(r8),parameter :: SHR_CONST_CO2_ALN       = exo_co2_aln         ! albedo of CO2 snow North 
   real(r8),parameter :: SHR_CONST_CO2_EMIS      = exo_co2_emis        ! Emissivity of CO2 ice (Hourdin 1993) 
   real(r8),parameter :: SHR_CONST_DENSITYCO2ICE = exo_densityco2ice   ! Density of CO2 ice (kg/m3)
   real(r8),parameter :: SHR_CONST_CPCO2ICE      = exo_cpco2ice        ! specific heat of co2 ice [J kg-1 K-1]   
   real(r8),parameter :: SHR_CONST_CO2_TI        = exo_co2ice_ti       ! thermal inertia of co2 ice
   real(r8),parameter :: SHR_CONST_CO2ICE_THICK  = exo_co2ice_thick    ! threshold for CO2 ice thickness to cover underlying layers [m]
 
   ! Mars Atmosphere Constants
   real(r8),parameter :: SHR_CONST_CO2   = exo_co2vmr         ! co2 fraction (volume)
 !real(SHR_KIND_R8),parameter :: SHR_CONST_N2   = .027_SHR_KIND_R8         ! nitrogen fraction (volume)
 !  real(SHR_KIND_R8),parameter :: SHR_CONST_AR   = .016_SHR_KIND_R8         ! argon fraction (volume)
   real(r8),parameter :: SHR_CONST_O2   = 0.0                       ! .0013_SHR_KIND_R8         ! oxygen fraction (volume)
 
  ! mars viscosity terms
   real(r8),parameter :: SHR_CONST_DVISC     = 1.45e-5_R8       ! Dynamic viscosity of air at 293K (Pa.s, Pang 2005)
   real(r8),parameter :: SHR_CONST_KVISC     = SHR_CONST_DVISC/SHR_CONST_RHODAIR ! Kinematic viscosity of air (m^2/s)
   real(r8),parameter :: SHR_CONST_CBL       = 7000.0_R8  !(8000 for mars?  orig is 1000) Characteristic convective boundary layer height (m) 
   real(r8),parameter :: SHR_CONST_TKAIR     = exo_tkair  ! thermal conductivity of air

   real(r8),parameter :: SHR_CONST_PPATRIP   = 611.73_r8      ! triple point pressure of fresh water ~ Pa    
   real(r8),parameter :: SHR_CONST_CO2_PTRIP = 518000.0_r8  ! triple point pressure of CO2 ~Pa

END MODULE shr_const_mod
