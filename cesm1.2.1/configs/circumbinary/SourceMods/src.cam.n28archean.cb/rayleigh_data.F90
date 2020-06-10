
module rayleigh_data

  use shr_kind_mod,      only: r8 => shr_kind_r8
  !
  ! Rayleigh scattering parameters (Vardavas and Carter, 1984)
  !
  real(r8), parameter :: delN2 = 0.0305     ! depolarization factor, N2
  real(r8), parameter :: delO2 = 0.054      ! depolarization factor, O2
  real(r8), parameter :: delCO2 = 0.0805    ! depolarization factor, CO2
  real(r8), parameter :: delH2O = 0.17      ! depolarization factor, H2O (Marshall & Smith, 1990)


  ! From Allen (1976) Astrophysical Quantities pg. 92
  real(r8), parameter :: raylA_N2 = 29.06
  real(r8), parameter :: raylB_N2 = 7.70
  real(r8), parameter :: raylA_O2 = 26.63
  real(r8), parameter :: raylB_O2 = 5.07
  real(r8), parameter :: raylA_CO2 = 43.90
  real(r8), parameter :: raylB_CO2 = 6.40

end module rayleigh_data

