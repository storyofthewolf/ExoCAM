

module kabs
! version.archean

  implicit none
  public

  ! directory
  character(len=256), parameter :: dirk = 'data/kdist/n28archean/'

  ! K coefficient file names
  character(len=256), parameter :: k01_file = 'n28_bin01_highco2.nc'
  character(len=256), parameter :: k02_file = 'n28_bin02_highco2.nc'
  character(len=256), parameter :: k03_file = 'n28_bin03_highco2.nc'
  character(len=256), parameter :: k04_file = 'n28_bin04_highco2.nc'
  character(len=256), parameter :: k05_file = 'n28_bin05_highco2.nc'
  character(len=256), parameter :: k06_file = 'n28_bin06_highco2.nc'
  character(len=256), parameter :: k07_file = 'n28_bin07_highco2.nc'
  character(len=256), parameter :: k08_file = 'n28_bin08_highco2.nc'
  character(len=256), parameter :: k09_file = 'n28_bin09_highco2.nc'
  character(len=256), parameter :: k10_file = 'n28_bin10_highco2.nc'
  character(len=256), parameter :: k11_file = 'n28_bin11_highco2.nc'
  character(len=256), parameter :: k12_file = 'n28_bin12_highco2.nc'
  character(len=256), parameter :: k13_file = 'n28_bin13_highco2.nc'
  character(len=256), parameter :: k14_file = 'n28_bin14_highco2.nc'
  character(len=256), parameter :: k15_file = 'n28_bin15_highco2.nc'
  character(len=256), parameter :: k16_file = 'n28_bin16_highco2.nc'
  character(len=256), parameter :: k17_file = 'n28_bin17_highco2.nc'
  character(len=256), parameter :: k18_file = 'n28_bin18_highco2.nc'
  character(len=256), parameter :: k19_file = 'n28_bin19_highco2.nc'
  character(len=256), parameter :: k20_file = 'n28_bin20_highco2.nc'
  character(len=256), parameter :: k21_file = 'n28_bin21_highco2.nc'
  character(len=256), parameter :: k22_file = 'n28_bin22_highco2.nc'
  character(len=256), parameter :: k23_file = 'n28_bin23_highco2.nc'
  character(len=256), parameter :: k24_file = 'n28_bin24_highco2.nc'
  character(len=256), parameter :: k25_file = 'n28_bin25_highco2.nc'

  ! K coefficients for continuum files
  ! directory
  character(len=256), parameter :: dirct = 'data/continuum/'
  character(len=256), parameter :: kh2oself_file = 'KH2OSELF_MTCKD2.5_28bin.nc'
  character(len=256), parameter :: kco2cont_file = 'KCO2CONT_n28.nc'

  character(len=256), parameter :: dirci = 'data/cia/'
  character(len=256), parameter :: kn2n2cia_file = 'N2-N2_cia_28bin.nc'
  character(len=256), parameter :: kh2n2cia_file = 'N2-H2_cia_28bin.nc'
  character(len=256), parameter :: kh2h2cia_file = 'H2-H2_cia_28bin.nc'

end module kabs
