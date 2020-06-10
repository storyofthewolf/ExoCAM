
module cloud

implicit none
public

  ! directory
  character(len=256), parameter :: dircld = 'data/cloud/'

  ! Cloud mie data
  character(len=256), parameter :: cldoptsL_file = 'cloudoptics_liquid_mie_n28.nc'
  character(len=256), parameter :: cldoptsI_file = 'cloudoptics_ice_mie_n28.nc'

end module cloud
