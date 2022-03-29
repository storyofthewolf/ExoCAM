module shr_condense_mod
! move this to src.cam/exo_condense_mod.F90

! module containing shared utilities for species condensation
! AUTHOR: Wolf, E.T. 11/18/2020

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: SHR_CONST_CO2_PTRIP


  implicit none
  private
  save

  integer, parameter :: igCO2 = 1

  public :: get_condense_temp




  !============================================================================  
  contains
  !============================================================================  
  subroutine get_condense_temp(gas_index, pressure, tcond)
    ! contains gas condense curves for gases other than CO2
    ! used for condensing clouds, and condensing onto the surface

    integer, intent(in) :: gas_index   ! which gas are we condensing
    real(r8), intent(in) :: pressure       ! pressure of gridcell in pascals
    real(r8), intent(out) :: tcond         ! condensation temperature, output

    ! do CO2 condensation
    if (gas_index .eq. igCO2) then
       !! Compute CO2 frost point. Taken from the Mars Book (Kieffer et al 1992)
       !! in chap 27 page 959 (chapter by James et al) :                 
       tcond = 3182.48 / (23.3494 -log(0.01_r8*pressure))
       
       !! Compute CO2 Frost point based on Fanale et al. (1982)
       !if (pressure .lt. SHR_CONST_CO2_PTRIP) then
       !   tcond = (-3167.8)/(alog(.01*pressure)-23.23) 
       !else 
       !   tcond = 684.2-92.3*alog(pressure)+4.32*alog(pressure)^2 
       !endif
    endif 
  
    ! Add other gas species condensation curves here as needed

  end subroutine get_condense_temp



end module shr_condense_mod
