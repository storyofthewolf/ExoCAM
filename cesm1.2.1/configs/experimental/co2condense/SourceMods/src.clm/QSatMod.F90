
module QSatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: QSatMod
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: QSat
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat
!
! !INTERFACE:
  subroutine QSat (T, p, es, esdT, qs, qsdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ, &
                             SHR_CONST_RGAS, SHR_CONST_LATSUB, & ! mars
                             SHR_CONST_EPSQS, SHR_CONST_PPATRIP

!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es       ! vapor pressure (pa)
    real(r8), intent(out) :: esdT     ! d(es)/d(T)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2
    real(r8) :: epsqs    !ratio of h2o to dry air molecular weights ! mars
    real(r8) :: omeps    ! mars


!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------


    epsqs = shr_const_epsqs ! mars
    omeps = 1. - epsqs      ! mars

    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if ((td >= 0.0_r8) .and. (p > SHR_CONST_PPATRIP)) then   ! mars
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = (exp(9.550426 - (5723.265/T) + 3.53068*log(T) - 0.00728332*T ))/100.  ! mars
       esdT = (SHR_CONST_LATSUB*es/(SHR_CONST_RGAS*T*T))/100.   ! mars
!       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &  ! mars
!            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))   ! mars
!       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &  ! mars
!            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))   ! mars
    endif

    es    = es    * 100._r8            ! pa
    esdT  = esdT  * 100._r8            ! pa/K

!    vp    = 1.0_r8   / (p - 0.378_r8*es)
!    vp1   = 0.622_r8 * vp
!    vp2   = vp1   * vp

!    qs    = es    * vp1             ! kg/kg
!    qsdT  = esdT  * vp2 * p         ! 1 / K

    vp    = 1.0   / (p - omeps*es)    ! mars 
    vp1   = epsqs * vp                ! mars
    vp2   = vp1   * vp

    qs    = min(1.,max(1.e-10,es * vp1))             ! kg/kg  nars
    qsdT  = esdT  * vp2 * p         ! 1 / K
    if (qs .eq. 1.e-10) then
      qs = 1.e-10
      qsdT = 0.
    endif



  end subroutine QSat

end module QSatMod
