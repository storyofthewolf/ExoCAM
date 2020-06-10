
module exo_init_ref

!---------------------------------------------------------------------       
! Purpose:                                                                   

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_const_mod,    only: SHR_CONST_PI
  use physconst,        only: scon
  use radgrid
  use spmd_utils,       only: masterproc
 
  implicit none
  public  

  ! Approximate smallest double precision floating point difference 
  real(r8), parameter :: SMALLd = 1.0d-12                             
  real(r8), parameter :: SMALLe = 1.0e-12
  !real(r8), parameter :: SMALLd = 1.0d-8
  !real(r8), parameter :: SMALLe = 1.0e-8

  real(r8), parameter :: sqrt3 = 1.732050808d0      ! square root of 3
  real(r8), parameter :: mb_to_atm = 9.869233e-4    ! convert pressure from Pa to atm

  !------------------------------------------------------------------------------  
  ! Radiative transfer model variable/array declarations                           
  !                                   !                                                                              
  ! Assign beginning and end wavelength range and point indices for each         
  !  wavelength group                                               

  ! Do full spectrum LW and SW calculation
  !integer, parameter  :: lw_iwbeg = 1     ! thermal band wvl integration limits                  
  !integer, parameter  :: lw_iwend = ntot_wavlnrng
  !integer, parameter  :: sw_iwbeg = 1     ! solar band wvl integration limits                    
  !integer, parameter  :: sw_iwend = ntot_wavlnrng
  !integer, parameter  :: lw_ipbeg = 1     ! thermal band gpt integration limits                  
  !integer, parameter  :: lw_ipend = ntot_gpt
  !integer, parameter  :: sw_ipbeg = 1     ! solar band gpt integration limits                    
  !integer, parameter  :: sw_ipend = ntot_gpt

  ! Limit LW and SW spectrum calculation for improved efficiency
  integer, parameter  :: lw_iwbeg = 1     ! thermal band wvl integration limits                  
  integer, parameter  :: lw_iwend = 17
  integer, parameter  :: sw_iwbeg = 6     ! solar band wvl integration limits                    
  integer, parameter  :: sw_iwend = 28
  integer, parameter  :: lw_ipbeg = 1     ! thermal band gpt integration limits                  
  integer, parameter  :: lw_ipend = 160
  integer, parameter  :: sw_ipbeg = 65     ! solar band gpt integration limits                    
  integer, parameter  :: sw_ipend = ntot_gpt

  !                                                                              
  ! set two-stream model coefficients                                            
  !                                           
  ! for solar stream (quadrature)    
  real(r8), parameter :: U1Isol  = sqrt3                    ! 2*PI / mu1 factors 
  real(r8), parameter :: U1I2sol = 0.5d0*sqrt3              ! 2*PI / mu1 factors 
  real(r8), parameter :: U1Ssol  = 2.0*SHR_CONST_PI/U1Isol  ! mu1 factors 
  ! for thermal stream (hemispsheric mean)      
  real(r8), parameter :: U1Iir   = 2.d0                     ! mu1 factors 
  real(r8), parameter :: U1I2ir  = 1.d0                     ! mu1 factors 
  real(r8), parameter :: U1Sir   = 2.0*SHR_CONST_PI/U1Iir   ! mu1 factors 


  real(r8), dimension(ntot_gpt) :: gw_solflux
  real(r8), dimension(ntot_wavlnrng) :: solflux

  ! Define constants for Gauss quadratrue calculations
  integer, parameter  :: ngangles = 3               ! # of Gauss angles to use

  ! Gauss weight vector and zenith angle
  real(r8), dimension(ntot_gpt) :: g_weight
  real(r8), dimension(ngangles) :: g_ang_weight
  real(r8), dimension(ngangles) :: g_angle

  ! Weights for Gaussian quadrature (for each probability interval) over a
  ! hemisphere subdivided into 'ngangles' equal angles  [none] (there are
  ! 'ngangles' of these):

  ! For ngangles = 3:
  data g_angle / 0.21234054, 0.59053314, 0.91141204 /
  data g_ang_weight / 0.06982698, 0.22924111, 0.20093191 /

!============================================================================
contains
!============================================================================


  subroutine init_ref
!------------------------------------------------------------------------      
! Purpose: Initial reference value arrays 
!------------------------------------------------------------------------ 
  implicit none
!------------------------------------------------------------------------ 
! Local Variables                                                         
!
integer :: iq
integer :: iw
integer :: ig
integer :: ip

!------------------------------------------------------------------------      
!                                                                              
! Start Code                                                                   
!                                                                              
    !write(*,*)   "INIT_REF: INITIALIZE GUASS POINT ARRAYS"
    ! Arrange g_weight(ntot_pt) array 
    iq = 0
    do iw=1,ntot_wavlnrng
      do ig=1, ngauss_pts(iw)
        iq = iq + 1
        if (ngauss_pts(iw) .eq. ngauss_8gpt) g_weight(iq) =  g_weight_8gpt(ig)
        if (ngauss_pts(iw) .eq. ngauss_16gpt) g_weight(iq) =  g_weight_16gpt(ig)
!!        if (iw .eq. 27) g_weight(iq) = solar_gweight_sp27(ig)   !!  decpreciated from old O2, O3 version         
!!        if (iw .eq. 28) g_weight(iq) = solar_gweight_sp28(ig)   !!  decpreciated from old O2, O3 version         
      enddo
    enddo

     ! Scale solar constant to namelist value
    solarflux(:) = solarflux(:)*scon/S0

    !
    ! Calculate the "average" wavenumber (1/wavelength) <wavenum()> for each    
    !  wavelength interval in each wavelength group [1/cm]:                     
    !                                                                           
    iq = 0
    ip = lw_ipbeg-1
    do iw=1,ntot_wavlnrng  ! "avg" wavenumber over each band                    
      do ig=1,ngauss_pts(iw)
        ip = ip+1
        ! Gauss-weighted solar flux in each probability interval:               
        gw_solflux(ip) = solarflux(iw)*g_weight(ip)
!write(*,*) "init_ref,g_weight",g_weight(ip)
      enddo
    enddo

    if (masterproc) then
      write(*,*) "INIT_REF: total solar irradiance scaled to ",scon, "W m-2"
      write(*,*) "solar flux [W m-2] in each spectral interval"
      do iw=lw_iwbeg, sw_iwend
        write(*,*) iw, solarflux(iw)
      enddo
      write(*,*) "TOTAL SOLAR FLUX:", SUM(solarflux)
    endif

    call map_co2cont_gpt

 end subroutine init_ref


!============================================================================ 

  subroutine map_co2cont_gpt

!------------------------------------------------------------------------
! PURPOSE:  CO2 Continuum data sets are on 8 point gauss interval bin.  Adjust
!           to match number of gauss intervals used for major absorbing gases 
!           Required to populated the CO2 continuum from file!!
!------------------------------------------------------------------------ 

    implicit none

!------------------------------------------------------------------------ 
!
! Local Variables
!
  integer :: iw
  integer :: ig
  integer :: itc

!------------------------------------------------------------------------ 
! Start Code
!

    ! Initialize reduced continuum k coefficient arrays 
    kco2cont(:) = 0.0
    !
    ! Gauss point adjustment for CO2 continuum 
    ! 
    itc = 0
    do iw=1, ntot_wavlnrng
      if (ngauss_pts(iw) .eq. 8) then ! no adjustment needed
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          kco2cont(itc) = kco2cont_8gpt(ig,iw)
        enddo
      endif
      if (ngauss_pts(iw) .eq. 16) then
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          kco2cont(itc) = kco2cont_8gpt(map8to16gpt(ig),iw)
        enddo
      endif
    enddo

  end subroutine map_co2cont_gpt

end module exo_init_ref
