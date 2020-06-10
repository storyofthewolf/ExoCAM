
module rad_interp_mod

! version highco2
! Contains subroutines used for interpolating K coefficient and 
! cloud optical data in radiation.F90

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radgrid

  implicit none
  private
  save

!------------------------------------------------------------------------
!
!  Public Interfaces
!

  public :: trilinear_interpK_lower_8gpt
  public :: trilinear_interpK_upper_8gpt
  public :: trilinear_interpK_lower_16gpt
  public :: trilinear_interpK_upper_16gpt
  public :: bilinear_interpK_lower_8gpt
  public :: bilinear_interpK_upper_8gpt
  public :: bilinear_interpK_lower_16gpt
  public :: bilinear_interpK_upper_16gpt
  public :: interpKself
  public :: interpKself2
  public :: interpolate_cld
  public :: interpH2H2cia
  public :: interpH2N2cia
  public :: interpN2N2cia

!========================================================================
contains
!========================================================================

  subroutine trilinear_interpK_lower_8gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_weight, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
! for lower atmosphere pressure grid                                       
!
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(8, kc_npress_lower, kc_ntemp_lower, kc_nweight)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_weight
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- trilinear_interpK_lower ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_lower) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))           
    else
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))  
    endif
   
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_lower(p_ref_index),log10pgrid_lower(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_lower) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    else
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid_lower(t_ref_index),tgrid_lower(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    if (w_ref_index .eq. kc_nweight) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    else
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    endif

    !write(*,*) "interp_weight", weight
 
    ! perform trilinear interpolation between P,T,W

    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    !write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight

    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight

    !write(*,*) "ans", ans

    return

  end subroutine trilinear_interpK_lower_8gpt


!========================================================================

  subroutine trilinear_interpK_upper_8gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_weight, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W 
! for upper atmosphere pressure grid
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(8, kc_npress_upper, kc_ntemp_upper, kc_nweight)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_weight
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- trilinear_interpK_upper ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_upper) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))            
    else
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))        
    endif
    
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_upper(p_ref_index),log10pgrid_upper(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_upper) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    else
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    endif

    !write(*,*) "interp_temp", temp

    if (w_ref_index .eq. kc_nweight) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    else
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    endif

    !write(*,*) "interp_weight", weight
 
    ! perform trilinear interpolation between P,T,W

    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    !write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight

    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight

    !write(*,*) "ans", ans
    return

  end subroutine trilinear_interpK_upper_8gpt

  
!============================================================================

  subroutine trilinear_interpK_lower_16gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_weight, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
! for lower atmosphere pressure grid                                       
!
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(16, kc_npress_lower, kc_ntemp_lower, kc_nweight)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_weight
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- trilinear_interpK_lower ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_lower) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))           
    else
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))  
    endif
   
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_lower(p_ref_index),log10pgrid_lower(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_lower) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    else
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    if (w_ref_index .eq. kc_nweight) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    else
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    endif

    !write(*,*) "interp_weight", weight
 
    ! perform trilinear interpolation between P,T,W

    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    !write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight

    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight

    !write(*,*) "ans", ans

    return

  end subroutine trilinear_interpK_lower_16gpt


!========================================================================

  subroutine trilinear_interpK_upper_16gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_weight, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W 
! for upper atmosphere pressure grid
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(16, kc_npress_upper, kc_ntemp_upper, kc_nweight)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_weight
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- trilinear_interpK_upper ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_upper) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))            
    else
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))        
    endif
    
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_upper(p_ref_index),log10pgrid_upper(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_upper) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    else
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    endif

    !write(*,*) "interp_temp", temp

    if (w_ref_index .eq. kc_nweight) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    else
      weight = (species_weight - wgrid(w_ref_index))/(wgrid(w_ref_indexp1) - wgrid(w_ref_index))
    endif

    !write(*,*) "interp_weight", weight
 
    ! perform trilinear interpolation between P,T,W

    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    !write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight

    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight

    !write(*,*) "ans", ans
    return

  end subroutine trilinear_interpK_upper_16gpt

  
!============================================================================

  subroutine bilinear_interpK_lower_8gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) ::Kdata(8, kc_npress_lower, kc_ntemp_lower)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- bilinear_interpK_lower ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_lower) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))      
    else
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))
    endif
   
    !write(*,*) "------------------------------------------------------" 
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_lower(p_ref_index),log10pgrid_lower(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_lower) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    else
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "tgrid(t_ref_index)", tgrid(t_ref_index)
    !write(*,*) "tgrid(t_ref_indexp1)", tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! perform bilinear interpolation between P,T

    vbi(1) = kdata(ig, p_ref_index,   t_ref_index)
    vbi(2) = kdata(ig, p_ref_indexp1, t_ref_index)
    vbi(3) = kdata(ig, p_ref_index,   t_ref_indexp1)
    vbi(4) = kdata(ig, p_ref_indexp1, t_ref_indexp1)

    !write(*,*) "V", vbi(:)
 
    onemp = 1. - press
    onemt = 1. - temp

    ans = vbi(1)*onemp*onemt &
        + vbi(2)*press*onemt &
        + vbi(3)*onemp*temp &
        + vbi(4)*press*temp

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine bilinear_interpK_lower_8gpt


!============================================================================

  subroutine bilinear_interpK_upper_8gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) ::Kdata(8, kc_npress_upper, kc_ntemp_upper)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- bilinear_interpK_upper ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_upper) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))
    else
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))
    endif
   
    !write(*,*) "------------------------------------------------------" 
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_upper(p_ref_index), log10pgrid_upper(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_upper) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    else
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "tgrid(t_ref_index)", tgrid(t_ref_index)
    !write(*,*) "tgrid(t_ref_indexp1)", tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! perform bilinear interpolation between P,T

    vbi(1) = kdata(ig, p_ref_index,   t_ref_index)
    vbi(2) = kdata(ig, p_ref_indexp1, t_ref_index)
    vbi(3) = kdata(ig, p_ref_index,   t_ref_indexp1)
    vbi(4) = kdata(ig, p_ref_indexp1, t_ref_indexp1)

    !write(*,*) "V", vbi(:)
 
    onemp = 1. - press
    onemt = 1. - temp

    ans = vbi(1)*onemp*onemt &
        + vbi(2)*press*onemt &
        + vbi(3)*onemp*temp &        
        + vbi(4)*press*temp

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine bilinear_interpK_upper_8gpt
  

!============================================================================

  subroutine bilinear_interpK_lower_16gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) ::Kdata(16, kc_npress_lower, kc_ntemp_lower)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- bilinear_interpK_lower ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_lower) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))      
    else
      press = (pressure - log10pgrid_lower(p_ref_index)) / &
              (log10pgrid_lower(p_ref_indexp1) - log10pgrid_lower(p_ref_index))
    endif
   
    !write(*,*) "------------------------------------------------------" 
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_lower(p_ref_index),log10pgrid_lower(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_lower) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    else
      temp = (temperature - tgrid_lower(t_ref_index))/(tgrid_lower(t_ref_indexp1) - tgrid_lower(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "tgrid(t_ref_index)", tgrid(t_ref_index)
    !write(*,*) "tgrid(t_ref_indexp1)", tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! perform bilinear interpolation between P,T

    vbi(1) = kdata(ig, p_ref_index,   t_ref_index)
    vbi(2) = kdata(ig, p_ref_indexp1, t_ref_index)
    vbi(3) = kdata(ig, p_ref_index,   t_ref_indexp1)
    vbi(4) = kdata(ig, p_ref_indexp1, t_ref_indexp1)

    !write(*,*) "V", vbi(:)
 
    onemp = 1. - press
    onemt = 1. - temp

    ans = vbi(1)*onemp*onemt &
        + vbi(2)*press*onemt &
        + vbi(3)*onemp*temp &
        + vbi(4)*press*temp

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine bilinear_interpK_lower_16gpt


!============================================================================

  subroutine bilinear_interpK_upper_16gpt(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) ::Kdata(8, kc_npress_upper, kc_ntemp_upper)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(4) :: vbi
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- bilinear_interpK_upper ----------"
    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1

    vbi(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. kc_npress_upper) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))
    else
      press = (pressure - log10pgrid_upper(p_ref_index)) / &
              (log10pgrid_upper(p_ref_indexp1) - log10pgrid_upper(p_ref_index))
    endif
   
    !write(*,*) "------------------------------------------------------" 
    !write(*,*) "pressure", pressure
    !write(*,*) "p_ref_index", p_ref_index
    !write(*,*) "reference", log10pgrid_upper(p_ref_index), log10pgrid_upper(p_ref_indexp1)
    !write(*,*) "interp_press", press
  
    if (t_ref_index .eq. kc_ntemp_upper) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    else
      temp = (temperature - tgrid_upper(t_ref_index))/(tgrid_upper(t_ref_indexp1) - tgrid_upper(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "tgrid(t_ref_index)", tgrid(t_ref_index)
    !write(*,*) "tgrid(t_ref_indexp1)", tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! perform bilinear interpolation between P,T

    vbi(1) = kdata(ig, p_ref_index,   t_ref_index)
    vbi(2) = kdata(ig, p_ref_indexp1, t_ref_index)
    vbi(3) = kdata(ig, p_ref_index,   t_ref_indexp1)
    vbi(4) = kdata(ig, p_ref_indexp1, t_ref_indexp1)

    !write(*,*) "V", vbi(:)
 
    onemp = 1. - press
    onemt = 1. - temp

    ans = vbi(1)*onemp*onemt &
        + vbi(2)*press*onemt &
        + vbi(3)*onemp*temp &        
        + vbi(4)*press*temp

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine bilinear_interpK_upper_16gpt


!========================================================================

  subroutine interpKself(Kdata, ig, pressure, p_ref_index, temperature, t_ref_index, &
                           species_weight, w_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(ntot_gpt, ks_npress, ks_ntemp, ks_nweight)
    integer, intent(in) :: ig
    integer, intent(inout) :: p_ref_index
    integer, intent(inout) :: t_ref_index
    integer, intent(inout) :: w_ref_index

    real(r8), intent(in) :: pressure  
    real(r8), intent(in) :: temperature
    real(r8), intent(in) :: species_weight
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!

    p_ref_indexp1 = p_ref_index + 1    
    t_ref_indexp1 = t_ref_index + 1
    w_ref_indexp1 = w_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    ! find the interpolating factor for the various parameters.  at the tops of grids, extrapolate
    if (p_ref_index .eq. ks_npress) then
      p_ref_index = p_ref_index - 1
      p_ref_indexp1 = p_ref_index + 1
      press = (pressure - log10pgrid_self(p_ref_index))/(log10pgrid_self(p_ref_indexp1) - log10pgrid_self(p_ref_index))
    else
      press = (pressure - log10pgrid_self(p_ref_index))/(log10pgrid_self(p_ref_indexp1) - log10pgrid_self(p_ref_index))
    endif

       
    write(*,*) "interp_press", pressure, press
  
    if (t_ref_index .eq. ks_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    else
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    endif

    write(*,*) "interp_temp", temperature, temp

    if (w_ref_index .eq. ks_nweight) then
      w_ref_index = w_ref_index - 1
      w_ref_indexp1 = w_ref_index + 1
      weight = (species_weight - wgrid_self(w_ref_index))/(wgrid_self(w_ref_indexp1) - wgrid_self(w_ref_index))
    else
      weight = (species_weight - wgrid_self(w_ref_index))/(wgrid_self(w_ref_indexp1) - wgrid_self(w_ref_index))
    endif

    write(*,*) "interp_weight", species_weight, weight
 
    !perform trilinear interpolation between P,T,W

    !write(*,*) kdata(ig,p_ref_index, t_ref_index, w_ref_index)

    vtri(1) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_index)
    vtri(2) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_index)
    vtri(3) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_index)
    vtri(4) = kdata(ig, p_ref_index,   t_ref_index,   w_ref_indexp1)
    vtri(5) = kdata(ig, p_ref_indexp1, t_ref_index,   w_ref_indexp1)
    vtri(6) = kdata(ig, p_ref_index,   t_ref_indexp1, w_ref_indexp1)
    vtri(7) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_index)
    vtri(8) = kdata(ig, p_ref_indexp1, t_ref_indexp1, w_ref_indexp1)

    write(*,*) "V", vtri(:)
 
    onemp = 1. - press
    onemt = 1. - temp
    onemw = 1. - weight 
    ans = vtri(1)*onemp*onemt*onemw &
        + vtri(2)*press*onemt*onemw &
        + vtri(3)*onemp*temp*onemw &
        + vtri(4)*onemp*onemt*weight &
        + vtri(5)*press*onemt*weight &
        + vtri(6)*onemp*temp*weight &
        + vtri(7)*press*temp*onemw &
        + vtri(8)*press*temp*weight
  
    ans = abs(ans)
    write(*,*) "ans",ans

    return

  end subroutine interpKself


!========================================================================

  subroutine interpKself2(Kdata, iw, temperature, t_ref_index, ans) 

!------------------------------------------------------------------------
!
! Purpose: Linearly interpolate k-coefficient from reference P,T,W
!                                        
!------------------------------------------------------------------------
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in) :: Kdata(ntot_wavlnrng, ks_ntemp)
    integer, intent(in) :: iw
    integer, intent(inout) :: t_ref_index
    real(r8), intent(in) :: temperature
    real(r8), intent(out) :: ans

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: p_ref_indexp1
    integer :: t_ref_indexp1
    integer :: w_ref_indexp1
    integer :: ik
   
    real(r8) :: u_col
    real(r8) :: press
    real(r8) :: temp
    real(r8) :: weight
    real(r8) :: onemp
    real(r8) :: onemt
    real(r8) :: onemw
    real(r8), dimension(8) :: vtri
    
    logical :: interpp
    logical :: interpt

!------------------------------------------------------------------------
!
! Start Code
!

    t_ref_indexp1 = t_ref_index + 1

    vtri(:) = 0.
    ans = 0.

    if (t_ref_index .eq. ks_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    else
      temp = (temperature - tgrid_self(t_ref_index))/(tgrid_self(t_ref_indexp1) - tgrid_self(t_ref_index))
    endif


    ans = kdata(iw,t_ref_index)*(1.d0-temp)+kdata(iw,t_ref_index+1)*temp 
!! Experimental
!!EW Klugde, return answer at bottom temperaturem
!    ans = kdata(iw,t_ref_index+1)

!    write(*,*) "--------------------------------------------------------------"
!    write(*,*) temperature, tgrid_self(t_ref_index), tgrid_self(t_ref_indexp1)
!    write(*,*) kdata(iw,t_ref_index)*(1.d0-temp)+kdata(iw,t_ref_index+1)*temp 
!    write(*,*) "interpKself2", kdata(iw,t_ref_index), kdata(iw,t_ref_index+1)
!    write(*,*) "interpKself2", tgrid_self(t_ref_index)/temperature, tgrid_self(t_ref_indexp1)/temperature
!    write(*,*) "interpKself2", kdata(iw,t_ref_index)*tgrid_self(t_ref_index)/temperature, kdata(iw,t_ref_index+1)*tgrid_self(t_ref_indexp1)/temperature

    return

  end subroutine interpKself2


!============================================================================

  subroutine interpolate_cld(cldgrp, iw_indx, Rcld, Qcld, Wcld, Gcld, & 
                             Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)

!------------------------------------------------------------------------
!
! Purpose: Interpolate cloud optical data to mode effective radius for ice
!          ice particles and liquid cloud drops.
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    integer, intent(in) :: cldgrp    ! cloud group index, 1 => liquid, 2 => ice
    integer, intent(in) :: iw_indx  ! spectral interval index
    real(r8), intent(in) :: Rcld      ! model output radius of cloud particle
    real(r8), intent(out) :: Qcld     ! interpolated value of mass extinction coefficient
    real(r8), intent(out) :: Wcld     ! interpolated value of single scattering albedo
    real(r8), intent(out) :: Gcld     ! interpolated value of asymmetry parameter
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Qcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Qcldice
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Wcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Wcldice
    real(r8), intent(in), dimension(nrel,ntot_wavlnrng) :: Gcldliq 
    real(r8), intent(in), dimension(nrei,ntot_wavlnrng) :: Gcldice
 
!------------------------------------------------------------------------
!
! Local Variables
! 

    integer :: Rcld_ref_index
    integer :: ir
    real(r8) :: fr

!------------------------------------------------------------------------
!
! Start Code
! 
    ir = int(Rcld)
    fr = dble(Rcld-real(ir))

    if (cldgrp .eq. 1) then   ! interpolate for liquid cloud droplets

      ! if Rcld less than minimum, force to be minimum grid value 
      if (Rcld .le. minval(rel_grid)) then 
        Qcld = Qcldliq(1,iw_indx)
        Wcld = Wcldliq(1,iw_indx)
        Gcld = gcldliq(1,iw_indx)
      endif

      ! if Rcld greater than maximum, force to be maximum grid value 
      if (Rcld .ge. maxval(rel_grid)) then 
        Qcld = Qcldliq(nrel,iw_indx)
        Wcld = Wcldliq(nrel,iw_indx)
        Gcld = gcldliq(nrel,iw_indx)
      endif

      ! interpolate Rcld to grid
      if (Rcld .gt. minval(rel_grid) .and. Rcld .lt. maxval(rel_grid)) then
        Qcld = Qcldliq(ir,iw_indx)*(1.d0-fr)+Qcldliq(ir+1,iw_indx)*fr 
        Wcld = Wcldliq(ir,iw_indx)*(1.d0-fr)+Wcldliq(ir+1,iw_indx)*fr 
        Gcld = Gcldliq(ir,iw_indx)*(1.d0-fr)+Gcldliq(ir+1,iw_indx)*fr 
      endif   
 
    endif

    if (cldgrp .eq. 2) then ! interpolate for ice cloud droplet
 
      ! if Rcld less than minimum, force to be minimum grid value 
      if (Rcld .le. minval(rei_grid)) then 
        Qcld = Qcldice(1,iw_indx)
        Wcld = Wcldice(1,iw_indx)
        Gcld = gcldice(1,iw_indx)
      endif

      ! if Rcld greater than maximum, force to be maximum grid value 
      if (Rcld .ge. maxval(rei_grid)) then 
        Qcld = Qcldice(nrei,iw_indx)
        Wcld = Wcldice(nrei,iw_indx)
        Gcld = gcldice(nrei,iw_indx)
      endif

      ! interpolate Rcld to grid
      if (Rcld .gt. minval(rei_grid) .and. Rcld .lt. maxval(rei_grid)) then
        Qcld = Qcldice(ir,iw_indx)*(1.d0-fr)+Qcldice(ir+1,iw_indx)*fr 
        Wcld = Wcldice(ir,iw_indx)*(1.d0-fr)+Wcldice(ir+1,iw_indx)*fr 
        Gcld = Gcldice(ir,iw_indx)*(1.d0-fr)+Gcldice(ir+1,iw_indx)*fr  
      endif  

    endif
   
  end subroutine interpolate_cld


!============================================================================

  subroutine interpH2H2cia(Kdata, iw_indx, temperature, t_ref_index, ans )


!------------------------------------------------------------------------
!
! Purpose: Interpolate H2-H2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kh2h2_ntemp)
    integer, intent(in) :: iw_indx
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpH2H2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans = 0.

    if (t_ref_index .eq. kh2h2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_h2h2(t_ref_index))/(tgrid_h2h2(t_ref_indexp1) - tgrid_h2h2(t_ref_index))
    else
      temp = (temperature - tgrid_h2h2(t_ref_index))/(tgrid_h2h2(t_ref_indexp1) - tgrid_h2h2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! linear interpolation T
    vli(1) = kdata(iw_indx,   t_ref_index)
    vli(2) = kdata(iw_indx,   t_ref_indexp1)

    !write(*,*) "V", vli(:)
 
    ydiff = vli(2) - vli(1)
    ans = vli(1) + ydiff*temp 

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpH2H2cia


!============================================================================

  subroutine interpH2N2cia(Kdata, iw_indx, temperature, t_ref_index, ans )


!------------------------------------------------------------------------
!
! Purpose: Interpolate H2-N2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kh2n2_ntemp)
    integer, intent(in) :: iw_indx
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpH2N2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans = 0.

    if (t_ref_index .eq. kh2n2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_h2n2(t_ref_index))/(tgrid_h2n2(t_ref_indexp1) - tgrid_h2n2(t_ref_index))
    else
      temp = (temperature - tgrid_h2n2(t_ref_index))/(tgrid_h2n2(t_ref_indexp1) - tgrid_h2n2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! linear interpolation T
    vli(1) = kdata(iw_indx,   t_ref_index)
    vli(2) = kdata(iw_indx,   t_ref_indexp1)

    !write(*,*) "V", vli(:)
 
    ydiff = vli(2) - vli(1)
    ans = vli(1) + ydiff*temp 

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpH2N2cia


!============================================================================

  subroutine interpN2N2cia(Kdata, iw_indx, temperature, t_ref_index, ans )


!------------------------------------------------------------------------
!
! Purpose: Interpolate N2-N2 collision induced absorption coefficients
!
!------------------------------------------------------------------------
 
    implicit none

!------------------------------------------------------------------------
!
! Arguments
! 
    real(r8), intent(in) :: Kdata(ntot_wavlnrng,kn2n2_ntemp)
    integer, intent(in) :: iw_indx
    real(r8), intent(in) :: temperature
    integer, intent(inout) :: t_ref_index
    real(r8), intent(out) :: ans  

!------------------------------------------------------------------------
!
! Local Variables
!

    integer :: t_ref_indexp1
    real(r8) :: temp
    real(r8), dimension(2) :: vli
    real(r8) :: ydiff

!------------------------------------------------------------------------
!
! Start Code
!
    !write(*,*) "-------- interpN2N2cia ----------"
    t_ref_indexp1 = t_ref_index + 1
    ans = 0.

    if (t_ref_index .eq. kh2n2_ntemp) then
      t_ref_index = t_ref_index - 1
      t_ref_indexp1 = t_ref_index + 1
      temp = (temperature - tgrid_n2n2(t_ref_index))/(tgrid_n2n2(t_ref_indexp1) - tgrid_n2n2(t_ref_index))
    else
      temp = (temperature - tgrid_n2n2(t_ref_index))/(tgrid_n2n2(t_ref_indexp1) - tgrid_n2n2(t_ref_index))
    endif

    !write(*,*) "temperature", temperature
    !write(*,*) "t_ref_index", t_ref_index
    !write(*,*) "reference", tgrid(t_ref_index),tgrid(t_ref_indexp1)
    !write(*,*) "interp_temp", temp

    ! linear interpolation T
    vli(1) = kdata(iw_indx,   t_ref_index)
    vli(2) = kdata(iw_indx,   t_ref_indexp1)

    !write(*,*) "V", vli(:)
 
    ydiff = vli(2) - vli(1)
    ans = vli(1) + ydiff*temp 

    !write(*,*) "ans", ans
    !write(*,*) "------------------------------------------------------" 

    return

  end subroutine interpN2N2cia

end module rad_interp_mod
