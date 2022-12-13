module exo_condense_mod
!===========================================================================================
! Module containing routines to condense exotic cloud types
! Author:  storyofthewolf (eric.wolf@colorado.edu)
! Date: November 2020
!
! CO2 Cloud condensation routine based on 
! Forget F. Hourdin, F. & Talagrand, O. (1998) 
!  "CO2 Snow Fall on Mars: Simulation with a General Circulation model." 
!  Icarus,   
!
!
! NOTES:
! from physics_update:  Update constituents, all schemes use time split q: no tendency kept  
! physics_update does not update the constitutent mixing ratios due to timesplitting
! condense cloud
! 
! See CO2 Module Options located below
! Describe

use shr_const_mod, only: SHR_CONST_G, SHR_CONST_CPDAIR, SHR_CONST_CO2LATSUB, &
                         SHR_CONST_PI, SHR_CONST_DENSITYCO2ICE, SHR_CONST_DENSITYCO2FR, &
                         SHR_CONST_SOILCP, SHR_CONST_SOILBD, SHR_CONST_RGAS, SHR_CONST_AVOGAD
use shr_kind_mod,  only: r8=>shr_kind_r8
use ppgrid,        only: pcols, pver, pverp
use constituents,  only: pcnst, cnst_get_ind
use exoplanet_mod, only: do_exo_condense_co2, exo_co2vmr
use physics_types, only: physics_state, physics_ptend, physics_ptend_init
use cam_history,   only: outfld
use physconst,     only: rair, gravit
use spmd_utils,    only: masterproc
use abortutils,    only: endrun

implicit none
private
save

public exo_condense_register
public exo_condense_init
public exo_condense_init_cnst
public exo_condense_implements_cnst
public exo_condense_diag_calc
public exo_condense_tend


! Global variables
integer :: ixcldice_co2
integer :: rei_co2_idx = 0
integer :: cicewp_co2_idx = 0
!logical :: cnst_reg = .false.  ! not used


!!--------------------------------------------------------------------------------
!!
!! Modifiable parameters and module options
!!

!! constituent names 
integer, parameter :: ncnst = 1
character(len=10), parameter :: cnst_names(ncnst) = (/'CLDICE_CO2'/)

!!
!! CO2 cloud condensation module options 
!!
integer,  parameter :: nsubsteps_co2cld = 20      ! number of substeps for CO2 cloud computation
real(r8), parameter :: co2_supersat = 1.35        ! CO2 supersaturation threshold for nucleation


!!-----  CO2 cloud radii options  -----
!! select one only

! use fixed value for CO2 ice cloud radii
logical,  parameter :: co2_reff_fixed = .false.       ! 
real(r8), parameter :: rei_co2_fixed = 50.0           ! microns

! use fixed number of CCN everywhere (the LMD method, Forget et al. 2013)
logical,  parameter :: co2_reff_fixedccn = .true.    ! use a fixed value for co2 ice cloud ccn (LMD method)
real(r8), parameter :: ncloud_co2 = 1.0e4            ! fixed cloud number concentration per kg of air

! use scale height dependent variation in CO2 cloud radii (the DRAMATIC method, Kuroda et al. 2013) 
logical,  parameter :: co2_reff_scaleheight = .false.         
real(r8), parameter :: r0 = 50.e-6  ! [m]             
real(r8), parameter :: h0 = 20.0e3  ! [m]               

! use a scale height dependent variation on CCN
logical,  parameter :: co2_reff_scaleheightccn = .false.      
real(r8), parameter :: h0cn = 20.0e3  ! [m]               


!! cloud condesation schemes  (eventually remove)
logical, parameter :: do_forget1998 =.true.
logical, parameter :: do_lmdg =.false.


!======================================================================= 
contains
!======================================================================= 


!=======================================================================
subroutine exo_condense_register
!---------------------------------------------------------------------- 
! Purpose: Register the constituents 
!
  use constituents, only: cnst_add
  use physconst,    only: mwdry, cpair
  use physics_buffer, only : pbuf_add_field, dtype_r8, pbuf_times

  !
  ! Start Code
  !

  ! Must be added to the physics buffer even if do_exo_condense_co2 = .false.
  call pbuf_add_field('REI_CO2',       'physpkg', dtype_r8, (/pcols,pver/), rei_co2_idx) 
  call pbuf_add_field('CICEWP_CO2', 'physpkg', dtype_r8,(/pcols,pver/), cicewp_co2_idx)

  if (.not. do_exo_condense_co2) return

  call cnst_add('CLDICE_CO2', mwdry, cpair, 0._r8, ixcldice_co2, &
                 longname='Grid box averaged CO2 ice cloud', is_convtran1=.true.)


end subroutine exo_condense_register


!=======================================================================
subroutine exo_condense_init
!----------------------------------------------------------------------
! Purpose: Initialize exo_condense module
!
  use cam_history, only: addfld, phys_decomp

  ! local variables
  character(len=16) :: sampling_seq
  integer:: option_count 

  !
  ! Start Code
  !

  if (.not. do_exo_condense_co2) return

  !! checking module options selections
  option_count = 0
  if (co2_reff_fixed) then 
    if(masterproc) write(*,*) "using fixed CO2 cloud radii ", rei_co2_fixed*1.0e6, " [microns]"
    option_count = option_count + 1
  endif
  if(co2_reff_fixedccn) then 
    if(masterproc) write(*,*) "using fixed CO2 particle number ", ncloud_co2, " [per kg air]"
    option_count = option_count + 1
  endif
  if(co2_reff_scaleheight) then
    if(masterproc) write(*,*) "using CO2 radii variation with scaleheight ", r0, " [microns] ", h0, " [m]"
    option_count = option_count + 1
  endif
  if(co2_reff_scaleheightccn) then 
    if(masterproc) write(*,*) "using CO2 particle number variation with scaleheight "
    option_count = option_count + 1
  endif
  if (option_count .gt. 1) then
    write(*,*) "Cannot select more than one method for CO2 radii handling"
    call endrun('ERROR: exo_condense_mod, CO2 reff option setting duplication')
  endif

  !
  ! add history variables for CO2 Cloud Condensation
  !
  call addfld('CLDICE_CO2', 'kg/kg   ', pver, 'A', 'CO2 ice cloud condensate', phys_decomp )

  call addfld('CLDICE_CO2_COL', 'g/m2   ', pver, 'A', 'CO2 ice cloud condensate', phys_decomp )
  call addfld('CLDICE_CO2_RHO', 'kg/m3   ', pver, 'A', 'CO2 ice cloud condensate', phys_decomp )

  call addfld('CLDICE_CO2_TEND', 'kg/kg/s   ', pver, 'A', 'CO2 ice cloud tendency, cond -heat -pot', phys_decomp )
  call addfld('CLDICE_CO2_COND_TEND', 'kg/kg/s   ', pver, 'A', 'CO2 ice cloud condensation tendency', phys_decomp )
  call addfld('CLDICE_CO2_HEAT_TEND', 'kg/kg/s   ', pver, 'A', 'CO2 ice cloud heat ice to layer temperature tendency', phys_decomp )
  call addfld('CLDICE_CO2_POT_TEND', 'kg/kg/s   ', pver, 'A', 'CO2 ice cloud potential energy from falling tendency', phys_decomp )

  call addfld('CLDICE_CO2_SED_TEND', 'kg/kg/s   ', pver, 'A', 'CO2 ice cloud sedimentation tendency', phys_decomp )
  call addfld('CLDICE_CO2_TEMP_TEND', 'K/s   ', pver, 'A', 'CO2 ice cloud temperature tendency', phys_decomp )

  call addfld('CLDICE_CO2_PVEL', 'Pa/s   ', pverp, 'A', 'CO2 ice cloud fall speeds pressure units', phys_decomp )
  call addfld('CLDICE_CO2_VFALL', 'm/s   ', pverp, 'A', 'CO2 ice cloud fall speeds', phys_decomp )

  call addfld ('REI_CO2',       'micron',   pver, 'A', 'CO2 ice cloud effective particle radius'  ,phys_decomp)
  call addfld ('TGCLDIWP_CO2', 'gram/m2' ,1,    'A','Total grid-box CO2 ice cloud path'   ,phys_decomp, sampling_seq=sampling_seq)

  call addfld ('CO2_SNOWFALL_RATE','kg/m2/s',1,    'A','CO2 snow fall rate'   ,phys_decomp, sampling_seq=sampling_seq)


end subroutine exo_condense_init


!======================================================================= 
function exo_condense_implements_cnst(name)
!----------------------------------------------------------------------
! Purpose:  Return true if specified constituent is implemented
!
  ! input arguments
  character(len=*), intent(in) :: name      ! constituent name
  
  ! local varaibles
  logical :: exo_condense_implements_cnst     ! return value 
  integer :: m

  !
  ! Start Code
  !
  exo_condense_implements_cnst = .false.

  do m = 1, ncnst
    if (name == cnst_names(m)) then
      exo_condense_implements_cnst = .true.
      return
     end if
  end do

end function exo_condense_implements_cnst


!======================================================================= 
subroutine exo_condense_init_cnst(name, q, gcid)
!---------------------------------------------------------------------- 
! Purpose: Initialize the cloud  mixing ratios, if 
!          they are not read from the initial file 
!

  ! input arguments
  character(len=*), intent(in)  :: name     ! constituent name
  real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
  integer,          intent(in)  :: gcid(:)  ! global column id 

  !
  ! Start Code
  !
  if ( name == cnst_names(1) ) then 
     q = 0.0_r8
  endif

end subroutine exo_condense_init_cnst


!=======================================================================
subroutine exo_condense_diag_calc(state, pbuf)
!-----------------------------------------------------------------------
! Purpose: Compute co2 ice path
!
  use physics_types, only: physics_state 
  use physics_buffer,only: physics_buffer_desc, pbuf_get_field

  ! input arguments
  type(physics_state), intent(in)    :: state        ! state variables
  type(physics_buffer_desc), pointer :: pbuf(:)

  real(r8), pointer :: cicewp_co2(:,:)    ! in-cloud cloud ice water path 
  real(r8) :: gicewp_co2(pcols,pver)      ! grid-box cloud ice water path
  real(r8) :: tgicewp_co2(pcols)          ! Vertically integrated ice water path
  
  ! local arguments
  integer :: i,k
  integer :: ncol, lchnk

  if (.not. do_exo_condense_co2) return

  ! Start code
  lchnk = state%lchnk
  ncol  = state%ncol

  if (cicewp_co2_idx>0) then
    call pbuf_get_field(pbuf, cicewp_co2_idx, cicewp_co2 )
  else
    allocate(cicewp_co2(pcols,pver))
  endif

  do k=1,pver
    do i = 1,ncol
      gicewp_co2(i,k) = state%q(i,k,ixcldice_co2)*state%pdel(i,k)/SHR_CONST_G*1000.0_r8 ! Grid box co2 ice condensate  [g/m2]
      cicewp_co2(i,k) = gicewp_co2(i,k)  !/ max(0.01_r8,cld(i,k))   ! in cloud condensate, for now, no cloud fraction treatement
!      write(*,*) "exo_condense_diag_calc",i,k, state%q(i,k,ixcldice_co2), cicewp_co2(i,k) 
    end do
  enddo

  ! calculate column amount of CO2 cice cloud
  tgicewp_co2(:ncol) = 0._r8
  do k=1,pver
    tgicewp_co2(:ncol)  = tgicewp_co2(:ncol) + gicewp_co2(:ncol,k)
  end do

  call outfld('TGCLDIWP_CO2',tgicewp_co2, pcols,lchnk)
  call outfld('CLDICE_CO2_COL' ,cicewp_co2    ,          pcols,lchnk)

end subroutine exo_condense_diag_calc


!=======================================================================
subroutine exo_condense_tend(state, pbuf, cam_in, cam_out, ptend)
!-----------------------------------------------------------------------
! Purpose:  Driver for condensing species other than H2O.
!           Called from tphysac. 
!           Condensation routines for different gases can be called from 
!           here in future revisions.
!           Presently only CO2 is included
!
  use physics_buffer,     only: physics_buffer_desc
  use camsrfexch,         only: cam_out_t, cam_in_t

  ! input arguments
  type(physics_state),        intent(inout) :: state     ! State variables
  type(physics_buffer_desc),  pointer       :: pbuf(:)
  type(cam_in_t),             intent(in)    :: cam_in
  type(cam_out_t),            intent(inout) :: cam_out
  type(physics_ptend),        intent(out)   :: ptend     ! Package tendencies

  ! local variabls
  type(physics_ptend)   :: ptend_loc 
  integer :: lchnk, ncol, i
  real(r8), dimension(pver) :: cld_mass_out, cld_reff_out
  real(r8) :: srf_mass_out
  real(r8), dimension(pcols,pver) :: cld_mass_temp, cld_reff_temp
  logical  :: lq(pcnst)


  !
  ! Start Code
  !

  lq(:) = .false.
  call cnst_get_ind('CLDICE_CO2', ixcldice_co2)
  lq(ixcldice_co2) = .true.
  call physics_ptend_init(ptend,state%psetcols, 'exo_condense_tend', ls=.true., lq=lq)

  ! CO2 condensation 
  call exo_condense_co2(state, pbuf, cam_in%ts, cam_in%co2dp, cam_out%co2srf_snow, ptend)

end subroutine exo_condense_tend


!! Private subroutines !!

!=======================================================================
subroutine exo_condense_co2(state, pbuf, ext_ts, ext_co2dp, ext_co2srf_snow, ptend)
!-----------------------------------------------------------------------
! Purpose: CO2 condensation routine.
!          Considers atmospheric condensation, sublimation and precipitation.
!          Condensation directory onto the surface, and snow accumulation is
!          handled within the surface models (land, ice)
!    
! Notes:  Based on Forget F. Hourdin, F. & Talagrand, O. (1998)
!         "CO2 Snow Fall on Mars: Simulation with a General Circulation model."
!          Icarus 131, 302-316

  use shr_condense_mod, only: get_condense_temp
  use time_manager,     only: get_step_size
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_get_index


  ! arguments 
  type(physics_state), intent(inout)    :: state
  type(physics_buffer_desc), pointer :: pbuf(:)
  real(r8), dimension(pcols),  intent(in)  :: ext_ts           ! [K] surface temperature radiative  (from cam_in%ts)     
  real(r8), dimension(pcols),  intent(in)  :: ext_co2dp        ! [kg/m2] co2 surface ice mass (from cam_in%co2dp)
  real(r8), dimension(pcols),  intent(out) :: ext_co2srf_snow  ! [kg/m2] co2 snowfall 
  type(physics_ptend), intent(inout)   :: ptend                ! Package tendency struct


  ! local variables
  integer :: ncol, lchnk
  real(r8), dimension(pcols, pver)  :: tmid_tmp   ! [K]       mid-layer temperature, temporary
  real(r8), dimension(pcols, pver)  :: cld_tmp    ! [kg/kg]   mid-layer cloud mass mixing ratio, temporary
  real(r8), pointer, dimension(:,:) :: cld_reff   ! cloud effective radius [m]      
  real(r8), dimension(pcols, pverp) :: cld_pvel   ! [Pa/s]    cloud fall speed, pressure units
  real(r8), dimension(pcols, pver)  :: cld_vfall  ! [m/2]     cloud fall speed, velocity units
  real(r8), dimension(pcols, pver)  :: cld_rho    ! [kg/m3]   cloud mass density, diagnostic
  real(r8), dimension(pcols, pver)  :: tcond      ! [K]       atmosphere condensation temperature grid        
  real(r8), dimension(pcols, pver)  :: tnuc       ! [K]       atmosphere nucleation temperature grid        
  real(r8), dimension(pcols, pver)  :: sed_tend   ! [kg/kg/s] sedimentation tendency
  real(r8), dimension(pcols, pver)  :: cld_tend   ! [kg/kg/s] total cloud tendency
  real(r8), dimension(pcols, pver)  :: cond_tend  ! [kg/kg/s] condensation tendency 
  real(r8), dimension(pcols, pver)  :: pot_tend   ! [kg/kg/s] tendency from potential energy release of falling ice  
  real(r8), dimension(pcols, pver)  :: heat_tend  ! [kg/kg/s] tendency from energy for heating falling mass
  real(r8), dimension(pcols, pver)  :: temp_tend  ! [K/s]     temperature tendency
  real(r8), dimension(pcols) :: sf_tend           ! [kg/m2/s] surface flux of cloud material
  real(r8), dimension(pcols) :: co2_snowfall      ! [kg/m2]   accumulated co2 surface ice mass during the timestep
  real(r8) :: tcond_tmp, tnuc_tmp                 ! [K]       condensation and nucleation temperatures,  temporary 
  real(r8) :: dtime, subtime                      ! [s]       timestep and substep
  real(r8) :: cpco2ice                            ! [J/kg/K]  specific heat of co2 ice, 349+4.8T
  integer :: nt, k, i, kstar
  logical  :: lq(pcnst)

  ! across substep diagnostics
  real(r8), dimension(pcols, pver)  :: sed_tend_diag   ! [kg/kg/s] sedimentation tendency, diagnostic 
  real(r8), dimension(pcols, pver)  :: cld_tend_diag   ! [kg/kg/s] total cloud tendency, diagnostic 
  real(r8), dimension(pcols, pver)  :: cond_tend_diag  ! [kg/kg/s] condensation tendency, diagnostic 
  real(r8), dimension(pcols, pver)  :: pot_tend_diag   ! [kg/kg/s] tendency from potential energy release of falling ice, diagnostic 
  real(r8), dimension(pcols, pver)  :: heat_tend_diag  ! [kg/kg/s] tendency from energy for heating falling mass, diagnostic 
  real(r8), dimension(pcols, pver)  :: temp_tend_diag  ! [K/s]     temperature tendency, diagnostic 


  ! initialize some things to zero
  cld_tmp(:,:)       = 0.0     
  cld_pvel(:,:)      = 0.0    
  cld_rho(:,:)       = 0.0    
  tmid_tmp(:,:)      = 0.0    
  sed_tend(:,:)      = 0.0    
  cond_tend(:,:)     = 0.0    
  temp_tend(:,:)     = 0.0    
  pot_tend(:,:)      = 0.0     
  heat_tend(:,:)     = 0.0    
  sf_tend(:)         = 0.0      
  co2_snowfall(:)    = 0.0
  ext_co2srf_snow(:) = 0.0

  ! across substep diagnostics
  cld_tend_diag(:,:)  = 0.0
  sed_tend_diag(:,:)  = 0.0 
  cond_tend_diag(:,:) = 0.0
  temp_tend_diag(:,:) = 0.0
  pot_tend_diag(:,:)  = 0.0
  heat_tend_diag(:,:) = 0.0

  dtime = get_step_size()           ! seconds
  subtime = dtime/nsubsteps_co2cld  ! seconds

  lchnk = state%lchnk
  ncol  = state%ncol 

  lq(:) = .false.
  call cnst_get_ind('CLDICE_CO2', ixcldice_co2)
  lq(ixcldice_co2) = .true.
  call physics_ptend_init(ptend,state%psetcols, 'exo_condense_co2_atm', ls=.true., lq=lq)
  call pbuf_get_field(pbuf, rei_co2_idx, cld_reff)

  !
  !------------------------------------------------------
  ! Start Code
  !

  ! set temporary arrays from the state structure
  ! before start of substep loop
  do k=1,pver
    do i=1, ncol
      tmid_tmp(i,k) = state%t(i,k)
      cld_tmp(i,k)  = state%q(i,k,ixcldice_co2)
    enddo
  enddo


  ! calculate condensation and nucleation temperatures
  ! save to arrays
  do k=1,pver
    do i=1,ncol
      call get_condense_temp(1, state%pmid(i,k)*exo_co2vmr, tcond_tmp)
      call get_condense_temp(1, state%pmid(i,k)*exo_co2vmr/co2_supersat, tnuc_tmp)
      tcond(i,k) = tcond_tmp
      tnuc(i,k) = tnuc_tmp
    enddo
  enddo
       
       
  !-----------------------------------------------------------
  ! Atmospheric condensation, sublimation, and sedimentation
  !-----------------------------------------------------------

  !! -----------------------
  !! Start substep loop here
  !! ------------------------
  do nt = 1,nsubsteps_co2cld
  
   ! add tendencies from external processes on substep loop increment
!    do k = 1,pver
!      do i = 1,ncol
!       cld_tmp(i,k) = cld_tmp(i,k) + ptend%q(i,k,ixcldice_co2) * subtime
!       tmid_tmp(i,k) = tmid_tmp(i,k) + ptend%s(i,k)/SHR_CONST_CPDAIR * subtime
!     enddo
!   enddo

   ! determine particle sizes
    call calc_co2_reff(ncol, state%zm, cld_tmp, cld_reff)

    ! calculate particle sedimentation velocities and tendency
    call exo_cloud_sediment_vel(ncol, tmid_tmp, state%pmid, state%pdel, cld_tmp, cld_reff, cld_vfall, cld_pvel)
    call exo_cloud_sediment_tend(ncol, subtime, tmid_tmp, state%pmid, state%pint, state%pdel, &
                                  cld_tmp, cld_pvel, tcond, sed_tend, sf_tend) 

   ! add sedimentation tendencies
   do k = 1,pver
     do i = 1,ncol
       cld_tmp(i,k) = cld_tmp(i,k) + sed_tend(i,k) * subtime
!!if(masterproc) write(*,*) "sed_tend", nt,i,k,sed_tend(i,k), cld_tmp(i,k)
     enddo
   enddo

   ! surface CO2 accumulations during timestep
   do i=1, ncol
     co2_snowfall(i) = co2_snowfall(i) + sf_tend(i)*subtime   
!!if(masterproc) write(*,*) "sf_tend", nt,i,k,sf_tend(i), co2_snowfall(i), subtime
   enddo


   ! do condensation and sublimation in the atmospheres

   !! FORGET 1998 method
if(do_forget1998) then

   ! top layer of model
   k=1
   do i=i,ncol
     if ( (tmid_tmp(i,k) .lt. tnuc(i,k) ) ) then
       cond_tend(i,k) = SHR_CONST_CPDAIR/SHR_CONST_CO2LATSUB * (tcond(i,k)-tmid_tmp(i,k)) / subtime
       temp_tend(i,k) = (tcond(i,k)-tmid_tmp(i,k)) / subtime
     endif
   enddo
   ! layers below
   do k = 2,pver
     do i = 1,ncol
       ! nucleation occurs or cloud already exists
       if ( (tmid_tmp(i,k) .lt. tnuc(i,k) ) .or. (cld_tmp(i,k) .gt. 1.e-20) ) then 
         kstar=max(1,k-1) ! layer above
         !! condensation tendnecy [kg/kg/s]
         cond_tend(i,k) = SHR_CONST_CPDAIR/SHR_CONST_CO2LATSUB * (tcond(i,k)-tmid_tmp(i,k)) / subtime  

         !! tendency from potential energy released by ice falling into the grid box
         pot_tend(i,k)  = 1./SHR_CONST_CO2LATSUB*(SHR_CONST_G*(state%zm(i,kstar) - state%zm(i,k))) * sed_tend(i,k) 

         !! tendency from energy to heat ice to Tc of lower level as it falls into grid box
         cpco2ice = 349.0 + 4.8*tmid_tmp(i,k)
         heat_tend(i,k) = 1./SHR_CONST_CO2LATSUB*(cpco2ice*(tmid_tmp(i,kstar) - tmid_tmp(i,k))) * sed_tend(i,k)

         !! add tendencies
         cld_tend(i,k) = cond_tend(i,k) - pot_tend(i,k) - heat_tend(i,k)
         temp_tend(i,k) = (tcond(i,k)-tmid_tmp(i,k))/subtime

         ! cloud sublimes entirely
         if ( (cld_tmp(i,k) .lt. -cond_tend(i,k)*subtime) .and. (cld_tmp(i,k) .gt. 0.0) ) then
           cld_tend(i,k) = -cld_tmp(i,k) / subtime
           temp_tend(i,k) = 1./ SHR_CONST_CPDAIR * &
                 (-SHR_CONST_CO2LATSUB + &
                   SHR_CONST_G*(state%zm(i,kstar) - state%zm(i,k)) + &
                   cpco2ice*(tmid_tmp(i,kstar) - tmid_tmp(i,k)))  * cld_tmp(i,k) / subtime
         endif

         cld_tmp(i,k)  = cld_tmp(i,k)  + cld_tend(i,k)*subtime 
         tmid_tmp(i,k) = tmid_tmp(i,k) + temp_tend(i,k)*subtime 
        endif
        if (cld_tmp(i,k) .lt. 0.0) cld_tmp(i,k) = 0.0   ! remove small negative values that sometimes occur
        ! diagnostics across substeps
        sed_tend_diag(i,k)  = sed_tend_diag(i,k)  + sed_tend(i,k)/nsubsteps_co2cld
        cond_tend_diag(i,k) = cond_tend_diag(i,k) + cond_tend(i,k)/nsubsteps_co2cld
        temp_tend_diag(i,k) = temp_tend_diag(i,k) + temp_tend(i,k)/nsubsteps_co2cld
        pot_tend_diag(i,k)  = pot_tend_diag(i,k)  + pot_tend(i,k)/nsubsteps_co2cld
        heat_tend_diag(i,k) = heat_tend_diag(i,k) + heat_tend(i,k)/nsubsteps_co2cld
      enddo
    enddo
endif
!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!
if(do_lmdg) then
  !! old way modeled after LMDG
   !! do condensation and sublimation in the atmospheres
   do k = 1,pver
     do i = 1,ncol

       ! if nucleation/condensation or is occuring, or if cloud is already present in grid
       if ( (tmid_tmp(i,k) .lt. tnuc(i,k) ) .or. (cld_tmp(i,k) .gt. 0.0) ) then 
!!if (masterproc) write(*,*) "condense CO2 cloud", nt,i,k, tnuc(i,k), tcond(i,k), tmid_tmp(i,k), cld_tmp(i,k), state%q(i,k,ixcldice_co2)
          cond_tend(i,k) = SHR_CONST_CPDAIR/SHR_CONST_CO2LATSUB * (tcond(i,k)-tmid_tmp(i,k)) / subtime        !!! condensation tendnecy [kg/kg/s]
          temp_tend(i,k) = (tcond(i,k)-tmid_tmp(i,k))/subtime
         ! cloud sublimes entirely
         if ( (cld_tmp(i,k) .lt. -cond_tend(i,k)*subtime) .and. (cld_tmp(i,k) .gt. 0.0) ) then
            cond_tend(i,k) = -cld_tmp(i,k)/subtime
            temp_tend(i,k) = -cld_tmp(i,k)/(SHR_CONST_CPDAIR/SHR_CONST_CO2LATSUB)/subtime
          endif
          cld_tmp(i,k)  = cld_tmp(i,k)  + cond_tend(i,k)*subtime 
          tmid_tmp(i,k) = tmid_tmp(i,k) + temp_tend(i,k)*subtime  ! sets atmosphere temperature equal to condensation temperature
!!!if (masterproc)  write(*,*), "cond_tend", nt,i,k, cond_tend(i,k), cld_tmp(i,k)
!!!if (masterproc)  write(*,*), "temp_tend", nt,i,k, temp_tend(i,k), tmid_tmp(i,k)
        endif
        if (cld_tmp(i,k) .lt. 0.0) cld_tmp(i,k) = 0.0   ! remove small negative values that sometimes occur
        !cld_rho(i,k) = state%pmid(i,k) / (rair * tmid_tmp(i,k))*cld_tmp(i,k)  ! diagnostic for cloud density
      enddo
    enddo
endif
!!!!!!!!!!!!!!!!!!!!

  enddo    !! End Substep Loop




  ! Compute diagnostic quantities and tendencies after substep loop
  do k = 1,pver
     do i = 1,ncol
       cld_rho(i,k) = state%pmid(i,k) / (rair * tmid_tmp(i,k))*cld_tmp(i,k)  ! diagnostic for cloud density
       cld_tend_diag(i,k) =  (cld_tmp(i,k) - state%q(i,k,ixcldice_co2)) / dtime 
!if (masterproc)  write(*,*),  cld_tmp(i,k), state%q(i,k,ixcldice_co2), ( cld_tmp(i,k) - state%q(i,k,ixcldice_co2) ) / dtime 
       ptend%q(i,k,ixcldice_co2) =  (cld_tmp(i,k) - state%q(i,k,ixcldice_co2)) / dtime 
!if (masterproc)  write(*,*), "after", ptend%q(i,k,ixcldice_co2)
       ptend%s(i,k) =   (tmid_tmp(i,k) - state%t(i,k)) * SHR_CONST_CPDAIR / dtime 
     enddo
  enddo

  ! set snowfall output
  do i=1,ncol
    ext_co2srf_snow(i) = co2_snowfall(i)   ![kg/m2]  
  enddo
  
  cld_reff(:,:) = cld_reff(:,:)*1.0e6  ! convert from meters to microns
  call outfld('REI_CO2', cld_reff, pcols, lchnk)
  call outfld('CLDICE_CO2_RHO' ,cld_rho    ,          pcols,lchnk)
  call outfld('CO2_SNOWFALL_RATE',co2_snowfall/dtime      , pcols, lchnk) 

  call outfld('CLDICE_CO2_TEND', cld_tend_diag, pcols, lchnk) 
  call outfld('CLDICE_CO2_COND_TEND', cond_tend_diag, pcols, lchnk) 
  call outfld('CLDICE_CO2_POT_TEND', pot_tend_diag, pcols, lchnk) 
  call outfld('CLDICE_CO2_HEAT_TEND', heat_tend_diag, pcols, lchnk) 
  call outfld('CLDICE_CO2_SED_TEND', sed_tend_diag, pcols, lchnk)
  call outfld('CLDICE_CO2_TEMP_TEND', temp_tend_diag, pcols, lchnk)

  call outfld('CLDICE_CO2_PVEL', cld_pvel, pcols, lchnk)
  call outfld('CLDICE_CO2_VFALL', cld_vfall, pcols, lchnk)


end subroutine exo_condense_co2


!=======================================================================
subroutine calc_co2_reff(ncol, zm, cld_mass, reff)
!-----------------------------------------------------------------------
! Purpose: Compute the effective radii of co2 ice particles
!
  ! input arguments
  integer, intent(in) :: ncol  ! columns
  real(r8), intent(in), dimension(pcols, pver)  :: zm       ! [m] geopotential height at midpoints
  real(r8), intent(in), dimension(pcols, pver)  :: cld_mass ! [kg/kg] tracer mixing ratios
  real(r8), intent(out), dimension(pcols, pver) :: reff     ! [m] co2 ice particles radii

  ! local variables
  integer :: i, k

  ! 
  ! Start Code
  !
  if (co2_reff_fixed) then                   ! used a fixed particle size
    reff(:,:) = rei_co2_fixed*1.0e-6         ! CO2 ice effective radii [m]
  elseif (co2_reff_fixedccn) then            ! use a fixed number of CDNC/kg air
    do i=1, ncol
      do k=1, pver
        reff(i,k) = (3. * cld_mass(i,k) /(4.*ncloud_co2*SHR_CONST_PI*SHR_CONST_DENSITYCO2ICE))**(1./3)
        reff(i,k) = min(max(reff(i,k),1.e-6),1000.e-6)
      enddo
    enddo
  elseif (co2_reff_scaleheight) then   ! set by scale height
    do i=1, ncol
      do k=1, pver
        reff(i,k) = r0*exp(-zm(i,k)/h0)
      enddo
    enddo 
  elseif (co2_reff_scaleheightccn) then   ! my way, more aerosols near ground, less aloft
    do i=1, ncol
      do k=1, pver
        reff(i,k) = (3. * cld_mass(i,k) /(4.*ncloud_co2*exp(-zm(i,k)/h0cn)*SHR_CONST_PI*SHR_CONST_DENSITYCO2ICE))**(1./3)
        reff(i,k) = min(max(reff(i,k),1.e-6),100.e-6) 
      enddo
    enddo
  endif

end subroutine calc_co2_reff
!==================================================================


!=======================================================================
subroutine exo_cloud_sediment_vel(ncol, ext_tmid, ext_pmid, ext_pdel, ext_cld, ext_reff, vfall_out, cld_pvel_out)
!-----------------------------------------------------------------------
! Calculate the stokes fall velocity
! Note: need to add Cunningham slip-correction factor
!
  ! input arguments
  integer, intent(in) :: ncol
  real(r8), intent(in), dimension(pcols,pver)  :: ext_tmid        ! temperature (K)
  real(r8), intent(in), dimension(pcols,pver)  :: ext_pmid     ! pressure of midpoint levels (Pa)
  real(r8), intent(in), dimension(pcols,pver)  :: ext_pdel         ! pressure diff across layer (Pa)
  real(r8), intent(in), dimension(pcols,pver)  :: ext_cld      ! cloud water, liquid (kg/kg) 
  real(r8), intent(in), dimension(pcols,pver)  :: ext_reff     ! effective radii [m]
  real(r8), intent(out), dimension(pcols,pver)  :: vfall_out     ! particle fall velocities [m/s]
  real(r8), intent(out), dimension(pcols,pverp) :: cld_pvel_out  ! vertical velocity of cloud liquid drops (Pa/s)

  ! local variables
  integer :: i,k
  real(r8) :: visc = 1.e-5   ! viscosity of air, CO2 (kg m / s)
  real(r8) :: molrad = 2.2e-10  ! effective molecular radius
  real(r8) :: vfall, a, b
  real(r8), dimension(pcols, pver) :: rho  

 ! parameters for Stokes velocity
  real (r8), parameter :: r40 =  40._r8              !  40 micron radius
  real (r8), parameter :: r400= 400._r8              ! 400 micron radius
  real (r8), parameter :: v400= 1.00_r8              ! fall velocity of 400 micron sphere (m/s)
  real (r8)            :: v40 ! = (2._r8/9._r8) * rhoh2o * gravit/eta * r40**2 * 1.e-12_r8
                                                     ! Stokes fall velocity of 40 micron sphere (m/s)
  real (r8)            :: vslope !  = (v400 - v40)/(r400 -r40) ! linear slope for large particles m/s/micron 

!----------------------------------------------------------------------- 
!------- initialize linear ramp variables for fall velocity ------------ 
!----------------------------------------------------------------------- 
  v40 = (2._r8/9._r8) * SHR_CONST_DENSITYCO2ICE * SHR_CONST_G/visc * r40**2 * 1.e-12_r8
  vslope = (v400 - v40)/(r400 -r40)

  ! initialization
  cld_pvel_out(:,:) = 0.0

  !
  ! Start Code
  !

  do k = 1,pver
    do i = 1,ncol
      rho(i,k) = ext_pmid(i,k) / (rair * ext_tmid(i,k))
      if (ext_cld(i,k) > 0._r8) then
!CAM way
!         if(ext_reff(i,k) .lt. 40.e-6) then  
!           vfall = 2._r8/9._r8 * SHR_CONST_DENSITYCO2ICE * SHR_CONST_G * ext_reff(i,k)**2 / visc   !!! * 1.e-12_r8  ! micons^2 -> m^2     
!            ! convert the fall speed to pressure units, but do not apply the traditional                     
!            ! negative convention for pvel.                                                                  
!         else
!             vfall = v40 + vslope * (ext_reff(i,k)-r40)      ! linear above 40 microns
!         endif
!LMD way
        b = (2._r8/9._r8) * SHR_CONST_DENSITYCO2ICE * SHR_CONST_G / visc
        a = 0.707*SHR_CONST_RGAS/(4.*SHR_CONST_PI*molrad**2*SHR_CONST_AVOGAD)
        vfall = b*ext_reff(i,k)**2 * (1. + 1.333 * (a*ext_tmid(i,k)/ext_pmid(i,k))/ext_reff(i,k)) 

        ! set fall velocity output array
        vfall_out(i,k) = vfall
        cld_pvel_out(i,k+1) = vfall * rho(i,k)*SHR_CONST_G        ! meters/sec to pascals/sec                        

      endif
    end do
  end do


end subroutine exo_cloud_sediment_vel


!=======================================================================

subroutine exo_cloud_sediment_tend(ncol, dtime, ext_tmid, ext_pmid, ext_pint, ext_pdel, ext_cld, &
                                   ext_cld_pvel, tcond, sed_tend, sf_tend)
!-----------------------------------------------------------------------  
! Purpose: Apply cloud particle sedimentation to condensate
! Adatpted from cld_sediment_tend in cld_pkg_sediment
!
  use pkg_cld_sediment, only: getflx

  ! input arguments
  integer,  intent(in)  :: ncol                                  ! number of colums to process
  real(r8), intent(in)  :: dtime                                 ! time step
  real(r8), intent(in), dimension(pcols,pver)   :: ext_tmid      ! temperature (K)
  real(r8), intent(in), dimension(pcols,pver)   :: ext_pmid      ! midpoint pressures (Pa)
  real(r8), intent(in), dimension(pcols,pverp)  :: ext_pint      ! interfaces pressure (Pa)
  real(r8), intent(in), dimension(pcols,pver)   :: ext_pdel      ! pressure diff across layer (Pa)
  real(r8), intent(in), dimension(pcols,pver)   :: ext_cld       ! cloud condensate (kg/kg)
  real(r8), intent(in), dimension(pcols,pverp)  :: ext_cld_pvel  ! vertical velocity of liquid drops  (Pa/s) 
  real(r8), intent(in), dimension(pcols, pver) :: tcond          ! atmosphere condensation temperature grid 
  ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1)) 
  real(r8), intent(out), dimension(pcols,pver) :: sed_tend       ! condensate tend [kg/kg/s]
  real(r8), intent(out), dimension(pcols) :: sf_tend             ! surface mass flux of condensate (snow/rain, kg/m2/s)

!  real(r8), intent(out) :: wvtend (pcols,pver)       ! water vapor tend
!  real(r8), intent(out), dimensions(pcols,pver) :: h_tend    ! heating rate

  ! Local variables
  real(r8), dimension(pcols,pverp) :: fxcld(pcols,pverp)     ! fluxes at the interfaces, liquid (positive = down)
  real(r8), dimension(pcols) :: cldab                        ! cloud in layer above
  real(r8) :: evapcld                                ! evaporation of cloud into environment
  real(r8) :: cldovrl                                ! cloud overlap factor
  real (r8), parameter :: mxsedfac   = 0.99_r8
  integer :: i,k, kstar

  fxcld(:ncol,:) = 0.0
  sed_tend(:ncol,:) = 0.0
  sf_tend(:ncol) = 0.0

  !
  ! Start Code
  ! 
!if(masterproc)write(*,*) "dtime sediment", dtime
  ! fluxes at interior points                                                                                      
  call getflx(ncol, ext_pint, ext_cld, ext_cld_pvel, dtime, fxcld)
!if(masterproc) write(*,*) "fxcld", fxcld

  ! calculate fluxes at boundaries                                                                                 
  do i = 1,ncol
    fxcld(i,1) = 0._r8
  ! surface flux by upstream scheme                                                                                
    fxcld(i,pverp) = ext_cld(i,pver) * ext_cld_pvel(i,pverp) * dtime
  end do

  ! filter out any negative fluxes from the getflx routine                                                           
  ! (typical fluxes are of order > 1e-3 when clouds are present)                                                   
  do k = 2,pver
    fxcld(:ncol,k) = max(0._r8, fxcld(:ncol,k))
  end do

  ! Limit the flux out of the bottom of each cell to the water content in each phase.
  ! Apply mxsedfac to prevent generating very small negative cloud water/ice
  ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
  ! ***Should we include the flux in the top, to allow for thin surface layers?
  ! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1,pver
       do i = 1,ncol
          fxcld(i,k+1) = min( fxcld(i,k+1), mxsedfac * ext_cld(i,k) * ext_pdel(i,k) )
!!          fxcld(i,k+1) = min( fxcld(i,k+1), mxsedfac * ext_cld(i,k) * ext_pdel(i,k) )

!!          fxcld(i,k+1) = min( fxcld(i,k+1), mxsedfac * ext_cld(i,k) * ext_pdel(i,k) + fxliq(i,k) )


!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))                               
!!$        fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))                               
!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )                              
!!$        fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )                              
       end do
    end do


! Now calculate the tendencies assuming that condensate evaporates when                                          
! it falls into environment, and does not when it falls into cloud.                                              
! All flux out of the layer comes from the cloudy part.                                                          
! Assume maximum overlap for stratiform clouds                                                                   
!  if cloud above < cloud,  all water falls into cloud below
!  if cloud above > cloud,  water split between cloud  and environment

! FOR CO2 ICE clouds.  As the ice particles sediment downward in the atmosphere, the must
! the release of potential energy by CO ice particles during their fall ??
! the heat consumption to warm the particles as the CO frost temperature increases with the local pressure)
! In theory this has already been accounted for else where in the model

   do k = 1,pver
!       cldab(:ncol) = 0._r8
       do i = 1,ncol

! cloud overlap cloud factor
!           cldovrl  = min( cloud(i,k) / (cldab(i)+.0001_r8), 1._r8 )
!          cldab(i) = cloud(i,k)
! evaporation into environment cause moistening and cooling
!   evapliq = fxliq(i,k) * (1._r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s) 
!          wvtend(i,k) = evapcld                         ! evaporation into environment (kg/kg/s)      
!          h_tend (i,k) = -latvap*evapliq -(latvap+latice)*evapice   ! evaporation (W/kg)                          

!  ah I need a co2 cloud fraction 

! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)                                        
!          cld_tend(i,k)  = (fxcld(i,k)*cldovrl - fxcld(i,k+1)) / (dtime * ext_pdel(i,k))

         sed_tend(i,k)  = (fxcld(i,k) - fxcld(i,k+1)) / (dtime * ext_pdel(i,k)) 


       end do
    end do

! convert flux out the bottom to mass units Pa -> kg/m2/s                                                        
    sf_tend(:ncol) = fxcld(:ncol,pverp) / (dtime*SHR_CONST_G)
!if(masterproc) write(*,*) "in sed_tend", sed_tend(1,1)
    return

end subroutine exo_cloud_sediment_tend

end module exo_condense_mod
