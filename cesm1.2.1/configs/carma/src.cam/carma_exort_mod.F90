!!  This module provides routines for coupling ExoRT with CARMA
!!  This module does not conform with the style of the CESM/CARMA 
!!  package, and its linkages to CAM3 RT, and rrtmg.  Frankly...
!!  I don't care.
!!
!!  @author  Eric T Wolf
!!  2018

module carma_exort_mod

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use carma_model_mod,  only: NELEM, NBIN
  use ppgrid,         only: pcols, pver
  use radgrid,          only: ntot_wavlnrng
  use physics_types, only: physics_state
 
#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none
  private
  save

  ! Public interfaces 
  public carma_exort_optics_init
  public carma_exort_get_mmr

  ! Public Data
  ! Optical constants for CARMA aerosols
  character(len=256), parameter :: dircarma = 'data/aerosol/'
  character(len=256), parameter :: cname = 'Khare.haze.fractal_n28_nbin40.nc'
  real(r8), public               :: qcarma(NELEM, NBIN, ntot_wavlnrng)   ! carma constituent mass extinction coefficient
  real(r8), public               :: wcarma(NELEM, NBIN, ntot_wavlnrng)   ! carma constituent single scattering albedo                             
  real(r8), public               :: gcarma(NELEM, NBIN, ntot_wavlnrng)   ! carma constituent asymmetry parameter    

  integer, parameter     :: nbins_c        = NBIN        ! Number of carma bins                       
  integer, parameter     :: nelems_c       = NELEM       ! Number of carma elements                 

  contains

!===============================================================================                                                                  
  subroutine carma_exort_optics_init
!-------------------------------------------------------------------------------                                                                  
! Purpose:                                                                                                                                        
!                                                                                                                                                 
! Read in CARMA binwise optical constants from initial file 'carmaoptics'                                                                         
! found on namelist.  Broadcast to all nodes. Adapted from subroutine                                                                             
! aer_optics_initialize found in aer_optics.F90.                                                                                                  
!                                                                                                                                                 
! NOTES: This subrountine is only set up to handle optics from 1 carma element.                                                                   
! I don't anticipate needing more than 1 element type.                                                                                            
!                                                                                                                                                 
! Author: Wolf, E.T. 3/15/2009                                                                                                                    
!--------------------------------------------------------------------------------                                                                 
! Notes:                                                                                                                                          
!                                                                                                                                                 
! Currently this subroutine is not set up to handle effects of relative humidity.                                                                 
! Haskenkopft said haze particles may absorb some water.                                                                                          
!-------------------------------------------------------------------------------                                                                  
!    use shr_kind_mod,     only: r8 => shr_kind_r8
    use spmd_utils,           only: masterproc
    use ioFileMod,        only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,  only: pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite,  &
              file_desc_t, var_desc_t
    use sys_rootdir
!--------------------------------------------------------------------------------                                                                 
    implicit none

    include 'netcdf.inc'

    integer :: nbnds, nbins, nelems,i 
    integer :: nc_id, bins_id,bnds_id, elems_id
    integer :: qcarma_id, wcarma_id, gcarma_id
    character(len=256)  :: locfn, filename 
    integer ::  q_id, w_id, g_id, ierr
    type(file_desc_t) :: ncid

!---------------------------------------------------------------------------------                                                                

    if ( masterproc ) then

      write (6, '(2x, a)') '-------------------------------------------------------'
      write (6, '(2x, a)') '======== initializing carma optical constants ========='
      write (6, '(2x, a)') '-------------------------------------------------------'

    end if   ! masterproc                                               
    
   ! Load K water cloud optics file                                                                                                              
    filename = trim(exort_rootdir)//trim(dircarma)//trim(cname)
    call getfil(filename, locfn, 0)
    call cam_pio_openfile(ncid, locfn, PIO_NOWRITE)
    ierr =  pio_inq_varid(ncid, 'qext_carma',   q_id)
    ierr =  pio_get_var(ncid, q_id, qcarma)
    ierr =  pio_inq_varid(ncid, 'ssa_carma',   w_id)
    ierr =  pio_get_var(ncid, w_id, wcarma)
    ierr =  pio_inq_varid(ncid, 'asm_carma',   g_id)
    ierr =  pio_get_var(ncid, g_id, gcarma)
    call pio_closefile(ncid)


#if ( defined SPMD )
    call mpibcast(qcarma, NELEM * NBIN * ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(wcarma, NELEM * NBIN * ntot_wavlnrng, mpir8, 0, mpicom)
    call mpibcast(gcarma, NELEM * NBIN * ntot_wavlnrng, mpir8, 0, mpicom)
#endif

  end subroutine carma_exort_optics_init



!========================================================================================                                                         
  subroutine carma_exort_get_mmr(state, carmamix)
!========================================================================================                                                         
! Purpose:                                                                                                                                        
!                                                                                                                                                 
! Return carma constituent mass mixing ratios, specifically for used for radiative                                                                
! transfer calculation.  Adapted from subroutine set_aerosol_from_prognostics located                                                             
! in aerosol_intr.F90, line 840.                                                                                                                  
!                                                                                                                                                 
! Author: Wolf, E.T. 3/30/2009                                                                                                                    
!--------------------------------------------------------------------------------                                                                 

    use constituents,   only: cnst_get_ind

    implicit none

    type(physics_state), intent(in)      :: state                   ! Physics state variables                                                     
    real(r8),            intent(inout)   :: carmamix(pcols, pver, nelems_c, nbins_c)  ! CARMA mmr                                          

!--------------------------------------------------------------------------------                                                                 
!   local variables                                                                                                                               

    character(len=8)    ::   carmaname, shortname
    integer i,j,m,ncol,cnst_idxCARMAfirst
    
!--------------------------------------------------------------------------------                                                                 
    !Kludge
    shortname = 'HAZE'
    ncol = state%ncol

    do i = 1, NELEM
       write(carmaname, '(A, I2.2)') trim(shortname), 1
       call cnst_get_ind (carmaname,cnst_idxCARMAfirst)
       do j = 1, NBIN
          carmamix(:ncol,:,i,j) = state%q(:ncol,:,cnst_idxCARMAfirst+j-1)
       end do
    end do

  end subroutine carma_exort_get_mmr



end module carma_exort_mod
