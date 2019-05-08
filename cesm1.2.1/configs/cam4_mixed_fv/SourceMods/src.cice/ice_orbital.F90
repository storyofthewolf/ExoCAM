!=======================================================================
!BOP
!
! !MODULE: ice_orbital - computes orbital parameters for solar zenith angle
!
! !DESCRIPTION:
!
! Orbital parameters computed from date
!
! !REVISION HISTORY:
!  SVN:$Id: ice_orbital.F90 37 2006-11-29 18:06:44Z eclare $
!
! author:  Bruce P. Briegleb, NCAR 
!
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2006: BPB  26 December Modified to compute diurnal mean coszen
!
! !INTERFACE:
!
      module ice_orbital
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_domain_size
      use ice_constants
      use shr_orb_mod
      use exoplanet_mod, only: exo_ndays
! 2 Jan07 BPB
      use ice_diagnostics
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: iyear_AD  ! Year to calculate orbit for
 
      real(kind=dbl_kind) :: eccen  !Earth's orbital eccentricity
      real(kind=dbl_kind) :: obliqr !Earth's obliquity in radians
      real(kind=dbl_kind) :: lambm0 !Mean longitude of perihelion at the
                                    !vernal equinox (radians)
      real(kind=dbl_kind) :: mvelpp !Earth's moving vernal equinox longitude
                                    !of perihelion + pi (radians)
      real(kind=dbl_kind) :: obliq  ! obliquity in degrees
      real(kind=dbl_kind) :: mvelp  ! moving vernal equinox long
      real(kind=dbl_kind) :: delta  ! solar declination angle in radians
      real(kind=dbl_kind) :: eccf   ! earth orbit eccentricity factor

      logical(kind=log_kind) :: log_print ! Flags print of status/error
 
!=======================================================================
 
      contains
 
!=======================================================================
!BOP
!
! !IROUTINE: init_orbit - initialize orbital parameters
!
! !INTERFACE:
!
      subroutine init_orbit
!
! !DESCRIPTION:
!
! Uses share routines to compute orbital parameters
! for the specified date.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
!
! !USES: none
!
! !INPUT/OUTPUT PARAMETERS: none
!
!EOP
!
      iyear_AD  = 1950
      log_print = .false.   ! if true, write out orbital parameters
!
      call shr_orb_params( iyear_AD , eccen  , obliq , mvelp     , &
                           obliqr   , lambm0 , mvelpp, log_print )
 
      end subroutine init_orbit
 
!=======================================================================
!BOP
!
! !IROUTINE: compute_coszen - computes cosine solar zenith angle
!
! !INTERFACE:
!
      subroutine compute_coszen (nx_block, ny_block, &
                                 icells,             &
                                 indxi,    indxj,    &
                                 tlat,     tlon,     &
                                 coszen,   dt)
!
! !DESCRIPTION:
!
! Uses orbital and lat/lon info to compute cosine solar zenith angle
! for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 
!
! !USES:
!
      use ice_calendar, only: yday, sec, secday, days_per_year, &
                              calendar_type, nextsw_cday, nyr
! 
! !INPUT/OUTPUT PARAMETERS: 
! 
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj          ! indices for ice-covered cells
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         tlat, tlon          ! latitude and longitude (radians)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         coszen              ! cosine solar zenith angle 
                             ! negative for sun below horizon
 
      real (kind=dbl_kind), intent(in) :: &
         dt                  ! thermodynamic time step
!
!EOP
!
      real (kind=dbl_kind) :: ydayp1 ! day of year plus one time step
 
      integer (kind=int_kind) :: &
         i   , & ! domain longitude index
         j   , & ! domain latitude index
         ij      ! horizontal index, combines i and j loops
 
 !WOLF
       real (kind=dbl_kind) :: frac_day, day_in_year
       integer (kind=int_kind) :: ncol

! Solar declination for next time step
 
#ifdef CCSMCOUPLED
      if (calendar_type == "GREGORIAN") then
         ydayp1 = min(nextsw_cday, real(days_per_year,kind=dbl_kind))
      else
         ydayp1 = nextsw_cday
      endif
#else
      ydayp1 = yday + sec/secday
#endif
 
      if (ydayp1 > -0.5_dbl_kind) then

      call shr_orb_decl(ydayp1, eccen, mvelpp, lambm0, &
                        obliqr, delta, eccf)

      coszen(:,:) = c0  ! sun at horizon

      !Wolf  
      ydayp1 = ydayp1 + 365.*(nyr-1)
      ydayp1=ydayp1-1
      frac_day = ydayp1/exo_ndays - FLOOR(ydayp1 / exo_ndays)
      ! day_in_year = yday - (365)*FLOOR(yday/(365))  ! scaling for different length years not operable
      !write(*,*) "coszen, frac_day, yday", frac_day, ydayp1
      ! \Wolf

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
        !WOLF 
        if (do_exo_synchronous) then    !synchronous rotation
          do ij = 1, icells
             i = indxi(ij)
             j = indxj(ij)
             coszen(i,j) = sin(tlat(i,j))*sin(delta) - &
                           cos(tlat(i,j))*cos(delta)   &
                           *cos(tlon(i,j))
          enddo
        else
          do ij = 1, icells
             i = indxi(ij)
             j = indxj(ij)
             coszen(i,j) = sin(tlat(i,j))*sin(delta) - &
                           cos(tlat(i,j))*cos(delta)   &
                           *cos(frac_day*c2*pi + tlon(i,j))
          enddo
        endif
      

      endif

      end subroutine compute_coszen
 
!=======================================================================
 
      end module ice_orbital
 
!=======================================================================
