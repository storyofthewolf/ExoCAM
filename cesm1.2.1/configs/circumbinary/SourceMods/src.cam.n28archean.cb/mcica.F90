
module mcica

!-----------------------------------------------------------------------
!  Adapted from AER RRTMG, Monte Carlo Independent Column Approximation
!  Modified For use with Early Earth Terrestrial Model  E.T.Wolf, Sept, 2010
!-----------------------------------------------------------------------
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.
! Two options are possible:      
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0)
!-----------------------------------------------------------------------



  use shr_kind_mod,    only: r8 => shr_kind_r8
!  use abortutils,      only: endrun
  use shr_const_mod,   only: SHR_CONST_G

  use radgrid
  use ppgrid  

  implicit none
  private
  save

!------------------------------------------------------------------------
!
! Public interfaces
!

  public :: mcica_subcol

!============================================================================
  contains
!============================================================================

!============================================================================
!
! Public subroutines
!
!============================================================================

  subroutine mcica_subcol(icld, permuteseed, pmid, &
                       cldfrac, ciwp, clwp, tauci, ssaci, asmci, taucl,  &
                       ssacl, asmcl, cldf_mcica, ciwp_mcica, clwp_mcica, &
                       tauci_mcica, ssaci_mcica, asmci_mcica, &
                       taucl_mcica, ssacl_mcica, asmcl_mcica)


!------------------------------------------------------------------------
!
! Purpose:
!
!------------------------------------------------------------------------
  implicit none
!------------------------------------------------------------------------
!
! Input Arguments
!
    integer, intent(in) :: icld                   ! clear/cloud, cloud overlap flag
    integer, intent(in) :: permuteseed            ! if the cloud generator is called multiple times,
                                                        ! permute the seed between each call;
                                                        ! between calls for LW and SW, recommended
                                                        ! permuteseed differs by 'ngpt'

    real(r8), intent(in), dimension(pverp) :: pmid                     ! mid layer pressures (mb) 
    real(r8), intent(in), dimension(pverp) :: cldfrac                  ! layer cloud fraction
    real(r8), intent(in), dimension(pverp) :: ciwp                     ! cloud ice water path
    real(r8), intent(in), dimension(pverp) :: clwp                     ! cloud liquid water path
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: tauci      ! ice cloud optical depth
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: ssaci      ! ice cloud single scattering albedo (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: asmci      ! ice cloud asymmetry parameter (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: taucl      ! liquid cloud optical depth
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: ssacl      ! liquid cloud single scattering albedo (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: asmcl      ! liquid cloud asymmetry parameter (non-delta scaled)

    real(r8), intent(out), dimension(ntot_gpt,pverp) :: cldf_mcica         ! cloud fraction [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ciwp_mcica         ! cloud ice water path [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: clwp_mcica         ! cloud liquid water path [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: tauci_mcica         ! ice cloud optical depth [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ssaci_mcica         ! ice cloud single scattering albedo [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: asmci_mcica         ! ice cloud asymmetry parameter [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: taucl_mcica         ! liquid cloud optical depth [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ssacl_mcica         ! liquid cloud single scattering albedo [mcica]
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: asmcl_mcica         ! liquid cloud asymmetry parameter [mcica]

!------------------------------------------------------------------------
!
! Local Variables
!

    ! Stochastic cloud generator variables [mcica]
    integer :: km                               ! loop indices
    integer :: im
    integer :: nm

    ! Return if clear sky; or stop if icld out of range
    if (icld.eq.0) return
    if (icld.lt.0.or.icld.gt.3) then 
      !call endrun('MCICA_SUBCOL: INVALID ICLD')
      write(*,*) "MCICA error, invalid icld"
    endif 

    !  Generate the stochastic subcolumns of cloud optical properties for the shortwave;
    call generate_stochastic_clouds (icld, pmid, cldfrac, clwp, ciwp, tauci, &
                                     ssaci, asmci, taucl, ssacl, asmcl, cldf_mcica, clwp_mcica, ciwp_mcica, &
                                     tauci_mcica, ssaci_mcica, asmci_mcica, taucl_mcica, ssacl_mcica, asmcl_mcica, &
                                     permuteseed )

  end subroutine mcica_subcol

!================================================================================================


  subroutine generate_stochastic_clouds(icld, pmid, cldf, clwp, ciwp, tauci, & 
                                   ssaci, asmci, taucl, ssacl, asmcl, cld_stoch, clwp_stoch, ciwp_stoch, tauci_stoch, &
                                   ssaci_stoch, asmci_stoch, taucl_stoch, ssacl_stoch, asmcl_stoch, changeSeed) 

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Original code: Based on Raisanen et al., QJRMS, 2004.
  !
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irnd' below . Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------
  
    ! from mcica_random_number.f90
    use mcica_random_numbers
    use MersenneTwister,   only: randomNumberSequence, new_RandomNumberSequence, getRandomReal

    implicit none

    type(randomNumberSequence) :: randomNumbers

!------------------------------------------------------------------------
!
! Input Arguments
!
    integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
    integer, optional, intent(in) :: changeSeed     ! allows permuting seed

    ! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
    real(r8), intent(in), dimension(pverp) :: pmid             ! mid layer pressure (Pa)
    real(r8), intent(in), dimension(pverp) :: cldf             ! cloud fraction 
    real(r8), intent(in), dimension(pverp) :: clwp             ! cloud liquid water path (g/m2)
    real(r8), intent(in), dimension(pverp) :: ciwp             ! cloud ice water path (g/m2)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: tauci     ! ice cloud optical depth (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: ssaci     ! ice cloud single scattering albedo (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: asmci     ! ice cloud asymmetry parameter (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: taucl     ! liquid cloud optical depth (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: ssacl     ! liquid cloud single scattering albedo (non-delta scaled)
    real(r8), intent(in), dimension(ntot_wavlnrng,pverp) :: asmcl     ! liquid cloud asymmetry parameter (non-delta scaled)

    real(r8), intent(out), dimension(ntot_gpt,pverp) :: cld_stoch     ! subcolumn cloud fraction 
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: clwp_stoch    ! subcolumn cloud liquid water path
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ciwp_stoch    ! subcolumn cloud ice water path
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: tauci_stoch    ! subcolumn ice cloud optical depth
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ssaci_stoch    ! subcolumn ice cloud single scattering albedo
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: asmci_stoch    ! subcolumn ice cloud asymmetry parameter
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: taucl_stoch    ! subcolumn liquid cloud optical depth
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: ssacl_stoch    ! subcolumn liquid cloud single scattering albedo
    real(r8), intent(out), dimension(ntot_gpt,pverp) :: asmcl_stoch    ! subcolumn liquid cloud asymmetry parameter

!------------------------------------------------------------------------
!
! Local Variables
!

    ! Set overlap
    integer :: overlap                     ! 1 = random overlap, 2 = maximum/random, 3 = maximum overlap, 

    ! Constants (min value for cloud fraction and cloud water and ice)
    real(r8), parameter :: cldmin = 1.0e-80_r8     ! min cloud fraction

    ! Variables related to random number and seed 
    integer :: irnd                         ! flag for random number generator, 0 = kissvec, 1 = Mersenne Twister

    ! for Kissvec random numbers
    real(r8), dimension(ntot_gpt,pverp) :: CDF       ! random numbers
    integer, dimension(1) :: seed1                     ! seed to create random number
    integer, dimension(1) :: seed2  
    integer, dimension(1) :: seed3  
    integer, dimension(1) :: seed4
    real(r8), dimension(1) :: rand_num                 ! random number (kissvec)

    ! for Mersenns Twister
    integer :: iseed                     ! seed to create random number (Mersenne Twister)
    real(r8) :: rand_num_mt              ! random number (Mersenne Twister)

    ! Flag to identify cloud fraction in subcolumns
    logical,  dimension(ntot_gpt, pverp) :: isCloudy   ! flag that says whether a gridbox is cloudy

    ! Indices
    integer :: ilev
    integer :: isubcol
    integer :: i
    integer :: n
    integer :: ngbm             

    real(r8), dimension(pverp) :: cldf_wrk             ! cloud fraction working array

!------------------------------------------------------------------------
!
! Start Code
!
  
    !Set randum number generator to use (0 = kissvec; 1 = mersennetwister)
    irnd = 1

    ! Pass input cloud overlap setting to local variable
    overlap = icld

    
    ! Ensure that cloud fractions are in bounds 
    cldf_wrk(:) = cldf(:)
    where (cldf(:) < cldmin)
      cldf_wrk(:) = 0._r8
    end where

    ! ----- Create seed  --------
   
    ! Advance randum number generator by changeseed values
    if (irnd.eq.0) then   

      ! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
      ! Must use pmid from bottom four layers.      
  
      if (pmid(pver) .lt. pmid(pver-1)) then
        !call endrun('MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.')
        write(*,*) "CRASH: MCICA seed generator"
      endif

      seed1(1) = (pmid(pverp) - int(pmid(pverp)))  * 1000000000
      seed2(1) = (pmid(pverp-1) - int(pmid(pverp-1)))  * 1000000000
      seed3(1) = (pmid(pverp-2) - int(pmid(pverp-2)))  * 1000000000
      seed4(1) = (pmid(pverp-3) - int(pmid(pverp-3)))  * 1000000000
 
      do i=1,changeSeed
        call kissvec(seed1, seed2, seed3, seed4, rand_num)
      enddo
    elseif (irnd.eq.1) then
      randomNumbers = new_RandomNumberSequence(seed = changeSeed)
    endif 

    ! ------ Apply overlap assumption --------

    ! generate the random numbers  

    select case (overlap)

      case(1) 
      ! Random overlap
      ! i) pick a random value at every level
  
        if (irnd.eq.0) then 
          do isubcol = 1,ntot_gpt
            do ilev = 1,pverp
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF(isubcol,ilev) = rand_num(1)
            enddo
          enddo
        elseif (irnd.eq.1) then
          do isubcol = 1, ntot_gpt
            do ilev = 1, pverp
              rand_num_mt = getRandomReal(randomNumbers)
              CDF(isubcol,ilev) = rand_num_mt
            enddo
          enddo
        endif

      case(2) 
      ! Maximum-Random overlap
      ! i) pick  a random number for top layer.
      ! ii) walk down the column: 
      !    - if the layer above is cloudy, we use the same random number than in the layer above
      !    - if the layer above is clear, we use a new random number 

        if (irnd.eq.0) then 
          do isubcol = 1,ntot_gpt
            do ilev = 1,pverp
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF(isubcol,ilev) = rand_num(1)
            enddo
          enddo
        elseif (irnd.eq.1) then
          do isubcol = 1, ntot_gpt
            do ilev = 1, pverp
              rand_num_mt = getRandomReal(randomNumbers)
              CDF(isubcol,ilev) = rand_num_mt
            enddo
          enddo
        endif

        do ilev = 2,pverp
          do isubcol = 1, ntot_gpt
            if (CDF(isubcol,ilev-1) > 1._r8 - cldf_wrk(ilev-1) ) then
              CDF(isubcol,ilev) = CDF(isubcol,ilev-1)       !cloudy, use same random number           
             else
              CDF(isubcol,ilev) = CDF(isubcol,ilev) * (1._r8 - cldf_wrk(ilev-1))  ! clear, use different
            endif    
          enddo
        enddo

      case(3) 
      ! Maximum overlap
      ! i) pick same random numebr at every level  

        if (irnd.eq.0) then 
          do isubcol = 1,ntot_gpt
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
            do ilev = 1,pverp
               CDF(isubcol,ilev) = rand_num(1)
             enddo
          enddo
        elseif (irnd.eq.1) then
          do isubcol = 1, ntot_gpt
            rand_num_mt = getRandomReal(randomNumbers)
            do ilev = 1, pverp
              CDF(isubcol,ilev) = rand_num_mt
            enddo       
          enddo
        endif

    end select
 
    ! -- generate subcolumns for homogeneous clouds -----
    do ilev = 1, pverp
      isCloudy(:,ilev) = (CDF(:,ilev) >= 1._r8 - spread(cldf_wrk(ilev), dim=1, nCopies=ntot_gpt) )
    enddo

    ! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
    ! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0

    do ilev = 1, pverp
      do isubcol = 1, ntot_gpt
        if (iscloudy(isubcol,ilev) ) then
          cld_stoch(isubcol,ilev) = 1._r8
        else
          cld_stoch(isubcol,ilev) = 0._r8
        endif
      enddo
    enddo

    ! where there is a cloud, set the subcolumn cloud properties;
    ! Incoming clwp, ciwp and tauc should be in-cloud quantites and not grid-averaged quantities

    do ilev = 1, pverp
      do isubcol = 1, ntot_gpt
        if ( iscloudy(isubcol,ilev) .and. (cldf(ilev) > 0._r8) ) then
          clwp_stoch(isubcol,ilev) = clwp(ilev)
          ciwp_stoch(isubcol,ilev) = ciwp(ilev)
        else
          clwp_stoch(isubcol,ilev) = 0._r8
          ciwp_stoch(isubcol,ilev) = 0._r8
        end if        
      enddo
    enddo
    do ilev = 1,pverp
      do isubcol = 1, ntot_gpt
        if ( iscloudy(isubcol,ilev) .and. (cldf(ilev) > 0._r8) ) then
          n = ngb(isubcol) 
          tauci_stoch(isubcol,ilev) = tauci(n,ilev)
          ssaci_stoch(isubcol,ilev) = ssaci(n,ilev)
          asmci_stoch(isubcol,ilev) = asmci(n,ilev)
          taucl_stoch(isubcol,ilev) = taucl(n,ilev)
          ssacl_stoch(isubcol,ilev) = ssacl(n,ilev)
          asmcl_stoch(isubcol,ilev) = asmcl(n,ilev)
        else
          tauci_stoch(isubcol,ilev) = 0._r8
          ssaci_stoch(isubcol,ilev) = 1._r8
          asmci_stoch(isubcol,ilev) = 1._r8
          taucl_stoch(isubcol,ilev) = 0._r8
          ssacl_stoch(isubcol,ilev) = 1._r8
          asmcl_stoch(isubcol,ilev) = 1._r8
        endif
      enddo
    enddo

  end subroutine generate_stochastic_clouds
  
!================================================================================================

  subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)

!------------------------------------------------------------------------
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
  
    implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!

    real(r8), intent(inout), dimension(:)  :: ran_arr
    integer, intent(inout), dimension(:) :: seed1
    integer, intent(inout), dimension(:) :: seed2
    integer, intent(inout), dimension(:) :: seed3
    integer, intent(inout), dimension(:) :: seed4
    integer :: i
    integer :: sz
    integer :: kiss
    integer :: m
    integer :: k
    integer :: n

!------------------------------------------------------------------------
!
! Input Arguments
!

    ! inline function 
    m(k, n) = ieor (k, ishft (k, n) )
    sz = size(ran_arr)
    do i = 1, sz
      seed1(i) = 69069 * seed1(i) + 1327217885
      seed2(i) = m (m (m (seed2(i), 13), - 17), 5)
      seed3(i) = 18000 * iand (seed3(i), 65535) + ishft (seed3(i), - 16)
      seed4(i) = 30903 * iand (seed4(i), 65535) + ishft (seed4(i), - 16)
      kiss = seed1(i) + seed2(i) + ishft (seed3(i), 16) + seed4(i)
      ran_arr(i) = kiss*2.328306e-10_r8 + 0.5_r8
    enddo
    
  end subroutine kissvec

end module mcica