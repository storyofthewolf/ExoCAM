module clmtypeInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
! ETW: Initialize them to zero.  Why the F would you initialize them NaN?
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use clmtype
  use clm_varpar    , only : maxpatch_pft, nlevsno, nlevgrnd, numrad, nlevlak, &
                             numpft, ndst, nlevurb, nlevsoi
  use clm_varctl  , only : use_c13, use_cn, use_cndv, use_crop
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
! Modified by Colette L. Heald (05/06) for VOC emission factors
! 3/17/08 David Lawrence, changed nlevsoi to nlevgrnd where appropriate
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_ecophys_constants
  private :: init_pft_DGVMecophys_constants
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
  private :: init_pft_pdgvstate_type
  private :: init_pft_vstate_type
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_cstate_type
  private :: init_pft_nstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_nflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_pft_depvd_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_nstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_column_cflux_type
  private :: init_column_nflux_type
  private :: init_landunit_pstate_type
  private :: init_landunit_eflux_type
  private :: init_gridcell_pstate_type
  private :: init_gridcell_efstate_type
  private :: init_gridcell_wflux_type
!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:
  subroutine initClmtype()
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
    use abortutils, only : endrun
    use decompMod , only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    character(len=32), parameter :: subname = "initClmtype"
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    call init_pft_type     (begp, endp, pft)
    call init_column_type  (begc, endc, col)
    call init_landunit_type(begl, endl, lun)
    call init_gridcell_type(begg, endg, grc)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    ! pft DGVM-specific ecophysiological constants

    if (use_cndv) then
       call init_pft_DGVMecophys_constants()
    end if

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, pebal)
    call init_energy_balance_type(begc, endc, cebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, pwbal)
    call init_water_balance_type(begc, endc, cwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, pcbal)
    call init_carbon_balance_type(begc, endc, ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, pnbal)
    call init_nitrogen_balance_type(begc, endc, cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, pps)
    call init_pft_pstate_type(begc, endc, pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, pepv)

    ! pft DGVM state variables at pft level 

    if (use_cndv) then
       call init_pft_pdgvstate_type(begp, endp, pdgvs)
    end if
    call init_pft_vstate_type(begp, endp, pvs)

    ! pft energy state variables at the pft level

    call init_pft_estate_type(begp, endp, pes)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, pws)
    call init_pft_wstate_type(begc, endc, pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, pcs)
    call init_pft_cstate_type(begc, endc, pcs_a)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_pft_cstate_type(begp, endp, pc13s)
       call init_pft_cstate_type(begc, endc, pc13s_a)
       if (use_crop) then
          call endrun( trim(subname)//" ERROR:: CROP and C13 can NOT be on at the same time" )
       end if
    endif

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, pns)
    call init_pft_nstate_type(begc, endc, pns_a)

    ! pft energy flux variables at pft level 

    call init_pft_eflux_type(begp, endp, pef)

    ! pft momentum flux variables at pft level 

    call init_pft_mflux_type(begp, endp, pmf)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, pwf)
    call init_pft_wflux_type(begc, endc, pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, pcf)
    call init_pft_cflux_type(begc, endc, pcf_a)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_pft_cflux_type(begp, endp, pc13f)
       call init_pft_cflux_type(begc, endc, pc13f_a)
    endif

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, pnf)
    call init_pft_nflux_type(begc, endc, pnf_a)

    ! pft VOC flux variables at pft level

    call init_pft_vflux_type(begp, endp, pvf)

    ! gridcell VOC emission factors (heald, 05/06)

    call init_gridcell_efstate_type(begg, endg, gve)

    ! pft dust flux variables at pft level 

    call init_pft_dflux_type(begp, endp, pdf)

    ! pft dry dep velocity variables at pft level 

    call init_pft_depvd_type(begp, endp, pdd)

    ! column physical state variables at column level 

    call init_column_pstate_type(begc, endc, cps)

    ! column energy state variables at column level 


    call init_column_estate_type(begc, endc, ces)

    ! column water state variables at column level 

    call init_column_wstate_type(begc, endc, cws)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, ccs)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_column_cstate_type(begc, endc, cc13s)
    endif

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, cns)

    ! column energy flux variables at column level 

    call init_column_eflux_type(begc, endc, cef)

    ! column water flux variables at column level 

    call init_column_wflux_type(begc, endc, cwf)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, ccf)
    if (use_c13) then
       ! 4/14/05: PET
       ! Adding isotope code
       call init_column_cflux_type(begc, endc, cc13f)
    endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, lps)

    ! land unit energy flux variables 

    call init_landunit_eflux_type(begl, endl, lef)

    ! gridcell DGVM variables

    if (use_cndv) then
       call init_gridcell_dgvstate_type(begg, endg, gdgvs)
    end if

    ! gridcell physical state variables


    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, gwf)

    ! gridcell: energy flux variables

    call init_gridcell_eflux_type(begg, endg, gef)

    ! gridcell: water state variables

    call init_gridcell_wstate_type(begg, endg, gws)

    ! gridcell: energy state variables

    call init_gridcell_estate_type(begg, endg, ges)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pft%gridcell(beg:end),&
             pft%wtgcell(beg:end), &
             pft%landunit(beg:end),&
             pft%wtlunit(beg:end), &
             pft%column(beg:end),  &
             pft%wtcol(beg:end),   &
             pft%itype(beg:end),   &
             pft%mxy(beg:end))

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(col%gridcell(beg:end),&
            col%wtgcell(beg:end), &
            col%landunit(beg:end),&
            col%wtlunit(beg:end), &
            col%pfti(beg:end),    &
            col%pftf(beg:end),    &
            col%npfts(beg:end),   &
            col%itype(beg:end))

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(lun%gridcell(beg:end), &
            lun%wtgcell(beg:end),  &
            lun%coli(beg:end),     &
            lun%colf(beg:end),     &
            lun%ncolumns(beg:end), &
            lun%pfti(beg:end),     &
            lun%pftf(beg:end),     &
            lun%npfts(beg:end),    & 
            lun%itype(beg:end),    &
            lun%ifspecial(beg:end),& 
            lun%lakpoi(beg:end),   &
            lun%urbpoi(beg:end),   &
            lun%glcmecpoi(beg:end))

   ! These should be moved to landunit physical state
   allocate(lun%canyon_hwr(beg:end),   &
            lun%wtroad_perv(beg:end),  &
            lun%ht_roof(beg:end),      &
            lun%wtlunit_roof(beg:end), &
            lun%wind_hgt_canyon(beg:end),&
            lun%z_0_town(beg:end),     &
            lun%z_d_town(beg:end))

   ival = 0.0
   lun%canyon_hwr(beg:end)     = ival !nan
   lun%wtroad_perv(beg:end)    = ival !nan
   lun%ht_roof(beg:end)        = ival !nan
   lun%wtlunit_roof(beg:end)   = ival !
   lun%wind_hgt_canyon(beg:end)= ival
   lun%z_0_town(beg:end)       = ival
   lun%z_d_town(beg:end)       = ival

   lun%glcmecpoi(beg:end) = .false.

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(grc%luni(beg:end),      &
            grc%lunf(beg:end),      &
            grc%nlandunits(beg:end),&
            grc%coli(beg:end),      &
            grc%colf(beg:end),      &
            grc%ncolumns(beg:end),  &
            grc%pfti(beg:end),      &
            grc%pftf(beg:end),      &
            grc%npfts(beg:end),     &
            grc%gindex(beg:end),    &
            grc%area(beg:end),      &
            grc%lat(beg:end),       &
            grc%lon(beg:end),       &
            grc%latdeg(beg:end),    &
            grc%londeg(beg:end),    &
            grc%gindex_a(beg:end),  &
            grc%lat_a(beg:end),     &
            grc%lon_a(beg:end),     &
            grc%latdeg_a(beg:end),  &
            grc%londeg_a(beg:end),  &
            grc%gris_mask(beg:end), &
            grc%gris_area(beg:end), &
            grc%aais_mask(beg:end), &
            grc%aais_area(beg:end))

   ival = 0.0
   grc%gris_mask(beg:end) = ival
   grc%gris_area(beg:end) = ival
   grc%aais_mask(beg:end) = ival
   grc%aais_area(beg:end) = ival

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(beg, end, ebal)
!
! !DESCRIPTION:
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    allocate(ebal%errseb(beg:end))
    allocate(ebal%errsol(beg:end))
    allocate(ebal%errlon(beg:end))

    ival = 0.0
    ebal%errsoi(beg:end) = ival
    ebal%errseb(beg:end) = ival
    ebal%errsol(beg:end) = ival
    ebal%errlon(beg:end) = ival

  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(beg, end, wbal)
!
! !DESCRIPTION:
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    allocate(wbal%endwb(beg:end))
    allocate(wbal%errh2o(beg:end))

    ival = 0.0
    wbal%begwb(beg:end) = ival
    wbal%endwb(beg:end) = ival
    wbal%errh2o(beg:end) = ival

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_balance_type
!
! !INTERFACE:
  subroutine init_carbon_balance_type(beg, end, cbal)
!
! !DESCRIPTION:
! Initialize carbon balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    allocate(cbal%endcb(beg:end))
    allocate(cbal%errcb(beg:end))

    ival = 0.0
    cbal%begcb(beg:end) = ival
    cbal%endcb(beg:end) = ival
    cbal%errcb(beg:end) = ival

  end subroutine init_carbon_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nitrogen_balance_type
!
! !INTERFACE:
  subroutine init_nitrogen_balance_type(beg, end, nbal)
!
! !DESCRIPTION:
! Initialize nitrogen balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    allocate(nbal%endnb(beg:end))
    allocate(nbal%errnb(beg:end))

    ival = 0.0
    nbal%begnb(beg:end) = ival
    nbal%endnb(beg:end) = ival
    nbal%errnb(beg:end) = ival

  end subroutine init_nitrogen_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%tree(0:numpft))
    allocate(pftcon%smpso(0:numpft)) 
    allocate(pftcon%smpsc(0:numpft)) 
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%sla(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%graincn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%resist(0:numpft))
    allocate(pftcon%dwood(0:numpft))

    ival = 0.0
    pftcon%noveg(:) = huge(1)
    pftcon%tree(:) = huge(1)
    pftcon%smpso(:) = ival
    pftcon%smpsc(:) = ival
    pftcon%fnitr(:) = ival
    pftcon%foln(:) = ival
    pftcon%dleaf(:) = ival
    pftcon%c3psn(:) = ival
    pftcon%mp(:) = ival
    pftcon%qe25(:) = ival
    pftcon%xl(:) = ival
    pftcon%rhol(:,:numrad) = ival
    pftcon%rhos(:,:numrad) = ival
    pftcon%taul(:,:numrad) = ival
    pftcon%taus(:,:numrad) = ival
    pftcon%z0mr(:) = ival
    pftcon%displar(:) = ival
    pftcon%roota_par(:) = ival
    pftcon%rootb_par(:) = ival
    pftcon%sla(:) = ival
    pftcon%slatop(:) = ival
    pftcon%dsladlai(:) = ival
    pftcon%leafcn(:) = ival
    pftcon%flnr(:) = ival
    pftcon%woody(:) = ival
    pftcon%lflitcn(:) = ival
    pftcon%frootcn(:) = ival
    pftcon%livewdcn(:) = ival
    pftcon%deadwdcn(:) = ival
    pftcon%graincn(:) = ival
    pftcon%froot_leaf(:) = ival
    pftcon%stem_leaf(:) = ival
    pftcon%croot_stem(:) = ival
    pftcon%flivewd(:) = ival
    pftcon%fcur(:) = ival
    pftcon%lf_flab(:) = ival
    pftcon%lf_fcel(:) = ival
    pftcon%lf_flig(:) = ival
    pftcon%fr_flab(:) = ival
    pftcon%fr_fcel(:) = ival
    pftcon%fr_flig(:) = ival
    pftcon%leaf_long(:) = ival
    pftcon%evergreen(:) = ival
    pftcon%stress_decid(:) = ival
    pftcon%season_decid(:) = ival
    pftcon%resist(:) = ival
    pftcon%dwood(:) = ival

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    ival = 0.0
    dgv_pftcon%crownarea_max(:) = ival
    dgv_pftcon%tcmin(:) = ival
    dgv_pftcon%tcmax(:) = ival
    dgv_pftcon%gddmin(:) = ival
    dgv_pftcon%twmax(:) = ival
    dgv_pftcon%reinickerp(:) = ival
    dgv_pftcon%allom1(:) = ival
    dgv_pftcon%allom2(:) = ival
    dgv_pftcon%allom3(:) = ival

  end subroutine init_pft_DGVMecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    ival = 0.0
    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%emv(beg:end))
    allocate(pps%z0mv(beg:end))
    allocate(pps%z0hv(beg:end))
    allocate(pps%z0qv(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevgrnd))
    allocate(pps%rootr(beg:end,1:nlevgrnd))
    allocate(pps%rresis(beg:end,1:nlevgrnd))
    allocate(pps%dewmx(beg:end))
    allocate(pps%rssun(beg:end))
    allocate(pps%rssha(beg:end))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%btran(beg:end))
    allocate(pps%fsun(beg:end))
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%fwet(beg:end))
    allocate(pps%fdry(beg:end))
    allocate(pps%dt_veg(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%z0m(beg:end))
    allocate(pps%displa(beg:end))
    allocate(pps%albd(beg:end,1:numrad))
    allocate(pps%albi(beg:end,1:numrad))
    allocate(pps%fabd(beg:end,1:numrad))
    allocate(pps%fabi(beg:end,1:numrad))
    allocate(pps%ftdd(beg:end,1:numrad))
    allocate(pps%ftid(beg:end,1:numrad))
    allocate(pps%ftii(beg:end,1:numrad))
    allocate(pps%u10(beg:end))
    allocate(pps%u10_clm(beg:end))
    allocate(pps%va(beg:end))
    allocate(pps%fv(beg:end))
    allocate(pps%ram1(beg:end))
    if ( crop_prog )then
       allocate(pps%hdidx(beg:end))
       allocate(pps%cumvd(beg:end))
       allocate(pps%htmx(beg:end))
       allocate(pps%vf(beg:end))
       allocate(pps%gddmaturity(beg:end))
       allocate(pps%gdd0(beg:end))
       allocate(pps%gdd8(beg:end))
       allocate(pps%gdd10(beg:end))
       allocate(pps%gdd020(beg:end))
       allocate(pps%gdd820(beg:end))
       allocate(pps%gdd1020(beg:end))
       allocate(pps%gddplant(beg:end))
       allocate(pps%gddtsoi(beg:end))
       allocate(pps%huileaf(beg:end))
       allocate(pps%huigrain(beg:end))
       allocate(pps%aleafi(beg:end))
       allocate(pps%astemi(beg:end))
       allocate(pps%aleaf(beg:end))
       allocate(pps%astem(beg:end))
       allocate(pps%croplive(beg:end))
       allocate(pps%cropplant(beg:end)) !,numpft)) ! make 2-D if using
       allocate(pps%harvdate(beg:end))  !,numpft)) ! crop rotation
       allocate(pps%idop(beg:end))
       allocate(pps%peaklai(beg:end))
    end if
    allocate(pps%vds(beg:end))
    allocate(pps%slasun(beg:end))
    allocate(pps%slasha(beg:end))
    allocate(pps%lncsun(beg:end))
    allocate(pps%lncsha(beg:end))
    allocate(pps%vcmxsun(beg:end))
    allocate(pps%vcmxsha(beg:end))
    allocate(pps%gdir(beg:end))
    allocate(pps%omega(beg:end,1:numrad))
    allocate(pps%eff_kid(beg:end,1:numrad))
    allocate(pps%eff_kii(beg:end,1:numrad))
    allocate(pps%sun_faid(beg:end,1:numrad))
    allocate(pps%sun_faii(beg:end,1:numrad))
    allocate(pps%sha_faid(beg:end,1:numrad))
    allocate(pps%sha_faii(beg:end,1:numrad))
    allocate(pps%forc_hgt_u_pft(beg:end))
    allocate(pps%forc_hgt_t_pft(beg:end))
    allocate(pps%forc_hgt_q_pft(beg:end))
    ! 4/14/05: PET
    ! Adding isotope code
    allocate(pps%cisun(beg:end))
    allocate(pps%cisha(beg:end))
    allocate(pps%alphapsnsun(beg:end))
    allocate(pps%alphapsnsha(beg:end))

    allocate(pps%sandfrac(beg:end))
    allocate(pps%clayfrac(beg:end))
    pps%sandfrac(beg:end) = ival
    pps%clayfrac(beg:end) = ival
    allocate(pps%mlaidiff(beg:end))
    allocate(pps%rb1(beg:end))
    allocate(pps%annlai(12,beg:end))
    pps%mlaidiff(beg:end) = ival
    pps%rb1(beg:end) = ival
    pps%annlai(:,:) = ival
    
    pps%frac_veg_nosno(beg:end) = huge(1)
    pps%frac_veg_nosno_alb(beg:end) = 0
    pps%emv(beg:end) = ival
    pps%z0mv(beg:end) = ival
    pps%z0hv(beg:end) = ival
    pps%z0qv(beg:end) = ival
    pps%rootfr(beg:end,:nlevgrnd) = spval
    pps%rootr (beg:end,:nlevgrnd) = spval
    pps%rresis(beg:end,:nlevgrnd) = spval
    pps%dewmx(beg:end) = ival
    pps%rssun(beg:end) = ival
    pps%rssha(beg:end) = ival
    pps%laisun(beg:end) = ival
    pps%laisha(beg:end) = ival
    pps%btran(beg:end) = spval
    pps%fsun(beg:end) = spval
    pps%tlai(beg:end) = 0._r8
    pps%tsai(beg:end) = 0._r8
    pps%elai(beg:end) = 0._r8
    pps%esai(beg:end) = 0._r8
    pps%fwet(beg:end) = ival
    pps%fdry(beg:end) = ival
    pps%dt_veg(beg:end) = ival
    pps%htop(beg:end) = 0._r8
    pps%hbot(beg:end) = 0._r8
    pps%z0m(beg:end) = ival
    pps%displa(beg:end) = ival
    pps%albd(beg:end,:numrad) = ival
    pps%albi(beg:end,:numrad) = ival
    pps%fabd(beg:end,:numrad) = ival
    pps%fabi(beg:end,:numrad) = ival
    pps%ftdd(beg:end,:numrad) = ival
    pps%ftid(beg:end,:numrad) = ival
    pps%ftii(beg:end,:numrad) = ival
    pps%u10(beg:end) = ival
    pps%u10_clm(beg:end) = ival
    pps%va(beg:end) = ival
    pps%fv(beg:end) = ival
    pps%ram1(beg:end) = ival
    if ( crop_prog )then
       pps%hdidx(beg:end)       = ival
       pps%cumvd(beg:end)       = ival
       pps%htmx(beg:end)        = 0.0_r8
       pps%vf(beg:end)          = 0.0_r8
       pps%gddmaturity(beg:end) = spval
       pps%gdd0(beg:end)        = spval
       pps%gdd8(beg:end)        = spval
       pps%gdd10(beg:end)       = spval
       pps%gdd020(beg:end)      = spval
       pps%gdd820(beg:end)      = spval
       pps%gdd1020(beg:end)     = spval
       pps%gddplant(beg:end)    = spval
       pps%gddtsoi(beg:end)     = spval
       pps%huileaf(beg:end)     = ival
       pps%huigrain(beg:end)    = ival
       pps%aleafi(beg:end)      = ival
       pps%astemi(beg:end)      = ival
       pps%aleaf(beg:end)       = ival
       pps%astem(beg:end)       = ival
       pps%croplive(beg:end)    = .false.
       pps%cropplant(beg:end)   = .false.
       pps%harvdate(beg:end)    = huge(1)
       pps%idop(beg:end)        = huge(1)
       pps%peaklai(beg:end)     = 0
    end if
    pps%vds(beg:end) = ival
    pps%slasun(beg:end) = ival
    pps%slasha(beg:end) = ival
    pps%lncsun(beg:end) = ival
    pps%lncsha(beg:end) = ival
    pps%vcmxsun(beg:end) = ival
    pps%vcmxsha(beg:end) = ival
    pps%gdir(beg:end) = ival
    pps%omega(beg:end,1:numrad) = ival
    pps%eff_kid(beg:end,1:numrad) = ival
    pps%eff_kii(beg:end,1:numrad) = ival
    pps%sun_faid(beg:end,1:numrad) = ival
    pps%sun_faii(beg:end,1:numrad) = ival
    pps%sha_faid(beg:end,1:numrad) = ival
    pps%sha_faii(beg:end,1:numrad) = ival
    pps%forc_hgt_u_pft(beg:end) = ival
    pps%forc_hgt_t_pft(beg:end) = ival
    pps%forc_hgt_q_pft(beg:end) = ival
    ! 4/14/05: PET
    ! Adding isotope code
    pps%cisun(beg:end) = ival
    pps%cisha(beg:end) = ival
    if (use_c13) then
       pps%alphapsnsun(beg:end) = ival
       pps%alphapsnsha(beg:end) = ival
    endif

  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_epv_type
!
! !INTERFACE:
  subroutine init_pft_epv_type(beg, end, pepv)
!
! !DESCRIPTION:
! Initialize pft ecophysiological variables
!
! !USES:
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pepv%dormant_flag(beg:end))
    allocate(pepv%days_active(beg:end))
    allocate(pepv%onset_flag(beg:end))
    allocate(pepv%onset_counter(beg:end))
    allocate(pepv%onset_gddflag(beg:end))
    allocate(pepv%onset_fdd(beg:end))
    allocate(pepv%onset_gdd(beg:end))
    allocate(pepv%onset_swi(beg:end))
    allocate(pepv%offset_flag(beg:end))
    allocate(pepv%offset_counter(beg:end))
    allocate(pepv%offset_fdd(beg:end))
    allocate(pepv%offset_swi(beg:end))
    allocate(pepv%lgsf(beg:end))
    allocate(pepv%bglfr(beg:end))
    allocate(pepv%bgtr(beg:end))
    allocate(pepv%dayl(beg:end))
    allocate(pepv%prev_dayl(beg:end))
    allocate(pepv%annavg_t2m(beg:end))
    allocate(pepv%tempavg_t2m(beg:end))
    allocate(pepv%gpp(beg:end))
    allocate(pepv%availc(beg:end))
    allocate(pepv%xsmrpool_recover(beg:end))
    allocate(pepv%xsmrpool_c13ratio(beg:end))
    allocate(pepv%alloc_pnow(beg:end))
    allocate(pepv%c_allometry(beg:end))
    allocate(pepv%n_allometry(beg:end))
    allocate(pepv%plant_ndemand(beg:end))
    allocate(pepv%tempsum_potential_gpp(beg:end))
    allocate(pepv%annsum_potential_gpp(beg:end))
    allocate(pepv%tempmax_retransn(beg:end))
    allocate(pepv%annmax_retransn(beg:end))
    allocate(pepv%avail_retransn(beg:end))
    allocate(pepv%plant_nalloc(beg:end))
    allocate(pepv%plant_calloc(beg:end))
    allocate(pepv%excess_cflux(beg:end))
    allocate(pepv%downreg(beg:end))
    allocate(pepv%prev_leafc_to_litter(beg:end))
    allocate(pepv%prev_frootc_to_litter(beg:end))
    allocate(pepv%tempsum_npp(beg:end))
    allocate(pepv%annsum_npp(beg:end))
    allocate(pepv%tempsum_litfall(beg:end))
    allocate(pepv%annsum_litfall(beg:end))
    ! 4/21/05, PET
    ! Adding isotope code
    allocate(pepv%rc13_canair(beg:end))
    allocate(pepv%rc13_psnsun(beg:end))
    allocate(pepv%rc13_psnsha(beg:end))

    ival = 0.0
    pepv%dormant_flag(beg:end) = ival
    pepv%days_active(beg:end) = ival
    pepv%onset_flag(beg:end) = ival
    pepv%onset_counter(beg:end) = ival
    pepv%onset_gddflag(beg:end) = ival
    pepv%onset_fdd(beg:end) = ival
    pepv%onset_gdd(beg:end) = ival
    pepv%onset_swi(beg:end) = ival
    pepv%offset_flag(beg:end) = ival
    pepv%offset_counter(beg:end) = ival
    pepv%offset_fdd(beg:end) = ival
    pepv%offset_swi(beg:end) = ival
    pepv%lgsf(beg:end) = ival
    pepv%bglfr(beg:end) = ival
    pepv%bgtr(beg:end) = ival
    pepv%dayl(beg:end) = ival
    pepv%prev_dayl(beg:end) = ival
    pepv%annavg_t2m(beg:end) = ival
    pepv%tempavg_t2m(beg:end) = ival
    pepv%gpp(beg:end) = ival
    pepv%availc(beg:end) = ival
    pepv%xsmrpool_recover(beg:end) = ival
    if (use_c13) then
       pepv%xsmrpool_c13ratio(beg:end) = ival
    endif
    pepv%alloc_pnow(beg:end) = ival
    pepv%c_allometry(beg:end) = ival
    pepv%n_allometry(beg:end) = ival
    pepv%plant_ndemand(beg:end) = ival
    pepv%tempsum_potential_gpp(beg:end) = ival
    pepv%annsum_potential_gpp(beg:end) = ival
    pepv%tempmax_retransn(beg:end) = ival
    pepv%annmax_retransn(beg:end) = ival
    pepv%avail_retransn(beg:end) = ival
    pepv%plant_nalloc(beg:end) = ival
    pepv%plant_calloc(beg:end) = ival
    pepv%excess_cflux(beg:end) = ival
    pepv%downreg(beg:end) = ival
    pepv%prev_leafc_to_litter(beg:end) = ival
    pepv%prev_frootc_to_litter(beg:end) = ival
    pepv%tempsum_npp(beg:end) = ival
    pepv%annsum_npp(beg:end) = ival
    pepv%tempsum_litfall(beg:end) = ival
    pepv%annsum_litfall(beg:end) = ival
    if (use_c13) then
       ! 4/21/05, PET
       ! Adding isotope code
       pepv%rc13_canair(beg:end) = ival
       pepv%rc13_psnsun(beg:end) = ival
       pepv%rc13_psnsha(beg:end) = ival
    endif
    
  end subroutine init_pft_epv_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !USES:
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%pftmayexist(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%fpcgridold(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%greffic(beg:end))
    allocate(pdgvs%heatstress(beg:end))

    ival = 0.0
    pdgvs%agddtw(beg:end)           = ival
    pdgvs%agdd(beg:end)             = ival
    pdgvs%t_mo(beg:end)             = ival
    pdgvs%t_mo_min(beg:end)         = ival
    pdgvs%prec365(beg:end)          = ival
    pdgvs%present(beg:end)          = .false.
    pdgvs%pftmayexist(beg:end)      = .true.
    pdgvs%nind(beg:end)             = ival
    pdgvs%lm_ind(beg:end)           = ival
    pdgvs%lai_ind(beg:end)          = ival
    pdgvs%fpcinc(beg:end)           = ival
    pdgvs%fpcgrid(beg:end)          = ival
    pdgvs%fpcgridold(beg:end)       = ival
    pdgvs%crownarea(beg:end)        = ival
    pdgvs%greffic(beg:end)          = ival
    pdgvs%heatstress(beg:end)       = ival

  end subroutine init_pft_pdgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vstate_type
!
! !INTERFACE:
  subroutine init_pft_vstate_type(beg, end, pvs)
!
! !DESCRIPTION:
! Initialize pft VOC variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vstate_type), intent(inout):: pvs
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!------------------------------------------------------------------------

    allocate(pvs%t_veg24 (beg:end))
    allocate(pvs%t_veg240(beg:end))
    allocate(pvs%fsd24   (beg:end))
    allocate(pvs%fsd240  (beg:end))
    allocate(pvs%fsi24   (beg:end))
    allocate(pvs%fsi240  (beg:end))
    allocate(pvs%fsun24  (beg:end))
    allocate(pvs%fsun240 (beg:end))
    allocate(pvs%elai_p  (beg:end))

    ival = 0.0
    pvs%t_veg24 (beg:end)   = spval
    pvs%t_veg240(beg:end)   = spval
    pvs%fsd24   (beg:end)   = spval
    pvs%fsd240  (beg:end)   = spval
    pvs%fsi24   (beg:end)   = spval
    pvs%fsi240  (beg:end)   = spval
    pvs%fsun24  (beg:end)   = spval
    pvs%fsun240 (beg:end)   = spval
    pvs%elai_p  (beg:end)   = spval
  end subroutine init_pft_vstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_estate_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod, only : crop_prog
! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    allocate(pes%t_ref1m(beg:end)) ! mars
    allocate(pes%t_ref2m_min(beg:end))
    allocate(pes%t_ref2m_max(beg:end))
    allocate(pes%t_ref2m_min_inst(beg:end))
    allocate(pes%t_ref2m_max_inst(beg:end))
    allocate(pes%q_ref2m(beg:end))
    allocate(pes%t_ref2m_u(beg:end))
    allocate(pes%t_ref2m_r(beg:end))
    allocate(pes%t_ref2m_min_u(beg:end))
    allocate(pes%t_ref2m_min_r(beg:end))
    allocate(pes%t_ref2m_max_u(beg:end))
    allocate(pes%t_ref2m_max_r(beg:end))
    allocate(pes%t_ref2m_min_inst_u(beg:end))
    allocate(pes%t_ref2m_min_inst_r(beg:end))
    allocate(pes%t_ref2m_max_inst_u(beg:end))
    allocate(pes%t_ref2m_max_inst_r(beg:end))
    allocate(pes%t10(beg:end))
    if ( crop_prog )then
       allocate(pes%a10tmin(beg:end))
       allocate(pes%a5tmin(beg:end))
    end if
    allocate(pes%rh_ref2m(beg:end))
    allocate(pes%rh_ref2m_u(beg:end))
    allocate(pes%rh_ref2m_r(beg:end))
    allocate(pes%t_veg(beg:end))
    allocate(pes%thm(beg:end))

    ival = 0.0
    pes%t_ref2m(beg:end) = ival
    pes%t_ref1m(beg:end) = ival  ! mars
    pes%t_ref2m_min(beg:end) = ival
    pes%t_ref2m_max(beg:end) = ival
    pes%t_ref2m_min_inst(beg:end) = ival
    pes%t_ref2m_max_inst(beg:end) = ival
    pes%q_ref2m(beg:end) = ival
    pes%t_ref2m_u(beg:end) = ival
    pes%t_ref2m_r(beg:end) = ival
    pes%t_ref2m_min_u(beg:end) = ival
    pes%t_ref2m_min_r(beg:end) = ival
    pes%t_ref2m_max_u(beg:end) = ival
    pes%t_ref2m_max_r(beg:end) = ival
    pes%t_ref2m_min_inst_u(beg:end) = ival
    pes%t_ref2m_min_inst_r(beg:end) = ival
    pes%t_ref2m_max_inst_u(beg:end) = ival
    pes%t_ref2m_max_inst_r(beg:end) = ival
    pes%t10(beg:end)                = spval
    if ( crop_prog )then
       pes%a10tmin(beg:end)     = spval
       pes%a5tmin(beg:end)      = spval
    end if
    pes%rh_ref2m(beg:end) = ival
    pes%rh_ref2m_u(beg:end) = ival
    pes%rh_ref2m_r(beg:end) = ival
    pes%t_veg(beg:end) = ival
    pes%thm(beg:end) = ival

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(beg, end, pws)
!
! !DESCRIPTION:
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    ival = 0.0
    allocate(pws%h2ocan(beg:end))
    pws%h2ocan(beg:end) = ival

  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cstate_type
!
! !INTERFACE:
  subroutine init_pft_cstate_type(beg, end, pcs)
!
! !DESCRIPTION:
! Initialize pft carbon state
!
! !USES:
    use surfrdMod, only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pcs%leafc(beg:end))
    allocate(pcs%leafc_storage(beg:end))
    allocate(pcs%leafc_xfer(beg:end))
    allocate(pcs%frootc(beg:end))
    allocate(pcs%frootc_storage(beg:end))
    allocate(pcs%frootc_xfer(beg:end))
    allocate(pcs%livestemc(beg:end))
    allocate(pcs%livestemc_storage(beg:end))
    allocate(pcs%livestemc_xfer(beg:end))
    allocate(pcs%deadstemc(beg:end))
    allocate(pcs%deadstemc_storage(beg:end))
    allocate(pcs%deadstemc_xfer(beg:end))
    allocate(pcs%livecrootc(beg:end))
    allocate(pcs%livecrootc_storage(beg:end))
    allocate(pcs%livecrootc_xfer(beg:end))
    allocate(pcs%deadcrootc(beg:end))
    allocate(pcs%deadcrootc_storage(beg:end))
    allocate(pcs%deadcrootc_xfer(beg:end))
    allocate(pcs%gresp_storage(beg:end))
    allocate(pcs%gresp_xfer(beg:end))
    allocate(pcs%cpool(beg:end))
    allocate(pcs%xsmrpool(beg:end))
    allocate(pcs%pft_ctrunc(beg:end))
    allocate(pcs%dispvegc(beg:end))
    allocate(pcs%storvegc(beg:end))
    allocate(pcs%totvegc(beg:end))
    allocate(pcs%totpftc(beg:end))
    allocate(pcs%leafcmax(beg:end))
    if ( crop_prog )then
       allocate(pcs%grainc(beg:end))
       allocate(pcs%grainc_storage(beg:end))
       allocate(pcs%grainc_xfer(beg:end))
    end if
    allocate(pcs%woodc(beg:end))

    ival = 0.0
    pcs%leafc(beg:end) = ival
    pcs%leafc_storage(beg:end) = ival
    pcs%leafc_xfer(beg:end) = ival
    pcs%frootc(beg:end) = ival
    pcs%frootc_storage(beg:end) = ival
    pcs%frootc_xfer(beg:end) = ival
    pcs%livestemc(beg:end) = ival
    pcs%livestemc_storage(beg:end) = ival
    pcs%livestemc_xfer(beg:end) = ival
    pcs%deadstemc(beg:end) = ival
    pcs%deadstemc_storage(beg:end) = ival
    pcs%deadstemc_xfer(beg:end) = ival
    pcs%livecrootc(beg:end) = ival
    pcs%livecrootc_storage(beg:end) = ival
    pcs%livecrootc_xfer(beg:end) = ival
    pcs%deadcrootc(beg:end) = ival
    pcs%deadcrootc_storage(beg:end) = ival
    pcs%deadcrootc_xfer(beg:end) = ival
    pcs%gresp_storage(beg:end) = ival
    pcs%gresp_xfer(beg:end) = ival
    pcs%cpool(beg:end) = ival
    pcs%xsmrpool(beg:end) = ival
    pcs%pft_ctrunc(beg:end) = ival
    pcs%dispvegc(beg:end) = ival
    pcs%storvegc(beg:end) = ival
    pcs%totvegc(beg:end) = ival
    pcs%totpftc(beg:end) = ival
    pcs%leafcmax(beg:end) = ival
    if ( crop_prog )then
       pcs%grainc(beg:end)         = ival
       pcs%grainc_storage(beg:end) = ival
       pcs%grainc_xfer(beg:end)    = ival
    end if
    pcs%woodc(beg:end) = ival

  end subroutine init_pft_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nstate_type
!
! !INTERFACE:
  subroutine init_pft_nstate_type(beg, end, pns)
!
! !DESCRIPTION:
! Initialize pft nitrogen state
!
! !USES:
    use surfrdMod, only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    if ( crop_prog )then
       allocate(pns%grainn(beg:end))
       allocate(pns%grainn_storage(beg:end))
       allocate(pns%grainn_xfer(beg:end))
    end if
    allocate(pns%leafn(beg:end))
    allocate(pns%leafn_storage(beg:end))
    allocate(pns%leafn_xfer(beg:end))
    allocate(pns%frootn(beg:end))
    allocate(pns%frootn_storage(beg:end))
    allocate(pns%frootn_xfer(beg:end))
    allocate(pns%livestemn(beg:end))
    allocate(pns%livestemn_storage(beg:end))
    allocate(pns%livestemn_xfer(beg:end))
    allocate(pns%deadstemn(beg:end))
    allocate(pns%deadstemn_storage(beg:end))
    allocate(pns%deadstemn_xfer(beg:end))
    allocate(pns%livecrootn(beg:end))
    allocate(pns%livecrootn_storage(beg:end))
    allocate(pns%livecrootn_xfer(beg:end))
    allocate(pns%deadcrootn(beg:end))
    allocate(pns%deadcrootn_storage(beg:end))
    allocate(pns%deadcrootn_xfer(beg:end))
    allocate(pns%retransn(beg:end))
    allocate(pns%npool(beg:end))
    allocate(pns%pft_ntrunc(beg:end))
    allocate(pns%dispvegn(beg:end))
    allocate(pns%storvegn(beg:end))
    allocate(pns%totvegn(beg:end))
    allocate(pns%totpftn(beg:end))

    ival = 0.0
    if ( crop_prog )then
       pns%grainn(beg:end)         = ival
       pns%grainn_storage(beg:end) = ival
       pns%grainn_xfer(beg:end)    = ival
    end if
    pns%leafn(beg:end) = ival
    pns%leafn_storage(beg:end) = ival
    pns%leafn_xfer(beg:end) = ival
    pns%frootn(beg:end) = ival
    pns%frootn_storage(beg:end) = ival
    pns%frootn_xfer(beg:end) = ival
    pns%livestemn(beg:end) = ival
    pns%livestemn_storage(beg:end) = ival
    pns%livestemn_xfer(beg:end) = ival
    pns%deadstemn(beg:end) = ival
    pns%deadstemn_storage(beg:end) = ival
    pns%deadstemn_xfer(beg:end) = ival
    pns%livecrootn(beg:end) = ival
    pns%livecrootn_storage(beg:end) = ival
    pns%livecrootn_xfer(beg:end) = ival
    pns%deadcrootn(beg:end) = ival
    pns%deadcrootn_storage(beg:end) = ival
    pns%deadcrootn_xfer(beg:end) = ival
    pns%retransn(beg:end) = ival
    pns%npool(beg:end) = ival
    pns%pft_ntrunc(beg:end) = ival
    pns%dispvegn(beg:end) = ival
    pns%storvegn(beg:end) = ival
    pns%totvegn(beg:end) = ival
    pns%totpftn(beg:end) = ival

  end subroutine init_pft_nstate_type
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(beg, end, pef)
!
! !DESCRIPTION:
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pef%sabg(beg:end))
    allocate(pef%sabv(beg:end))
    allocate(pef%fsa(beg:end))
    allocate(pef%fsa_u(beg:end))
    allocate(pef%fsa_r(beg:end))
    allocate(pef%fsr(beg:end))
    allocate(pef%parsun(beg:end))
    allocate(pef%parsha(beg:end))
    allocate(pef%dlrad(beg:end))
    allocate(pef%ulrad(beg:end))
    allocate(pef%co2ice_lh(beg:end))  ! mars
    allocate(pef%eflx_lh_co2(beg:end)) ! Mars
    allocate(pef%eflx_lh_tot(beg:end))
    allocate(pef%eflx_lh_tot_u(beg:end))
    allocate(pef%eflx_lh_tot_r(beg:end))
    allocate(pef%eflx_lh_grnd(beg:end))
    allocate(pef%eflx_soil_grnd(beg:end))
    allocate(pef%eflx_soil_grnd_u(beg:end))
    allocate(pef%eflx_soil_grnd_r(beg:end))
    allocate(pef%eflx_sh_tot(beg:end))
    allocate(pef%eflx_sh_tot_u(beg:end))
    allocate(pef%eflx_sh_tot_r(beg:end))
    allocate(pef%eflx_sh_grnd(beg:end))
    allocate(pef%eflx_sh_veg(beg:end))
    allocate(pef%eflx_lh_vege(beg:end))
    allocate(pef%eflx_lh_vegt(beg:end))
    allocate(pef%eflx_wasteheat_pft(beg:end))
    allocate(pef%eflx_heat_from_ac_pft(beg:end))
    allocate(pef%eflx_traffic_pft(beg:end))
    allocate(pef%eflx_anthro(beg:end))
    allocate(pef%cgrnd(beg:end))
    allocate(pef%cgrndl(beg:end))
    allocate(pef%cgrnds(beg:end))
    allocate(pef%eflx_gnet(beg:end))
    allocate(pef%dgnetdT(beg:end))
    allocate(pef%eflx_lwrad_out(beg:end))
    allocate(pef%eflx_lwrad_net(beg:end))
    allocate(pef%eflx_lwrad_net_u(beg:end))
    allocate(pef%eflx_lwrad_net_r(beg:end))
    allocate(pef%netrad(beg:end))
    allocate(pef%fsds_vis_d(beg:end))
    allocate(pef%fsds_nir_d(beg:end))
    allocate(pef%fsds_vis_i(beg:end))
    allocate(pef%fsds_nir_i(beg:end))
    allocate(pef%fsr_vis_d(beg:end))
    allocate(pef%fsr_nir_d(beg:end))
    allocate(pef%fsr_vis_i(beg:end))
    allocate(pef%fsr_nir_i(beg:end))
    allocate(pef%fsds_vis_d_ln(beg:end))
    allocate(pef%fsds_nir_d_ln(beg:end))
    allocate(pef%fsr_vis_d_ln(beg:end))
    allocate(pef%fsr_nir_d_ln(beg:end))
    allocate(pef%sun_add(beg:end,1:numrad))
    allocate(pef%tot_aid(beg:end,1:numrad))
    allocate(pef%sun_aid(beg:end,1:numrad))
    allocate(pef%sun_aii(beg:end,1:numrad))
    allocate(pef%sha_aid(beg:end,1:numrad))
    allocate(pef%sha_aii(beg:end,1:numrad))
    allocate(pef%sun_atot(beg:end,1:numrad))
    allocate(pef%sha_atot(beg:end,1:numrad))
    allocate(pef%sun_alf(beg:end,1:numrad))
    allocate(pef%sha_alf(beg:end,1:numrad))
    allocate(pef%sun_aperlai(beg:end,1:numrad))
    allocate(pef%sha_aperlai(beg:end,1:numrad))
    allocate(pef%sabg_lyr(beg:end,-nlevsno+1:1))
    allocate(pef%sfc_frc_aer(beg:end))
    allocate(pef%sfc_frc_bc(beg:end))
    allocate(pef%sfc_frc_oc(beg:end))
    allocate(pef%sfc_frc_dst(beg:end))
    allocate(pef%sfc_frc_aer_sno(beg:end))
    allocate(pef%sfc_frc_bc_sno(beg:end))
    allocate(pef%sfc_frc_oc_sno(beg:end))
    allocate(pef%sfc_frc_dst_sno(beg:end))
    allocate(pef%fsr_sno_vd(beg:end))
    allocate(pef%fsr_sno_nd(beg:end))
    allocate(pef%fsr_sno_vi(beg:end))
    allocate(pef%fsr_sno_ni(beg:end))
    allocate(pef%fsds_sno_vd(beg:end))
    allocate(pef%fsds_sno_nd(beg:end))
    allocate(pef%fsds_sno_vi(beg:end))
    allocate(pef%fsds_sno_ni(beg:end))

    ival = 0.0
    pef%sabg(beg:end) = ival
    pef%sabv(beg:end) = ival
    pef%fsa(beg:end) = ival
    pef%fsa_u(beg:end) = ival
    pef%fsa_r(beg:end) = ival
    pef%fsr(beg:end) = ival
    pef%parsun(beg:end) = ival
    pef%parsha(beg:end) = ival
    pef%dlrad(beg:end) = ival
    pef%ulrad(beg:end) = ival
    pef%co2ice_lh(beg:end) = ival   ! mars
    pef%eflx_lh_co2(beg:end) = ival    ! mars
    pef%eflx_lh_tot(beg:end) = ival
    pef%eflx_lh_tot_u(beg:end) = ival
    pef%eflx_lh_tot_r(beg:end) = ival
    pef%eflx_lh_grnd(beg:end) = ival
    pef%eflx_soil_grnd(beg:end) = ival
    pef%eflx_soil_grnd_u(beg:end) = ival
    pef%eflx_soil_grnd_r(beg:end) = ival
    pef%eflx_sh_tot(beg:end) = ival
    pef%eflx_sh_tot_u(beg:end) = ival
    pef%eflx_sh_tot_r(beg:end) = ival
    pef%eflx_sh_grnd(beg:end) = ival
    pef%eflx_sh_veg(beg:end) = ival
    pef%eflx_lh_vege(beg:end) = ival
    pef%eflx_lh_vegt(beg:end) = ival
    pef%eflx_wasteheat_pft(beg:end) = ival
    pef%eflx_heat_from_ac_pft(beg:end) = ival
    pef%eflx_traffic_pft(beg:end) = ival
    pef%eflx_anthro(beg:end) = ival
    pef%cgrnd(beg:end) = ival
    pef%cgrndl(beg:end) = ival
    pef%cgrnds(beg:end) = ival
    pef%eflx_gnet(beg:end) = ival
    pef%dgnetdT(beg:end) = ival
    pef%eflx_lwrad_out(beg:end) = ival
    pef%eflx_lwrad_net(beg:end) = ival
    pef%eflx_lwrad_net_u(beg:end) = ival
    pef%eflx_lwrad_net_r(beg:end) = ival
    pef%netrad(beg:end) = ival
    pef%fsds_vis_d(beg:end) = ival
    pef%fsds_nir_d(beg:end) = ival
    pef%fsds_vis_i(beg:end) = ival
    pef%fsds_nir_i(beg:end) = ival
    pef%fsr_vis_d(beg:end) = ival
    pef%fsr_nir_d(beg:end) = ival
    pef%fsr_vis_i(beg:end) = ival
    pef%fsr_nir_i(beg:end) = ival
    pef%fsds_vis_d_ln(beg:end) = ival
    pef%fsds_nir_d_ln(beg:end) = ival
    pef%fsr_vis_d_ln(beg:end) = ival
    pef%fsr_nir_d_ln(beg:end) = ival
    pef%sun_add(beg:end,1:numrad) = ival
    pef%tot_aid(beg:end,1:numrad) = ival
    pef%sun_aid(beg:end,1:numrad) = ival
    pef%sun_aii(beg:end,1:numrad) = ival
    pef%sha_aid(beg:end,1:numrad) = ival
    pef%sha_aii(beg:end,1:numrad) = ival
    pef%sun_atot(beg:end,1:numrad) = ival
    pef%sha_atot(beg:end,1:numrad) = ival
    pef%sun_alf(beg:end,1:numrad) = ival
    pef%sha_alf(beg:end,1:numrad) = ival
    pef%sun_aperlai(beg:end,1:numrad) = ival
    pef%sha_aperlai(beg:end,1:numrad) = ival
    pef%sabg_lyr(beg:end,-nlevsno+1:1) = ival
    pef%sfc_frc_aer(beg:end) = ival
    pef%sfc_frc_bc(beg:end) = ival
    pef%sfc_frc_oc(beg:end) = ival
    pef%sfc_frc_dst(beg:end) = ival
    pef%sfc_frc_aer_sno(beg:end) = ival
    pef%sfc_frc_bc_sno(beg:end) = ival
    pef%sfc_frc_oc_sno(beg:end) = ival
    pef%sfc_frc_dst_sno(beg:end) = ival
    pef%fsr_sno_vd(beg:end) = ival
    pef%fsr_sno_nd(beg:end) = ival
    pef%fsr_sno_vi(beg:end) = ival
    pef%fsr_sno_ni(beg:end) = ival
    pef%fsds_sno_vd(beg:end) = ival
    pef%fsds_sno_nd(beg:end) = ival
    pef%fsds_sno_vi(beg:end) = ival
    pef%fsds_sno_ni(beg:end) = ival
  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(beg, end, pmf)
!
! !DESCRIPTION:
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pmf%taux(beg:end))
    allocate(pmf%tauy(beg:end))

    ival = 0.0
    pmf%taux(beg:end) = ival
    pmf%tauy(beg:end) = ival

  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(beg, end, pwf)
!
! !DESCRIPTION:
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
    real(r8) :: ival ! mars
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pwf%qflx_prec_intr(beg:end))
    allocate(pwf%qflx_prec_grnd(beg:end))
    allocate(pwf%qflx_rain_grnd(beg:end))
    allocate(pwf%qflx_snow_grnd(beg:end))
    allocate(pwf%qflx_snwcp_liq(beg:end))
    allocate(pwf%qflx_snwcp_ice(beg:end))
    allocate(pwf%qflx_evap_veg(beg:end))
    allocate(pwf%qflx_tran_veg(beg:end))
    allocate(pwf%qflx_evap_can(beg:end))
    allocate(pwf%qflx_evap_soi(beg:end))
    allocate(pwf%qflx_evap_tot(beg:end))
    allocate(pwf%qflx_evap_grnd(beg:end))
    allocate(pwf%qflx_dew_grnd(beg:end))
    allocate(pwf%qflx_sub_snow(beg:end))
    allocate(pwf%qflx_dew_snow(beg:end))

    ival = 0.0
    pwf%qflx_prec_intr(beg:end) = ival
    pwf%qflx_prec_grnd(beg:end) = ival
    pwf%qflx_rain_grnd(beg:end) = ival
    pwf%qflx_snow_grnd(beg:end) = ival
    pwf%qflx_snwcp_liq(beg:end) = ival
    pwf%qflx_snwcp_ice(beg:end) = ival
    pwf%qflx_evap_veg(beg:end) = ival
    pwf%qflx_tran_veg(beg:end) = ival
    pwf%qflx_evap_can(beg:end) = ival
    pwf%qflx_evap_soi(beg:end) = ival
    pwf%qflx_evap_tot(beg:end) = ival
    pwf%qflx_evap_grnd(beg:end) = ival
    pwf%qflx_dew_grnd(beg:end) = ival
    pwf%qflx_sub_snow(beg:end) = ival
    pwf%qflx_dew_snow(beg:end) = ival

  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !USES:
    use clm_varcon, only : spval
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    allocate(pcf%psnsha(beg:end))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%fco2(beg:end))

    allocate(pcf%m_leafc_to_litter(beg:end))
    allocate(pcf%m_frootc_to_litter(beg:end))
    allocate(pcf%m_leafc_storage_to_litter(beg:end))
    allocate(pcf%m_frootc_storage_to_litter(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_to_litter(beg:end))
    allocate(pcf%m_deadstemc_to_litter(beg:end))
    allocate(pcf%m_livecrootc_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_to_litter(beg:end))
    allocate(pcf%m_gresp_storage_to_litter(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter(beg:end))
    allocate(pcf%hrv_leafc_to_litter(beg:end))             
    allocate(pcf%hrv_leafc_storage_to_litter(beg:end))     
    allocate(pcf%hrv_leafc_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_frootc_to_litter(beg:end))            
    allocate(pcf%hrv_frootc_storage_to_litter(beg:end))    
    allocate(pcf%hrv_frootc_xfer_to_litter(beg:end))       
    allocate(pcf%hrv_livestemc_to_litter(beg:end))         
    allocate(pcf%hrv_livestemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_livestemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(pcf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(pcf%hrv_deadstemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_deadstemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_livecrootc_to_litter(beg:end))        
    allocate(pcf%hrv_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_livecrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_deadcrootc_to_litter(beg:end))        
    allocate(pcf%hrv_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_gresp_storage_to_litter(beg:end))     
    allocate(pcf%hrv_gresp_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_xsmrpool_to_atm(beg:end))                 
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%leafc_xfer_to_leafc(beg:end))
    allocate(pcf%frootc_xfer_to_frootc(beg:end))
    allocate(pcf%livestemc_xfer_to_livestemc(beg:end))
    allocate(pcf%deadstemc_xfer_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_xfer_to_livecrootc(beg:end))
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(beg:end))
    allocate(pcf%leafc_to_litter(beg:end))
    allocate(pcf%frootc_to_litter(beg:end))
    allocate(pcf%leaf_mr(beg:end))
    allocate(pcf%froot_mr(beg:end))
    allocate(pcf%livestem_mr(beg:end))
    allocate(pcf%livecroot_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%psnsun_to_cpool(beg:end))
    allocate(pcf%psnshade_to_cpool(beg:end))
    allocate(pcf%cpool_to_xsmrpool(beg:end))
    allocate(pcf%cpool_to_leafc(beg:end))
    allocate(pcf%cpool_to_leafc_storage(beg:end))
    allocate(pcf%cpool_to_frootc(beg:end))
    allocate(pcf%cpool_to_frootc_storage(beg:end))
    allocate(pcf%cpool_to_livestemc(beg:end))
    allocate(pcf%cpool_to_livestemc_storage(beg:end))
    allocate(pcf%cpool_to_deadstemc(beg:end))
    allocate(pcf%cpool_to_deadstemc_storage(beg:end))
    allocate(pcf%cpool_to_livecrootc(beg:end))
    allocate(pcf%cpool_to_livecrootc_storage(beg:end))
    allocate(pcf%cpool_to_deadcrootc(beg:end))
    allocate(pcf%cpool_to_deadcrootc_storage(beg:end))
    allocate(pcf%cpool_to_gresp_storage(beg:end))
    allocate(pcf%cpool_leaf_gr(beg:end))
    allocate(pcf%cpool_leaf_storage_gr(beg:end))
    allocate(pcf%transfer_leaf_gr(beg:end))
    allocate(pcf%cpool_froot_gr(beg:end))
    allocate(pcf%cpool_froot_storage_gr(beg:end))
    allocate(pcf%transfer_froot_gr(beg:end))
    allocate(pcf%cpool_livestem_gr(beg:end))
    allocate(pcf%cpool_livestem_storage_gr(beg:end))
    allocate(pcf%transfer_livestem_gr(beg:end))
    allocate(pcf%cpool_deadstem_gr(beg:end))
    allocate(pcf%cpool_deadstem_storage_gr(beg:end))
    allocate(pcf%transfer_deadstem_gr(beg:end))
    allocate(pcf%cpool_livecroot_gr(beg:end))
    allocate(pcf%cpool_livecroot_storage_gr(beg:end))
    allocate(pcf%transfer_livecroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_storage_gr(beg:end))
    allocate(pcf%transfer_deadcroot_gr(beg:end))
    allocate(pcf%leafc_storage_to_xfer(beg:end))
    allocate(pcf%frootc_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_storage_to_xfer(beg:end))
    allocate(pcf%deadstemc_storage_to_xfer(beg:end))
    allocate(pcf%livecrootc_storage_to_xfer(beg:end))
    allocate(pcf%deadcrootc_storage_to_xfer(beg:end))
    allocate(pcf%gresp_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_to_deadcrootc(beg:end))
    allocate(pcf%gpp(beg:end))
    allocate(pcf%mr(beg:end))
    allocate(pcf%current_gr(beg:end))
    allocate(pcf%transfer_gr(beg:end))
    allocate(pcf%storage_gr(beg:end))
    allocate(pcf%gr(beg:end))
    allocate(pcf%ar(beg:end))
    allocate(pcf%rr(beg:end))
    allocate(pcf%npp(beg:end))
    allocate(pcf%agnpp(beg:end))
    allocate(pcf%bgnpp(beg:end))
    allocate(pcf%litfall(beg:end))
    allocate(pcf%vegfire(beg:end))
    allocate(pcf%wood_harvestc(beg:end))
    allocate(pcf%pft_cinputs(beg:end))
    allocate(pcf%pft_coutputs(beg:end))
    allocate(pcf%pft_fire_closs(beg:end))
    if ( crop_prog )then
       allocate(pcf%xsmrpool_to_atm(beg:end))
       allocate(pcf%grainc_xfer_to_grainc(beg:end))
       allocate(pcf%livestemc_to_litter(beg:end))
       allocate(pcf%grainc_to_food(beg:end))
       allocate(pcf%cpool_to_grainc(beg:end))
       allocate(pcf%cpool_to_grainc_storage(beg:end))
       allocate(pcf%cpool_grain_gr(beg:end))
       allocate(pcf%cpool_grain_storage_gr(beg:end))
       allocate(pcf%transfer_grain_gr(beg:end))
       allocate(pcf%grainc_storage_to_xfer(beg:end))
    end if
    if (use_cn) then
       allocate(pcf%frootc_alloc(beg:end))
       allocate(pcf%frootc_loss(beg:end))
       allocate(pcf%leafc_alloc(beg:end))
       allocate(pcf%leafc_loss(beg:end))
       allocate(pcf%woodc_alloc(beg:end))
       allocate(pcf%woodc_loss(beg:end))
    end if

    ival = 0.0
    pcf%psnsun(beg:end) = ival
    pcf%psnsha(beg:end) = ival
    pcf%fpsn(beg:end) = spval
    pcf%fco2(beg:end) = ival !0.0_r8

    pcf%m_leafc_to_litter(beg:end) = ival
    pcf%m_frootc_to_litter(beg:end) = ival
    pcf%m_leafc_storage_to_litter(beg:end) = ival
    pcf%m_frootc_storage_to_litter(beg:end) = ival
    pcf%m_livestemc_storage_to_litter(beg:end) = ival
    pcf%m_deadstemc_storage_to_litter(beg:end) = ival
    pcf%m_livecrootc_storage_to_litter(beg:end) = ival
    pcf%m_deadcrootc_storage_to_litter(beg:end) = ival
    pcf%m_leafc_xfer_to_litter(beg:end) = ival
    pcf%m_frootc_xfer_to_litter(beg:end) = ival
    pcf%m_livestemc_xfer_to_litter(beg:end) = ival
    pcf%m_deadstemc_xfer_to_litter(beg:end) = ival
    pcf%m_livecrootc_xfer_to_litter(beg:end) = ival
    pcf%m_deadcrootc_xfer_to_litter(beg:end) = ival
    pcf%m_livestemc_to_litter(beg:end) = ival
    pcf%m_deadstemc_to_litter(beg:end) = ival
    pcf%m_livecrootc_to_litter(beg:end) = ival
    pcf%m_deadcrootc_to_litter(beg:end) = ival
    pcf%m_gresp_storage_to_litter(beg:end) = ival
    pcf%m_gresp_xfer_to_litter(beg:end) = ival
    pcf%hrv_leafc_to_litter(beg:end) = ival             
    pcf%hrv_leafc_storage_to_litter(beg:end) = ival     
    pcf%hrv_leafc_xfer_to_litter(beg:end) = ival        
    pcf%hrv_frootc_to_litter(beg:end) = ival            
    pcf%hrv_frootc_storage_to_litter(beg:end) = ival    
    pcf%hrv_frootc_xfer_to_litter(beg:end) = ival       
    pcf%hrv_livestemc_to_litter(beg:end) = ival         
    pcf%hrv_livestemc_storage_to_litter(beg:end) = ival 
    pcf%hrv_livestemc_xfer_to_litter(beg:end) = ival    
    pcf%hrv_deadstemc_to_prod10c(beg:end) = ival        
    pcf%hrv_deadstemc_to_prod100c(beg:end) = ival       
    pcf%hrv_deadstemc_storage_to_litter(beg:end) = ival 
    pcf%hrv_deadstemc_xfer_to_litter(beg:end) = ival    
    pcf%hrv_livecrootc_to_litter(beg:end) = ival        
    pcf%hrv_livecrootc_storage_to_litter(beg:end) = ival
    pcf%hrv_livecrootc_xfer_to_litter(beg:end) = ival   
    pcf%hrv_deadcrootc_to_litter(beg:end) = ival        
    pcf%hrv_deadcrootc_storage_to_litter(beg:end) = ival
    pcf%hrv_deadcrootc_xfer_to_litter(beg:end) = ival   
    pcf%hrv_gresp_storage_to_litter(beg:end) = ival     
    pcf%hrv_gresp_xfer_to_litter(beg:end) = ival        
    pcf%hrv_xsmrpool_to_atm(beg:end) = ival                 
    pcf%m_leafc_to_fire(beg:end) = ival
    pcf%m_frootc_to_fire(beg:end) = ival
    pcf%m_leafc_storage_to_fire(beg:end) = ival
    pcf%m_frootc_storage_to_fire(beg:end) = ival
    pcf%m_livestemc_storage_to_fire(beg:end) = ival
    pcf%m_deadstemc_storage_to_fire(beg:end) = ival
    pcf%m_livecrootc_storage_to_fire(beg:end) = ival
    pcf%m_deadcrootc_storage_to_fire(beg:end) = ival
    pcf%m_leafc_xfer_to_fire(beg:end) = ival
    pcf%m_frootc_xfer_to_fire(beg:end) = ival
    pcf%m_livestemc_xfer_to_fire(beg:end) = ival
    pcf%m_deadstemc_xfer_to_fire(beg:end) = ival
    pcf%m_livecrootc_xfer_to_fire(beg:end) = ival
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = ival
    pcf%m_livestemc_to_fire(beg:end) = ival
    pcf%m_deadstemc_to_fire(beg:end) = ival
    pcf%m_deadstemc_to_litter_fire(beg:end) = ival
    pcf%m_livecrootc_to_fire(beg:end) = ival
    pcf%m_deadcrootc_to_fire(beg:end) = ival
    pcf%m_deadcrootc_to_litter_fire(beg:end) = ival
    pcf%m_gresp_storage_to_fire(beg:end) = ival
    pcf%m_gresp_xfer_to_fire(beg:end) = ival
    pcf%leafc_xfer_to_leafc(beg:end) = ival
    pcf%frootc_xfer_to_frootc(beg:end) = ival
    pcf%livestemc_xfer_to_livestemc(beg:end) = ival
    pcf%deadstemc_xfer_to_deadstemc(beg:end) = ival
    pcf%livecrootc_xfer_to_livecrootc(beg:end) = ival
    pcf%deadcrootc_xfer_to_deadcrootc(beg:end) = ival
    pcf%leafc_to_litter(beg:end) = ival
    pcf%frootc_to_litter(beg:end) = ival
    pcf%leaf_mr(beg:end) = ival
    pcf%froot_mr(beg:end) = ival
    pcf%livestem_mr(beg:end) = ival
    pcf%livecroot_mr(beg:end) = ival
    pcf%leaf_curmr(beg:end) = ival
    pcf%froot_curmr(beg:end) = ival
    pcf%livestem_curmr(beg:end) = ival
    pcf%livecroot_curmr(beg:end) = ival
    pcf%leaf_xsmr(beg:end) = ival
    pcf%froot_xsmr(beg:end) = ival
    pcf%livestem_xsmr(beg:end) = ival
    pcf%livecroot_xsmr(beg:end) = ival
    pcf%psnsun_to_cpool(beg:end) = ival
    pcf%psnshade_to_cpool(beg:end) = ival
    pcf%cpool_to_xsmrpool(beg:end) = ival
    pcf%cpool_to_leafc(beg:end) = ival
    pcf%cpool_to_leafc_storage(beg:end) = ival
    pcf%cpool_to_frootc(beg:end) = ival
    pcf%cpool_to_frootc_storage(beg:end) = ival
    pcf%cpool_to_livestemc(beg:end) = ival
    pcf%cpool_to_livestemc_storage(beg:end) = ival
    pcf%cpool_to_deadstemc(beg:end) = ival
    pcf%cpool_to_deadstemc_storage(beg:end) = ival
    pcf%cpool_to_livecrootc(beg:end) = ival
    pcf%cpool_to_livecrootc_storage(beg:end) = ival
    pcf%cpool_to_deadcrootc(beg:end) = ival
    pcf%cpool_to_deadcrootc_storage(beg:end) = ival
    pcf%cpool_to_gresp_storage(beg:end) = ival
    pcf%cpool_leaf_gr(beg:end) = ival
    pcf%cpool_leaf_storage_gr(beg:end) = ival
    pcf%transfer_leaf_gr(beg:end) = ival
    pcf%cpool_froot_gr(beg:end) = ival
    pcf%cpool_froot_storage_gr(beg:end) = ival
    pcf%transfer_froot_gr(beg:end) = ival
    pcf%cpool_livestem_gr(beg:end) = ival
    pcf%cpool_livestem_storage_gr(beg:end) = ival
    pcf%transfer_livestem_gr(beg:end) = ival
    pcf%cpool_deadstem_gr(beg:end) = ival
    pcf%cpool_deadstem_storage_gr(beg:end) = ival
    pcf%transfer_deadstem_gr(beg:end) = ival
    pcf%cpool_livecroot_gr(beg:end) = ival
    pcf%cpool_livecroot_storage_gr(beg:end) = ival
    pcf%transfer_livecroot_gr(beg:end) = ival
    pcf%cpool_deadcroot_gr(beg:end) = ival
    pcf%cpool_deadcroot_storage_gr(beg:end) = ival
    pcf%transfer_deadcroot_gr(beg:end) = ival
    pcf%leafc_storage_to_xfer(beg:end) = ival
    pcf%frootc_storage_to_xfer(beg:end) = ival
    pcf%livestemc_storage_to_xfer(beg:end) = ival
    pcf%deadstemc_storage_to_xfer(beg:end) = ival
    pcf%livecrootc_storage_to_xfer(beg:end) = ival
    pcf%deadcrootc_storage_to_xfer(beg:end) = ival
    pcf%gresp_storage_to_xfer(beg:end) = ival
    pcf%livestemc_to_deadstemc(beg:end) = ival
    pcf%livecrootc_to_deadcrootc(beg:end) = ival
    pcf%gpp(beg:end) = ival
    pcf%mr(beg:end) = ival
    pcf%current_gr(beg:end) = ival
    pcf%transfer_gr(beg:end) = ival
    pcf%storage_gr(beg:end) = ival
    pcf%gr(beg:end) = ival
    pcf%ar(beg:end) = ival
    pcf%rr(beg:end) = ival
    pcf%npp(beg:end) = ival
    pcf%agnpp(beg:end) = ival
    pcf%bgnpp(beg:end) = ival
    pcf%litfall(beg:end) = ival
    pcf%vegfire(beg:end) = ival
    pcf%wood_harvestc(beg:end) = ival
    pcf%pft_cinputs(beg:end) = ival
    pcf%pft_coutputs(beg:end) = ival
    pcf%pft_fire_closs(beg:end) = ival
    if ( crop_prog )then
       pcf%xsmrpool_to_atm(beg:end)         = ival
       pcf%grainc_xfer_to_grainc(beg:end)   = ival
       pcf%livestemc_to_litter(beg:end)     = ival
       pcf%grainc_to_food(beg:end)          = ival
       pcf%cpool_to_grainc(beg:end)         = ival
       pcf%cpool_to_grainc_storage(beg:end) = ival
       pcf%cpool_grain_gr(beg:end)          = ival
       pcf%cpool_grain_storage_gr(beg:end)  = ival
       pcf%transfer_grain_gr(beg:end)       = ival
       pcf%grainc_storage_to_xfer(beg:end)  = ival
    end if
    if (use_cn) then
       pcf%frootc_alloc(beg:end) = ival
       pcf%frootc_loss(beg:end) = ival
       pcf%leafc_alloc(beg:end) = ival
       pcf%leafc_loss(beg:end) = ival
       pcf%woodc_alloc(beg:end) = ival
       pcf%woodc_loss(beg:end) = ival
    end if

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nflux_type
!
! !INTERFACE:
  subroutine init_pft_nflux_type(beg, end, pnf)
!
! !DESCRIPTION:
! Initialize pft nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pnf%m_leafn_to_litter(beg:end))
    allocate(pnf%m_frootn_to_litter(beg:end))
    allocate(pnf%m_leafn_storage_to_litter(beg:end))
    allocate(pnf%m_frootn_storage_to_litter(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_to_litter(beg:end))
    allocate(pnf%m_deadstemn_to_litter(beg:end))
    allocate(pnf%m_livecrootn_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_to_litter(beg:end))
    allocate(pnf%m_retransn_to_litter(beg:end))
    allocate(pnf%hrv_leafn_to_litter(beg:end))             
    allocate(pnf%hrv_frootn_to_litter(beg:end))            
    allocate(pnf%hrv_leafn_storage_to_litter(beg:end))     
    allocate(pnf%hrv_frootn_storage_to_litter(beg:end))    
    allocate(pnf%hrv_livestemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_deadstemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_leafn_xfer_to_litter(beg:end))        
    allocate(pnf%hrv_frootn_xfer_to_litter(beg:end))       
    allocate(pnf%hrv_livestemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_deadstemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_livecrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_livestemn_to_litter(beg:end))         
    allocate(pnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(pnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(pnf%hrv_livecrootn_to_litter(beg:end))        
    allocate(pnf%hrv_deadcrootn_to_litter(beg:end))        
    allocate(pnf%hrv_retransn_to_litter(beg:end))              
    allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))
    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_litter(beg:end))
    allocate(pnf%retransn_to_npool(beg:end))
    allocate(pnf%sminn_to_npool(beg:end))
    allocate(pnf%npool_to_leafn(beg:end))
    allocate(pnf%npool_to_leafn_storage(beg:end))
    allocate(pnf%npool_to_frootn(beg:end))
    allocate(pnf%npool_to_frootn_storage(beg:end))
    allocate(pnf%npool_to_livestemn(beg:end))
    allocate(pnf%npool_to_livestemn_storage(beg:end))
    allocate(pnf%npool_to_deadstemn(beg:end))
    allocate(pnf%npool_to_deadstemn_storage(beg:end))
    allocate(pnf%npool_to_livecrootn(beg:end))
    allocate(pnf%npool_to_livecrootn_storage(beg:end))
    allocate(pnf%npool_to_deadcrootn(beg:end))
    allocate(pnf%npool_to_deadcrootn_storage(beg:end))
    allocate(pnf%leafn_storage_to_xfer(beg:end))
    allocate(pnf%frootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_storage_to_xfer(beg:end))
    allocate(pnf%deadstemn_storage_to_xfer(beg:end))
    allocate(pnf%livecrootn_storage_to_xfer(beg:end))
    allocate(pnf%deadcrootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_to_deadstemn(beg:end))
    allocate(pnf%livestemn_to_retransn(beg:end))
    allocate(pnf%livecrootn_to_deadcrootn(beg:end))
    allocate(pnf%livecrootn_to_retransn(beg:end))
    allocate(pnf%ndeploy(beg:end))
    allocate(pnf%pft_ninputs(beg:end))
    allocate(pnf%pft_noutputs(beg:end))
    allocate(pnf%wood_harvestn(beg:end))
    allocate(pnf%pft_fire_nloss(beg:end))
    if ( crop_prog )then
       allocate(pnf%grainn_xfer_to_grainn(beg:end))
       allocate(pnf%livestemn_to_litter(beg:end))
       allocate(pnf%grainn_to_food(beg:end))
       allocate(pnf%npool_to_grainn(beg:end))
       allocate(pnf%npool_to_grainn_storage(beg:end))
       allocate(pnf%grainn_storage_to_xfer(beg:end))
    end if

    ival = 0.0
    pnf%m_leafn_to_litter(beg:end) = ival
    pnf%m_frootn_to_litter(beg:end) = ival
    pnf%m_leafn_storage_to_litter(beg:end) = ival
    pnf%m_frootn_storage_to_litter(beg:end) = ival
    pnf%m_livestemn_storage_to_litter(beg:end) = ival
    pnf%m_deadstemn_storage_to_litter(beg:end) = ival
    pnf%m_livecrootn_storage_to_litter(beg:end) = ival
    pnf%m_deadcrootn_storage_to_litter(beg:end) = ival
    pnf%m_leafn_xfer_to_litter(beg:end) = ival
    pnf%m_frootn_xfer_to_litter(beg:end) = ival
    pnf%m_livestemn_xfer_to_litter(beg:end) = ival
    pnf%m_deadstemn_xfer_to_litter(beg:end) = ival
    pnf%m_livecrootn_xfer_to_litter(beg:end) = ival
    pnf%m_deadcrootn_xfer_to_litter(beg:end) = ival
    pnf%m_livestemn_to_litter(beg:end) = ival
    pnf%m_deadstemn_to_litter(beg:end) = ival
    pnf%m_livecrootn_to_litter(beg:end) = ival
    pnf%m_deadcrootn_to_litter(beg:end) = ival
    pnf%m_retransn_to_litter(beg:end) = ival
    pnf%hrv_leafn_to_litter(beg:end) = ival             
    pnf%hrv_frootn_to_litter(beg:end) = ival            
    pnf%hrv_leafn_storage_to_litter(beg:end) = ival     
    pnf%hrv_frootn_storage_to_litter(beg:end) = ival    
    pnf%hrv_livestemn_storage_to_litter(beg:end) = ival 
    pnf%hrv_deadstemn_storage_to_litter(beg:end) = ival 
    pnf%hrv_livecrootn_storage_to_litter(beg:end) = ival
    pnf%hrv_deadcrootn_storage_to_litter(beg:end) = ival
    pnf%hrv_leafn_xfer_to_litter(beg:end) = ival        
    pnf%hrv_frootn_xfer_to_litter(beg:end) = ival       
    pnf%hrv_livestemn_xfer_to_litter(beg:end) = ival    
    pnf%hrv_deadstemn_xfer_to_litter(beg:end) = ival    
    pnf%hrv_livecrootn_xfer_to_litter(beg:end) = ival   
    pnf%hrv_deadcrootn_xfer_to_litter(beg:end) = ival   
    pnf%hrv_livestemn_to_litter(beg:end) = ival         
    pnf%hrv_deadstemn_to_prod10n(beg:end) = ival        
    pnf%hrv_deadstemn_to_prod100n(beg:end) = ival       
    pnf%hrv_livecrootn_to_litter(beg:end) = ival        
    pnf%hrv_deadcrootn_to_litter(beg:end) = ival        
    pnf%hrv_retransn_to_litter(beg:end) = ival           
    pnf%m_leafn_to_fire(beg:end) = ival
    pnf%m_frootn_to_fire(beg:end) = ival
    pnf%m_leafn_storage_to_fire(beg:end) = ival
    pnf%m_frootn_storage_to_fire(beg:end) = ival
    pnf%m_livestemn_storage_to_fire(beg:end) = ival
    pnf%m_deadstemn_storage_to_fire(beg:end) = ival
    pnf%m_livecrootn_storage_to_fire(beg:end) = ival
    pnf%m_deadcrootn_storage_to_fire(beg:end) = ival
    pnf%m_leafn_xfer_to_fire(beg:end) = ival
    pnf%m_frootn_xfer_to_fire(beg:end) = ival
    pnf%m_livestemn_xfer_to_fire(beg:end) = ival
    pnf%m_deadstemn_xfer_to_fire(beg:end) = ival
    pnf%m_livecrootn_xfer_to_fire(beg:end) = ival
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = ival
    pnf%m_livestemn_to_fire(beg:end) = ival
    pnf%m_deadstemn_to_fire(beg:end) = ival
    pnf%m_deadstemn_to_litter_fire(beg:end) = ival
    pnf%m_livecrootn_to_fire(beg:end) = ival
    pnf%m_deadcrootn_to_fire(beg:end) = ival
    pnf%m_deadcrootn_to_litter_fire(beg:end) = ival
    pnf%m_retransn_to_fire(beg:end) = ival
    pnf%leafn_xfer_to_leafn(beg:end) = ival
    pnf%frootn_xfer_to_frootn(beg:end) = ival
    pnf%livestemn_xfer_to_livestemn(beg:end) = ival
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = ival
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = ival
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = ival
    pnf%leafn_to_litter(beg:end) = ival
    pnf%leafn_to_retransn(beg:end) = ival
    pnf%frootn_to_litter(beg:end) = ival
    pnf%retransn_to_npool(beg:end) = ival
    pnf%sminn_to_npool(beg:end) = ival
    pnf%npool_to_leafn(beg:end) = ival
    pnf%npool_to_leafn_storage(beg:end) = ival
    pnf%npool_to_frootn(beg:end) = ival
    pnf%npool_to_frootn_storage(beg:end) = ival
    pnf%npool_to_livestemn(beg:end) = ival
    pnf%npool_to_livestemn_storage(beg:end) = ival
    pnf%npool_to_deadstemn(beg:end) = ival
    pnf%npool_to_deadstemn_storage(beg:end) = ival
    pnf%npool_to_livecrootn(beg:end) = ival
    pnf%npool_to_livecrootn_storage(beg:end) = ival
    pnf%npool_to_deadcrootn(beg:end) = ival
    pnf%npool_to_deadcrootn_storage(beg:end) = ival
    pnf%leafn_storage_to_xfer(beg:end) = ival
    pnf%frootn_storage_to_xfer(beg:end) = ival
    pnf%livestemn_storage_to_xfer(beg:end) = ival
    pnf%deadstemn_storage_to_xfer(beg:end) = ival
    pnf%livecrootn_storage_to_xfer(beg:end) = ival
    pnf%deadcrootn_storage_to_xfer(beg:end) = ival
    pnf%livestemn_to_deadstemn(beg:end) = ival
    pnf%livestemn_to_retransn(beg:end) = ival
    pnf%livecrootn_to_deadcrootn(beg:end) = ival
    pnf%livecrootn_to_retransn(beg:end) = ival
    pnf%ndeploy(beg:end) = ival
    pnf%pft_ninputs(beg:end) = ival
    pnf%pft_noutputs(beg:end) = ival
    pnf%wood_harvestn(beg:end) = ival
    pnf%pft_fire_nloss(beg:end) = ival
    if ( crop_prog )then
       pnf%grainn_xfer_to_grainn(beg:end)   = ival
       pnf%livestemn_to_litter(beg:end)     = ival
       pnf%grainn_to_food(beg:end)          = ival
       pnf%npool_to_grainn(beg:end)         = ival
       pnf%npool_to_grainn_storage(beg:end) = ival
       pnf%grainn_storage_to_xfer(beg:end)  = ival
    end if

  end subroutine init_pft_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(beg, end, pvf)
!
! !DESCRIPTION:
! Initialize pft VOC flux variables
!
    use clm_varcon, only : spval
    use shr_megan_mod, only: shr_megan_megcomps_n
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf
    real(r8) :: ival

    integer :: i
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! (heald, 08/06)
!
!EOP
!------------------------------------------------------------------------

    if (shr_megan_megcomps_n<1) return

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:shr_megan_megcomps_n))
    allocate(pvf%Eopt_out(beg:end))
    allocate(pvf%topt_out(beg:end))
    allocate(pvf%alpha_out(beg:end))
    allocate(pvf%cp_out(beg:end))
    allocate(pvf%para_out(beg:end))
    allocate(pvf%par24a_out(beg:end))
    allocate(pvf%par240a_out(beg:end))
    allocate(pvf%paru_out(beg:end))
    allocate(pvf%par24u_out(beg:end))
    allocate(pvf%par240u_out(beg:end))
    allocate(pvf%gamma_out(beg:end))
    allocate(pvf%gammaL_out(beg:end))
    allocate(pvf%gammaT_out(beg:end))
    allocate(pvf%gammaP_out(beg:end))
    allocate(pvf%gammaA_out(beg:end))
    allocate(pvf%gammaS_out(beg:end))
    allocate(pvf%gammaC_out(beg:end))

    ival = 0.0
    pvf%vocflx_tot(beg:end) = ival
    pvf%vocflx(beg:end,1:shr_megan_megcomps_n) = ival
    pvf%Eopt_out(beg:end) = ival
    pvf%topt_out(beg:end) = ival
    pvf%alpha_out(beg:end) = ival
    pvf%cp_out(beg:end) = ival
    pvf%para_out(beg:end) = ival
    pvf%par24a_out(beg:end) = ival
    pvf%par240a_out(beg:end) = ival
    pvf%paru_out(beg:end) = ival
    pvf%par24u_out(beg:end) = ival
    pvf%par240u_out(beg:end) = ival
    pvf%gamma_out(beg:end) = ival
    pvf%gammaL_out(beg:end) = ival
    pvf%gammaT_out(beg:end) = ival
    pvf%gammaP_out(beg:end) = ival
    pvf%gammaA_out(beg:end) = ival
    pvf%gammaS_out(beg:end) = ival
    pvf%gammaC_out(beg:end) = ival

    allocate(pvf%meg(shr_megan_megcomps_n))

    do i=1,shr_megan_megcomps_n
       allocate(pvf%meg(i)%flux_out(beg:end))
       pvf%meg(i)%flux_out(beg:end) = ival
    enddo

  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(beg, end, pdf)
!
! !DESCRIPTION:
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    allocate(pdf%vlc_trb_1(beg:end))
    allocate(pdf%vlc_trb_2(beg:end))
    allocate(pdf%vlc_trb_3(beg:end))
    allocate(pdf%vlc_trb_4(beg:end))

    ival = 0.0
    pdf%flx_mss_vrt_dst(beg:end,1:ndst) = ival
    pdf%flx_mss_vrt_dst_tot(beg:end) = ival
    pdf%vlc_trb(beg:end,1:ndst) = ival
    pdf%vlc_trb_1(beg:end) = ival
    pdf%vlc_trb_2(beg:end) = ival
    pdf%vlc_trb_3(beg:end) = ival
    pdf%vlc_trb_4(beg:end) = ival

  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_depvd_type
!
! !INTERFACE:
  subroutine init_pft_depvd_type(beg, end, pdd)

    use seq_drydep_mod, only:  n_drydep, drydep_method, DD_XLND
!
! !DESCRIPTION:
! Initialize pft dep velocity variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_depvd_type), intent(inout):: pdd
    integer :: i
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by James Sulzman 541-929-6183
!
!EOP
!------------------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(pdd%drydepvel(beg:end,n_drydep))
       ival = 0.0
       pdd%drydepvel = ival
    end if

  end subroutine init_pft_depvd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cps%snl(beg:end))      !* cannot be averaged up
    allocate(cps%isoicol(beg:end))  !* cannot be averaged up
    allocate(cps%ithincol(beg:end))  !* cannot be averaged up ! mars
    allocate(cps%ithincol_i(beg:end))  !* cannot be averaged up ! mars
    allocate(cps%bsw(beg:end,nlevgrnd))
    allocate(cps%watsat(beg:end,nlevgrnd))
    allocate(cps%watfc(beg:end,nlevgrnd))
    allocate(cps%watdry(beg:end,nlevgrnd)) 
    allocate(cps%watopt(beg:end,nlevgrnd)) 
    allocate(cps%hksat(beg:end,nlevgrnd))
    allocate(cps%sucsat(beg:end,nlevgrnd))
    allocate(cps%csol(beg:end,nlevgrnd))
    allocate(cps%tkmg(beg:end,nlevgrnd))
    allocate(cps%tkdry(beg:end,nlevgrnd))
    allocate(cps%tksatu(beg:end,nlevgrnd))
    allocate(cps%smpmin(beg:end))
    allocate(cps%hkdepth(beg:end))
    allocate(cps%wtfact(beg:end))
    allocate(cps%fracice(beg:end,nlevgrnd))
    allocate(cps%gwc_thr(beg:end))
    allocate(cps%mss_frc_cly_vld(beg:end))
    allocate(cps%mbl_bsn_fct(beg:end))
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snowdp(beg:end))
    allocate(cps%srf_emiss(beg:end))   ! mars 
    allocate(cps%co2dp(beg:end))   ! mars
    allocate(cps%co2_mass_change(beg:end))   ! mars
    allocate(cps%frac_co2(beg:end))   ! mars
    allocate(cps%frac_sno (beg:end))
    allocate(cps%zi(beg:end,-nlevsno+0:nlevgrnd))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%z (beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%imelt(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%eff_porosity(beg:end,nlevgrnd))
    allocate(cps%emg(beg:end))
    allocate(cps%z0m_ini(beg:end))  !mars
    allocate(cps%z0mg(beg:end))
    allocate(cps%z0hg(beg:end))
    allocate(cps%z0qg(beg:end))
    allocate(cps%htvp(beg:end))
    allocate(cps%beta(beg:end))
    allocate(cps%zii(beg:end))
    allocate(cps%albgrd(beg:end,numrad))
    allocate(cps%albgri(beg:end,numrad))
    allocate(cps%rootr_column(beg:end,nlevgrnd))
    allocate(cps%rootfr_road_perv(beg:end,nlevgrnd))
    allocate(cps%rootr_road_perv(beg:end,nlevgrnd))
    allocate(cps%wf(beg:end))
!   allocate(cps%xirrig(beg:end))
    allocate(cps%max_dayl(beg:end))
    allocate(cps%bsw2(beg:end,nlevgrnd))
    allocate(cps%psisat(beg:end,nlevgrnd))
    allocate(cps%vwcsat(beg:end,nlevgrnd))
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    allocate(cps%decl(beg:end))
    allocate(cps%coszen(beg:end))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%cannavg_t2m(beg:end))
    allocate(cps%me(beg:end))
    allocate(cps%fire_prob(beg:end))
    allocate(cps%mean_fire_prob(beg:end))
    allocate(cps%fireseasonl(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%ann_farea_burned(beg:end))
    allocate(cps%albsnd_hst(beg:end,numrad))
    allocate(cps%albsni_hst(beg:end,numrad))
    allocate(cps%albsod(beg:end,numrad))
    allocate(cps%albsoi(beg:end,numrad))
    allocate(cps%flx_absdv(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absdn(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absiv(beg:end,-nlevsno+1:1))
    allocate(cps%flx_absin(beg:end,-nlevsno+1:1))
    allocate(cps%snw_rds(beg:end,-nlevsno+1:0))
    allocate(cps%snw_rds_top(beg:end))
    allocate(cps%sno_liq_top(beg:end))
    allocate(cps%mss_bcpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bcphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bctot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_bc_col(beg:end))
    allocate(cps%mss_bc_top(beg:end))
    allocate(cps%mss_ocpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_ocphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_octot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_oc_col(beg:end))
    allocate(cps%mss_oc_top(beg:end))
    allocate(cps%mss_dst1(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst2(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst3(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst4(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dsttot(beg:end,-nlevsno+1:0))
    allocate(cps%mss_dst_col(beg:end))
    allocate(cps%mss_dst_top(beg:end))
    allocate(cps%h2osno_top(beg:end))
    allocate(cps%mss_cnc_bcphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_bcpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocphi(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_ocpho(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst1(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst2(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst3(beg:end,-nlevsno+1:0))
    allocate(cps%mss_cnc_dst4(beg:end,-nlevsno+1:0))
    allocate(cps%albgrd_pur(beg:end,numrad))
    allocate(cps%albgri_pur(beg:end,numrad))
    allocate(cps%albgrd_bc(beg:end,numrad))
    allocate(cps%albgri_bc(beg:end,numrad))
    allocate(cps%albgrd_oc(beg:end,numrad))
    allocate(cps%albgri_oc(beg:end,numrad))
    allocate(cps%albgrd_dst(beg:end,numrad))
    allocate(cps%albgri_dst(beg:end,numrad))
    allocate(cps%dTdz_top(beg:end))
    allocate(cps%snot_top(beg:end))
    allocate(cps%irrig_rate(beg:end))
    allocate(cps%n_irrig_steps_left(beg:end))
    allocate(cps%forc_pbot(beg:end))
    allocate(cps%forc_rho(beg:end))
    allocate(cps%glc_topo(beg:end))

    ival = 0.0
    cps%isoicol(beg:end) = huge(1)
    cps%ithincol(beg:end) = ival  ! mars
    cps%ithincol_i(beg:end) = ival ! mars
    cps%bsw(beg:end,1:nlevgrnd) = ival
    cps%watsat(beg:end,1:nlevgrnd) = ival
    cps%watfc(beg:end,1:nlevgrnd) = ival
    cps%watdry(beg:end,1:nlevgrnd) = ival
    cps%watopt(beg:end,1:nlevgrnd) = ival
    cps%hksat(beg:end,1:nlevgrnd) = ival
    cps%sucsat(beg:end,1:nlevgrnd) = ival
    cps%csol(beg:end,1:nlevgrnd) = ival
    cps%tkmg(beg:end,1:nlevgrnd) = ival
    cps%tkdry(beg:end,1:nlevgrnd) = ival
    cps%tksatu(beg:end,1:nlevgrnd) = ival
    cps%smpmin(beg:end) = ival
    cps%hkdepth(beg:end) = ival
    cps%wtfact(beg:end) = ival
    cps%fracice(beg:end,1:nlevgrnd) = ival
    cps%gwc_thr(beg:end) = ival
    cps%mss_frc_cly_vld(beg:end) = ival
    cps%mbl_bsn_fct(beg:end) = ival
    cps%do_capsnow (beg:end)= .false.
    cps%snowdp(beg:end) = ival
    cps%srf_emiss(beg:end) = ival   ! mars
    cps%co2dp(beg:end) = ival   ! mars
    cps%co2_mass_change(beg:end) = ival  ! mars
    cps%frac_co2(beg:end) = ival   ! mars
    cps%frac_sno(beg:end) = ival
    cps%zi(beg:end,-nlevsno+0:nlevgrnd) = ival
    cps%dz(beg:end,-nlevsno+1:nlevgrnd) = ival
    cps%z (beg:end,-nlevsno+1:nlevgrnd) = ival
    cps%frac_iceold(beg:end,-nlevsno+1:nlevgrnd) = spval
    cps%imelt(beg:end,-nlevsno+1:nlevgrnd) = huge(1)
    cps%eff_porosity(beg:end,1:nlevgrnd) = spval
    cps%emg(beg:end) = ival
    cps%z0m_ini(beg:end) = ival ! mars
    cps%z0mg(beg:end) = ival
    cps%z0hg(beg:end) = ival
    cps%z0qg(beg:end) = ival
    cps%htvp(beg:end) = ival
    cps%beta(beg:end) = ival
    cps%zii(beg:end) = ival
    cps%albgrd(beg:end,:numrad) = ival
    cps%albgri(beg:end,:numrad) = ival
    cps%rootr_column(beg:end,1:nlevgrnd) = spval
    cps%rootfr_road_perv(beg:end,1:nlevurb) = ival
    cps%rootr_road_perv(beg:end,1:nlevurb) = ival
    cps%wf(beg:end) = ival
!   cps%xirrig(beg:end) = 0._r8
    cps%bsw2(beg:end,1:nlevgrnd) = ival
    cps%psisat(beg:end,1:nlevgrnd) = ival
    cps%vwcsat(beg:end,1:nlevgrnd) = ival
    cps%soilpsi(beg:end,1:nlevgrnd) = spval
    cps%decl(beg:end) = ival
    cps%coszen(beg:end) = ival
    cps%fpi(beg:end) = ival
    cps%fpg(beg:end) = ival
    cps%annsum_counter(beg:end) = ival
    cps%cannsum_npp(beg:end) = ival
    cps%cannavg_t2m(beg:end) = ival
    cps%me(beg:end) = ival
    cps%fire_prob(beg:end) = ival
    cps%mean_fire_prob(beg:end) = ival
    cps%fireseasonl(beg:end) = ival
    cps%farea_burned(beg:end) = ival
    cps%ann_farea_burned(beg:end) = ival
    cps%albsnd_hst(beg:end,:numrad) = spval
    cps%albsni_hst(beg:end,:numrad) = spval
    cps%albsod(beg:end,:numrad) = ival
    cps%albsoi(beg:end,:numrad) = ival
    cps%flx_absdv(beg:end,-nlevsno+1:1) = spval
    cps%flx_absdn(beg:end,-nlevsno+1:1) = spval
    cps%flx_absiv(beg:end,-nlevsno+1:1) = spval
    cps%flx_absin(beg:end,-nlevsno+1:1) = spval
    cps%snw_rds(beg:end,-nlevsno+1:0) = ival
    cps%snw_rds_top(beg:end) = ival
    cps%sno_liq_top(beg:end) = ival
    cps%mss_bcpho(beg:end,-nlevsno+1:0) = ival
    cps%mss_bcphi(beg:end,-nlevsno+1:0) = ival
    cps%mss_bctot(beg:end,-nlevsno+1:0) = ival
    cps%mss_bc_col(beg:end) = ival
    cps%mss_bc_top(beg:end) = ival
    cps%mss_ocpho(beg:end,-nlevsno+1:0) = ival
    cps%mss_ocphi(beg:end,-nlevsno+1:0) = ival
    cps%mss_octot(beg:end,-nlevsno+1:0) = ival
    cps%mss_oc_col(beg:end) = ival
    cps%mss_oc_top(beg:end) = ival
    cps%mss_dst1(beg:end,-nlevsno+1:0) = ival
    cps%mss_dst2(beg:end,-nlevsno+1:0) = ival
    cps%mss_dst3(beg:end,-nlevsno+1:0) = ival
    cps%mss_dst4(beg:end,-nlevsno+1:0) = ival
    cps%mss_dsttot(beg:end,-nlevsno+1:0) = ival
    cps%mss_dst_col(beg:end) = ival
    cps%mss_dst_top(beg:end) = ival
    cps%h2osno_top(beg:end) = ival
    cps%mss_cnc_bcphi(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_bcpho(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_ocphi(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_ocpho(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_dst1(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_dst2(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_dst3(beg:end,-nlevsno+1:0) = ival
    cps%mss_cnc_dst4(beg:end,-nlevsno+1:0) = ival
    cps%albgrd_pur(beg:end,:numrad) = ival
    cps%albgri_pur(beg:end,:numrad) = ival
    cps%albgrd_bc(beg:end,:numrad) = ival
    cps%albgri_bc(beg:end,:numrad) = ival
    cps%albgrd_oc(beg:end,:numrad) = ival
    cps%albgri_oc(beg:end,:numrad) = ival 
    cps%albgrd_dst(beg:end,:numrad) = ival
    cps%albgri_dst(beg:end,:numrad) = ival
    cps%dTdz_top(beg:end) = ival
    cps%snot_top(beg:end) = ival
    cps%irrig_rate(beg:end) = ival
    cps%n_irrig_steps_left(beg:end) = 0
    cps%forc_pbot(beg:end) = ival
    cps%forc_rho(beg:end) = ival
    cps%glc_topo(beg:end) = ival

  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(ces%t_grnd(beg:end))
    allocate(ces%t_grnd_u(beg:end))
    allocate(ces%t_grnd_r(beg:end))
    allocate(ces%dt_grnd(beg:end))
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd))
    allocate(ces%t_soi_10cm(beg:end))
    allocate(ces%t_lake(beg:end,1:nlevlak))
    allocate(ces%tssbef(beg:end,-nlevsno+1:nlevgrnd))
    allocate(ces%thv(beg:end))
    allocate(ces%hc_soi(beg:end))
    allocate(ces%hc_soisno(beg:end))
    allocate(ces%forc_t(beg:end))
    allocate(ces%forc_th(beg:end))

    ival = 0.0
    ces%t_grnd(beg:end)    = ival
    ces%t_grnd_u(beg:end)  = ival
    ces%t_grnd_r(beg:end)  = ival
    ces%dt_grnd(beg:end)   = ival
    ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd) = spval
    ces%t_soi_10cm(beg:end) = spval
    ces%t_lake(beg:end,1:nlevlak)            = ival
    ces%tssbef(beg:end,-nlevsno+1:nlevgrnd)   = ival
    ces%thv(beg:end)       = ival
    ces%hc_soi(beg:end)    = ival
    ces%hc_soisno(beg:end) = ival
    ces%forc_t(beg:end) = ival
    ces%forc_th(beg:end) = ival

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cws%h2osno(beg:end))
    allocate(cws%errh2osno(beg:end))
    allocate(cws%snow_sources(beg:end))
    allocate(cws%snow_sinks(beg:end))
    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_ice(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cws%h2osoi_liqice_10cm(beg:end))
    allocate(cws%h2osoi_vol(beg:end,1:nlevgrnd))
    allocate(cws%h2osno_old(beg:end))
    allocate(cws%qg(beg:end))
    allocate(cws%dqgdT(beg:end))
    allocate(cws%snowice(beg:end))
    allocate(cws%snowliq(beg:end))
    allocate(cws%soilalpha(beg:end))
    allocate(cws%soilbeta(beg:end))
    allocate(cws%soilalpha_u(beg:end))
    allocate(cws%zwt(beg:end))
    allocate(cws%fcov(beg:end))
    allocate(cws%fsat(beg:end))
    allocate(cws%wa(beg:end))
    allocate(cws%wt(beg:end))
    allocate(cws%qcharge(beg:end))
    allocate(cws%smp_l(beg:end,1:nlevgrnd))
    allocate(cws%hk_l(beg:end,1:nlevgrnd))
    allocate(cws%forc_q(beg:end))

    ival = 0.0
    cws%h2osno(beg:end) = ival
    cws%errh2osno(beg:end) = ival
    cws%snow_sources(beg:end) = ival
    cws%snow_sinks(beg:end) = ival
    cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd)= spval
    cws%h2osoi_ice(beg:end,-nlevsno+1:nlevgrnd) = spval
    cws%h2osoi_liqice_10cm(beg:end) = spval
    cws%h2osoi_vol(beg:end,1:nlevgrnd) = spval
    cws%h2osno_old(beg:end) = ival
    cws%qg(beg:end) = ival
    cws%dqgdT(beg:end) = ival
    cws%snowice(beg:end) = ival
    cws%snowliq(beg:end) = ival
    cws%soilalpha(beg:end) = ival
    cws%soilbeta(beg:end) = ival
    cws%soilalpha_u(beg:end) = ival
    cws%zwt(beg:end) = ival
    cws%fcov(beg:end) = ival
    cws%fsat(beg:end) = ival
    cws%wa(beg:end) = ival
    cws%wt(beg:end) = ival
    cws%qcharge(beg:end) = ival
    cws%smp_l(beg:end,1:nlevgrnd) = spval
    cws%hk_l(beg:end,1:nlevgrnd) = spval
    cws%forc_q(beg:end) = ival

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ccs%soilc(beg:end))
    allocate(ccs%cwdc(beg:end))
    allocate(ccs%litr1c(beg:end))
    allocate(ccs%litr2c(beg:end))
    allocate(ccs%litr3c(beg:end))
    allocate(ccs%soil1c(beg:end))
    allocate(ccs%soil2c(beg:end))
    allocate(ccs%soil3c(beg:end))
    allocate(ccs%soil4c(beg:end))
    allocate(ccs%seedc(beg:end))
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%prod10c(beg:end))
    allocate(ccs%prod100c(beg:end))
    allocate(ccs%totprodc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

    ival = 0.0
    ccs%soilc(beg:end) = ival
    ccs%cwdc(beg:end) = ival
    ccs%litr1c(beg:end) = ival
    ccs%litr2c(beg:end) = ival
    ccs%litr3c(beg:end) = ival
    ccs%soil1c(beg:end) = ival
    ccs%soil2c(beg:end) = ival
    ccs%soil3c(beg:end) = ival
    ccs%soil4c(beg:end) = ival
    ccs%seedc(beg:end) = ival
    ccs%col_ctrunc(beg:end) = ival
    ccs%prod10c(beg:end) = ival
    ccs%prod100c(beg:end) = ival
    ccs%totprodc(beg:end) = ival
    ccs%totlitc(beg:end) = ival
    ccs%totsomc(beg:end) = ival
    ccs%totecosysc(beg:end) = ival
    ccs%totcolc(beg:end) = ival

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nstate_type
!
! !INTERFACE:
  subroutine init_column_nstate_type(beg, end, cns)
!
! !DESCRIPTION:
! Initialize column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cns%cwdn(beg:end))
    allocate(cns%litr1n(beg:end))
    allocate(cns%litr2n(beg:end))
    allocate(cns%litr3n(beg:end))
    allocate(cns%soil1n(beg:end))
    allocate(cns%soil2n(beg:end))
    allocate(cns%soil3n(beg:end))
    allocate(cns%soil4n(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%seedn(beg:end))
    allocate(cns%prod10n(beg:end))
    allocate(cns%prod100n(beg:end))
    allocate(cns%totprodn(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    ival = 0.0
    cns%cwdn(beg:end) = ival
    cns%litr1n(beg:end) = ival
    cns%litr2n(beg:end) = ival
    cns%litr3n(beg:end) = ival
    cns%soil1n(beg:end) = ival
    cns%soil2n(beg:end) = ival
    cns%soil3n(beg:end) = ival
    cns%soil4n(beg:end) = ival
    cns%sminn(beg:end) = ival
    cns%col_ntrunc(beg:end) = ival
    cns%seedn(beg:end) = ival
    cns%prod10n(beg:end) = ival
    cns%prod100n(beg:end) = ival
    cns%totprodn(beg:end) = ival
    cns%totlitn(beg:end) = ival
    cns%totsomn(beg:end) = ival
    cns%totecosysn(beg:end) = ival
    cns%totcoln(beg:end) = ival

  end subroutine init_column_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(beg, end, cef)
!
! !DESCRIPTION:
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cef%eflx_snomelt(beg:end))
    allocate(cef%eflx_snomelt_u(beg:end))
    allocate(cef%eflx_snomelt_r(beg:end))
    allocate(cef%eflx_impsoil(beg:end))
    allocate(cef%eflx_fgr12(beg:end))
    allocate(cef%eflx_building_heat(beg:end))
    allocate(cef%eflx_urban_ac(beg:end))
    allocate(cef%eflx_urban_heat(beg:end))
    allocate(cef%eflx_bot(beg:end))

    ival = 0.0
    cef%eflx_snomelt(beg:end)       = ival
    cef%eflx_snomelt_u(beg:end)       = ival
    cef%eflx_snomelt_r(beg:end)       = ival
    cef%eflx_impsoil(beg:end)       = ival
    cef%eflx_fgr12(beg:end)         = ival
    cef%eflx_building_heat(beg:end) = ival
    cef%eflx_urban_ac(beg:end) = ival
    cef%eflx_urban_heat(beg:end) = ival
    cef%eflx_bot(beg:end) = ival

  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cwf%qflx_infl(beg:end))
    allocate(cwf%qflx_surf(beg:end))
    allocate(cwf%qflx_drain(beg:end))
    allocate(cwf%qflx_top_soil(beg:end))
    allocate(cwf%qflx_sl_top_soil(beg:end))
    allocate(cwf%qflx_snomelt(beg:end))
    allocate(cwf%qflx_qrgwl(beg:end))
    allocate(cwf%qflx_runoff(beg:end))
    allocate(cwf%qflx_runoff_u(beg:end))
    allocate(cwf%qflx_runoff_r(beg:end))
    allocate(cwf%qmelt(beg:end))
    allocate(cwf%h2ocan_loss(beg:end))
    allocate(cwf%qflx_rsub_sat(beg:end))
    allocate(cwf%flx_bc_dep_dry(beg:end))
    allocate(cwf%flx_bc_dep_wet(beg:end))
    allocate(cwf%flx_bc_dep_pho(beg:end))
    allocate(cwf%flx_bc_dep_phi(beg:end))
    allocate(cwf%flx_bc_dep(beg:end))
    allocate(cwf%flx_oc_dep_dry(beg:end))
    allocate(cwf%flx_oc_dep_wet(beg:end))
    allocate(cwf%flx_oc_dep_pho(beg:end))
    allocate(cwf%flx_oc_dep_phi(beg:end))
    allocate(cwf%flx_oc_dep(beg:end))
    allocate(cwf%flx_dst_dep_dry1(beg:end))
    allocate(cwf%flx_dst_dep_wet1(beg:end))
    allocate(cwf%flx_dst_dep_dry2(beg:end))
    allocate(cwf%flx_dst_dep_wet2(beg:end))
    allocate(cwf%flx_dst_dep_dry3(beg:end))
    allocate(cwf%flx_dst_dep_wet3(beg:end))
    allocate(cwf%flx_dst_dep_dry4(beg:end))
    allocate(cwf%flx_dst_dep_wet4(beg:end))
    allocate(cwf%flx_dst_dep(beg:end))
    allocate(cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0))
    allocate(cwf%qflx_snofrz_col(beg:end))
    allocate(cwf%qflx_irrig(beg:end))
    allocate(cwf%qflx_glcice(beg:end))
    allocate(cwf%qflx_glcice_frz(beg:end))
    allocate(cwf%qflx_glcice_melt(beg:end))
    allocate(cwf%glc_rofi(beg:end))
    allocate(cwf%glc_rofl(beg:end))
    allocate(cwf%qflx_floodc(beg:end))
    allocate(cwf%qflx_snow_melt(beg:end))

    ival = 0.0
    cwf%qflx_infl(beg:end) = ival
    cwf%qflx_surf(beg:end) = ival
    cwf%qflx_drain(beg:end) = ival
    cwf%qflx_top_soil(beg:end) = spval
    cwf%qflx_sl_top_soil(beg:end) = ival
    cwf%qflx_snomelt(beg:end) = ival
    cwf%qflx_qrgwl(beg:end) = ival
    cwf%qflx_runoff(beg:end) = ival
    cwf%qflx_runoff_u(beg:end) = ival
    cwf%qflx_runoff_r(beg:end) = ival
    cwf%qmelt(beg:end) = ival
    cwf%h2ocan_loss(beg:end) = ival
    cwf%qflx_rsub_sat(beg:end) = ival
    cwf%flx_bc_dep_dry(beg:end) = ival
    cwf%flx_bc_dep_wet(beg:end) = ival
    cwf%flx_bc_dep_pho(beg:end) = ival
    cwf%flx_bc_dep_phi(beg:end) = ival
    cwf%flx_bc_dep(beg:end) = ival
    cwf%flx_oc_dep_dry(beg:end) = ival
    cwf%flx_oc_dep_wet(beg:end) = ival
    cwf%flx_oc_dep_pho(beg:end) = ival
    cwf%flx_oc_dep_phi(beg:end) = ival
    cwf%flx_oc_dep(beg:end) = ival
    cwf%flx_dst_dep_dry1(beg:end) = ival
    cwf%flx_dst_dep_wet1(beg:end) = ival
    cwf%flx_dst_dep_dry2(beg:end) = ival
    cwf%flx_dst_dep_wet2(beg:end) = ival
    cwf%flx_dst_dep_dry3(beg:end) = ival
    cwf%flx_dst_dep_wet3(beg:end) = ival
    cwf%flx_dst_dep_dry4(beg:end) = ival
    cwf%flx_dst_dep_wet4(beg:end) = ival
    cwf%flx_dst_dep(beg:end) = ival
    cwf%qflx_snofrz_lyr(beg:end,-nlevsno+1:0) = spval
    cwf%qflx_snofrz_col(beg:end) = ival
    cwf%qflx_irrig(beg:end)  = ival
    cwf%qflx_glcice(beg:end) = ival
    cwf%qflx_glcice_frz(beg:end) = ival
    cwf%qflx_glcice_melt(beg:end) = ival
    cwf%glc_rofi(beg:end)    = ival
    cwf%glc_rofl(beg:end)    = ival
    cwf%qflx_floodc(beg:end) = spval
    cwf%qflx_snow_melt(beg:end) = spval

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cflux_type
!
! !INTERFACE:
  subroutine init_column_cflux_type(beg, end, ccf)
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(ccf%m_leafc_to_litr1c(beg:end))
    allocate(ccf%m_leafc_to_litr2c(beg:end))
    allocate(ccf%m_leafc_to_litr3c(beg:end))
    allocate(ccf%m_frootc_to_litr1c(beg:end))
    allocate(ccf%m_frootc_to_litr2c(beg:end))
    allocate(ccf%m_frootc_to_litr3c(beg:end))
    allocate(ccf%m_leafc_storage_to_litr1c(beg:end))
    allocate(ccf%m_frootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_leafc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_frootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_to_cwdc(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc(beg:end))
    allocate(ccf%m_livecrootc_to_cwdc(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc(beg:end))
    allocate(ccf%m_gresp_storage_to_litr1c(beg:end))
    allocate(ccf%m_gresp_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc_fire(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc_fire(beg:end))
    allocate(ccf%hrv_leafc_to_litr1c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr2c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr3c(beg:end))             
    allocate(ccf%hrv_frootc_to_litr1c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr2c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr3c(beg:end))            
    allocate(ccf%hrv_livestemc_to_cwdc(beg:end))           
    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(ccf%hrv_livecrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_deadcrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_leafc_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_frootc_storage_to_litr1c(beg:end))    
    allocate(ccf%hrv_livestemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_deadstemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_gresp_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_leafc_xfer_to_litr1c(beg:end))        
    allocate(ccf%hrv_frootc_xfer_to_litr1c(beg:end))       
    allocate(ccf%hrv_livestemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_deadstemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_livecrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_gresp_xfer_to_litr1c(beg:end))        
    allocate(ccf%m_litr1c_to_fire(beg:end))
    allocate(ccf%m_litr2c_to_fire(beg:end))
    allocate(ccf%m_litr3c_to_fire(beg:end))
    allocate(ccf%m_cwdc_to_fire(beg:end))
    if ( crop_prog )then
       allocate(ccf%grainc_to_litr1c(beg:end))
       allocate(ccf%grainc_to_litr2c(beg:end))
       allocate(ccf%grainc_to_litr3c(beg:end))
       allocate(ccf%livestemc_to_litr1c(beg:end))
       allocate(ccf%livestemc_to_litr2c(beg:end))
       allocate(ccf%livestemc_to_litr3c(beg:end))
    end if
    allocate(ccf%leafc_to_litr1c(beg:end))
    allocate(ccf%leafc_to_litr2c(beg:end))
    allocate(ccf%leafc_to_litr3c(beg:end))
    allocate(ccf%frootc_to_litr1c(beg:end))
    allocate(ccf%frootc_to_litr2c(beg:end))
    allocate(ccf%frootc_to_litr3c(beg:end))
    allocate(ccf%cwdc_to_litr2c(beg:end))
    allocate(ccf%cwdc_to_litr3c(beg:end))
    allocate(ccf%litr1_hr(beg:end))
    allocate(ccf%litr1c_to_soil1c(beg:end))
    allocate(ccf%litr2_hr(beg:end))
    allocate(ccf%litr2c_to_soil2c(beg:end))
    allocate(ccf%litr3_hr(beg:end))
    allocate(ccf%litr3c_to_soil3c(beg:end))
    allocate(ccf%soil1_hr(beg:end))
    allocate(ccf%soil1c_to_soil2c(beg:end))
    allocate(ccf%soil2_hr(beg:end))
    allocate(ccf%soil2c_to_soil3c(beg:end))
    allocate(ccf%soil3_hr(beg:end))
    allocate(ccf%soil3c_to_soil4c(beg:end))
    allocate(ccf%soil4_hr(beg:end))
    if (use_cn) then
       allocate(ccf%dwt_seedc_to_leaf(beg:end))
       allocate(ccf%dwt_seedc_to_deadstem(beg:end))
       allocate(ccf%dwt_conv_cflux(beg:end))
       allocate(ccf%dwt_prod10c_gain(beg:end))
       allocate(ccf%dwt_prod100c_gain(beg:end))
       allocate(ccf%dwt_frootc_to_litr1c(beg:end))
       allocate(ccf%dwt_frootc_to_litr2c(beg:end))
       allocate(ccf%dwt_frootc_to_litr3c(beg:end))
       allocate(ccf%dwt_livecrootc_to_cwdc(beg:end))
       allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end))
       allocate(ccf%dwt_closs(beg:end))
       allocate(ccf%landuseflux(beg:end))
       allocate(ccf%landuptake(beg:end))
       allocate(ccf%prod10c_loss(beg:end))
       allocate(ccf%prod100c_loss(beg:end))
       allocate(ccf%product_closs(beg:end))
    end if
    allocate(ccf%lithr(beg:end))
    allocate(ccf%somhr(beg:end))
    allocate(ccf%hr(beg:end))
    allocate(ccf%sr(beg:end))
    allocate(ccf%er(beg:end))
    allocate(ccf%litfire(beg:end))
    allocate(ccf%somfire(beg:end))
    allocate(ccf%totfire(beg:end))
    allocate(ccf%nep(beg:end))
    allocate(ccf%nbp(beg:end))
    allocate(ccf%nee(beg:end))
    allocate(ccf%col_cinputs(beg:end))
    allocate(ccf%col_coutputs(beg:end))
    allocate(ccf%col_fire_closs(beg:end))

    if (use_cn) then
       allocate(ccf%cwdc_hr(beg:end))
       allocate(ccf%cwdc_loss(beg:end))
       allocate(ccf%litterc_loss(beg:end))
    end if

    ival = 0.0
    ccf%m_leafc_to_litr1c(beg:end)                = ival
    ccf%m_leafc_to_litr2c(beg:end)                = ival
    ccf%m_leafc_to_litr3c(beg:end)                = ival
    ccf%m_frootc_to_litr1c(beg:end)               = ival
    ccf%m_frootc_to_litr2c(beg:end)               = ival
    ccf%m_frootc_to_litr3c(beg:end)               = ival
    ccf%m_leafc_storage_to_litr1c(beg:end)        = ival
    ccf%m_frootc_storage_to_litr1c(beg:end)       = ival
    ccf%m_livestemc_storage_to_litr1c(beg:end)    = ival
    ccf%m_deadstemc_storage_to_litr1c(beg:end)    = ival
    ccf%m_livecrootc_storage_to_litr1c(beg:end)   = ival
    ccf%m_deadcrootc_storage_to_litr1c(beg:end)   = ival
    ccf%m_leafc_xfer_to_litr1c(beg:end)           = ival
    ccf%m_frootc_xfer_to_litr1c(beg:end)          = ival
    ccf%m_livestemc_xfer_to_litr1c(beg:end)       = ival
    ccf%m_deadstemc_xfer_to_litr1c(beg:end)       = ival
    ccf%m_livecrootc_xfer_to_litr1c(beg:end)      = ival
    ccf%m_deadcrootc_xfer_to_litr1c(beg:end)      = ival
    ccf%m_livestemc_to_cwdc(beg:end)              = ival
    ccf%m_deadstemc_to_cwdc(beg:end)              = ival
    ccf%m_livecrootc_to_cwdc(beg:end)             = ival
    ccf%m_deadcrootc_to_cwdc(beg:end)             = ival
    ccf%m_gresp_storage_to_litr1c(beg:end)        = ival
    ccf%m_gresp_xfer_to_litr1c(beg:end)           = ival
    ccf%m_deadstemc_to_cwdc_fire(beg:end)         = ival
    ccf%m_deadcrootc_to_cwdc_fire(beg:end)        = ival
    ccf%hrv_leafc_to_litr1c(beg:end)              = ival             
    ccf%hrv_leafc_to_litr2c(beg:end)              = ival             
    ccf%hrv_leafc_to_litr3c(beg:end)              = ival             
    ccf%hrv_frootc_to_litr1c(beg:end)             = ival            
    ccf%hrv_frootc_to_litr2c(beg:end)             = ival            
    ccf%hrv_frootc_to_litr3c(beg:end)             = ival            
    ccf%hrv_livestemc_to_cwdc(beg:end)            = ival           
    ccf%hrv_deadstemc_to_prod10c(beg:end)         = ival        
    ccf%hrv_deadstemc_to_prod100c(beg:end)        = ival       
    ccf%hrv_livecrootc_to_cwdc(beg:end)           = ival          
    ccf%hrv_deadcrootc_to_cwdc(beg:end)           = ival          
    ccf%hrv_leafc_storage_to_litr1c(beg:end)      = ival     
    ccf%hrv_frootc_storage_to_litr1c(beg:end)     = ival    
    ccf%hrv_livestemc_storage_to_litr1c(beg:end)  = ival 
    ccf%hrv_deadstemc_storage_to_litr1c(beg:end)  = ival 
    ccf%hrv_livecrootc_storage_to_litr1c(beg:end) = ival
    ccf%hrv_deadcrootc_storage_to_litr1c(beg:end) = ival
    if ( crop_prog )then
       ccf%grainc_to_litr1c(beg:end)    = ival
       ccf%grainc_to_litr2c(beg:end)    = ival
       ccf%grainc_to_litr3c(beg:end)    = ival
       ccf%livestemc_to_litr1c(beg:end) = ival
       ccf%livestemc_to_litr2c(beg:end) = ival
       ccf%livestemc_to_litr3c(beg:end) = ival
    end if
    ccf%hrv_gresp_storage_to_litr1c(beg:end)      = ival     
    ccf%hrv_leafc_xfer_to_litr1c(beg:end)         = ival        
    ccf%hrv_frootc_xfer_to_litr1c(beg:end)        = ival       
    ccf%hrv_livestemc_xfer_to_litr1c(beg:end)     = ival    
    ccf%hrv_deadstemc_xfer_to_litr1c(beg:end)     = ival    
    ccf%hrv_livecrootc_xfer_to_litr1c(beg:end)    = ival   
    ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end)    = ival   
    ccf%hrv_gresp_xfer_to_litr1c(beg:end)         = ival        
    ccf%m_litr1c_to_fire(beg:end)                 = ival
    ccf%m_litr2c_to_fire(beg:end)                 = ival
    ccf%m_litr3c_to_fire(beg:end)                 = ival
    ccf%m_cwdc_to_fire(beg:end)                   = ival
    ccf%leafc_to_litr1c(beg:end)                  = ival
    ccf%leafc_to_litr2c(beg:end)                  = ival
    ccf%leafc_to_litr3c(beg:end)                  = ival
    ccf%frootc_to_litr1c(beg:end)                 = ival
    ccf%frootc_to_litr2c(beg:end)                 = ival
    ccf%frootc_to_litr3c(beg:end)                 = ival
    ccf%cwdc_to_litr2c(beg:end)                   = ival
    ccf%cwdc_to_litr3c(beg:end)                   = ival
    ccf%litr1_hr(beg:end)                         = ival
    ccf%litr1c_to_soil1c(beg:end)                 = ival
    ccf%litr2_hr(beg:end)                         = ival
    ccf%litr2c_to_soil2c(beg:end)                 = ival
    ccf%litr3_hr(beg:end)                         = ival
    ccf%litr3c_to_soil3c(beg:end)                 = ival
    ccf%soil1_hr(beg:end)                         = ival
    ccf%soil1c_to_soil2c(beg:end)                 = ival
    ccf%soil2_hr(beg:end)                         = ival
    ccf%soil2c_to_soil3c(beg:end)                 = ival
    ccf%soil3_hr(beg:end)                         = ival
    ccf%soil3c_to_soil4c(beg:end)                 = ival
    ccf%soil4_hr(beg:end)                         = ival
    if (use_cn) then
       ccf%dwt_seedc_to_leaf(beg:end)                = ival
       ccf%dwt_seedc_to_deadstem(beg:end)            = ival
       ccf%dwt_conv_cflux(beg:end)                   = ival
       ccf%dwt_prod10c_gain(beg:end)                 = ival
       ccf%dwt_prod100c_gain(beg:end)                = ival
       ccf%dwt_frootc_to_litr1c(beg:end)             = ival
       ccf%dwt_frootc_to_litr2c(beg:end)             = ival
       ccf%dwt_frootc_to_litr3c(beg:end)             = ival
       ccf%dwt_livecrootc_to_cwdc(beg:end)           = ival
       ccf%dwt_deadcrootc_to_cwdc(beg:end)           = ival
       ccf%dwt_closs(beg:end)                        = ival
       ccf%landuseflux(beg:end)                      = ival
       ccf%landuptake(beg:end)                       = ival
       ccf%prod10c_loss(beg:end)                     = ival
       ccf%prod100c_loss(beg:end)                    = ival
       ccf%product_closs(beg:end)                    = ival
    end if
    ccf%lithr(beg:end)                            = ival
    ccf%somhr(beg:end)                            = ival
    ccf%hr(beg:end)                               = ival
    ccf%sr(beg:end)                               = ival
    ccf%er(beg:end)                               = ival
    ccf%litfire(beg:end)                          = ival
    ccf%somfire(beg:end)                          = ival
    ccf%totfire(beg:end)                          = ival
    ccf%nep(beg:end)                              = ival
    ccf%nbp(beg:end)                              = ival
    ccf%nee(beg:end)                              = ival
    ccf%col_cinputs(beg:end)                      = ival
    ccf%col_coutputs(beg:end)                     = ival
    ccf%col_fire_closs(beg:end)                   = ival

    if (use_cn) then
       ccf%cwdc_hr(beg:end)                          = ival
       ccf%cwdc_loss(beg:end)                        = ival
       ccf%litterc_loss(beg:end)                     = ival
    end if

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !USES:
    use surfrdMod , only : crop_prog
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cnf%ndep_to_sminn(beg:end))
    allocate(cnf%nfix_to_sminn(beg:end))
    allocate(cnf%m_leafn_to_litr1n(beg:end))
    allocate(cnf%m_leafn_to_litr2n(beg:end))
    allocate(cnf%m_leafn_to_litr3n(beg:end))
    allocate(cnf%m_frootn_to_litr1n(beg:end))
    allocate(cnf%m_frootn_to_litr2n(beg:end))
    allocate(cnf%m_frootn_to_litr3n(beg:end))
    allocate(cnf%m_leafn_storage_to_litr1n(beg:end))
    allocate(cnf%m_frootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_leafn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_frootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_to_cwdn(beg:end))
    allocate(cnf%m_deadstemn_to_cwdn(beg:end))
    allocate(cnf%m_livecrootn_to_cwdn(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%m_retransn_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_to_litr1n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr2n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr3n(beg:end))             
    allocate(cnf%hrv_frootn_to_litr1n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr2n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr3n(beg:end))            
    allocate(cnf%hrv_livestemn_to_cwdn(beg:end))           
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(cnf%hrv_livecrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_deadcrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_retransn_to_litr1n(beg:end))          
    allocate(cnf%hrv_leafn_storage_to_litr1n(beg:end))     
    allocate(cnf%hrv_frootn_storage_to_litr1n(beg:end))    
    allocate(cnf%hrv_livestemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_deadstemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_xfer_to_litr1n(beg:end))        
    allocate(cnf%hrv_frootn_xfer_to_litr1n(beg:end))       
    allocate(cnf%hrv_livestemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_deadstemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_livecrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%m_deadstemn_to_cwdn_fire(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn_fire(beg:end))
    allocate(cnf%m_litr1n_to_fire(beg:end))
    allocate(cnf%m_litr2n_to_fire(beg:end))
    allocate(cnf%m_litr3n_to_fire(beg:end))
    allocate(cnf%m_cwdn_to_fire(beg:end))
    if ( crop_prog )then
       allocate(cnf%grainn_to_litr1n(beg:end))
       allocate(cnf%grainn_to_litr2n(beg:end))
       allocate(cnf%grainn_to_litr3n(beg:end))
       allocate(cnf%livestemn_to_litr1n(beg:end))
       allocate(cnf%livestemn_to_litr2n(beg:end))
       allocate(cnf%livestemn_to_litr3n(beg:end))
    end if
    allocate(cnf%leafn_to_litr1n(beg:end))
    allocate(cnf%leafn_to_litr2n(beg:end))
    allocate(cnf%leafn_to_litr3n(beg:end))
    allocate(cnf%frootn_to_litr1n(beg:end))
    allocate(cnf%frootn_to_litr2n(beg:end))
    allocate(cnf%frootn_to_litr3n(beg:end))
    allocate(cnf%cwdn_to_litr2n(beg:end))
    allocate(cnf%cwdn_to_litr3n(beg:end))
    allocate(cnf%litr1n_to_soil1n(beg:end))
    allocate(cnf%sminn_to_soil1n_l1(beg:end))
    allocate(cnf%litr2n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_l2(beg:end))
    allocate(cnf%litr3n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_l3(beg:end))
    allocate(cnf%soil1n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_s1(beg:end))
    allocate(cnf%soil2n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_s2(beg:end))
    allocate(cnf%soil3n_to_soil4n(beg:end))
    allocate(cnf%sminn_to_soil4n_s3(beg:end))
    allocate(cnf%soil4n_to_sminn(beg:end))
    allocate(cnf%sminn_to_denit_l1s1(beg:end))
    allocate(cnf%sminn_to_denit_l2s2(beg:end))
    allocate(cnf%sminn_to_denit_l3s3(beg:end))
    allocate(cnf%sminn_to_denit_s1s2(beg:end))
    allocate(cnf%sminn_to_denit_s2s3(beg:end))
    allocate(cnf%sminn_to_denit_s3s4(beg:end))
    allocate(cnf%sminn_to_denit_s4(beg:end))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached(beg:end))
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    allocate(cnf%dwt_conv_nflux(beg:end))
    allocate(cnf%dwt_prod10n_gain(beg:end))
    allocate(cnf%dwt_prod100n_gain(beg:end))
    allocate(cnf%dwt_frootn_to_litr1n(beg:end))
    allocate(cnf%dwt_frootn_to_litr2n(beg:end))
    allocate(cnf%dwt_frootn_to_litr3n(beg:end))
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_nloss(beg:end))
    allocate(cnf%prod10n_loss(beg:end))
    allocate(cnf%prod100n_loss(beg:end))
    allocate(cnf%product_nloss(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))

    ival = 0.0
    cnf%ndep_to_sminn(beg:end) = ival
    cnf%nfix_to_sminn(beg:end) = ival
    cnf%m_leafn_to_litr1n(beg:end) = ival
    cnf%m_leafn_to_litr2n(beg:end) = ival
    cnf%m_leafn_to_litr3n(beg:end) = ival
    cnf%m_frootn_to_litr1n(beg:end) = ival
    cnf%m_frootn_to_litr2n(beg:end) = ival
    cnf%m_frootn_to_litr3n(beg:end) = ival
    cnf%m_leafn_storage_to_litr1n(beg:end) = ival
    cnf%m_frootn_storage_to_litr1n(beg:end) = ival
    cnf%m_livestemn_storage_to_litr1n(beg:end) = ival
    cnf%m_deadstemn_storage_to_litr1n(beg:end) = ival
    cnf%m_livecrootn_storage_to_litr1n(beg:end) = ival
    cnf%m_deadcrootn_storage_to_litr1n(beg:end) = ival
    cnf%m_leafn_xfer_to_litr1n(beg:end) = ival
    cnf%m_frootn_xfer_to_litr1n(beg:end) = ival
    cnf%m_livestemn_xfer_to_litr1n(beg:end) = ival
    cnf%m_deadstemn_xfer_to_litr1n(beg:end) = ival
    cnf%m_livecrootn_xfer_to_litr1n(beg:end) = ival
    cnf%m_deadcrootn_xfer_to_litr1n(beg:end) = ival
    cnf%m_livestemn_to_cwdn(beg:end) = ival
    cnf%m_deadstemn_to_cwdn(beg:end) = ival
    cnf%m_livecrootn_to_cwdn(beg:end) = ival
    cnf%m_deadcrootn_to_cwdn(beg:end) = ival
    cnf%m_retransn_to_litr1n(beg:end) = ival
    cnf%hrv_leafn_to_litr1n(beg:end) = ival             
    cnf%hrv_leafn_to_litr2n(beg:end) = ival             
    cnf%hrv_leafn_to_litr3n(beg:end) = ival             
    cnf%hrv_frootn_to_litr1n(beg:end) = ival            
    cnf%hrv_frootn_to_litr2n(beg:end) = ival            
    cnf%hrv_frootn_to_litr3n(beg:end) = ival            
    cnf%hrv_livestemn_to_cwdn(beg:end) = ival           
    cnf%hrv_deadstemn_to_prod10n(beg:end) = ival        
    cnf%hrv_deadstemn_to_prod100n(beg:end) = ival       
    cnf%hrv_livecrootn_to_cwdn(beg:end) = ival          
    cnf%hrv_deadcrootn_to_cwdn(beg:end) = ival          
    cnf%hrv_retransn_to_litr1n(beg:end) = ival          
    cnf%hrv_leafn_storage_to_litr1n(beg:end) = ival     
    cnf%hrv_frootn_storage_to_litr1n(beg:end) = ival    
    cnf%hrv_livestemn_storage_to_litr1n(beg:end) = ival 
    cnf%hrv_deadstemn_storage_to_litr1n(beg:end) = ival 
    cnf%hrv_livecrootn_storage_to_litr1n(beg:end) = ival
    cnf%hrv_deadcrootn_storage_to_litr1n(beg:end) = ival
    cnf%hrv_leafn_xfer_to_litr1n(beg:end) = ival        
    cnf%hrv_frootn_xfer_to_litr1n(beg:end) = ival       
    cnf%hrv_livestemn_xfer_to_litr1n(beg:end) = ival    
    cnf%hrv_deadstemn_xfer_to_litr1n(beg:end) = ival    
    cnf%hrv_livecrootn_xfer_to_litr1n(beg:end) = ival   
    cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end) = ival   
    cnf%m_deadstemn_to_cwdn_fire(beg:end) = ival
    cnf%m_deadcrootn_to_cwdn_fire(beg:end) = ival
    cnf%m_litr1n_to_fire(beg:end) = ival
    cnf%m_litr2n_to_fire(beg:end) = ival
    cnf%m_litr3n_to_fire(beg:end) = ival
    cnf%m_cwdn_to_fire(beg:end) = ival
    if ( crop_prog )then
       cnf%grainn_to_litr1n(beg:end)    = ival
       cnf%grainn_to_litr2n(beg:end)    = ival
       cnf%grainn_to_litr3n(beg:end)    = ival
       cnf%livestemn_to_litr1n(beg:end) = ival
       cnf%livestemn_to_litr2n(beg:end) = ival
       cnf%livestemn_to_litr3n(beg:end) = ival
    end if
    cnf%leafn_to_litr1n(beg:end) = ival
    cnf%leafn_to_litr2n(beg:end) = ival
    cnf%leafn_to_litr3n(beg:end) = ival
    cnf%frootn_to_litr1n(beg:end) = ival
    cnf%frootn_to_litr2n(beg:end) = ival
    cnf%frootn_to_litr3n(beg:end) = ival
    cnf%cwdn_to_litr2n(beg:end) = ival
    cnf%cwdn_to_litr3n(beg:end) = ival
    cnf%litr1n_to_soil1n(beg:end) = ival
    cnf%sminn_to_soil1n_l1(beg:end) = ival
    cnf%litr2n_to_soil2n(beg:end) = ival
    cnf%sminn_to_soil2n_l2(beg:end) = ival
    cnf%litr3n_to_soil3n(beg:end) = ival
    cnf%sminn_to_soil3n_l3(beg:end) = ival
    cnf%soil1n_to_soil2n(beg:end) = ival
    cnf%sminn_to_soil2n_s1(beg:end) = ival
    cnf%soil2n_to_soil3n(beg:end) = ival
    cnf%sminn_to_soil3n_s2(beg:end) = ival
    cnf%soil3n_to_soil4n(beg:end) = ival
    cnf%sminn_to_soil4n_s3(beg:end) = ival
    cnf%soil4n_to_sminn(beg:end) = ival
    cnf%sminn_to_denit_l1s1(beg:end) = ival
    cnf%sminn_to_denit_l2s2(beg:end) = ival
    cnf%sminn_to_denit_l3s3(beg:end) = ival
    cnf%sminn_to_denit_s1s2(beg:end) = ival
    cnf%sminn_to_denit_s2s3(beg:end) = ival
    cnf%sminn_to_denit_s3s4(beg:end) = ival
    cnf%sminn_to_denit_s4(beg:end) = ival
    cnf%sminn_to_denit_excess(beg:end) = ival
    cnf%sminn_leached(beg:end) = ival
    cnf%dwt_seedn_to_leaf(beg:end) = ival
    cnf%dwt_seedn_to_deadstem(beg:end) = ival
    cnf%dwt_conv_nflux(beg:end) = ival
    cnf%dwt_prod10n_gain(beg:end) = ival
    cnf%dwt_prod100n_gain(beg:end) = ival
    cnf%dwt_frootn_to_litr1n(beg:end) = ival
    cnf%dwt_frootn_to_litr2n(beg:end) = ival
    cnf%dwt_frootn_to_litr3n(beg:end) = ival
    cnf%dwt_livecrootn_to_cwdn(beg:end) = ival
    cnf%dwt_deadcrootn_to_cwdn(beg:end) = ival
    cnf%dwt_nloss(beg:end) = ival
    cnf%prod10n_loss(beg:end) = ival
    cnf%prod100n_loss(beg:end) = ival
    cnf%product_nloss(beg:end) = ival
    cnf%potential_immob(beg:end) = ival
    cnf%actual_immob(beg:end) = ival
    cnf%sminn_to_plant(beg:end) = ival
    cnf%supplement_to_sminn(beg:end) = ival
    cnf%gross_nmin(beg:end) = ival
    cnf%net_nmin(beg:end) = ival
    cnf%denit(beg:end) = ival
    cnf%col_ninputs(beg:end) = ival
    cnf%col_noutputs(beg:end) = ival
    cnf%col_fire_nloss(beg:end) = ival

  end subroutine init_column_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(beg, end, lps)
!
! !DESCRIPTION:
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(lps%t_building(beg:end))
    allocate(lps%t_building_max(beg:end))
    allocate(lps%t_building_min(beg:end))
    allocate(lps%tk_wall(beg:end,nlevurb))
    allocate(lps%tk_roof(beg:end,nlevurb))
    allocate(lps%tk_improad(beg:end,nlevgrnd))
    allocate(lps%cv_wall(beg:end,nlevurb))
    allocate(lps%cv_roof(beg:end,nlevurb))
    allocate(lps%cv_improad(beg:end,nlevgrnd))
    allocate(lps%thick_wall(beg:end))
    allocate(lps%thick_roof(beg:end))
    allocate(lps%nlev_improad(beg:end))
    allocate(lps%vf_sr(beg:end))
    allocate(lps%vf_wr(beg:end))
    allocate(lps%vf_sw(beg:end))
    allocate(lps%vf_rw(beg:end))
    allocate(lps%vf_ww(beg:end))
    allocate(lps%taf(beg:end))
    allocate(lps%qaf(beg:end))
    allocate(lps%sabs_roof_dir(beg:end,1:numrad))
    allocate(lps%sabs_roof_dif(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dir(beg:end,1:numrad))
    allocate(lps%sabs_sunwall_dif(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dir(beg:end,1:numrad))
    allocate(lps%sabs_shadewall_dif(beg:end,1:numrad))
    allocate(lps%sabs_improad_dir(beg:end,1:numrad))
    allocate(lps%sabs_improad_dif(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dir(beg:end,1:numrad))
    allocate(lps%sabs_perroad_dif(beg:end,1:numrad))

    ival = 0.0
    lps%t_building(beg:end) = ival
    lps%t_building_max(beg:end) = ival
    lps%t_building_min(beg:end) = ival
    lps%tk_wall(beg:end,1:nlevurb) = ival
    lps%tk_roof(beg:end,1:nlevurb) = ival
    lps%tk_improad(beg:end,1:nlevgrnd) = ival
    lps%cv_wall(beg:end,1:nlevurb) = ival
    lps%cv_roof(beg:end,1:nlevurb) = ival
    lps%cv_improad(beg:end,1:nlevgrnd) = ival
    lps%cv_improad(beg:end,1:5) = ival
    lps%thick_wall(beg:end) = ival
    lps%thick_roof(beg:end) = ival
    lps%nlev_improad(beg:end) = huge(1)
    lps%vf_sr(beg:end) = ival
    lps%vf_wr(beg:end) = ival
    lps%vf_sw(beg:end) = ival
    lps%vf_rw(beg:end) = ival
    lps%vf_ww(beg:end) = ival
    lps%taf(beg:end) = ival
    lps%qaf(beg:end) = ival
    lps%sabs_roof_dir(beg:end,1:numrad) = ival
    lps%sabs_roof_dif(beg:end,1:numrad) = ival
    lps%sabs_sunwall_dir(beg:end,1:numrad) = ival
    lps%sabs_sunwall_dif(beg:end,1:numrad) = ival
    lps%sabs_shadewall_dir(beg:end,1:numrad) = ival
    lps%sabs_shadewall_dif(beg:end,1:numrad) = ival
    lps%sabs_improad_dir(beg:end,1:numrad) = ival
    lps%sabs_improad_dif(beg:end,1:numrad) = ival
    lps%sabs_perroad_dir(beg:end,1:numrad) = ival
    lps%sabs_perroad_dif(beg:end,1:numrad) = ival

  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_eflux_type
!
! !INTERFACE:
  subroutine init_landunit_eflux_type(beg, end, lef)
!
! !DESCRIPTION: 
! Initialize landunit energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end 
    type (landunit_eflux_type), intent(inout):: lef 
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Keith Oleson
!
!EOP
!------------------------------------------------------------------------

    allocate(lef%eflx_traffic(beg:end))
    allocate(lef%eflx_traffic_factor(beg:end))
    allocate(lef%eflx_wasteheat(beg:end))
    allocate(lef%eflx_heat_from_ac(beg:end))

    ival = 0.0
    lef%eflx_traffic(beg:end) = ival
    lef%eflx_traffic_factor(beg:end) = ival
    lef%eflx_wasteheat(beg:end) = ival
    lef%eflx_heat_from_ac(beg:end) = ival

  end subroutine init_landunit_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%agdd20(beg:end))
    allocate(gps%tmomin20(beg:end))
    allocate(gps%t10min(beg:end))

    ival = 0.0
    gps%agdd20(beg:end) = ival
    gps%tmomin20(beg:end) = ival
    gps%t10min(beg:end) = ival

  end subroutine init_gridcell_dgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    
    
    !allocate(gps%bcphiwet2t(beg:end,1:2))
    !allocate(gps%bcphidry2t(beg:end,1:2))
    !allocate(gps%bcphodry2t(beg:end,1:2))
    !allocate(gps%ocphiwet2t(beg:end,1:2))
    !allocate(gps%ocphidry2t(beg:end,1:2))
    !allocate(gps%ocphodry2t(beg:end,1:2))
    !allocate(gps%dstx01wd2t(beg:end,1:2))
    !allocate(gps%dstx01dd2t(beg:end,1:2))
    !allocate(gps%dstx02wd2t(beg:end,1:2))
    !allocate(gps%dstx02dd2t(beg:end,1:2))
    !allocate(gps%dstx03wd2t(beg:end,1:2))
    !allocate(gps%dstx03dd2t(beg:end,1:2))
    !allocate(gps%dstx04wd2t(beg:end,1:2))
    !allocate(gps%dstx04dd2t(beg:end,1:2))
    
    !ival = 0.0
    !gps%bcphiwet2t(beg:end,1:2) = ival
    !gps%bcphidry2t(beg:end,1:2) = ival
    !gps%bcphodry2t(beg:end,1:2) = ival
    !gps%ocphiwet2t(beg:end,1:2) = ival
    !gps%ocphidry2t(beg:end,1:2) = ival
    !gps%ocphodry2t(beg:end,1:2) = ival
    !gps%dstx01wd2t(beg:end,1:2) = ival
    !gps%dstx01dd2t(beg:end,1:2) = ival
    !gps%dstx02wd2t(beg:end,1:2) = ival
    !gps%dstx02dd2t(beg:end,1:2) = ival
    !gps%dstx03wd2t(beg:end,1:2) = ival
    !gps%dstx03dd2t(beg:end,1:2) = ival
    !gps%dstx04wd2t(beg:end,1:2) = ival
    !gps%dstx04dd2t(beg:end,1:2) = ival

  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_efstate_type
!
! !INTERFACE:
  subroutine init_gridcell_efstate_type(beg, end, gve)
!
! !DESCRIPTION:
! Initialize gridcell isoprene emission factor variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_efstate_type), intent(inout) :: gve
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein (heald)
!
!EOP
!------------------------------------------------------------------------

    allocate(gve%efisop(6,beg:end))
    ival = 0.0
    gve%efisop(:,beg:end) = ival

  end subroutine init_gridcell_efstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(beg, end, gwf)
!
! !DESCRIPTION:
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gwf%qflx_runoffg(beg:end))
    allocate(gwf%qflx_snwcp_iceg(beg:end))
    allocate(gwf%qflx_liq_dynbal(beg:end))
    allocate(gwf%qflx_ice_dynbal(beg:end))
    allocate(gwf%qflx_floodg(beg:end))

    ival = 0.0
    gwf%qflx_runoffg(beg:end) = 0._r8
    gwf%qflx_snwcp_iceg(beg:end) = 0._r8
    gwf%qflx_liq_dynbal(beg:end) = ival
    gwf%qflx_ice_dynbal(beg:end) = ival
    gwf%qflx_floodg(beg:end) = 0._r8 !rtm_flood: initialize to zero for 1st time step instead of ival

  end subroutine init_gridcell_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_eflux_type
!
! !INTERFACE:
  subroutine init_gridcell_eflux_type(beg, end, gef)
!
! !DESCRIPTION:
! Initialize gridcell energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_eflux_type), intent(inout):: gef
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gef%eflx_sh_totg(beg:end))
    allocate(gef%eflx_dynbal(beg:end))

    ival = 0.0
    gef%eflx_sh_totg(beg:end) = ival
    gef%eflx_dynbal(beg:end) = ival

  end subroutine init_gridcell_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wstate_type
!
! !INTERFACE:
  subroutine init_gridcell_wstate_type(beg, end, gws)
!
! !DESCRIPTION:
! Initialize gridcell water state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wstate_type), intent(inout):: gws
    real(r8) :: ival
!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(gws%gc_liq1(beg:end))
    allocate(gws%gc_liq2(beg:end))
    allocate(gws%gc_ice1(beg:end))     
    allocate(gws%gc_ice2(beg:end))    

    ival = 0.0
    gws%gc_liq1(beg:end) = ival
    gws%gc_liq2(beg:end) = ival
    gws%gc_ice1(beg:end) = ival     
    gws%gc_ice2(beg:end) = ival    

  end subroutine init_gridcell_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_estate_type
!
! !INTERFACE:
  subroutine init_gridcell_estate_type(beg, end, ges)
!
! !DESCRIPTION:
! Initialize gridcell energy state variables     
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_estate_type), intent(inout):: ges
    real(r8) :: ival

!
! !REVISION HISTORY:
! Created by David Lawrence
!
!EOP
!------------------------------------------------------------------------
    allocate(ges%gc_heat1(beg:end))     
    allocate(ges%gc_heat2(beg:end))    

    ival = 0.0
    ges%gc_heat1(beg:end) = ival     
    ges%gc_heat2(beg:end) = ival    

  end subroutine init_gridcell_estate_type

end module clmtypeInitMod
