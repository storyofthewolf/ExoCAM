
module exo_clm_condense

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  exo_clm_condense
!
! !DESCRIPTION:
!  collection of CO2 surface ice, CO2 surface ice flux calculated from atmosphere
!  model and  passed to land model
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: exo_clm_co2srfice
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2Frost
!
! !INTERFACE:
  subroutine exo_clm_co2srfice(lbc, ubc, lbp, ubp, num_nolakec, filter_nolakec, &
                         num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Calculation of
!
! !USES:

    use spmdMod         ,only : masterproc
    use shr_kind_mod , only : r8 => SHR_KIND_R8
    use shr_const_mod, only : SHR_CONST_CPICE, SHR_CONST_DENSITYCO2FR, SHR_CONST_CO2LATSUB, &
             		      SHR_CONST_PI, SHR_CONST_STEBOL, SHR_CONST_CDAY, SHR_CONST_SOILBD, &
    			      SHR_CONST_CO2_TK, SHR_CONST_CO2_EMIS,SHR_CONST_SOILCP, SHR_CONST_G, SHR_CONST_SOILTK
    use shr_condense_mod,     only : get_condense_temp
    use clm_atmlnd,           only : clm_a2l
    use clmtype
    use clm_varcon,           only : tfrz, istice, istwet, istsoil
    use clm_varpar,           only : nlevsno, nlevsoi, nlevgrnd
    use clm_time_manager,     only : get_step_size
    use subgridAveMod,        only : p2c
    use SoilTemperatureMod,   only : SoilThermProp
    !use nanMod
    use CO2substepMod
    use abortutils,           only: endrun
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: lbc, ubc                     ! column bounds
    integer, intent(in) :: num_nolakec                  ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)    ! column filter for non-lake points
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
    !
    ! !CALLED FROM:
    ! subroutine biogeophys2

    ! !LOCAL VARIABLES:
    !
    ! local pointers to original implicit in arrays
    !
    integer , pointer :: cgridcell(:)      ! columns's gridcell
    integer , pointer :: clandunit(:)      ! columns's landunit
    integer , pointer :: pgridcell(:)      ! pft's gridcell
    integer , pointer :: plandunit(:)      ! pft's landunit
    integer , pointer :: pcolumn(:)        ! pft's column
    integer , pointer :: npfts(:)          ! number of pfts in column
    integer , pointer :: pfti(:)           ! column's beginning pft index
    integer , pointer :: itype(:)          ! landunit type
    real(r8), pointer :: wtcol(:)          ! weight relative to column
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_pbot(:)      ! atmospheric surface pressure  (Pa)
    real(r8), pointer :: forc_lwrad(:)     ! downward infrared (longwave) radiation (W/m**2)
    integer , pointer :: snl(:)            ! number of snow layers
    real(r8), pointer :: emg(:)            ! ground emissivity
    real(r8), pointer :: eflx_sh_tot(:)    ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)    ! total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_co2(:)    ! co2 latent heat flux (W/m**2)  [+ to atm] ! Mars
    real(r8), pointer :: eflx_soil_grnd(:) ! soil heat flux (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:) ! emitted IR radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:) ! net emitted IR radiation (W/m**2)
    real(r8), pointer :: co2ice_lh(:)      ! co2 latent heat flux energy change (W/m**2) [+ to grnd]
    real(r8), pointer :: z(:,:)            ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)    
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)  
    real(r8), pointer :: sabg(:)           ! solar radiation absorbed by ground (W/m**2)
    !
    ! local pointers to original implicit inout arrays
    !
    real(r8), pointer :: snowdp(:)           ! snow height (m)
    real(r8), pointer :: co2dp(:)            ! CO2 frost mass (kg/m2)
    real(r8), pointer :: co2_mass_change(:)  !  CO2 frost mass change (kg/m2)

    real(r8), pointer :: forc_co2srf_snow(:)         ! [kg/m2]
    real(r8), pointer :: forc_co2srf_cond(:)         ! [kg/m3/s] atmospheric temperature (Kelvin)  {not used}

    !
    ! local pointers to original implicit out arrays
    !
    !
    !EOP
    !
    ! !OTHER LOCAL VARIABLES:
    !
    integer  :: f                            ! filter index
    integer  :: c                            ! column index
    integer  :: l                            ! landunit index
    integer  :: g                            ! gridcell index
    integer  :: p                            ! pft index
    integer  :: num							 ! work variable
    real(r8) :: dtime                        ! land model time step (sec)
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: T_co2frostpoint, tcond       ! The temperature of the CO2 frost point (K) 
    real(r8) :: co2dp_ip	             ! co2 ice mass from previous step + precip
    real(r8) :: co2dp_i	             ! co2 ice mass from previous step
    real(r8) :: t_grnd_i	             ! co2 ice mass from previous step
    real(r8) :: co2z						! co2 snow depth
    real(r8) :: co2_temp					! temporary co2 mass change
    real(r8) :: cv (lbc:ubc,-nlevsno+1:nlevgrnd)  !heat capacity [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsno+1:nlevgrnd)  !thermal conductivity [W/(m K)]
    !real(r8), pointer :: latdeg(:)	  ! column latitude in degrees
    real(r8) :: layer1_thick

    real(r8) :: surf_cond, temp_change

    logical :: do_urata, do_forget

    !-----------------------------------------------------------------------
    !  choose method
    ! 
    do_urata = .true.
    do_forget = .false.

    ! Assign local pointers to derived type members (gridcell-level)
    forc_t         => clm_a2l%forc_t
    forc_pbot	   => clm_a2l%forc_pbot  ! or possibly we should use pbot?
    forc_lwrad     => clm_a2l%forc_lwrad

    forc_co2srf_snow     => clm_a2l%forc_co2srf_snow
    forc_co2srf_cond     => clm_a2l%forc_co2srf_cond

    !forc_t        => cps%forc_t
    !forc_pbot	   => cps%forc_pbot  ! or possibly we should use pbot?
    !forc_lwrad    => cps%forc_lwrad

    ! Assign local pointers to derived type members (landunit-level)
    clandunit          => col%landunit
    itype              => col%itype

    ! Assign local pointers to derived type members (column-level)
    cgridcell          => col%gridcell
    pfti               => col%pfti
    npfts              => col%npfts
    emg      	       => cps%emg
    t_grnd             => ces%t_grnd
    snl                => cps%snl
    snowdp             => cps%snowdp
    co2dp              => cps%co2dp
    co2_mass_change    => cps%co2_mass_change
    z                  => cps%z
    t_soisno           => ces%t_soisno
    !latdeg  => clm3%g%l%c%latdeg
   
    ! Assign local pointers to derived type members (pft-level)
    pgridcell          => pft%gridcell
    plandunit          => pft%landunit
    pcolumn            => pft%column
    wtcol              => pft%wtcol
    eflx_sh_tot        => pef%eflx_sh_tot
    eflx_lh_tot        => pef%eflx_lh_tot
    eflx_lh_co2        => pef%eflx_lh_co2  ! Mars     
    eflx_soil_grnd     => pef%eflx_soil_grnd
    eflx_lwrad_out     => pef%eflx_lwrad_out
    eflx_lwrad_net     => pef%eflx_lwrad_net
    co2ice_lh          => pef%co2ice_lh
    sabg               => pef%sabg
	
    ! Compute time step
    dtime = get_step_size()


! SHR_CONST_DENSITYCO2FR
! SHR_CONST_DENSITYCO2ICE kg/m3
! SHR_CONST_CPCO2ICE J kg-1 K-1

    !
    ! Adapted from Richard Urata
!     
if(do_urata) then 
!if(masterproc)write(*,*) "urata surface cond"
! presently this method yields ~100 Kg m2 at the notrh pole
! high terrain points have ever growing CO2 ice
! over high terrains this is due to the CO2 snowfall rates being much larger than the

! even with truncated mountain tops, still experiencing too CO2 over high terrains


    ! Thermal conductivity and Heat capacity
    ! returns thermal conductivity and heat capacity of soil
    call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv)
	
    do f = 1, num_nolakep  ! no lake? what if its a frozen lake cell
      p = filter_nolakep(f)
      g = pgridcell(p)
      c = pcolumn(p)
      l = plandunit(p)

      ! set grid CO2 ice quantities before calculation
      co2dp_i = co2dp(c)                        ! save initial co2 ice mass before accumulation of snow
      co2dp_ip = co2dp(c) + forc_co2srf_snow(c) ! save initial co2 ice mass + co2 snow fall   
      co2dp(c) = co2dp_ip
   
      t_grnd_i = t_grnd(c)                      ! save initial ground temperature
      co2_mass_change(c) = 0.d0                 ! surface ice mass change
      co2ice_lh(p) = 0.d0                       ! co2 latent heat release
      co2z = co2dp(c) / SHR_CONST_DENSITYCO2FR  ! calculate CO2 depth [m]

      if (co2dp(c) .eq. 0.d0) then
        layer1_thick = (abs(z(c,snl(c)+2)-z(c,snl(c)+1)))
      else if (co2dp(c) .gt. 0.d0) then
        layer1_thick = (abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z)
      endif

     ! get condensation temperature at surface
      call get_condense_temp(1, forc_pbot(g), tcond)

if(forc_pbot(g) .lt. 40000.) write(*,*) "---------------------------------------"
if(forc_pbot(g) .lt. 40000.) write(*,*) "forc_pbot(g), tcond", forc_pbot(g), tcond, co2dp_i, co2dp_ip, forc_co2srf_snow(c)

      if (t_grnd_i .lt. tcond) then  ! CO2 condensation logic
if(forc_pbot(g) .lt. 40000.) write(*,*) "start condensation", t_grnd_i
        if (co2dp_ip .eq. 0.d0) then  ! CO2 condensing onto bare ground, no snow fall accumulated
! old way
!           co2ice_lh(p) = SHR_CONST_SOILBD*SHR_CONST_SOILCP*(tcond-t_grnd(c))*layer1_thick/dtime
!
!           write(*,*) "co2dp_ip eq 0, 1", co2ice_lh(p)

! why don't we calculate using this way
!           co2ice_lh(p) = (-sabg(p) & 
!                           + eflx_sh_tot(p) + (SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
!                           - tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

           ! use soil thermal conductivity 
           co2ice_lh(p) = (-sabg(p) & 
                           + eflx_sh_tot(p) + (SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                           - SHR_CONST_SOILTK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

!write(*,*) "tk1", tk(c,snl(c)+1), SHR_CONST_SOILTK, SHR_CONST_SOILTK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick
           write(*,*) "co2dp_ip eq 0, 2", co2ice_lh(p)
!
           co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
           co2dp(c) = co2_mass_change(c)
           t_grnd(c) = tcond
           t_soisno(c,snl(c)+1) = tcond
         endif
         if (co2dp_ip .gt. 0.d0) then ! CO2 condensing onto existing CO2 frost

!           co2ice_lh(p) = (-sabg(p) & 
!                           + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
!                           - tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

    ! user CO2 ice thermal conductivity
           co2ice_lh(p) = (-sabg(p) & 
                           + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                           - SHR_CONST_CO2_TK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

!write(*,*) "tk 2", tk(c,snl(c)+1), SHR_CONST_CO2_TK, SHR_CONST_CO2_TK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick

if(forc_pbot(g) .lt. 40000.)write(*,*) "terms", -sabg(p), eflx_sh_tot(p), SHR_CONST_CO2_EMIS*SHR_CONST_STEBOL*tcond**4,forc_lwrad(g),tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick, tk(c,snl(c)+1), (t_soisno(c,snl(c)+1)-tcond), layer1_thick

           co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
           co2dp(c) = co2dp(c) + co2_mass_change(c)
           t_grnd(c) = tcond
           t_soisno(c,snl(c)+1) = t_grnd(c)
         endif
if(forc_pbot(g) .lt. 40000.) write(*,*) "co2dp1", co2dp_i, co2dp(c), forc_co2srf_snow(c), co2z, co2_mass_change(c)
       endif  ! end condensation

       if ((t_grnd_i .gt. tcond) .and. (co2dp_ip .gt. 0.d0)) then  ! CO2 sublimation logic
         ! sublimation occurs
if(forc_pbot(g) .lt. 40000.)write(*,*) "start sublimation"
         ! i don't think this logic works for snow on warm surfaces
         ! I think that teh sign of eflx_sh_tot should be reversed
!         co2ice_lh(p) = (-sabg(p) &
!                         + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
!                         - tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

         co2ice_lh(p) = (-sabg(p) &
                         + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                         - SHR_CONST_CO2_TK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick)

!write(*,*) "tk 3", tk(c,snl(c)+1), SHR_CONST_CO2_TK, SHR_CONST_CO2_TK*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick

if(forc_pbot(g) .lt. 40000.)write(*,*) "terms", -sabg(p), eflx_sh_tot(p), SHR_CONST_CO2_EMIS*SHR_CONST_STEBOL*tcond**4,forc_lwrad(g),tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond)/layer1_thick, tk(c,snl(c)+1), (t_soisno(c,snl(c)+1)-tcond), layer1_thick

         co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)

         co2dp(c) = co2dp(c) + co2_mass_change(c)

if(forc_pbot(g) .lt. 35000.) write(*,*) "co2dp2", co2dp_i, co2dp_ip, co2dp(c), forc_co2srf_snow(c), co2z, t_grnd(c), co2_mass_change(c)

         if (co2dp(c) .lt. 0.d0) then  ! more CO2 sublimation occured than ice was available

! this equation just blows up
!           t_grnd(c) = tcond - co2dp(c)*SHR_CONST_CO2LATSUB*sqrt((4.*SHR_CONST_PI/SHR_CONST_CDAY) &
!                       /(tk(c,snl(c)+1)*SHR_CONST_SOILCP*SHR_CONST_SOILBD))

! kludge logic
           co2_mass_change(c) = co2dp_ip
           temp_change = SHR_CONST_CO2LATSUB/(SHR_CONST_SOILBD*SHR_CONST_SOILCP*layer1_thick)*co2dp_i
           t_grnd(c) = t_grnd(c) - temp_change
           co2ice_lh(p) = SHR_CONST_CO2LATSUB/dtime*co2_mass_change(c)
!kludge logic

           t_soisno(c,snl(c)+1) = t_grnd(c)
           co2dp(c) = 0.d0
           co2_mass_change(c) = -co2dp_ip
if(forc_pbot(g) .lt. 35000.) write(*,*) "co2dp2 .lt. zero", t_grnd(c), co2dp(c), co2_mass_change(c)
         else      ! CO2 ice remains on ground after sublimation occurs
           t_grnd(c) = tcond
           t_soisno(c,snl(c)+1) = tcond
if(forc_pbot(g) .lt. 40000.) write(*,*) "ice remains after sublime", t_grnd(c), co2dp(c), co2_mass_change(c)
         endif
if(forc_pbot(g) .lt. 40000.) write(*,*) "co2dp2", co2dp_i, co2dp(c), forc_co2srf_snow(c), co2z, t_grnd(c)
      else   ! no sublimation takes place, nothing changes?
        t_grnd(c) = t_grnd(c)
        t_soisno(c,snl(c)+1) = t_grnd(c)
      endif
if(forc_pbot(g) .lt. 40000.) write(*,*) "end step", t_grnd_i, t_grnd(c), co2ice_lh(p), co2dp(c), co2dp_i, co2dp_ip

      ! error trap small negative values which can sometimes occur
      if (co2dp(c) .lt. 0.d0) co2dp(c) = 0.d0
    enddo

endif



    !
    ! simplified implementation and logic based on Forget
    !
if(do_forget) then 
if(masterproc)write(*,*) "forget surface cond"
! Presently the high terrains have low, acceptable amounts of CO2 ice and are stable
! However, the north polar ice exopentially unphysically grows
! This is caused by scaling the layer1_thick with co2z.  For cells as co2z grows, the condensation/sublimation rate enhances

! trying... limit the layer depth to a certain value for thick co2 ice, because it is unrealistic to think that the surface heat capacity extends 10s to 100s of meters in depth
! surface heat capacity J m-2 K-1  Soil bulk density * soil cp * thickness  ... but the layer thickness is somewhat arbitrary? no?
! this doesn't seem to work either, still have amplifying effect

!***** 
! with an arbitrary soil thickness, and the three high points truncated to 10000, I get good answers.


    do f = 1, num_nolakep     !nolakec?
       p = filter_nolakep(f)
       g = pgridcell(p)
       c = pcolumn(p)
       l = plandunit(p)

       ! set initial co2 snow mass [kg/m2], and depth [m]
       co2dp_i = co2dp(c)  
       t_grnd_i = t_grnd(c)
       co2z = co2dp(c) / SHR_CONST_DENSITYCO2FR  

!      if (co2dp(c) .eq. 0.d0) then
!         layer1_thick = (abs(z(c,snl(c)+2)-z(c,snl(c)+1)))
         layer1_thick = 0.2
!      else if (co2dp(c) .gt. 0.d0) then
!         layer1_thick = (abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z)
!      else if (co2dp(c) .gt. 1.d0) then
!         layer1_thick = (abs(z(c,snl(c)+2)-z(c,snl(c)+1))+1.0)  ! kludge to limit effective layer depth relevant to surface specific heat
!       endif

       ! get condensation temperature at surface 
       call get_condense_temp(1, forc_pbot(g), tcond) 

       if ( (t_grnd(c) .lt. tcond) .or.   &                                                ! condensation
         (forc_co2srf_snow(c) .gt. 1.0e-20)  .or.   &                                        ! co2 snow fall on dry surface
         ( (t_grnd(c) .gt. tcond) .and. (co2dp(c)+forc_co2srf_snow(c) .gt. 0.0)) ) then  ! sublimation

if((forc_pbot(g) .lt. 40000.))write(*,*) "--------------------------------- ", c, g
if((forc_pbot(g) .lt. 40000.))write(*,*) "co2dp, snow, p, tc, tg", co2dp(c), forc_co2srf_snow(c), forc_pbot(g), tcond, t_grnd(c)

         ! Condensation of partial sublimation of CO2 ice
         ! if t_grnd > tcond then surfcond_rate becomes negative
         surf_cond = SHR_CONST_SOILBD*SHR_CONST_SOILCP/SHR_CONST_CO2LATSUB*layer1_thick*(tcond-t_grnd(c)) ! [kg/m2]
         temp_change = (tcond-t_grnd(c))  ! 

if((forc_pbot(g) .lt. 40000.))write(*,*) "layer1_thick", layer1_thick
if((forc_pbot(g) .lt. 40000.))write(*,*) "abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z", z(c,snl(c)+2), z(c,snl(c)+1), snl(c),co2z
if((forc_pbot(g) .lt. 40000.))write(*,*) "initial calc,  surf cond, temp_change", surf_cond, temp_change

       ! If entire CO2 ice layer sublimes
         if ((co2dp(c)+forc_co2srf_snow(c)) .le. -surf_cond) then
!        if ((co2dp(c)) .le. -surf_cond) then
           surf_cond = -co2dp(c) - forc_co2srf_snow(c)
!           surf_cond = -co2dp(c) !- forc_co2srf_snow(c)
           temp_change = SHR_CONST_CO2LATSUB/(SHR_CONST_SOILBD*SHR_CONST_SOILCP*layer1_thick)*surf_cond
if((forc_pbot(g) .lt. 40000.))write(*,*) "everything sublimes", co2dp(c), forc_co2srf_snow(c), surf_cond, temp_change
         endif

       ! change CO2 ice amount
       ! CO2 depth equals the surface cond/subl rate + the snow fall rate
       co2dp(c) = co2dp(c) + forc_co2srf_snow(c) + surf_cond
       !co2dp(c) = co2dp(c) + surf_cond
       t_grnd(c) = t_grnd(c) + temp_change
       t_soisno(c,snl(c)+1) = t_grnd(c)
       co2_mass_change(c) = co2dp(c) - co2dp_i
       co2ice_lh(p) = SHR_CONST_CO2LATSUB/dtime*co2_mass_change(c)

if((forc_pbot(g) .lt. 40000.))write(*,*)"end logic", co2dp_i, co2dp(c), forc_co2srf_snow(c), co2_mass_change(c), co2ice_lh(p), t_grnd(c)

      endif
    enddo
endif


  end subroutine exo_clm_co2srfice

end module exo_clm_condense
