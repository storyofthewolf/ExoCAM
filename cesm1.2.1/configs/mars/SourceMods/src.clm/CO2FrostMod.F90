
module CO2FrostMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  CO2FrostMod
!
! !DESCRIPTION:
! Calculation of
! (1) CO2 Condensation and sublimation at the surface
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CO2Frost
!
! !REVISION HISTORY:
! Created by Richard Urata
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
  subroutine CO2Frost(lbc, ubc, lbp, ubp, num_nolakec, filter_nolakec, &
                         num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Calculation of
!
! !USES:

    use spmdMod         ,only : masterproc
    use shr_kind_mod , only : r8 => SHR_KIND_R8
    use SHR_CONST_mod, only : SHR_CONST_CPICE, SHR_CONST_DENSITYCO2FR, SHR_CONST_CO2LATSUB, &
             		      SHR_CONST_PI, SHR_CONST_STEBOL, SHR_CONST_CDAY, SHR_CONST_SOILBD, &
    			      SHR_CONST_CO2_TK, SHR_CONST_CO2_EMIS,SHR_CONST_SOILCP, SHR_CONST_G
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
    real(r8), pointer :: wtcol(:)         ! weight relative to column
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_pbot(:)      ! atmospheric surface pressure  (Pa)
    real(r8), pointer :: forc_lwrad(:)     !downward infrared (longwave) radiation (W/m**2)
    integer , pointer :: snl(:)             !number of snow layers
    real(r8), pointer :: emg(:)            !ground emissivity
    real(r8), pointer :: eflx_sh_tot(:)    !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)    ! total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_co2(:)    ! total latent heat flux (W/m**2)  [+ to atm]  ! Mars
    real(r8), pointer :: eflx_soil_grnd(:) ! soil heat flux (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:) ! emitted IR radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:) ! net emitted IR radiation (W/m**2)
    real(r8), pointer :: co2ice_lh(:)  !co2 latent heat flux energy change (W/m**2) [+ to grnd]
    real(r8), pointer :: z(:,:)            !layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)    
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)  
    real(r8), pointer :: sabg(:)           ! solar radiation absorbed by ground (W/m**2)
    !
    ! local pointers to original implicit inout arrays
    !
    real(r8), pointer :: snowdp(:)           ! snow height (m)
    real(r8), pointer :: co2dp(:)            ! CO2 frost mass (kg/m2)
    real(r8), pointer :: co2_mass_change(:)  !  CO2 frost mass change (kg/m2)

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
    real(r8) :: tcond                        ! The temperature of the CO2 frost point (K) 
!    real(r8) :: Cp_co2_frost		     ! specific heat for solid CO2
    real(r8) :: old_co2_mass				! temporary co2 mass from last time step 
    real(r8) :: co2z						! co2 snow depth
    real(r8) :: co2_temp					! temporary co2 mass change
    real(r8) :: cv (lbc:ubc,-nlevsno+1:nlevgrnd)  !heat capacity [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsno+1:nlevgrnd)  !thermal conductivity [W/(m K)]
    !real(r8), pointer :: latdeg(:)	  ! column latitude in degrees

    !-----------------------------------------------------------------------
    ! Assign local pointers to derived type members (gridcell-level)
    forc_t         => clm_a2l%forc_t
    forc_pbot	   => clm_a2l%forc_pbot  ! or possibly we should use pbot?
    forc_lwrad     => clm_a2l%forc_lwrad

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
     eflx_lh_co2        => pef%eflx_lh_co2  ! Mars (not currently used)
     eflx_soil_grnd     => pef%eflx_soil_grnd
     eflx_lwrad_out     => pef%eflx_lwrad_out
     eflx_lwrad_net     => pef%eflx_lwrad_net
     co2ice_lh      => pef%co2ice_lh
     sabg               => pef%sabg
	
     ! Compute time step
     dtime = get_step_size()

     ! Thermal conductivity and Heat capacity
     ! returns thermal conductivity and heat capacity of soil
     call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv)
	
     !-------------------Start CO2 surface condensation loop---------------------
     do f = 1, num_nolakep	!nolakec?	
       p = filter_nolakep(f)
       g = pgridcell(p)
       c = pcolumn(p)			
       l = plandunit(p)

       ! Compute CO2 frost point. Taken from the Mars Book (Kieffer et al 1992) 
       ! in chap 27 page 959 (chapter by James et al) : 

       call get_condense_temp(1, forc_pbot(g), tcond)  ! new subroutin

       if (tcond .le. 100.) then
         write(6,*)'co2 frost temp < 100 K',tcond,'p=',forc_pbot(g) !,'lat=',latdeg(c)
       endif

       ! Compute Cp for solid CO2:  (From Mars Book (Kieffer et al 1992) pg 958)
!       Cp_co2_frost= 349.0d0 + 4.8d0 * T_co2frostpoint
	
       ! reset negative co2dp values if they exist
       if (co2dp(c) < 0.d0) co2dp(c) = 0.d0
		
       ! reset mass change values	
       co2_mass_change(c) = 0.d0
       old_co2_mass = 0.d0
		
       ! reset energy_change values (i.e. default is 0)	
       co2ice_lh(p) = 0.d0
       
       ! If the ground T is below the frost point for CO2
       ! then condense CO2 onto the surface:
       if (t_grnd(c) .lt. tcond) then

         if (co2dp(c) .eq. 0.d0) then 	! case if there was no frost on the ground

           ! the units are wrong here.  this comes out to Jm-2   where as Wm-2 = Js-1m-2
           co2ice_lh(p) = SHR_CONST_SOILBD*SHR_CONST_SOILCP*(tcond-t_grnd(c))*(abs(z(c,snl(c)+2)-z(c,snl(c)+1)))/dtime

           co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)

!write(*,*) "wtcol", wtcol(p)
! commented out CO2Substep loop
!           if (co2_mass_change(c) > (forc_pbot(g)/(50.*SHR_CONST_g))) then
!             num=0
!             do while ((co2_mass_change(c) > (forc_pbot(g)/(50.*SHR_CONST_G))) .and. (num .le. 5))
!               num=num+1
!               call CO2substep(dtime,forc_pbot(g),t_grnd(c),t_soisno(c,snl(c)+1), &
!	                       sabg(p),eflx_sh_tot(p),forc_lwrad(g),tk(c,snl(c)+1), &
!			       z(c,snl(c)+1),z(c,snl(c)+2),wtcol(p),co2dp(c))
!			       co2_mass_change(c) = co2dp(c)
!             end do				
!             write(6,*) 'condensed more than sfc_p of CO2. called CO2 substep',num,'times'
!           else
!\ commented out
             co2dp(c) = co2_mass_change(c)
             t_grnd(c) = tcond
!           endif  !CO2substep
           t_soisno(c,snl(c)+1) = t_grnd(c) 
	 else		! if condensing CO2 onto CO2 frost
           co2z = co2dp(c) / SHR_CONST_DENSITYCO2FR
           old_co2_mass = co2dp(c)
           ! I am confused sbout this "energy_change"

           ! - absorbed solar by ground + sensible heat tp atm + thermal emissions - downwelling LW - thermal conductivity * soil temp delta / depth
           co2ice_lh(p) = (-sabg(p) + & ! correct units Wm-2
                               eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - &
                               forc_lwrad(g)) - tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond) &
                               /(abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z))
				! dz is the thickness of the first layer, and then added to the co2 thickness
                                !co2 thickness added to the thickness of the first layer.  This would work well in suppressing CO2 condensing onto CO2

                                 
	   co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
! commented CO2substep loop
!           if (co2_mass_change(c) > (forc_pbot(g)/(50.*SHR_CONST_G))) then
!             num=0
!             do while ((co2_mass_change(c) > (forc_pbot(g)/(50.*SHR_CONST_G))) .and. (num .le. 5))
!               num=num+1
!               call CO2substep(dtime,forc_pbot(g),t_grnd(c),t_soisno(c,snl(c)+1), &
!                               sabg(p),eflx_sh_tot(p),forc_lwrad(g),tk(c,snl(c)+1), &
!                               z(c,snl(c)+1),z(c,snl(c)+2),wtcol(p),co2dp(c))
!                               co2_mass_change(c) = co2dp(c) - old_co2_mass
!             end do
!             write(6,*) 'condensed more than sfc_p of CO2. called CO2 substep',num,'times'
!           else
!\ commented
            co2dp(c) = co2dp(c) + co2_mass_change(c)
!           endif  !c

	   if (co2dp(c) < 0.d0) then  ! if co2 ice depth is negative???/ 
write(*,*) "co2dp negative" , co2dp(c)
           ! how can this ever be reached?  condition one is condensation temperature limit.  t_soisno > tcond, but tgrnd < tcond
             t_grnd(c) = tcond - co2dp(c)*SHR_CONST_CO2LATSUB*sqrt((4.*SHR_CONST_PI/SHR_CONST_CDAY) &
                         /(tk(c,snl(c)+1)*SHR_CONST_SOILCP*SHR_CONST_SOILBD))
             t_soisno(c,snl(c)+1) = t_grnd(c)
	     co2dp(c) = 0.d0
             co2_mass_change(c) = -old_co2_mass
           else
             t_grnd(c) = tcond
             t_soisno(c,snl(c)+1) = t_grnd(c)
           endif
	 endif


       ! If the ground T is above the frost point for CO2
       ! then sublimate co2 off the ground	
       else if ((t_grnd(c) > tcond) .and. (co2dp(c) > 0.d0)) then

         co2z = co2dp(c) / SHR_CONST_DENSITYCO2FR
         co2ice_lh(p) = (-sabg(p) + & 
                             eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - &
                             forc_lwrad(g)) - tk(c,snl(c)+1)*(t_soisno(c,snl(c)+1)-tcond) &
                             /(abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z))
         co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
	 old_co2_mass = co2dp(c)
! comment CO2 Substep loop
!	 if (co2_mass_change(c) > (forc_pbot(g)/(5.*SHR_CONST_G))) then
!	   write(6,*) 'condensing more than .2*sfc_p. sfc_p=',forc_pbot(g),'condense=',co2_mass_change(c)*SHR_CONST_G
!           call CO2substep(dtime,forc_pbot(g),t_grnd(c),t_soisno(c,snl(c)+1), &
!                           sabg(p),eflx_sh_tot(p),forc_lwrad(g),tk(c,snl(c)+1), &
!                           z(c,snl(c)+1),z(c,snl(c)+2),wtcol(p),co2dp(c))
!         else
! \comment
           co2dp(c) = co2dp(c) + co2_mass_change(c)
!         endif !c
         if (co2dp(c) < 0.d0) then  ! if co2 ice depth is negative ???
write(*,*) "co2dp negative 2", co2dp(c)
           t_grnd(c) = tcond - co2dp(c)*SHR_CONST_CO2LATSUB*sqrt((4.*SHR_CONST_pi/SHR_CONST_CDAY) &
                       /(tk(c,snl(c)+1)*SHR_CONST_SOILCP*SHR_CONST_SOILBD))
           t_soisno(c,snl(c)+1) = t_grnd(c)
           co2dp(c) = 0.d0
           co2_mass_change(c) = -old_co2_mass
         else
           t_grnd(c) = tcond
           t_soisno(c,snl(c)+1) = t_grnd(c)
         endif	
       else  ! Not sublimating CO2, Not condensing CO2. 
         t_grnd(c) = t_grnd(c)
         t_soisno(c,snl(c)+1) = t_grnd(c)
       endif
    end do

  end subroutine CO2Frost

end module CO2FrostMod
