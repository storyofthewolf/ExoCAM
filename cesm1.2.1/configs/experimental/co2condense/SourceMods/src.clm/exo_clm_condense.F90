
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
  subroutine exo_clm_co2srfice(lbc, ubc, lbp, ubp, &
                               num_nolakec, filter_nolakec, &
                               num_condensep, filter_condensep)
!
! !DESCRIPTION:
! Calculation of
!
! !USES:

    use spmdMod,              only : masterproc
    use shr_kind_mod ,        only : r8 => SHR_KIND_R8
    use shr_const_mod,        only : SHR_CONST_CPICE, SHR_CONST_DENSITYCO2FR, SHR_CONST_CO2LATSUB, &
             		             SHR_CONST_PI, SHR_CONST_STEBOL, SHR_CONST_CDAY, SHR_CONST_SOILBD, &
         			     SHR_CONST_CO2_EMIS,SHR_CONST_SOILCP, SHR_CONST_G
    use shr_condense_mod,     only : get_condense_temp
    use clm_atmlnd,           only : clm_a2l
    use clmtype
    use clm_varcon,           only : tfrz, tkwat, tkice, &
                                     istice, istwet, istsoil, istdlak, istslak
    use clm_varpar,           only : nlevsno, nlevsoi, nlevgrnd
    use clm_time_manager,     only : get_step_size
    use subgridAveMod,        only : p2c
    use SoilTemperatureMod,   only : SoilThermProp
    use abortutils,           only : endrun

    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: lbc, ubc                     ! column bounds
    integer, intent(in) :: num_nolakec                  ! number of column non-lake points in column filter                                    
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)    ! column filter for non-lake points                      
    integer, intent(in) :: num_condensep                ! number of pfts in condense filter 
    integer, intent(in) :: filter_condensep(ubp-lbp+1)  ! pft filter for condense points

    !integer, intent(in) :: num_nolakec                  ! number of column non-lake points in column filter
    !integer, intent(in) :: filter_nolakec(ubc-lbc+1)    ! column filter for non-lake points
    !integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    !integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
    !integer, intent(in) :: num_nourbanc                ! number of columns in non-urban filter
    !integer, intent(in) :: filter_nourbanc(ubc-lbc+1)  ! column filter for non-urban points
    !integer, intent(in) :: num_nourbanp                ! number of pfts in non-urban filter 
    !integer, intent(in) :: filter_nourbanp(ubp-lbp+1)  ! pft filter for non-urban points

    !
    ! !CALLED FROM:
    ! subroutine clm_driver

    ! !LOCAL VARIABLES:
    !
    ! local pointers to original implicit in arrays
    !
    integer , pointer :: ltype(:)           ! landunit type 
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
    real(r8), pointer :: eflx_lh_co2(:)    ! co2 ice latent heat flux (W/m**2)  [+ to atm]  !Mars  
    real(r8), pointer :: eflx_soil_grnd(:) ! soil heat flux (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:) ! emitted IR radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:) ! net emitted IR radiation (W/m**2)
    real(r8), pointer :: z(:,:)            ! layer thickness (m)
    real(r8), pointer :: sabg(:)           ! solar radiation absorbed by ground (W/m**2)
    !
    ! local pointers to original, implicit inout arrays
    !
    real(r8), pointer :: snowdp(:)           ! snow height (m)
    real(r8), pointer :: co2dp(:)            ! CO2 frost mass (kg/m2)
    real(r8), pointer :: co2_mass_change(:)  ! CO2 frost mass change (kg/m2)
    real(r8), pointer :: co2_temp_change(:)  ! CO2 surface temperature change (K)

    real(r8), pointer :: forc_co2srf_snow(:)         ! [kg/m2]
    real(r8), pointer :: forc_co2srf_cond(:)         ! [kg/m3/s] atmospheric temperature (Kelvin)  {not used}

    real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)    
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)  
    real(r8), pointer :: t_lake(:,:)       ! lake temperature (Kelvin)
    real(r8), pointer :: co2ice_lh(:)      ! co2 latent heat flux energy change (W/m**2) [+ to grnd]

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
    integer  :: num				  ! work variable
    real(r8) :: dtime                             ! land model time step (sec)
    real(r8) :: tcond                             ! The temperature of the CO2 frost point (K) 
    real(r8) :: bifall                            ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: co2dp_ip	                  ! co2 ice mass from previous step + precip
    real(r8) :: co2dp_i	                          ! co2 ice mass from previous step
    real(r8) :: t_grnd_i	                  ! co2 ice mass from previous step
    real(r8) :: co2z			          ! co2 snow depth
    real(r8) :: co2_temp		          ! temporary co2 mass change
    real(r8) :: cv (lbc:ubc,-nlevsno+1:nlevgrnd)  !heat capacity [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsno+1:nlevgrnd)  !thermal conductivity [W/(m K)]
    real(r8) :: thick_layer1
    real(r8) :: t_layer1
    real(r8) :: tk_layer1
    real(r8) :: surf_cond, temp_change
    real(r8) :: tk_local(lbc:ubc)



    !-----------------------------------------------------------------------


    ! Assign local pointers to derived type members (gridcell-level)
    forc_t         => clm_a2l%forc_t
    forc_pbot	   => clm_a2l%forc_pbot  ! or possibly we should use pbot?
    forc_lwrad     => clm_a2l%forc_lwrad

    forc_co2srf_snow     => clm_a2l%forc_co2srf_snow
    forc_co2srf_cond     => clm_a2l%forc_co2srf_cond

    ! Assign local pointers to derived type members (landunit-level)
    clandunit          => col%landunit
    itype              => col%itype

    ! Assign local pointers to derived subtypes components (landunit-level)                                                         
    ltype      => lun%itype

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
    co2_temp_change    => ces%co2_temp_change
    z                  => cps%z
    t_soisno           => ces%t_soisno
    t_lake         => ces%t_lake
    
    ! Assign local pointers to derived type members (pft-level)
    pgridcell          => pft%gridcell
    plandunit          => pft%landunit
    pcolumn            => pft%column
    wtcol              => pft%wtcol
    eflx_sh_tot        => pef%eflx_sh_tot
    eflx_lh_tot        => pef%eflx_lh_tot
    eflx_lh_co2        => pef%eflx_lh_co2
    eflx_soil_grnd     => pef%eflx_soil_grnd
    eflx_lwrad_out     => pef%eflx_lwrad_out
    eflx_lwrad_net     => pef%eflx_lwrad_net
    co2ice_lh          => pef%co2ice_lh
    sabg               => pef%sabg
	
    ! Compute time step
    dtime = get_step_size()


    !
    ! Adapted from Richard Urata
    !

    ! Latent heat is the energy involved in a phase change at constant temperature, 
    ! and goes into breaking or recombining molecular bonds (of CO2 in this case).
    ! In the real world, as a phase change occurs the temperature of the system stays constant.
    ! However, in the discrete modeling world, the temperature changes between timesteps
    ! and then is adjusted back to the condensation temperature depending on condensation of sublimation.
    ! Latent heat exchange for the surface is accounted for here, as the temperature is directly changed.    


    ! Thermal conductivity and Heat capacity
    ! returns thermal conductivity and heat capacity of soil
    tk(:,:) = 0.0
    cv(:,:) = 0.0
    call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv)

    do f = 1,num_condensep
      p = filter_condensep(f)
      g = pgridcell(p)
      c = pcolumn(p) 
      l = plandunit(p)
       
      if (ltype(l) == istdlak .or. ltype(l) == istslak) then
        ! lake cells
        t_layer1 = t_lake(c,1)
        if (t_layer1 .gt. tfrz) then 
          tk_layer1 = tkwat  ! this is never used
        else
          tk_layer1 = tkice
        endif
      else
        ! soil/glacier cells
        t_layer1 = t_soisno(c,snl(c)+1)
        tk_layer1 = tk(c,snl(c)+1)
      endif

      ! calculate layer one thickness
      co2z = co2dp(c) / SHR_CONST_DENSITYCO2FR  ! calculate CO2 depth [m]
      if (co2dp(c) .eq. 0.d0) then
        thick_layer1= (abs(z(c,snl(c)+2)-z(c,snl(c)+1)))
      else if (co2dp(c) .gt. 0.d0) then
        thick_layer1 = (abs(z(c,snl(c)+2)-z(c,snl(c)+1))+co2z)
      endif
             

! is the lake layer1 thickness correct here
!      write(*,*) "ltype", ltype, thick_layer1, t_layer1, tk_layer1

      ! set grid CO2 ice quantities before calculation
      co2dp_i = co2dp(c)                        ! save initial co2 ice mass before accumulation of snow
      co2dp_ip = co2dp(c) + forc_co2srf_snow(c) ! save initial co2 ice mass + co2 snow fall   
      co2dp(c) = co2dp_ip
   
      t_grnd_i = t_grnd(c)                      ! save initial ground temperature
      co2_mass_change(c) = 0.d0                 ! surface ice mass change
      co2_temp_change(c) = 0.d0                 ! surface temperature change from CO2 cond/subl
      co2ice_lh(p) = 0.d0                       ! co2 latent heat release


     ! get condensation temperature at surface
      call get_condense_temp(1, forc_pbot(g), tcond)

      if (t_grnd_i .lt. tcond) then  ! CO2 condensation logic
        if (co2dp_ip .eq. 0.d0) then  ! CO2 condensing onto bare ground, no snow fall accumulated
          ! use soil thermal conductivity 
           co2ice_lh(p) = (-sabg(p) & 
                           + eflx_sh_tot(p) + (SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                           - tk_layer1*(t_layer1-tcond)/thick_layer1)
           co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
           co2dp(c) = co2_mass_change(c)
           t_grnd(c) = tcond
           co2_temp_change(c) = t_grnd(c) - t_grnd_i
           if (ltype(l) == istdlak .or. ltype(l) == istslak) then 
             t_lake(c,1) = tcond
           else 
             t_soisno(c,snl(c)+1) = tcond     
           endif
         endif
         if (co2dp_ip .gt. 0.d0) then ! CO2 condensing onto existing CO2 frost
           co2ice_lh(p) = (-sabg(p) &         
                           + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                           - tk_layer1*(t_layer1-tcond)/thick_layer1)
           co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
           co2dp(c) = co2dp(c) + co2_mass_change(c)
           t_grnd(c) = tcond
           co2_temp_change(c) = t_grnd(c) - t_grnd_i
           if (ltype(l) == istdlak .or. ltype(l) == istslak) then 
             t_lake(c,1) = tcond
           else
             t_soisno(c,snl(c)+1) = tcond
           endif
         endif
       endif  ! end condensation

      if ((t_grnd_i .gt. tcond) .and. (co2dp_ip .gt. 0.d0)) then  ! CO2 sublimation logic
         co2ice_lh(p) = (-sabg(p) &
                         + eflx_sh_tot(p) + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*tcond**4 - forc_lwrad(g)) &
                         - tk_layer1*(t_layer1-tcond)/thick_layer1)
         co2_mass_change(c) = (dtime/SHR_CONST_CO2LATSUB)*co2ice_lh(p)*wtcol(p)
         co2dp(c) = co2dp(c) + co2_mass_change(c)
         if (co2dp(c) .lt. 0.d0) then      ! more CO2 sublimation occured than ice was available
           co2_mass_change(c) = -co2dp_ip   ! set to CO2 existant mass
           co2ice_lh(p) = SHR_CONST_CO2LATSUB/dtime*co2_mass_change(c)
           temp_change = SHR_CONST_CO2LATSUB/(SHR_CONST_SOILBD*SHR_CONST_SOILCP*thick_layer1)*co2dp_ip
           t_grnd(c) = t_grnd(c) - temp_change
           co2_temp_change(c) = t_grnd(c) - t_grnd_i
           co2dp(c) = 0.d0
           if (ltype(l) == istdlak .or. ltype(l) == istslak) then 
              t_lake(c,1) = t_grnd(c)
           else 
             t_soisno(c,snl(c)+1) = t_grnd(c)
           endif
         else      ! CO2 ice remains on ground after sublimation occurs
           t_grnd(c) = tcond
           co2_temp_change(c) = t_grnd(c) - t_grnd_i
           if (ltype(l) == istdlak .or. ltype(l) == istslak) then 
              t_lake(c,1) = tcond
           else 
             t_soisno(c,snl(c)+1) = tcond
           endif
         endif
      else   ! no sublimation takes place, nothing changes?
        t_grnd(c) = t_grnd(c)
        if (ltype(l) == istdlak .or. ltype(l) == istslak) then
          t_lake(c,1) = t_grnd(c)
        else
          t_soisno(c,snl(c)+1) = t_grnd(c)
        endif
        co2_temp_change(c) = t_grnd(c) - t_grnd_i
      endif
      ! error trap small negative values which can sometimes occur
      if (co2dp(c) .lt. 0.d0) co2dp(c) = 0.d0
      ! set CO2 latent flux coupler array
      eflx_lh_co2(p) = -co2ice_lh(p)
    enddo

  end subroutine exo_clm_co2srfice

end module exo_clm_condense
