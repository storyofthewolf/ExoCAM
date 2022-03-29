!#include <misc.h>
!#include <preproc.h>

module CO2substepMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  CO2substep
!
! !DESCRIPTION:
! Substep calculation of
! (1) CO2 Condensation and sublimation at the surface
! output is the co2 depth and surface temperature after the substeps
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CO2substep
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
  subroutine CO2substep(dtime,sfcp,sfct,tsoi,solin,sh,lrad,tk,z1,z2,wt,co2d)
!
! !DESCRIPTION:
!
! !USES:

    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_const_mod, only : SHR_CONST_CPICE, SHR_CONST_DENSITYCO2FR, SHR_CONST_CO2LATSUB, &
                              SHR_CONST_PI, SHR_CONST_STEBOL, SHR_CONST_PI, SHR_CONST_CDAY, SHR_CONST_SOILBD, &
                              SHR_CONST_CO2_TK, SHR_CONST_CO2_EMIS, SHR_CONST_SOILCP, SHR_CONST_G 
    use shr_condense_mod, only: get_condense_temp
    use clmtype
    use clm_varcon   , only : tfrz, istice, istwet, istsoil
    use clm_varpar    , only : nlevsno, nlevsoi
    use clm_time_manager , only : get_step_size
    use subgridAveMod, only : p2c
    use SoilTemperatureMod, only : SoilThermProp
    !use nanMod
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: dtime                        ! time step length
    real(r8), intent(in) :: sfcp						 ! old sfc pressure (pa)
    real(r8), intent(inout) :: sfct						 ! old sfc temperature (k)
    real(r8), intent(inout) :: tsoi						 ! snow temperature
    real(r8), intent(in) :: solin						 ! solar flux (w/m2)
    real(r8), intent(in) :: sh							 ! sensible heat (w/m2)
    real(r8), intent(in) :: lrad						 ! incoming lw rad (w/m2)
    real(r8), intent(in) :: tk							 ! soil thermal conductivity
    real(r8), intent(in) :: z1							 ! soil level depth layer 1
    real(r8), intent(in) :: z2							 ! soil level depth layer 2
    real(r8), intent(in) :: wt							 ! column weight
    real(r8), intent(inout) :: co2d							 ! old co2 depth
    !
    ! !CALLED FROM:
    ! subroutine CO2FrostMod


    ! !LOCAL VARIABLES:
    !
    ! local pointers to original implicit in arrays
    !
    ! local pointers to original implicit inout arrays
    !
    !
    ! local pointers to original implicit out arrays
    !
    !
    !EOP
    !
    ! !OTHER LOCAL VARIABLES:
    !
    integer  :: n							 ! number of substeps
    integer  :: p							 ! loop index
    real(r8) :: dttmp						 ! substep length (dtime/n)
    real(r8) :: ptemp						 ! temporary sfc pressure
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: T_co2frostpoint, tcond       ! The temperature of the CO2 frost point (K) 
!    real(r8) :: Cp_co2_frost		     ! specific heat for solid CO2
    real(r8) :: old_co2_mass				! temporary co2 mass from last time step 
    real(r8) :: co2z						! co2 snow depth work var
    real(r8) :: change						! mass difference (kg/m2)
    real(r8) :: echangetmp						! substep heat difference (w/m2)
    !-----------------------------------------------------------------------

    ! do n number of substeps
    n=10
    ! find substep time length
    dttmp = dtime/float(n)
    ptemp=sfcp
		
    !-------------------Start CO2 stuff---------------------
    do p=1,n
      ! Compute CO2 frost point. Taken from the Mars Book (Kieffer et al 1992) 
      ! in chap 27 page 959 (chapter by James et al) : 
      T_co2frostpoint = 3182.48 / (23.3494 - log (0.01d0*ptemp))
      call get_condense_temp(1, sfcp, tcond)  ! new subroutine   

      ! Compute Cp for solid CO2:  (From Mars Book (Kieffer et al 1992) pg 958)
!      Cp_co2_frost= 349.0d0 + 4.8d0 * T_co2frostpoint
	
      ! reset negative co2dp values if they exist
      if (co2d < 0.d0) co2d = 0.d0
		
      ! reset mass change values	
      change = 0.d0
      old_co2_mass = 0.d0
		
      ! reset energy_change values (i.e. default is 0)	
      echangetmp = 0.d0
       
      ! If the ground T is below the frost point for CO2
      ! then condense CO2 onto the surface:
      if (sfct <= T_co2frostpoint) then
        if (co2d == 0.d0) then 	! case if there was no frost on the ground
	  echangetmp = SHR_CONST_SOILBD*SHR_CONST_SOILCP*(T_co2frostpoint-sfct)*(abs(z2-z1))
	  change = (dttmp/SHR_CONST_CO2LATSUB)*echangetmp*wt
          sfct = T_co2frostpoint
	  co2d = change
          tsoi = sfct
          ptemp=ptemp-change*shr_const_g
	else		! condensing co2 onto co2 frost
	  co2z = co2d / SHR_CONST_DENSITYCO2FR 
	  old_co2_mass = co2d
	  echangetmp = (-solin + & 
                        sh + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*T_co2frostpoint**4 - &
			lrad) - tk*(tsoi-T_co2frostpoint) &
			/(abs(z1-z2)+co2z))
          ! dz is the thickness of the first layer, and then added to the co2 thickness
	  change = (dttmp/SHR_CONST_CO2LATSUB)*echangetmp*wt
          co2d = co2d + change
          if (co2d < 0.d0) then
	    sfct = T_co2frostpoint - co2d*SHR_CONST_CO2LATSUB*sqrt((4.*SHR_CONST_PI/SHR_CONST_CDAY) &
	           /(tk*SHR_CONST_SOILCP*SHR_CONST_SOILBD))
            tsoi = sfct
	    co2d = 0.d0
	    change = -old_co2_mass
            ptemp=ptemp-change*SHR_CONST_G
          else
            sfct = T_co2frostpoint
  	    tsoi = sfct
 	    ptemp=ptemp-change*SHR_CONST_G
	  endif
	endif
      else if ((sfct > T_co2frostpoint) .and. (co2d > 0.d0)) then
        ! If the ground T is above the frost point for CO2
        ! then sublimate co2 off the ground	
	co2z = co2d / SHR_CONST_DENSITYCO2FR
	echangetmp = (-solin + & 
            	      sh + SHR_CONST_CO2_EMIS*(SHR_CONST_STEBOL*T_co2frostpoint**4 - &
		      lrad) - tk*(tsoi-T_co2frostpoint) &
		      /(abs(z1-z2)+co2z))
	change = (dttmp/SHR_CONST_CO2LATSUB)*echangetmp*wt
	old_co2_mass = co2d
	co2d = co2d + change
	if (co2d < 0.d0) then
	  sfct = T_co2frostpoint - co2d*SHR_CONST_CO2LATSUB*sqrt((4.*SHR_CONST_PI/SHR_CONST_CDAY) &
	         /(tk*SHR_CONST_SOILCP*SHR_CONST_SOILBD))
          tsoi = sfct
	  co2d = 0.d0
	  change = -old_co2_mass
	  ptemp=ptemp-change*SHR_CONST_G
	else
	  sfct = T_co2frostpoint
	  tsoi = sfct
	  ptemp=ptemp-change*SHR_CONST_G
	endif
      else
	sfct = sfct
	tsoi = sfct
      endif

    end do

  end subroutine CO2substep

end module CO2substepMod
