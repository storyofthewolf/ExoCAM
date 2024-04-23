!! This CARMA model is for Titan/early Earth-like fractal organic hazes 
!! based upon Wolf & Toon (2010)
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!
!! @version Jan-2011
!! @author  Chuck Bardeen 
!! @modified by Eric T Wolf 2024
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
  use abortutils,     only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc
  use phys_grid,      only: get_rlat_all_p, get_rlon_all_p

#if ( defined SPMD )
  use mpishorthand
#endif  

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition

  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 1               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 40              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 1               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_ORGHAZE   = 1       !! organic hazes

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_ORGHAZE     = 1  !! organic hazes

  integer, public, parameter      :: I_ELEM_ORGHAZE    = 1     !! organic hazes

  ! These variables are all set during initialization and are used to calculate
  ! emission tendencies.
  integer                             :: carma_emis_nLevs        ! number of emission levels
  real(r8), allocatable, dimension(:) :: carma_emis_lev          ! emission levels (Pa)
  integer                             :: carma_emis_ilev_min     ! index of minimum level in table 
  integer                             :: carma_emis_ilev_max     ! index of maximum level in table 
  integer                             :: carma_emis_ilev_incr    ! index increment to increase level 
  real(r8)                            :: carma_emis_expected     ! Expected emission rate per column (kg/m2/s)
 
  integer                             :: carma_emis_nZens        ! number of zenith angles 
  real(r8), allocatable, dimension(:) :: carma_emis_zen          ! zenith angles (degrees)
  integer                             :: carma_emis_izen_min     ! index of minimum level in table 
  integer                             :: carma_emis_izen_max     ! index of maximum level in table 
  integer                             :: carma_emis_izen_incr    ! index increment to increase level 

  real(r8), allocatable, dimension(:,:) :: carma_emis_rate         ! emission rate lookup table (g cm-3 s-1), zenith angle dependence

  real(r8), parameter  :: illumFrac = 1.9250898                  ! illumination fraction, used for scaling photochemically driven aerosol mass
                                                                 ! value of 1.9250898 scales mass when source only on daylight hemisphere
                                                                 ! value is =/= 2.0 due to grid boxes at 90 deg zenith angles being removed

contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!  @modified Eric T Wolf 2024
  subroutine CARMA_DefineModel(carma, rc)
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: RHO_ORGHAZE = 0.64_f ! density of photochemical haze particles (g/cm3 ) Trainer+ (2016)
    real(kind=f), parameter            :: rmin        = 1e-7_f    ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat       = 2.5_f     ! volume/mass ratio
    real(kind=f), parameter            :: rmon_in     = 50.e-7_f  ! monomer size
    real(kind=f), parameter            :: falpha_in   = 1.0      ! fractal packing coefficient
    real(kind=f), allocatable          :: df_in(:)
    real(kind=f), allocatable          :: nmon(:)

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    
    ! Default return code.
    rc = RC_OK
    
    allocate(df_in(NBIN))     ! 40 bin, from Wolf & Toon (2010)
    df_in = (/ &              ! bin dependent fractal dimensions 
      3.00000,      3.00000,      3.00000,      3.00000,      3.00000, &
      3.00000,      3.00000,      3.00000,      3.00000,      3.00000, &
      3.00000,      3.00000,      3.00000,      1.50214,      1.50535, &
      1.51331,      1.53291,      1.58003,      1.68694,      1.89714, &
      2.18998,      2.37633,      2.39990,      2.40000,      2.40000, &
      2.40000,      2.40000,      2.40000,      2.40000,      2.40000, &
      2.40000,      2.40000,      2.40000,      2.40000,      2.40000, &
      2.40000,      2.40000,      2.40000,      2.40000,      2.40000 /)

    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_do_escale  = ', carma_do_escale
      if (do_print) write(LUNOPRT,*) '  carma_emis_total = ', carma_emis_total
      if (do_print) write(LUNOPRT,*) '  carma_emis_file  = ', trim(carma_emis_file)
      if (do_print) write(LUNOPRT,*) '  carma_escale_file= ', trim(carma_escale_file)
    end if
    
    
    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    call CARMAGROUP_Create(carma, I_GRP_ORGHAZE, "fractal haze", rmin, vmrat, I_SPHERE, 1._f, .false., &
                          rc, shortname="HAZE", is_fractal=.true., do_wetdep=.true., do_drydep=.true.,  &
                          rmon=rmon_in, df=df_in, falpha=falpha_in)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_ORGHAZE, I_GRP_ORGHAZE, "photochemical haze", RHO_ORGHAZE, I_INVOLATILE, I_ORGHAZE, rc, shortname="HAZE")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
    
    
    ! Define the Solutes
    
    
    ! Define the Gases
    
    
    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_ORGHAZE, I_GRP_ORGHAZE, I_GRP_ORGHAZE, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')
    
    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
     tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none

    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)

    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none

    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins


  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version Jan-2011
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use camsrfexch,       only: cam_in_t
    use time_manager,  only: get_curr_calday, is_perpetual, get_perp_date, get_curr_date
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    integer                            :: ilat                  ! latitude index 
    integer                            :: iltime                ! local time index
    integer                            :: lchnk                 ! chunk identifier
    integer                            :: ncol                  ! number of columns in chunk
    integer                            :: icol                  ! column index
    integer                            :: igroup                ! the index of the carma aerosol group
    integer                            :: k                     ! vertical index
    integer                            :: ilev                  ! level index in emissions data
    integer                            :: izen                  ! zenith index in emissions data
    character(len=32)                  :: shortname             ! the shortname of the group
    real(r8)                           :: r(NBIN)               ! bin center
    real(r8)                           :: dr(NBIN)              ! bin width
    real(r8)                           :: rmass(NBIN)           ! bin mass
    real(r8)                           :: pressure              ! pressure (Pa)
    real(r8)                           :: cosz                  ! cosine of the zenith angle
    real(r8)                           :: thickness             ! layer thickness (m)
    real(r8)                           :: rate                  ! emission rate (#/cm-3/s)
    real(r8)                           :: massflux              ! emission mass flux (kg/m2/s)
    real(r8)                           :: columnMass            ! mass of the total column (kg/m2/s)
    real(r8)                           :: scale                 ! scaling factor to conserve the expected mass

    real(r8)                           :: calday                ! current calendar day
    integer                            :: yr, mon, day, ncsec, doy
    integer                            :: ncdate
    real(r8)                           :: ltime                 ! local time
    real(r8)                           :: clat(pcols)           ! current latitudes(radians)
    real(r8)                           :: clon(pcols)           ! current longitudes(radians)
    real(r8)                           :: coszrs(pcols)         ! Cosine solar zenith angle
    


    ! variables for bilinear interpolation
    integer :: p_ref_index, p_ref_indexp1  ! pressure indices
    integer :: z_ref_index, z_ref_indexp1  ! zenith angle indices
    real(r8) :: press, onemp
    real(r8) :: zen,   onemz
    real(r8), dimension(4) :: vbi


    ! Default return code.
    rc = RC_OK

    ! Get the current date and time.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)

    ! NOTE: The global relative flux file is based upon a noleap calendar,
    ! so don't let the doy get too big.
    doy = min(365, doy)

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    ! Cosine solar zenith angle for current time step
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8


    ! Only do emission for the first bin of the meteor smoke group.
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    ! For photochemical organic hazes, the source from the haze only goes into the
    ! smallest bin, being less than or equal to the monomer size
    if ((shortname .eq. "HAZE") .and. (ibin .eq. 1)) then

      ! Set tendencies for any sources or sinks in the atmosphere.
      do k = 1, pver
        do icol = 1, ncol


          pressure = state%pmid(icol, k)
          cosz = coszrs(icol)
 
          ! NOTE: Based upon US Standard Atmosphere 1976.
          if ((pressure >= carma_emis_lev(carma_emis_ilev_min)) .and. &
              (pressure <= carma_emis_lev(carma_emis_ilev_max)) .and. &
              (cosz > 0.0)) then

            ! The rates in are in terms of g cm-3 s-1, 
            !
            ! The values are in a lookup table, so find the two numbers
            ! surrounding the pressure and do a linear interpolation on the
            ! rate. This linear search is kind of expensive, particularly if
            ! there are a lot of points.
            !
            ! We turn this into a bilinear interpolation wtih zenith angle added
            ! 
            ! NOTE: The tendency is on a mass mixing ratio (kg/kg/s)

             do ilev = carma_emis_ilev_min, (carma_emis_ilev_max - carma_emis_ilev_incr), carma_emis_ilev_incr 
               if ((pressure >= carma_emis_lev(ilev)) .and. (pressure <= carma_emis_lev(ilev+carma_emis_ilev_incr))) then 
                 p_ref_index = ilev 
                 p_ref_indexp1 = p_ref_index + 1
                 exit
               endif
             end do
             if (p_ref_index .eq. carma_emis_ilev_max) then
                p_ref_index = p_ref_index - 1
                p_ref_indexp1 = p_ref_index + 1
                press = (pressure - carma_emis_lev(p_ref_index)) / &
                        (carma_emis_lev(p_ref_indexp1) - carma_emis_lev(p_ref_index))
             else
               press = (pressure - carma_emis_lev(p_ref_index)) / &
                       (carma_emis_lev(p_ref_indexp1) - carma_emis_lev(p_ref_index))
             endif

             if (cosz < carma_emis_zen(carma_emis_izen_min)) cosz = carma_emis_zen(carma_emis_izen_min)
             if (cosz > carma_emis_zen(carma_emis_izen_max)) cosz = carma_emis_zen(carma_emis_izen_max)
             do izen = carma_emis_izen_min, (carma_emis_izen_max - carma_emis_izen_incr), carma_emis_izen_incr
               if ((cosz >= carma_emis_zen(izen)) .and. (cosz <= carma_emis_zen(izen+carma_emis_izen_incr))) then 
                 z_ref_index = izen
                 z_ref_indexp1 = z_ref_index + 1
                 exit
               endif
             end do
             if (z_ref_index .eq. carma_emis_izen_max) then
                z_ref_index = z_ref_index - 1
                z_ref_indexp1 = z_ref_index + 1
                zen = (cosz - carma_emis_zen(z_ref_index)) / &
                        (carma_emis_zen(z_ref_indexp1) - carma_emis_zen(z_ref_index))
             else
                zen = (cosz - carma_emis_zen(z_ref_index)) / &
                       (carma_emis_zen(z_ref_indexp1) - carma_emis_zen(z_ref_index))
             endif

             vbi(1) = carma_emis_rate(p_ref_index  ,  z_ref_index)
             vbi(2) = carma_emis_rate(p_ref_index+1,  z_ref_index)
             vbi(3) = carma_emis_rate(p_ref_index  ,  z_ref_index+1)
             vbi(4) = carma_emis_rate(p_ref_index+1,  z_ref_index+1)

             onemp = 1. - press
             onemz = 1. - zen

             rate    = vbi(1)*onemp*onemz &
                     + vbi(2)*press*onemz &
                     + vbi(3)*onemp*zen &
                     + vbi(4)*press*zen

            ! Calculate the mass flux in terms of kg/m3/s
            massflux = rate * 1.0e-3_r8 * 1.0e6_r8  ! g cm-3 s-1 to kg/m3/s
             
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            thickness = state%zi(icol, k) - state%zi(icol, k+1)
            tendency(icol, k) = (massflux * thickness) / (state%pdel(icol, k) / gravit)        
          end if ! if particle emission presssure and zentih criteria met
        enddo
      enddo

      ! Scale the columns to conserve mass emission rate, carma_emis_total
      !
      ! Presently this scales the haze mass to carma_emis_total by fixing
      ! each column integrated emission rate to a constant value [kg/m2/s].
      ! Naturally, grid cells near the equator have larger spatial areas 
      ! and yield more total haze emission.  However, a future implementation should 
      ! permit column densities to differ across the grid
      do icol = 1, ncol
        cosz = coszrs(icol)
        if ( (cosz > 0.00001)) then
          ! haze prodcution only on the night side
          columnMass = sum(tendency(icol, :) * (state%pdel(icol, :) / gravit))
          scale = carma_emis_expected / columnMass     
          tendency(icol, :) = tendency(icol, :) * scale 
        else
           ! no photochemical haze production on the nightside  
           tendency(icol, :) = 0.0 
           columnMass = 0.0
        endif
      
      end do
    end if  ! if "HAZE" type
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use ioFileMod,    only: getfil
    use constituents, only: pcnst
    use wrap_nf

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ilev                  ! level index
    integer                            :: fid                   ! file id
    integer                            :: lev_did               ! level dimension id
    integer                            :: zen_did               ! zenith dimension id
    integer                            :: lev_vid               ! level variable id
    integer                            :: zen_vid               ! zenith variable id
    integer                            :: rate_vid              ! rate variable
    integer                            :: tmp
    integer                            :: lat_did               ! latitude dimension id
    integer                            :: ltime_did             ! local time dimension id
    integer                            :: time_did              ! time
    integer                            :: lat_vid               ! latitude variable id
    integer                            :: lrf_vid               ! local relative flux variable id
    integer                            :: grf_vid               ! global relative flux variable id
    integer                            :: ltime_vid             ! local time variable id
    character(len=256)                 :: efile                 ! emission file name

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
    ! Initialize the emissions rate table.
    if (carma_do_emission) then 
      if (masterproc) then

        ! Open the netcdf file (read only)
        call getfil(carma_emis_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading particle emission rates from ', efile

        call wrap_open(efile, 0, fid)

        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "lev", lev_did)
        call wrap_inq_dimlen(fid, lev_did, carma_emis_nLevs)

        call wrap_inq_dimid(fid, "zenith", zen_did)
        call wrap_inq_dimlen(fid, zen_did, carma_emis_nZens)

      endif
    
#if ( defined SPMD )
      call mpibcast(carma_emis_nLevs, 1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_nZens, 1, mpiint, 0, mpicom)
#endif

      allocate(carma_emis_lev(carma_emis_nLevs))
      allocate(carma_emis_zen(carma_emis_nZens))
      allocate(carma_emis_rate(carma_emis_nLevs, carma_emis_nZens))

      if (masterproc) then
        call wrap_inq_varid(fid, 'lev', lev_vid)
        call wrap_get_var_realx(fid, lev_vid, carma_emis_lev)

        call wrap_inq_varid(fid, 'cosz', zen_vid)
        call wrap_get_var_realx(fid, zen_vid, carma_emis_zen)

        write(*,*) "CARMA emis levs", carma_emis_lev
        write(*,*) "CARMA emis zenith", carma_emis_zen

        ! Read in the tables.
        call wrap_inq_varid(fid, 'MHAZE', rate_vid)
        !call wrap_get_var_realx(fid, rate_vid, carma_emis_rate)
        tmp = nf90_get_var (fid, rate_vid, carma_emis_rate)
        if (tmp/=NF90_NOERR) then
           write(iulog,*) 'CARMA_InitializeModel: error reading varid =', rate_vid
           call handle_error (tmp)
        end if

        write(*,*) "CARMA emis rate", carma_emis_rate

        ! Close the file.
        call wrap_close(fid)

        ! Find out where the bounds of the table are and in what order
        ! the pressures levels are in.
        carma_emis_ilev_min = 1
        carma_emis_ilev_max = carma_emis_nLevs

        if (carma_emis_lev(carma_emis_ilev_min) < carma_emis_lev(carma_emis_ilev_max)) then
          carma_emis_ilev_incr = 1
        else
          carma_emis_ilev_incr = -1
          tmp = carma_emis_ilev_min
          carma_emis_ilev_min = carma_emis_ilev_max
          carma_emis_iLev_max = tmp 
        endif

        ! Find out where the bounds of the table are and in what order
        ! the zenith angles are in.
        carma_emis_izen_min = 1
        carma_emis_izen_max = carma_emis_nZens

        if (carma_emis_zen(carma_emis_izen_min) < carma_emis_zen(carma_emis_izen_max)) then
          carma_emis_izen_incr = 1
        else
          carma_emis_izen_incr = -1
          tmp = carma_emis_izen_min
          carma_emis_izen_min = carma_emis_izen_max
          carma_emis_izen_max = tmp 
        endif

        if (do_print) write(LUNOPRT,*) ' '
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_nLevs     = ', carma_emis_nLevs
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_min  = ', carma_emis_ilev_min 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_max  = ', carma_emis_ilev_max 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_incr = ', carma_emis_ilev_incr 

        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_nZens     = ', carma_emis_nZens
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_izen_min  = ', carma_emis_izen_min 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_izen_max  = ', carma_emis_izen_max 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_izen_incr = ', carma_emis_izen_incr 
        if (do_print) write(LUNOPRT,*) ''
        
        if (do_print) write(LUNOPRT,*) 'level, pressure (Pa), emission rate (g cm-3 sec-1)'
        do ilev = carma_emis_ilev_min, carma_emis_ilev_max, carma_emis_ilev_incr
          if (do_print) write(LUNOPRT,*) ilev, carma_emis_lev(ilev), carma_emis_rate(ilev, carma_emis_izen_min)
        enddo
        
        if (do_print) write(LUNOPRT, *) 'carma_init(): Total Emission = ', carma_emis_total, ' (kt/yr)'
        ! expected smission rate per column (kg/m2/s)
        carma_emis_expected = ((carma_emis_total * 1e6_r8) / (3600.0_r8 * 24.0_r8 * 365.0_r8)) / (4.0_r8 * PI * ((REARTH / 100._r8) ** 2)) * illumFrac
        if (do_print) write(LUNOPRT,*) 'carma_init(): Done with emission table.'

      endif

#if ( defined SPMD )
      call mpibcast(carma_emis_lev,  carma_emis_nLevs, mpir8, 0, mpicom)
      call mpibcast(carma_emis_zen,  carma_emis_nZens, mpir8, 0, mpicom)
      call mpibcast(carma_emis_rate, carma_emis_nLevs*carma_emis_nZens, mpir8, 0, mpicom)
      call mpibcast(carma_emis_expected,  1, mpir8,  0, mpicom)
      call mpibcast(carma_emis_ilev_min,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_max,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_incr, 1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_izen_min,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_izen_max,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_izen_incr, 1, mpiint, 0, mpicom)
#endif

    endif
    
    
    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, q, gcid, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    real(r8), intent(inout)            :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer, intent(in)                :: gcid(:)  ! global column id
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

    return
  end subroutine CARMA_InitializeParticle
  
  
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_WetDeposition 

end module
