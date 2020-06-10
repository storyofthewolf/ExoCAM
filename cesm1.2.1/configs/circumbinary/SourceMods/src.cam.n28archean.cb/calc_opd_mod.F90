
module calc_opd_mod
! version n28archean

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use rayleigh_data
  use shr_const_mod,    only: SHR_CONST_PI,SHR_CONST_PI, SHR_CONST_G, &
                              SHR_CONST_RGAS, SHR_CONST_AVOGAD, &
                              SHR_CONST_STEBOL, &
                              SHR_CONST_BOLTZ, &
                              SHR_CONST_RHOFW, SHR_CONST_RHOICE, &
                              SHR_CONST_LOSCHMIDT
  use physconst,        only: mwn2, mwco2, mwch4, mwh2o, mwo2, mwh2, mwo3, mwdry, cpair
  use radgrid
  use rad_interp_mod
  use ppgrid

  implicit none
  private
  save

  public :: calc_gasopd
  public :: calc_cldopd


!============================================================================
contains
!============================================================================

!============================================================================

  subroutine calc_gasopd(tmid, pmid, pdel, coldens, coldens_dry, qh2o, qco2, qch4, qO2, qO3, qH2, qN2, pathlength, &
                         tau_gas, tau_ray)

!------------------------------------------------------------------------
!
! Purpose: Calculate the optical depths of gases
!          Optical depths stored in tau_gas(g_value,wavelength_band,vertical_level)
!
!------------------------------------------------------------------------
!
                 
    !use time_manager,   only: get_nstep
   
    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in), dimension(pverp) :: tmid         ! temperatures at mid layers [K]
    real(r8), intent(in), dimension(pverp) :: pmid         ! pressure midlayers [mb]  
    real(r8), intent(in), dimension(pver) :: pdel          ! layer thickness [mb]  
    real(r8), intent(in), dimension(pverp) :: coldens      ! Wet Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: coldens_dry  ! Dry Column density profile [molec m-2]
    real(r8), intent(in), dimension(pverp) :: qh2o         ! mass mixing ratio h2o profile [kg/kg] wet
    real(r8), intent(in), dimension(pverp) :: qco2         ! mass mixing ratio co2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qch4         ! mass mixing ratio ch4 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qo2          ! mass mixing ratio o2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qo3          ! mass mixing ratio o3 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qh2          ! mass mixing ratio h2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: qn2          ! mass mixing ratio n2 profile [kg/kg] dry
    real(r8), intent(in), dimension(pverp) :: pathlength   ! thickness of layer [cm]   
    real(r8), intent(out), dimension(ntot_gpt,pverp) ::  tau_gas 
    real(r8), intent(out), dimension(ntot_wavlnrng,pverp) ::  tau_ray   


!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8) :: species_weight_h2o_co2
    real(r8) :: species_weight_h2o_ch4
    real(r8) :: species_weight_o3_co2
    real(r8) :: species_weight_h2o_o3
    real(r8) :: species_weight_h2o_o2
    real(r8) :: h2ovmr
    real(r8) :: co2vmr
    real(r8) :: ch4vmr
    real(r8) :: o2vmr
    real(r8) :: o3vmr
    real(r8) :: h2vmr
    real(r8) :: n2vmr
    real(r8) :: pressure
    real(r8) :: pressure_s
    real(r8) :: temperature
    integer :: p_ref_index  
    integer :: t_ref_index_lower
    integer :: t_ref_index_upper
    integer :: w_ref_index_h2o_co2 
    integer :: w_ref_index_h2o_ch4
    integer :: w_ref_index_o3_co2
    integer :: w_ref_index_h2o_o3
    integer :: w_ref_index_h2o_o2
    integer :: t_ref_index_s
    integer :: w_ref_index_s
    integer :: t_ref_index_h2h2
    integer :: t_ref_index_h2n2
    integer :: t_ref_index_n2n2
    integer :: ik
    integer :: iw
    integer :: ig
    integer :: itu
    integer :: itl
    integer :: itc
    !integer :: nstep
    integer  :: iwbeg
    integer  :: iwend
    real(r8) :: ans
    real(r8) :: u_col
    real(r8) :: u_h2o        ! column amount of H2O
    real(r8) :: u_co2        ! column amount of CO2
    real(r8) :: u_ch4        ! column amount of CH4
    real(r8) :: u_o2         ! column amount of O2
    real(r8) :: u_o3         ! column amount of O3
    real(r8) :: u_h2         ! column amount of H2
    real(r8) :: u_n2         ! column amount of N2
    real(r8) :: kh2o
    real(r8) :: ph2o
    real(r8) :: cv
    real(r8) :: ct
    real(r8) :: wm
    real(r8) :: wl, wla
    real(r8) :: r
    real(r8) :: ns
    real(r8) :: pdelik
    real(r8) :: depolN2      ! for rayleigh scattering calc, depolarization
    real(r8) :: depolCO2
    real(r8) :: depolH2O      
    real(r8) :: allenN2      ! for rayleigh scattering calc, Allen (1976) coefficients
    real(r8) :: allenCO2
    real(r8) :: sigmaRaylCO2 ! rayleigh scattering cross sections [cm2 molecule-1]
    real(r8) :: sigmaRaylN2
    real(r8) :: sigmaRaylH2
    real(r8) :: sigmaRaylH2O
    real(r8) :: sigmaRayl
    real(r8) :: kg_sw_minval  !! minimum value to check sw_abs error
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_n2n2cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2n2cia
    real(r8), dimension(ntot_wavlnrng,pverp) ::  tau_h2h2cia


!------------------------------------------------------------------------
!
! Start Code
!

    !
    kg_sw_minval =1.0e-99

    ! calculate optical depth over all spectral intervals
    iwbeg = 1
    iwend = ntot_wavlnrng

    !nstep = get_nstep()
 
    !open (unit = 10, file = "gas_opd.txt", status = "replace", position = "append", iostat = openstatus)

    do ik = 1,pverp 

      ! qco2, qh2o, qch4 are in mass mixing ratio, multiply by mwdry/mw*** to get volume mixing ratio
      ! factor of 10^-4 converts molecules m-2 to molecules cm-2

      !qh2o is defined relative to wet atm mass
      !all others are defined relative to dry atm mass
      !set all vmr quantities relative to dry mass

      h2ovmr = qh2o(ik)/(1.0-qh2o(ik))*mwdry/mwh2o         ! H2O volume mixing ratio
      co2vmr = qco2(ik)*mwdry/mwco2         ! CO2 volume mixing ratio
      ch4vmr = qch4(ik)*mwdry/mwch4         ! CH4 volume mixing ratio
      o2vmr = qo2(ik)*mwdry/mwo2            ! O2 volume mixing ratio
      o3vmr = qo3(ik)*mwdry/mwo3            ! O3 volume mixing ratio
      h2vmr = qh2(ik)*mwdry/mwh2            ! H2 volume mixing ratio
      n2vmr = qn2(ik)*mwdry/mwn2            ! N2 volume mixing ratio

      u_h2o = h2ovmr*coldens_dry(ik)/10000.     !   water column amount [ molecules cm-2 ]
      u_co2 = co2vmr*coldens_dry(ik)/10000.     !   co2 column amount [ molecules cm-2 ]
      u_ch4 = ch4vmr*coldens_dry(ik)/10000.     !   ch4 column amount [ molecules cm-2 ]
      u_o2 = o2vmr*coldens_dry(ik)/10000.       !   o2 column amount [ molecules cm-2 ]
      u_o3 = o3vmr*coldens_dry(ik)/10000.       !   o3 column amount [ molecules cm-2 ]
      u_h2 = h2vmr*coldens_dry(ik)/10000.       !   h2 column amount [ molecules cm-2 ]
      u_n2 = n2vmr*coldens_dry(ik)/10000.       !   n2 column amount [ molecules cm-2 ]

      pressure = log10(pmid(ik))       ! log pressure
  
      p_ref_index = kc_npress          ! index of reference pressure, default is max
!      p_ref_index_s = ks_npress

      ! find the reference pressure value, exit if pressure less than minimum of grid     
      do  ! for K coefficient data sets        
        if (p_ref_index .le. 1) exit
        if (log10pgrid(p_ref_index) .le. pressure ) exit        
        p_ref_index = p_ref_index - 1
      enddo
      ! if pressure less than minimum of grid, force reference to minimum grid value
      ! force reference pressure for interpolation to minimum pressure in pgrid
      !if (pressure .le. log10pgrid(1)) then
      !  p_ref_index = 1
      !  pressure = log10pgrid(p_ref_index)        
      !endif      

      ! if pressure greater than the maximum of the grid, force reference to maximum grid value
      ! force reference pressure for interpolation to maximum pressure in pgrid
      if (pressure .ge. log10pgrid(kc_npress)) then
        p_ref_index = kc_npress
	pressure = log10pgrid(kc_npress)
      endif


      if (pressure .ge. 2.0) then
        p_ref_index = p_ref_index - (kc_npress_upper-1)
      endif

      temperature = tmid(ik)              ! actual temperature [K]
      t_ref_index_lower = kc_ntemp_lower  ! index of reference temperature lower
      t_ref_index_upper = kc_ntemp_upper  ! index of reference temperature upper
      t_ref_index_s = ks_ntemp
      t_ref_index_h2h2 = kh2h2_ntemp
      t_ref_index_h2n2 = kh2n2_ntemp
      t_ref_index_n2n2 = kn2n2_ntemp

      ! For gas absorption k coefficients
      ! find the reference lower temperature value, exit if temperature less than minimum of grid
      do  ! for k coefficient data sets
        if (t_ref_index_lower .le. 1) exit
        if ((tgrid_lower(t_ref_index_lower) .le. temperature)) exit 
        t_ref_index_lower = t_ref_index_lower - 1
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_lower .lt. 1) then
        t_ref_index_lower = 1
        temperature = tgrid_lower(t_ref_index_lower)
      endif

      ! find the reference upper temperature value, exit if temperature less than minimum of grid
      do  ! for k coefficient data sets
        if (t_ref_index_upper .le. 1) exit
        if ((tgrid_upper(t_ref_index_upper) .le. temperature)) exit 
        t_ref_index_upper = t_ref_index_upper - 1
      enddo

      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_upper .lt. 1) then
        t_ref_index_upper = 1
        temperature = tgrid_upper(t_ref_index_upper)
      endif     

      ! For water vapor self continuum, find temperature index
      !if (pressure .ge. 2.0) then

        ! find the reference temperature value, exit if temperature less than minimum of grid
        do  ! for k coefficient data sets
          if (t_ref_index_s .le. 1) exit
          if ((tgrid_self(t_ref_index_s) .le. temperature)) exit 
          t_ref_index_s = t_ref_index_s - 1
        enddo
        ! if temperature less than minimum of grid, force reference to minimum grid value
        ! force reference temperature for interpolation to minimum temperature in tgrid
        if (t_ref_index_s .lt. 1) then
          t_ref_index_s = 1
          temperature = tgrid_self(t_ref_index_s)
        endif

      !endif ! H2O self temperature index

      ! For H2-N2 CIA, find temperature index
      !! find the reference temperature value, exit if temperature less than minimum of grid
      do  
        if (t_ref_index_h2n2 .le. 1) exit
        if ((tgrid_h2n2(t_ref_index_h2n2) .le. temperature)) exit
        t_ref_index_h2n2 = t_ref_index_h2n2 - 1
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_h2n2 .lt. 1) then
        t_ref_index_h2n2 = 1
        temperature = tgrid_h2n2(t_ref_index_h2n2)
      endif

      ! For H2-H2 CIA, find temperature index
      ! find the reference temperature value, exit if temperature less than minimum of grid
      do  
        if (t_ref_index_h2h2 .le. 1) exit
        if ((tgrid_h2h2(t_ref_index_h2h2) .le. temperature)) exit
        t_ref_index_h2h2 = t_ref_index_h2h2 - 1
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_h2h2 .lt. 1) then
        t_ref_index_h2h2 = 1
        temperature = tgrid_h2h2(t_ref_index_h2h2)
      endif

      ! For N2-N2 CIA, find temperature index
      ! find the reference temperature value, exit if temperature less than minimum of grid
      do  
        if (t_ref_index_n2n2 .le. 1) exit
        if ((tgrid_n2n2(t_ref_index_n2n2) .le. temperature)) exit
        t_ref_index_n2n2 = t_ref_index_n2n2 - 1
      enddo
      ! if temperature less than minimum of grid, force reference to minimum grid value
      ! force reference temperature for interpolation to minimum temperature in tgrid
      if (t_ref_index_n2n2 .lt. 1) then
        t_ref_index_n2n2 = 1
        temperature = tgrid_n2n2(t_ref_index_n2n2)
      endif

      ! normalized species weight
      ! add index finding loops to these if statement calls
      if ((u_h2o .gt. 0.0) .or. (u_co2 .gt. 0.0)) species_weight_h2o_co2 = u_h2o/(u_co2+u_h2o)
      if ((u_h2o .gt. 0.0) .or. (u_ch4 .gt. 0.0)) species_weight_h2o_ch4 = u_h2o/(u_ch4+u_h2o)
      if ((u_o3 .gt. 0.0) .or. (u_co2 .gt. 0.0))  species_weight_o3_co2 = u_o3/(u_co2+u_o3)
      if ((u_h2o .gt. 0.0) .or. (u_o3 .gt. 0.0))  species_weight_h2o_o3 = u_h2o/(u_o3+u_h2o)
      if ((u_h2o .gt. 0.0) .or. (u_o2 .gt. 0.0))  species_weight_h2o_o2 = u_h2o/(u_o2+u_h2o)

      w_ref_index_h2o_co2 = kc_nweight
      w_ref_index_h2o_ch4 = kc_nweight
      w_ref_index_o3_co2 = kc_nweight
      w_ref_index_h2o_o3 = kc_nweight
      w_ref_index_h2o_o2 = kc_nweight
      w_ref_index_s = ks_nweight

      ! find the reference species weight H2O/CO2, exit if weight less than minimum of grid
      do
        if (w_ref_index_h2o_co2 .le. 1) exit
        if ((wgrid(w_ref_index_h2o_co2) .le. species_weight_h2o_co2)) exit
        w_ref_index_h2o_co2 =  w_ref_index_h2o_co2 - 1
      enddo

      ! if species weight less than minimum of grid, force reference to minimum grid value
      if (w_ref_index_h2o_co2 .lt. 1) then
        w_ref_index_h2o_co2 = 1
      endif

      ! find the reference species weight H2O/CH4, exit if weight less than minimum of grid
      do
        if (w_ref_index_h2o_ch4 .le. 1) exit
        if ((wgrid(w_ref_index_h2o_ch4) .le. species_weight_h2o_ch4)) exit
        w_ref_index_h2o_ch4 =  w_ref_index_h2o_ch4 - 1
      enddo

      ! if species weight less than minimum of grid, force reference to minimum grid value
      if (w_ref_index_h2o_ch4 .lt. 1) then
        w_ref_index_h2o_ch4 = 1
      endif

      ! find the reference species weight O3/CO2, exit if weight less than minimum of grid
      do
        if (w_ref_index_o3_co2 .le. 1) exit
        if ((wgrid(w_ref_index_o3_co2) .le. species_weight_o3_co2)) exit
        w_ref_index_o3_co2 =  w_ref_index_o3_co2 - 1
      enddo

      ! if species weight less than minimum of grid, force reference to minimum grid value
      if (w_ref_index_o3_co2 .lt. 1) then
        w_ref_index_o3_co2 = 1
      endif

      ! find the reference species weight H2O/O3, exit if weight less than minimum of grid
      do
        if (w_ref_index_h2o_o3 .le. 1) exit
        if ((wgrid(w_ref_index_h2o_o3) .le. species_weight_h2o_o3)) exit
        w_ref_index_h2o_o3 =  w_ref_index_h2o_o3 - 1
      enddo

      ! if species weight less than minimum of grid, force reference to minimum grid value
      if (w_ref_index_h2o_o3 .lt. 1) then
        w_ref_index_h2o_o3 = 1
      endif

      ! find the reference species weight H2O/O2, exit if weight less than minimum of grid
      do
        if (w_ref_index_h2o_o2 .le. 1) exit
        if ((wgrid(w_ref_index_h2o_o2) .le. species_weight_h2o_o2)) exit
        w_ref_index_h2o_o2 =  w_ref_index_h2o_o2 - 1
      enddo

      ! if species weight less than minimum of grid, force reference to minimum grid value
      if (w_ref_index_h2o_o2 .lt. 1) then
        w_ref_index_h2o_o2 = 1
      endif

!NOT USED?
      ! H2O weight(vmr) reference, find index
      !if (pressure .ge. 2.0) then
!
!        do
!          if (w_ref_index_s .le. 1) exit
!          if ((wgrid_self(w_ref_index_s) .le. h2ovmr)) exit
!          w_ref_index_s =  w_ref_index_s - 1
!        enddo
!
!        ! if species weight less than minimum of grid, force reference to minimum grid value
!        if (w_ref_index_s .lt. 1) then
!          w_ref_index_s = 1
!        endif

      !endif

      !====== interpolate to find k(p,t,w), calculate opd ======! 
      itl = 0
      itu = 0
  
      !=====  interval 1: 10 - 350 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o      
        do ig = 1,ngauss_pts(1)      
          itl = itl + 1         
          call bilinear_interpK_lower_8gpt(k01_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o      
        do ig = 1,ngauss_pts(1)         
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k01_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !===== interval 2: 350 - 500 cm-1 =====!
      
      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o
        do ig = 1,ngauss_pts(2)    
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k02_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o      
        do ig = 1,ngauss_pts(2)      
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k02_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif 

      !===== interval 3: 500 - 630 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(3)   
          itl = itl + 1
          call trilinear_interpK_lower_16gpt(k03_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species,  H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(3)   
          itu = itu + 1
          call trilinear_interpK_upper_16gpt(k03_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo  
      endif

      !===== interval 4: 630 - 700 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(4)  
          itl = itl + 1
          call trilinear_interpK_lower_16gpt(k04_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower , &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(4)  
          itu = itu + 1
          call trilinear_interpK_upper_16gpt(k04_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper , &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col
        enddo   
      endif

      !===== interval 5: 700 - 820 cm-1 =====!

      if (pressure .ge. log10Pcut) then   
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(5)   
          itl = itl + 1
          call trilinear_interpK_lower_16gpt(k05_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
         enddo   
      else
        ! Upper Atmoshere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(5)   
          itu = itu + 1
          call trilinear_interpK_upper_16gpt(k05_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !===== interval 6: 820 - 980 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(6)    
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k06_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(6)   
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k06_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !===== interval 7: 980 - 1100 cm-1 =====!
      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosophere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(7)		  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k07_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(7)   
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k07_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !=====interval 8: 1100 - 1180 cm-1 =====!
      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4
        do ig = 1,ngauss_pts(8)             
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k08_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, & 
                                       species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4  
        do ig = 1,ngauss_pts(8)             
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k08_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, & 
                                       species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo
      endif

      !===== interval 9: 1180 - 1390 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(9)   
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k09_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col 
        enddo
      else   
        ! Upper Atmosphere:  major species,H2O, CH4
        u_col = u_h2o + u_ch4
        do ig = 1,ngauss_pts(9)   
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k09_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itu,ik) = ans*u_col 
        enddo   
      endif

      !===== interval 10:  1390 - 1480 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4
        do ig = 1,ngauss_pts(10)             
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k10_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                      species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4    
        do ig = 1,ngauss_pts(10)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k10_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                      species_weight_h2o_ch4, w_ref_index_h2o_ch4,  ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !===== interval 11:  1480 - 1800 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(11)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k11_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                       species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo
      else
        ! Upper Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(11)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k11_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, & 
                                       species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo
      endif

      !==== interval 12: 1800 - 2080 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(12)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k12_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                               species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo      
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(12)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k12_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo      
      endif

      !===== interval 13: 2080 - 2200 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(13)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k13_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, & 
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(13) 
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k13_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !==== interval 14: 2200 - 2380 cm-1, CO2 ====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, CO2
        u_col = u_co2        
        do ig = 1,ngauss_pts(14)   
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k14_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, CO2
        u_col = u_co2        
        do ig = 1,ngauss_pts(14)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k14_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !==== interval 15: 2380 - 2600 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, CO2
        u_col = u_co2        
        do ig = 1,ngauss_pts(15)  
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k15_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, CO2
        u_col = u_co2        
        do ig = 1,ngauss_pts(15)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k15_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
     endif 

      !===== interval 16: 2600 - 3250 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(16)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k16_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                       species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(16)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k16_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !===== interval 17:  3250 - 4000 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(17)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k17_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(17)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k17_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo   
      endif

      !===== interval 18: 4000 - 4650 cm-1 =====!
   
      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4       
        do ig = 1,ngauss_pts(18)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k18_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo      
      else
        ! Upper Atmosphere:  major species, H2O, CH4
        u_col = u_h2o + u_ch4        
        do ig = 1,ngauss_pts(18)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k18_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                 species_weight_h2o_ch4, w_ref_index_h2o_ch4, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo      
      endif

      !===== interval 19: 4650 - 5150 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(19)  
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k19_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                 species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo      	    
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(19)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k19_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo      	    
      endif

      !===== interval 20: 5150 - 6150 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o
        do ig = 1,ngauss_pts(20)  
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k20_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o
        do ig = 1,ngauss_pts(20)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k20_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !===== interval 21: 6150 - 7700 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2        
        do ig = 1,ngauss_pts(21)
          itl = itl + 1
          call trilinear_interpK_lower_8gpt(k21_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, &
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itl,ik) = ans*u_col          
        enddo      
      else
        ! Upper Atmosphere:  major species, H2O, CO2
        u_col = u_h2o + u_co2
        do ig = 1,ngauss_pts(21)  
          itu = itu + 1
          call trilinear_interpK_upper_8gpt(k21_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, &
                                       species_weight_h2o_co2, w_ref_index_h2o_co2, ans)
          tau_gas(itu,ik) = ans*u_col          
        enddo      
      endif

      !===== interval 22: 7700 - 8050 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o        
        do ig = 1,ngauss_pts(22)  
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k22_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o
        do ig = 1,ngauss_pts(22)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k22_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !===== interval 23: 8050 - 12850 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o        
        do ig = 1,ngauss_pts(23)  
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k23_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o        
        do ig = 1,ngauss_pts(23)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k23_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !===== interval 24: 12850 - 16000 cm-1 =====!

      if (pressure .ge. log10Pcut) then 
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o       
        do ig = 1,ngauss_pts(24)  
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k24_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
     else
        ! Upper Atmosphere: major species, H2O
        u_col = u_h2o        
        do ig = 1,ngauss_pts(24)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k24_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      !===== interval 25: 16000 - 22650 cm-1 =====!

      if (pressure .ge. log10Pcut) then       
        ! Lower Atmosphere:  major species, H2O
        u_col = u_h2o        
        do ig = 1,ngauss_pts(25)   
          itl = itl + 1
          call bilinear_interpK_lower_8gpt(k25_lower, ig, pressure, p_ref_index, temperature, t_ref_index_lower, ans)
          tau_gas(itl,ik) = ans*u_col           
        enddo   
      else
        ! Upper Atmosphere:  major species, H2O
        u_col = u_h2o       
        do ig = 1,ngauss_pts(25)  
          itu = itu + 1
          call bilinear_interpK_upper_8gpt(k25_upper, ig, pressure, p_ref_index, temperature, t_ref_index_upper, ans)
          tau_gas(itu,ik) = ans*u_col           
        enddo   
      endif

      ! Used for both Water Vapor Self Continuum and Rayleigh Scattering
      if (ik .eq. 1) then
        pdelik = pmid(1)
      else
        pdelik = pdel(ik-1)
      endif
      
      !
      !  Water Vapor Self Continuum 
      !  MT_CKD
      itc=0
      do iw = iwbeg,iwend     ! loop over bands      
        call interpKself2(kh2oself, iw, temperature, t_ref_index_s, ans)
        do ig=1, ngauss_pts(iw)
          itc=itc+1
          tau_gas(itc,ik) = tau_gas(itc,ik) +ans*u_h2o*(pmid(ik)/1013.)*h2ovmr          
        enddo
      enddo    ! close band loop

      
      !! KASTING 1984
      !!  CO2 Continuum 
      !!  Results fit to LBLRTM/MT_CKD v2.5
      itc = 0
      do iw = iwbeg,iwend     ! loop over bands      
        do ig=1, ngauss_pts(iw)
       
          itc = itc + 1
          tau_gas(itc,ik) = tau_gas(itc,ik) + kco2cont(itc)*u_co2 &
                          * (pmid(ik)/pref_co2)*(1.0+0.3*co2vmr)*(co2vmr/vref_co2) &
                          * (296.0/tmid(ik))
           !write(*,*) "CO2 cont", iw, kco2cont(itc)
         enddo
      enddo    ! close band loop

      !
      !  Calculate N2-N2 collision induced absorption from HITRAN
      !
      do iw=iwbeg,iwend      ! loop over bands
        ! N2-N2 CIA coefficients from Borysow et al. (1986)
        call interpN2N2cia(kn2n2, iw, temperature, t_ref_index_n2n2, ans)
        tau_n2n2cia(iw,ik) = ans * pmid(ik)*100./(SHR_CONST_BOLTZ*tmid(ik)) / SHR_CONST_LOSCHMIDT &
                               * u_n2 / (SHR_CONST_LOSCHMIDT*1.0e-6) &
                               * n2vmr
        !write(*,*) "N2-N2 CIA",iw, ans, n2vmr !, tau_n2n2cia(iw,ik)
      enddo
     

      !
      !  Calculate H2-H2 collision induced absorption
      !
      do iw=iwbeg,iwend      ! loop over bands
        ! H2-H2 CIA coefficients from Borysow et al. (1986)
        call interpH2H2cia(kh2h2, iw, temperature, t_ref_index_h2h2, ans)        
        tau_h2h2cia(iw,ik) = ans * pmid(ik)*100./(SHR_CONST_BOLTZ*tmid(ik)) / SHR_CONST_LOSCHMIDT &
                               * u_h2 / (SHR_CONST_LOSCHMIDT*1.0e-6) &
                               * h2vmr 
        !write(*,*) "H2-H2 CIA",iw, ans, h2vmr !, tau_h2h2cia(iw,ik)
      enddo
        
      !
      !  Calculate H2-N2 collision induced absorption
      !    
      do iw=iwbeg,iwend      ! loop over bands
        ! H2-N2 CIA coefficients from Borysow et al. (1986)
        call interpH2N2cia(kh2n2, iw, temperature, t_ref_index_h2n2, ans)        
        tau_h2n2cia(iw,ik) = ans * pmid(ik)*100./(SHR_CONST_BOLTZ*tmid(ik)) / SHR_CONST_LOSCHMIDT &
                               * u_h2 / (SHR_CONST_LOSCHMIDT*1.0e-6) &
                               * (1.-h2vmr)
        !write(*,*) "H2-N2 CIA",iw, ans, h2vmr !, tau_h2n2cia(iw,ik)
      enddo

      !
      ! Add CIA optical depths to total optical depth
      !
      itc = 0
      do iw=iwbeg,iwend      ! loop over bands
        do ig=1, ngauss_pts(iw)
          itc = itc + 1
          tau_gas(itc, ik) = tau_gas(itc, ik) + tau_n2n2cia(iw,ik) + tau_h2n2cia(iw,ik) + tau_h2h2cia(iw,ik)
        enddo
      enddo

      !
      !  Calculate Rayleigh scattering optical depth
      !
      do iw = iwbeg,iwend     ! loop over bands

        wm = (wavenum_edge(iw) + wavenum_edge(iw+1))/2.0    ! wavenumber at bin midpoint
        wl = 1.e4/wm  ! wavelength in microns
        wla = wl*1.0e4 ! wavelength in angstroms

        !
        ! Vardavas and Carter (1984), Allen (1976) N2, CO2 Rayleigh scattering
        !
        ! Rayleigh scattering for CO2, N2
        depolCO2 = (6+3*delCO2)/(6-7*delCO2)
        depolN2  = (6+3*delN2)/(6-7*delN2)
        allenCO2 = (1.0E-5*raylA_CO2*(1.0+1.0E-3*raylB_CO2/wl**2))**2
        allenN2 = (1.0E-5*raylA_N2*(1.0+1.0E-3*raylB_N2/wl**2))**2
        sigmaRaylCO2 = 4.577E-21/wl**4*depolCO2*allenCO2
        sigmaRaylN2 = 4.577E-21/wl**4*depolN2*allenN2         
        !  Rayleigh scattering from H2O 
        ns = (5791817./(238.0185-(1./wl)**2) + 167909./(57.362-(1./wl)**2))/1.0E8  ! Bucholtz (1995)  !new
        r = 0.85*ns 
        depolH2O = (6+3*delH2O)/(6-7*delH2O)
        sigmaRaylH2O = 4.577e-21*depolH2O*(r**2)/(wl**4)  !new
        ! Rayleigh scattering from H2, Dalgarno & Williams 1962, ApJ, 136, 690D
        sigmaRaylH2 = 8.14e-13/(wla**4) + 1.28e-6/(wla**6) + 1.61/(wla**8)
        ! Total Rayleigh scattering
        tau_ray(iw,ik) = sigmaRaylCO2*u_co2 + sigmaRaylN2*u_n2 + sigmaRaylH2O*u_h2o + sigmaRaylH2*u_h2

      enddo  ! close band loop
      
    enddo  ! close level loop    

!    close(10)
    return
  
  end subroutine calc_gasopd  
  


!============================================================================

  subroutine calc_cldopd(ext_pmid, ext_cICE, ext_cLIQ, ext_REI, ext_REL, ext_cFRC, & 
                         tau_cld_mcica, singscat_cld_mcica, asym_cld_mcica, cFRC_mcica, &
                         cICE_mcica, cLIQ_mcica )

!------------------------------------------------------------------------
!
!  Purpose: Calculate the optical depth and related properties of clouds.
!  Uses MIE optical properties for both ice clouds and liquid water clouds.
!  Uses Monte Carlo Independent Column Approximation (MCICA) to treat cloud
!    overlap. Pincus (2003), Raisanen (2004)  
!------------------------------------------------------------------------
   
   use mcica,            only: mcica_subcol      
   use time_manager,     only: get_nstep   

   implicit none

!------------------------------------------------------------------------
!
! Input Arguments
!
    real(r8), intent(in), dimension(pverp)  ::  ext_pmid      ! [Pa] interface pressure
    real(r8), intent(in), dimension(pverp)  ::  ext_cICE   ! [g/m2]  ice 
    real(r8), intent(in), dimension(pverp)  ::  ext_cLIQ   ! [g/m2]  liquid
    real(r8), intent(in), dimension(pverp)  ::  ext_REI    ! [microns]  ice radii
    real(r8), intent(in), dimension(pverp)  ::  ext_REL    ! [microns]  liquid
    real(r8), intent(in), dimension(pverp)  ::  ext_cFRC   ! cloud fraction

    ! output bulk cloud optical properties after MCICA
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  tau_cld_mcica       ! cloud optical depth
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  singscat_cld_mcica  ! cloud single scattering albedo
    real(r8), intent(out), dimension(ncld_grp, ntot_gpt, pverp)  ::  asym_cld_mcica      ! cloud asymmetry parameter

    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cFRC_mcica  ! stochastic cloud fraction
    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cICE_mcica  ! stochastic ice cloud water content
    real(r8), intent(out), dimension(ntot_gpt, pverp) :: cLIQ_mcica  ! stochastic liquid cloud water content

!------------------------------------------------------------------------
!
! Local Variables
!
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  tau_cld_temp       ! cloud optical depth temporary array
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  singscat_cld_temp  ! cloud single scattering albedo temporary array
    real(r8), dimension(ncld_grp, ntot_wavlnrng, pverp)  ::  asym_cld_temp      ! cloud asymmetry parameter temporary array
    real(r8), dimension(ntot_gpt, pverp) :: tau_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: ssa_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: asym_mcica_ice
    real(r8), dimension(ntot_gpt, pverp) :: tau_mcica_liq
    real(r8), dimension(ntot_gpt, pverp) :: ssa_mcica_liq
    real(r8), dimension(ntot_gpt, pverp) :: asym_mcica_liq

    real(r8) :: Qliq
    real(r8) :: Wliq
    real(r8) :: Gliq

    real(r8) :: Qice
    real(r8) :: Wice
    real(r8) :: Gice

    real(r8) :: rho_liq
    real(r8) :: rho_ice
    real(r8) :: r_liq
    real(r8) :: r_ice

    integer :: iwend
    integer :: iwbeg
    integer :: iw
    integer :: ik
    integer :: ig

    integer :: icldovr 
    integer :: permuteseed
    integer :: liqcld
    integer :: icecld    

    integer :: nstep

!------------------------------------------------------------------------
!
! Start Code
!

    rho_liq = SHR_CONST_RHOFW*1000.0   ! density of fresh water [g m-3]
    rho_ice = SHR_CONST_RHOICE*1000.0  ! density of ice water [g m-3]

    ! Calculate gas optical depths over all spectral intervals
    iwbeg = 1
    iwend = ntot_wavlnrng

    liqcld = 1
    icecld = 2

    !initialize cloud optics array
    tau_cld_mcica(:,:,:) = 0.0
    singscat_cld_mcica(:,:,:) = 0.0
    asym_cld_mcica(:,:,:) = 0.0
    tau_cld_temp(:,:,:) = 0.0
    singscat_cld_temp(:,:,:) = 0.0
    asym_cld_temp(:,:,:) = 0.0

    !
    ! Water clouds 
    ! 
    do ik=1, pverp
      do iw=iwbeg, iwend

      r_liq = ext_REL(ik) * 1.0e-6   ! liquid cloud drop size [m] 

        if (ext_cLIQ(ik) .le. cldmin) then
          tau_cld_temp(liqcld,iw,ik) = 0.0
          singscat_cld_temp(liqcld,iw,ik) = 0.0
          asym_cld_temp(liqcld,iw,ik) = 0.0
        else
          call interpolate_cld(liqcld, iw, r_liq*1.0e6, Qliq, Wliq, Gliq, Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)
          !write(*,*) "interpolate_cld, liquid: r,q,w,g ",ext_REL(ik), Qliq, Wliq, Gliq, ext_cLIQ(ik)
          tau_cld_temp(liqcld,iw,ik) = 3.*Qliq  / (4.*rho_liq*r_liq) * ext_cLIQ(ik)
          singscat_cld_temp(liqcld,iw,ik) = Wliq 
          asym_cld_temp(liqcld,iw,ik) = Gliq 
        endif
    
      enddo
    enddo   

    !
    ! Ice clouds
    !
    do ik=1, pverp
      do iw=iwbeg, iwend

        r_ice = ext_REI(ik) * 1.0e-6   ! ice cloud particle size [m] 

        if (ext_cICE(ik) .le. cldmin) then
          tau_cld_temp(icecld,iw,ik) = 0.0
          singscat_cld_temp(icecld,iw,ik) = 0.0
          asym_cld_temp(icecld,iw,ik) = 0.0
        else          
          call interpolate_cld(icecld, iw, ext_REI(ik), Qice, Wice, Gice, Qcldliq, Qcldice, Wcldliq, Wcldice, Gcldliq, Gcldice)
          !write(*,*) "interpolate_cld, ice: r,q,w,g ",ext_REI(ik), Qice, Wice, Gice, ext_cICE(ik)
          tau_cld_temp(icecld,iw,ik) = 3.*Qice  / (4.*rho_ice*r_ice) * ext_cICE(ik)
          singscat_cld_temp(icecld,iw,ik) = Wice 
          asym_cld_temp(icecld,iw,ik) = Gice 
        endif
 
      enddo
    enddo

    ! Call MCICA sub-column generator (Monte Carlo Independent Column Approximation for cloud overlap)
   
    ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
    icldovr = 2
    nstep = get_nstep()
    permuteseed = nstep

    call mcica_subcol(icldovr, permuteseed, ext_pmid, ext_cFRC, ext_cICE, ext_cLIQ, &
                      tau_cld_temp(icecld,:,:), singscat_cld_temp(icecld,:,:), asym_cld_temp(icecld,:,:), &
                      tau_cld_temp(liqcld,:,:), singscat_cld_temp(liqcld,:,:), asym_cld_temp(liqcld,:,:), &
                      cFRC_mcica, cICE_mcica, cLIQ_mcica, & 
                      tau_mcica_ice, ssa_mcica_ice, asym_mcica_ice, &
                      tau_mcica_liq, ssa_mcica_liq, asym_mcica_liq ) 

    tau_cld_mcica(icecld,:,:) = tau_mcica_ice
    singscat_cld_mcica(icecld,:,:) = ssa_mcica_ice
    asym_cld_mcica(icecld,:,:) = asym_mcica_ice
    tau_cld_mcica(liqcld,:,:) = tau_mcica_liq
    singscat_cld_mcica(liqcld,:,:) = ssa_mcica_liq
    asym_cld_mcica(liqcld,:,:) = asym_mcica_liq
        
    return

  end subroutine calc_cldopd

!============================================================================

  subroutine calc_aeropd( )

!------------------------------------------------------------------------
!
! Purpose: Calculate the optical depths of aerosols
!          Optical depths stored in calculates the current 'tau_aer(1:ip,wavelength_band,vertical_level)',
!          'singscat_albd(1:ip,wavelength_band,vertical_level)', 
!          'asym_fact(1:ip,wavelength_band,vertical_level)' [i.e., layer opacity,
!          single-scattering albedo, and asymmetry factor of aerosol specie 'ip'].
!          
!------------------------------------------------------------------------
! NOTES: Incomplete, hook up with CARMA haze aerosols   

    implicit none
   
!------------------------------------------------------------------------
!
! Input Arguments
!
!
! Local Variables
!

!------------------------------------------------------------------------
!
! Start Code
!
   !loop over wavelengths
   !loop oer levels find optical coefficients
   !find optical depths over particles sizes
   ! find single scatter albedo

    return
  end subroutine calc_aeropd



end module calc_opd_mod
