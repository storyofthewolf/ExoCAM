module circumbinary_mod

!---------------------------------------------------------------------
! Purpose:
!
! Provides interfacing routines for treating circumbinary systems
! Designed to interface with twostars_m.F90 module by S. Eggl
!
! Revisions history
! June 2019, E.T. Wolf
!
! Notes:
!---------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod
  use exoplanet_mod  
  use ppgrid
 
 
  implicit none
  private
  save

  integer, parameter :: dp = selected_real_kind(15, 307)    !precision parameter  
  real(r8), dimension(1:2), public :: semia_adj

!------------------------------------------------------------------------
!
! Public interfaces
!
  public :: circumbinary_addfld_stdoutput
  public :: circumbinary_init_orbit
  public :: circumbinary_set_stellar


!============================================================================
contains
!============================================================================

!============================================================================
!
! Public subroutines
!
!============================================================================

!============================================================================

  subroutine circumbinary_addfld_stdoutput()
!-----------------------------------------------------------------------
!
! Purpose: 
!-----------------------------------------------------------------------

    use cam_history,     only: addfld, add_default, phys_decomp
!------------------------------------------------------------------------
!
! Local Variables
!

!------------------------------------------------------------------------
!
! Start Code
!

    call addfld ('FDSTOA_cb1  ','W/m2    ',1,'A','Shortwave downward flux TOA, star 1',phys_decomp)
    call addfld ('FUS_cb1     ','W/m2    ',pverp,'A','Shortwave upward flux, star 1',phys_decomp)
    call addfld ('FDS_cb1     ','W/m2    ',pverp,'A','Shortwave downward flux, star 1',phys_decomp)
    call addfld ('FUSC_cb1    ','W/m2    ',pverp,'A','Shortwave clear-sky upward flux, star 1',phys_decomp)
    call addfld ('FDSC_cb1    ','W/m2    ',pverp,'A','Shortwave clear-sky downward flux, star 1',phys_decomp)
    call addfld ('FUL_cb1     ','W/m2    ',pverp,'A','Longwave upward flux, star 2',phys_decomp)
    call addfld ('FDL_cb1     ','W/m2    ',pverp,'A','Longwave downward flux, star 2',phys_decomp)
    call addfld ('FULC_cb1    ','W/m2    ',pverp,'A','Longwave clear-sky upward flux, star 2',phys_decomp)
    call addfld ('FDLC_cb1    ','W/m2    ',pverp,'A','Longwave clear-sky downward flux, star 2',phys_decomp)
    call addfld ('QRS_cb1     ','K/day   ',pver, 'A','Solar heating rate, star1',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRL_cb1     ','K/day   ',pver, 'A','Longwave heating rate, star 1',phys_decomp, sampling_seq='rad_lwsw')

    call addfld ('FDSTOA_cb2  ','W/m2    ',1,'A','Shortwave downward flux TOA, star 2',phys_decomp)
    call addfld ('FUS_cb2     ','W/m2    ',pverp,'A','Shortwave upward flux, star 2',phys_decomp)
    call addfld ('FDS_cb2     ','W/m2    ',pverp,'A','Shortwave downward flux, star 2',phys_decomp)
    call addfld ('FUSC_cb2    ','W/m2    ',pverp,'A','Shortwave clear-sky upward flux, star 2',phys_decomp)
    call addfld ('FDSC_cb2    ','W/m2    ',pverp,'A','Shortwave clear-sky downward flux, star 2',phys_decomp)
    call addfld ('FUL_cb2     ','W/m2    ',pverp,'A','Longwave upward flux, star 2',phys_decomp)
    call addfld ('FDL_cb2     ','W/m2    ',pverp,'A','Longwave downward flux, star 2',phys_decomp)
    call addfld ('FULC_cb2    ','W/m2    ',pverp,'A','Longwave clear-sky upward flux, star 2',phys_decomp)
    call addfld ('FDLC_cb2    ','W/m2    ',pverp,'A','Longwave clear-sky downward flux, star 2',phys_decomp)
    call addfld ('QRS_cb2     ','K/day   ',pver, 'A','Solar heating rate, star 2',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRL_cb2     ','K/day   ',pver, 'A','Longwave heating rate, star 1',phys_decomp, sampling_seq='rad_lwsw')

    call addfld ('QRS_test     ','K/day   ',pver, 'A','Solar heating rate, star 1',phys_decomp, sampling_seq='rad_lwsw')

  end subroutine circumbinary_addfld_stdoutput



!============================================================================
!
  subroutine circumbinary_init_orbit()
!-----------------------------------------------------------------------
!
! Purpose: Determine orbital system parameters for circumbinary planet, 
!          fitting for specified mean insolation, binary masses, 
!          and binary separation.  Output system specs

!-----------------------------------------------------------------------

!    use exoplanet_mod,    only: cb_mass, cb_teff, cb_radius, cb_spin, cb_eccen, &
!                                exo_scon
!    use shr_const_mod,    only: SHR_CONST_MSUN, SHR_CONST_LSUN, &
!                                SHR_CONST_AU, &
!                                SHR_CONST_PI, SHR_CONST_GRAVCON

    use cam_logfile,      only: iulog
    use spmd_utils,       only: masterproc
    use twostars_m,       only: insolation3, checkconf
    use error_messages, only: alloc_err
!------------------------------------------------------------------------
!
! Local Variables
!
!
!  real(kind=dp), parameter :: flux_error_spec = 1.0_R8
  real(kind=dp), parameter :: flux_error_spec = 0.1_R8
  real(kind=dp)::e1                                !eccentricity of the stellar binary  
  real(kind=dp),dimension(1:2)::teff,rs            !stellar effective temperatures, stellar radii
  real(kind=dp),dimension(1:2)::semia,man          !semimajor axes of binary and planetary orbits, initial mean longitudes
  real(kind=dp),dimension(1:3)::mass               !system    masses: m(1:2) stars, m(3) planet  
  real(kind=dp),dimension(1:4)::spin               !spin    quantities: (1) obliquity [deg]
                                              !                    (2) angle of precession [deg] 
                                              !                    (3) rotation    period [D]
                                              !                    (4) initial hour angle [deg] (=initial position of continents wrt x-axis) 
  real(kind=dp), dimension(1:2)::lum               ! stellar luminosities
  real(kind=dp)::time,tend,dt                      ! time, end time, output timestep
 
  !results    
  integer::ok                                                     !configuration check parameter: 0-> everything is ok, else...                          
  real(kind=dp), dimension(1:2)::decl,phi                          !declination and hour angle wrt each star                                              
  real(kind=dp), dimension(1:2)::d2,tf                            !distance squared factor, and transit factor                                           
  real(kind=dp), dimension(1:2)::TSI
  real(kind=dp), dimension(1:2)::semia_out
  real(kind=dp) :: TSI_total, TSI_mean
  real(kind=dp) :: error, estepm, estepp
  real(kind=dp) :: rcm, fac, fac2, pbin, pplanet, tend_p
  real(kind=dp) :: e2  ! planet eccentricity 
  real(kind=dp) :: semia_guess, norb
  real(kind=dp) :: maxv, minv, ampl
  real(kind=dp) :: omegaplanet, omegabinary, padj

  real(kind=dp), allocatable, dimension(:) :: TSI_1_vec     ! time vector of stellar insolation from star 1
  real(kind=dp), allocatable, dimension(:) :: TSI_2_vec     ! time vector of stellar insolation from star 2
  real(kind=dp), allocatable, dimension(:) :: TSI_total_vec ! time vector of stellar insolation from star 1 + star 2
  real(kind=dp), allocatable, dimension(:,:) :: decl_vec    ! time vector of stellar declinations
  real(kind=dp), allocatable, dimension(:,:) :: phi_vec     ! time vector of stellar hour angles

  integer::num,ti, eloop_count, istat


!------------------------------------------------------------------------
!
! Start Code
!
  !---------------------------------------------------------------
  !set local variables for orbit iteration                                               
  !---------------------------------------------------------------
  lum(:) = cb_lum(:)
  mass(:) = cb_mass(:)
  mass(3) = mass(3) * SHR_CONST_MEARTH/SHR_CONST_MSUN
  semia(:) = cb_semia(:)
  man(:) = cb_man(:)
  spin(:) = cb_spin(:)
  teff(:) = cb_teff(:)
  rs(:) = cb_radius(:)
  e1 = cb_eccen

  ! Make initial guess at semimajor axis
  semia_guess = sqrt ( (lum(1)+lum(2))*SHR_CONST_LSUN / (4.*SHR_CONST_PI*exo_scon))
  semia_out(1) = semia(1)
  if (masterproc) write(*,*) "semi-major axis initial guess: ", semia_guess/SHR_CONST_AU
  semia_out(2) = semia_guess/SHR_CONST_AU

  pplanet = sqrt( 4.*SHR_CONST_PI**2 / (SHR_CONST_GRAVCON * (SHR_CONST_MSUN*(mass(1)+mass(2)) + SHR_CONST_MEARTH*mass(3))) &
                  * (SHR_CONST_AU*semia_out(2))**3.) / SHR_CONST_CDAY

  pbin = sqrt( 4.*SHR_CONST_PI**2 / (SHR_CONST_GRAVCON* (SHR_CONST_MSUN*(mass(1)+mass(2)))) &
                  * (SHR_CONST_AU*semia_out(1))**3.) / SHR_CONST_CDAY

  if (masterproc) then
     write(iulog,*) '************************************************************'
     write(iulog,*) '************************************************************'
     write(iulog,*) '***               CIRCUMBINARY_ORBIT_INIT                ***'
     write(iulog,*) '************************************************************'
     write(iulog,*) '************************************************************'
     write(iulog,*) 'STELLAR MASS 1: (Msun) ', cb_mass(1), cb_mass(1)*SHR_CONST_MSUN
     write(iulog,*) 'STELLAR MASS 2: (Msun) ', cb_mass(2), cb_mass(2)*SHR_CONST_MSUN
     write(iulog,*) 'PLANET MASS: (Mearth) ', cb_mass(3), cb_mass(3)*SHR_CONST_MEARTH
     write(iulog,*) 'BINARY SEPARATION: (AU) ', semia_adj(1)
     write(iulog,*) 'BINARY PERIOD (days) ', Pbin
     write(iulog,*) 'BINARY ECCENTRICITY: ', cb_eccen
     write(iulog,*) 'LUMINOSITY 1: (Lsun) ', cb_lum(1), cb_lum(1)*SHR_CONST_LSUN
     write(iulog,*) 'LUMINOSITY 2: (Lsun) ', cb_lum(2), cb_lum(2)*SHR_CONST_LSUN
     write(iulog,*) 'TEMPERATURE 1: (K) ', cb_teff(1)
     write(iulog,*) 'TEMPERATURE 2: (K) ', cb_teff(2)
     write(iulog,*) 'RADIUS 1: (Rsun) ', cb_radius(1)
     write(iulog,*) 'RADIUS 2: (Rsun) ', cb_radius(2)
     write(iulog,*) '***********************************'
     write(iulog,*) '  orbit fitting ... '
     write(iulog,*) '***********************************'
  endif

  ! Start iteration loop
  error = 100000.
  eloop_count=0
  time=0.0
  dt=0.1
  do while (error > flux_error_spec)
 
    call checkconf(mass,semia_out,e1,ok)
    if (masterproc) write(*,*)'configuration ok?',ok

    !if(ok.gt.0.and.ok.le.4) 
    if (ok .eq. 0) then  ! orbital configuration valid

      ! determine number planet orbits to averaged over
      !norb = 3
      norb=1
      do while (mod(norb*pplanet,pbin) .ge. 1) 
        norb = norb + 1
        !write(*,*) norb*pplanet,pbin,mod(norb*pplanet,pbin)
      enddo                                        
!      if (norb .gt. 100) norb=100  ! limit max orbits to integrate

      tend_p = int(pplanet*norb)
      num = int(tend_p)/dt
!      if(masterproc) write(*,*) pplanet*norb, int(pplanet*norb), tend_p/dt, int(tend_p/dt) 
!      if(masterproc) write(*,*) "norb, num: ", norb, num
      allocate(TSI_1_vec(num))
      allocate(TSI_2_vec(num))
      allocate(TSI_total_vec(num))
      allocate(decl_vec(2,num))
      allocate(phi_vec(2,num))

      pplanet = sqrt( 4.*SHR_CONST_PI**2 / (SHR_CONST_GRAVCON * (SHR_CONST_MSUN*(mass(1)+mass(2)) + SHR_CONST_MEARTH*mass(3))) &
                        * (SHR_CONST_AU*semia_out(2))**3.) / SHR_CONST_CDAY

      pbin = sqrt( 4.*SHR_CONST_PI**2 / (SHR_CONST_GRAVCON* (SHR_CONST_MSUN*(mass(1)+mass(2)))) &
                       * (SHR_CONST_AU*semia_out(1))**3.) / SHR_CONST_CDAY

      omegabinary = sqrt(SHR_CONST_GRAVCON*(SHR_CONST_MSUN*(mass(1)+mass(2))) / (semia_out(1)*SHR_CONST_AU)**3)
      omegaplanet = sqrt(SHR_CONST_GRAVCON*(SHR_CONST_MSUN*(mass(1)+mass(2)) + SHR_CONST_MEARTH*mass(3)) / (semia_out(2)*SHR_CONST_AU)**3)
      padj = 2.*SHR_CONST_PI/(omegabinary-omegaplanet) / (24.*60.*60.)

      ! loop twostars calculation over norbits    
       ti = 1
       time=0.0
 !     do while(time.lt.tend_p)
       do while(ti.le.num)


        call insolation3(mass,teff,rs,semia_out,e1,man,spin,time, &
                         d2,tf,decl,phi,e2)

! transits turned off for orbit fitting
        if (cb_transit .eqv. .false.) then
          tf(1) = 1.0
          tf(2) = 1.0
        end if

        TSI(1) = lum(1) * SHR_CONST_LSUN / (4.*SHR_CONST_PI*d2(1) * SHR_CONST_AU**2 ) * tf(1)
        TSI(2) = lum(2) * SHR_CONST_LSUN / (4.*SHR_CONST_PI*d2(2) * SHR_CONST_AU**2 ) * tf(2)
        TSI_total = TSI(1) + TSI(2)
        TSI_1_vec(ti) = TSI(1)

        TSI_2_vec(ti) = TSI(2)
        TSI_total_vec(ti) = TSI(1) + TSI(2)
        decl_vec(:,ti) = decl
        phi_vec(:,ti) = phi

        time=time+dt
        ti=ti+1
!        if(masterproc) write(*,*) "time: ", time, ti
      end do
 
     ! compute TSI statistics 
      TSI_mean = SUM(TSI_total_vec)/size(TSI_total_vec)
      maxv = maxval(TSI_total_vec(:))
      minv = minval(TSI_total_vec(:))
      ampl = maxv-minv
      error = abs(TSI_mean-exo_scon)


      if (masterproc) write(iulog,*) eloop_count,  norb, error, TSI_mean, maxv, minv, ampl,pbin, pplanet, padj
      if (error > flux_error_spec) then
        ! step semimajor axis
        fac=1.
        if (abs(exo_scon - TSI_mean) .le. 100.) fac=5.
        if (abs(exo_scon - TSI_mean) .le. 10.) fac=50.
        if  (abs(exo_scon - TSI_mean) .le. 1.) fac=100.
        if ((exo_scon - TSI_mean) .gt. 0)  fac2 = -1.0
        if ((exo_scon - TSI_mean) .le. 0)  fac2 = 1.0
        estepm=1.0-0.01/fac
        estepp=1.0+0.01/fac
        rcm = semia_out(2)
        if (TSI_mean .le. exo_scon) rcm = rcm*estepm
        if (TSI_mean .ge. exo_scon) rcm = rcm*estepp
        semia_out(2) = rcm
!        semia_out(2) = semia_out(2) + fac*fac2
        eloop_count = eloop_count+1
!        if (masterproc) write(*,*) "semia_out", semia_out
      end if      

      deallocate(TSI_1_vec)
      deallocate(TSI_2_vec)
      deallocate(TSI_total_vec)
      deallocate(decl_vec)
      deallocate(phi_vec)    
    end if

  end do
  

  semia_adj(:) = semia_out(:)




  if (masterproc) then
     write(iulog,*) '***********************************'
     write(iulog,*)  'Target Flux: (Wm-2) ', exo_scon
     write(iulog,*)  'Planet Semimajor Axis (AU) ', semia_adj(2)
     write(iulog,*)  'Planet Period (days) ', pplanet
     write(iulog,*)  'Planet Eccentricity ', e2
     write(iulog,*)  'Mean Flux: (Wm-2) ', TSI_mean
     write(iulog,*)  '  max TSI: (Wm-2) ', maxv
     write(iulog,*)  '  min TSI: (Wm-2) ', minv
     write(iulog,*)  '  amplitude TSI: (Wm-2, %) ', ampl, ampl/TSI_mean*100. 
     write(iulog,*)  'Planet Omega (rad/s) ', omegaplanet
     write(iulog,*)  'Binary Omega (rad/s) ', omegabinary
     write(iulog,*)  'Adjusted Binary Period: ', padj
     write(iulog,*) '******************************************************'
     write(iulog,*) '******************************************************'
  endif

  end subroutine circumbinary_init_orbit


!============================================================================

  subroutine circumbinary_set_stellar(whichstar, d2, tf)
!-----------------------------------------------------------------------
!
! Purpose: Sets per timestep properties for the binary system.
!          Called separately for each star, before aerad_driver.
!
!-----------------------------------------------------------------------

    use exo_init_ref !, only: gw_solflux
    use radgrid !,      only: S0, S0_2, solarflux, solarflux_2, ntot_wavlnrng
    use spmd_utils,      only: masterproc    

!-----------------------------------------------------------------------
! Input Arguments
!
    integer, intent(in) :: whichstar ! toggle star 1 or star 2
    real(kind=dp), dimension(1:2), intent(in) :: d2, tf        ! toggle star 1 or star 2

!------------------------------------------------------------------------
!
! Local Variables
!
    integer :: iq
    integer :: iw
    integer :: ig
    integer :: ip
    real(kind=dp) :: TSI_single
    real(kind=dp), dimension(ntot_wavlnrng) :: solarflux_working
!------------------------------------------------------------------------
!
! Start Code
!
   TSI_single = 0.0
   if (whichstar == 1) then
     TSI_single = cb_lum(whichstar) * SHR_CONST_LSUN / (4.*SHR_CONST_PI*d2(whichstar) * SHR_CONST_AU**2 ) * tf(whichstar)
!     if(masterproc) write(*,*), "TSI_single 1?: ", cb_lum(whichstar), d2(whichstar)
!     if(masterproc) write(*,*), "TSI_single 1: ", TSI_single
     solarflux_working(:) = solarflux(:)*TSI_single/S0 ! TSI in column is weighted for the star
     !solarflux(:) = solarflux(:)*scon/S0
   endif
   if (whichstar == 2) then
     TSI_single = cb_lum(whichstar) * SHR_CONST_LSUN / (4.*SHR_CONST_PI*d2(whichstar) * SHR_CONST_AU**2 ) * tf(whichstar)
!     if(masterproc) write(*,*), "TSI_single 2?: ", cb_lum(whichstar), d2(whichstar)
!     if(masterproc) write(*,*), "TSI_single 2: ", TSI_single
     solarflux_working(:) = solarflux_2(:)*TSI_single/S0_2 ! TSI in column is weighted for the star
     !solarflux(:) = solarflux(:)*scon/S0
   endif

   !if (which-star != 2 and whichstar != 1) then begin
   !  !call ENDRUN: Critical error in circumbinary module.
   !endif

  !
  ! Calculate the "average" wavenumber (1/wavelength) <wavenum()> for each
  !  wavelength interval in each wavelength group [1/cm]:
  !
  iq = 0
  ip = lw_ipbeg-1
  do iw=1,ntot_wavlnrng  ! "avg" wavenumber over each band
    do ig=1,ngauss_pts(iw) 
      ip = ip+1
      ! Gauss-weighted solar flux in each probability interval:
      gw_solflux(ip) = solarflux_working(iw)*g_weight(ip)  
      !write(*,*) "init_ref,g_weight",g_weight(ip)
    enddo
  enddo

!   if (masterproc) then
!      write(*,*) "whichstar: ", whichstar
!      write(*,*) " circumbinary_set_stellar: total solar irradiance"
!      write(*,*) "solar flux [W m-2] in each spectral interval"
!      do iw=lw_iwbeg, sw_iwend
!        write(*,*) iw, solarflux_working(iw)
!      enddo
!      write(*,*) "TOTAL SOLAR FLUX:", SUM(solarflux_working)
!    endif


  end subroutine circumbinary_set_stellar

end module circumbinary_mod
