module twostars_m
!#############################################################################################################################################
!Routines to simulate the insolation and orbital evolution of a coplanar circumbinary planet. 
!
! constructed by Siegfried Eggl 	20140609
!	
! last modified 20150304
!
! modification history
! update of secular constant K3 to include mutual perturbation and Post Newtonian Correction
! inconsistencies in real and integer variable declarations resolved 
!
!
!##############################################################################################################################################

implicit none
!##############################################################################################################################################
!			Module Variables
!##############################################################################################################################################
PRIVATE
!REAL PRECISION (single,double,quadruple) 
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)  
  integer, parameter :: qp = selected_real_kind(33, 4931)
!INTEGER PRECISION 
  integer, parameter :: isp= selected_int_kind(4)
  integer, parameter :: idp= selected_int_kind(8)

!FIXED CONSTANTS
  real(kind=dp),parameter::pi=3.14159265358979323846264338328_dp	![] PI
  real(kind=dp),parameter::pix2=6.28318530717958647692528676656_dp	![] PI x 2
  real(kind=dp),parameter::pib2=1.57079632679489661923132169164_dp	![] PI / 2
  real(kind=dp),parameter::pi3b2=4.71238898038468985769396507492_dp     ![] 3 PI / 2
  real(kind=dp),parameter::deg2rad=0.0174532925199432957692369076849_dp ![1/deg] PI / 180
  real(kind=dp),parameter::k=0.01720209895_dp				![au^(3/2) Msun^(-1/2) D^(-1)] Gaussian gravitational constant
  real(kind=dp),parameter::clight=173.1446326743_dp 			![au/D] vacuum speed of light
  real(kind=dp),parameter::rsun=696342.E3     					![m] mean radius of the Sun (SOHO)
  real(kind=dp),parameter::rsunau=0.00465475876589465387357281416172_dp   	![au] mean radius of the Sun (SOHO) in au
  real(kind=dp),parameter::au=149597870700._dp 					![m] astronomical unit 
  real(kind=dp),parameter::teffsun=5777._dp    					![K] effective temperature of the Sun
  real(kind=dp),parameter::zsun=0.02_dp        					![]  metallicity of the sun
 
!##############################################################################################################################################
!			Public Subroutines
!##############################################################################################################################################


!--------------------------------------------------------------------------------------------------------------	       
!	       !CONFIGURATION CHECK (DYNAMICAL STABILITY, ETC.)
!--------------------------------------------------------------------------------------------------------------	


public::checkconf !check whether the binary planet configuration is suitable for the current propagation


!--------------------------------------------------------------------------------------------------------------	       
!	       !INSOLATION SUBROUTINES
!--------------------------------------------------------------------------------------------------------------	


public::insolation1   !Evolve the binary + planet and determine the top-of the atmosphere insolation a circumbinary (P-type) planet receives in specified wave-length bands

public::insolation2  !Same as insolation 1 but with direct input of stellar effective temperatures and radii

public::insolation3  !Same as insolation 1 but with direct input of stellar effective temperatures and radii, for ExoCAM linkage



!--------------------------------------------------------------------------------------------------------------	       
!	       !DISTANCE, POSITION AND ECCENTRICITY EVOLUTION SUBROUTINES
!--------------------------------------------------------------------------------------------------------------	

public::pconst !generates constants necessary for describing the outer (planetary) eccentricity evolution in P-type systems
               !HAS TO BE CALLED BEFORE "PDIST"

public::pdist  !distance vector, eccentricity and argument of pericenter evolution of a P-type circumbinary system


!--------------------------------------------------------------------------------------------------------------	       
!	       !MAXIMUM AND AVERAGE SQUARED ECCENTRICITY SUBROUTINES
!--------------------------------------------------------------------------------------------------------------	
 
public::e2av   !<e^2> of the outer orbit in a hierarchical triple system
	       !assumptions: inner orbit = const (only valid for planetary masses!), system coplanar

public::e2max  ! maximum eccentricity of the outer orbit in a hierarchical triple system
	       !assumptions: inner orbit = const (only valid for planetary masses!), system coplanar


 contains

!####################################################################################################################################### 
subroutine checkconf(mass,semia,e1,ok)
!-----------------------------------------------------------------------------------------------------------------------------------------
! Check whether the suggested binary - planet configuration makes sense and the dynamical theory is indeed applicable  
!
! written by Siegfried Eggl 20140609
!
! dependencies: hw99p
!-----------------------------------------------------------------------------------------------------------------------------------------
! Input
!
! real::
! mass(1:3)      ... [Msun] 	m0, m1,m2: 	masses of the primary star, secondary star, and the planet
! semia(1:2)     ... [au]   	a1, a2: 	semimajor axes of the binary (inner) and planetary (outer) orbit. a1<a2  
! e1             ... [] 			eccentricity of binary star's orbit  0 <= e1 < 1 / the planet will always start on a circular orbit      
!
!-----------------------------------------------------------------------------------------------------------------------------------------
! Output
!
! integer::
! ok            ...[]   0: all fine, 1: error due to mass arrangement, 2: error due to choice of semimajor axis, 3: error due to closeness of binary
!------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 !input
 real(kind=dp),intent(in)::e1
 real(kind=dp),dimension(1:2),intent(in)::semia
 real(kind=dp),dimension(1:3),intent(in)::mass

 !output
 integer::ok
 !local
 integer::i
 real(kind=dp)::ac

 !everything is fine
 ok=0
 
 !check whether the analytical propagation is applicable
 if(mass(3).gt.sum(mass(1:2)).or.mass(3).gt.mass(2).or.mass(3).gt.mass(1)) then
!  write(*,*)'The system configuration you have chosen is no binary star with a circumbinary planet.', &
!  'Please alter the masses of the system so that the outer mass: mass(3) is smaller than the stellar masses: mass(1:2).'
  ok=1
 end if
  
 !check whether this configuration is dynamically stable
 call hw99p(mass(1),mass(2),semia(1),e1,ac)
 if(semia(2).le.ac) then
!  write(*,*)'The orbital configuration you suggested is dynamically unstable.', & 
!  'Please move the planet to a semimajor axis of at least ',ac,' [au].'
 ok=2
 end if
 
 !check whether this is a tidally dominated or contact binary
 if(semia(1)*(1._dp-e1).lt.0.01_dp) then 
!   write(*,*)'The orbital theory used in this module does not support compact binaries.', & 
!  'The pericenter distance q=a*(1-e) of the binary should not be smaller than 0.01 [au].'
 ok=3
 end if
 
 
 !check for negative quantities 
 if(minval(mass(:)).lt.0._dp.or.minval(semia(:)).lt.0._dp.or.e1.lt.0.or.e1.ge.1._dp) then 
 !  write(*,*)'All masses and semimajor axes must be >=0! The binary stars eccentricity must be in [0, 1['
 ok=4
 end if 
 
   !check whether stellar masses are in reasonable region
 do i=1,2  
 if(mass(i).lt.0.1_dp.or.mass(i).gt.100._dp) then 
 !  write(*,*)'Warning: effective temperatures and stellar radii can not be interpolated reasonably for stars with masses', &
 !            'outside [0.1-100] solar masses. Please use the routine insolation 2 and specify effective temperatures and ',&
 !            'stellar radii explicitly!'
 ok=5
 end if 
 end do
 
return
end subroutine
 
!*****************************************************************************************************************************************
subroutine insolation1(mass,met,semia,e1,man,spin,t,nbands,bands,bndflux,decl,phi) 
!-----------------------------------------------------------------------------------------------------------------------------------------
! Determine the top-of the atmosphere insolation a circumbinary (P-type) planet receives in 
! specified wave-length bands
!
!*)stellar effective temperature and radii are derived via observational fits by Tout et al. 1996
! 
!*)planetary and stellar dynamics are solved analytically using relations from Georgakarakos & Eggl 2015
!
!*)spectral fluxes are integrated following Widger & Woodall 1976, Planck distribution assumed
!
! written by Siegfried Eggl 20140609
!
! dependencies: pconst, pdist, getspin, planetorient, zamsmrlt, intplanckwn, transitfactor
!-----------------------------------------------------------------------------------------------------------------------------------------
! Input
!
! real::
! mass(1:3)      		... [Msun] 	m0, m1,m2: 	masses of the primary star, secondary star, and the planet
! met(1:2)       		... []		Z1,Z2		metallicity of each star (0.0001-0.03); best stick to solar value 0.02
! semia(1:2)     		... [au]   	a1, a2: 	semimajor axes of the binary (inner) and planetary (outer) orbit. 0<a1<a2  
! e1             		... [] 			        eccentricity of binary star's orbit  0 <= e1 < 1 / the planet will always start on a circular orbit      
! man(1:2)       		... [deg]	ma10, ma20:	initial positions of the binary/planet on its orbit (initial mean longitudes) range [0,360[ deg
! spin(1:4)      		... [mixed] 	(1) obl [deg]	obliquity (ref: z-axis = normal to the orbital plane), range [0,180[ deg
! 						(2) pre [deg]	angle of precession (ref: x-axis = along initial pericenter of binary)  range [0,360[ deg
! 						(3) pp  [D]     rotation period of the planet
!						(4) eta0 [deg]  initial hour angle of zero-meridian of planet wrt x-axis, range [0,360[ deg
! t		 		... [D]				current time in (Gaussian) days 
! bands(1:2,1:nbands,1:2)     	... [cm^-1]			wave numbers of insolation band borders (first index: star, second index: band number, third index: lower and upper band limit wave numbers)
!
! integer::
! nbands			... []				number of spectral bands
!
!
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
! Output
!
! real::
! bndflux(1:2,1:nbands)		[W m^-2 sr^-1]	...top of the atmosphere flux on the planet for each of the two stars in each band in Watt per square meter per solid angle
! decl(1:2)                    	[rad]          	... declination of the planet wrt to each star (angle between equator and direction vector to each star)
! phi (1:2)      		[rad]		... hour angle (zero-meridian vector to direction vector of the respective star)
! 
!----------------------------------------------------------------------------------------------------------------------------------------------------------------- 
 implicit none
 !input
 integer,intent(in)::nbands
 real(kind=dp),intent(in)::e1,t
 real(kind=dp),dimension(1:2),intent(in)::met,semia,man
 real(kind=dp),dimension(1:3),intent(in)::mass
 real(kind=dp),dimension(1:4),intent(in)::spin
 real(kind=dp),dimension(1:2,1:nbands,1:2),intent(in)::bands
 
 !output
 real(kind=dp),dimension(1:2),intent(out)::decl,phi
 real(kind=dp),dimension(1:2,1:nbands),intent(out)::bndflux
 
 !local
 integer::i,j
 real(kind=dp),dimension(1:2)::r0,r1,r2 	!stellar and planetary positions at time t [au]
 real(kind=dp),dimension(1:2)::rs,teff,lum      !stellar radii [Rsun], effective temperature[K] and luminosity[Lsun]
 real(kind=dp)::cs(1:6),c(1:14) 		!dynamical constants
 real(kind=dp)::e2,w1,w2        		!eccentricity of the outer orbit,argument of pericenter of the inner and outer orbit [rad]
 real(kind=dp)::koos(1:2,1:3),koop(1:3)   	!vectors for transit factor determination
 real(kind=dp),dimension(1:2)::d2,dr1,dr2  	!distances between the planet and the stars [au], relative vectors
 real(kind=dp),dimension(1:4)::spinrad          !initial spin vector with angles in [rad] 
 real(kind=dp),dimension(1:2)::tf,df 		!transit and distance insolation diminishing factors for each star
 
 !---------------------------------------------------------------------
 !dynamical evolution of the planet and the binary
 !---------------------------------------------------------------------
 call pconst(semia(1),semia(2),e1,man(1)*deg2rad,man(2)*deg2rad,mass(1),mass(2),mass(3),c,cs)
 call pdist(semia(1),semia(2),e1,mass(1),mass(2),mass(3),c,cs,t,.true.,0._dp,r0,r1,r2,e2,w1,w2)
 
 !get the 3D planetary spin vector in Cartesian coordinates from obliquity, angle of precession and rotation period
 spinrad(1:2)=spin(1:2)*deg2rad !obliquity and angle of precession [deg] -> [rad]
 spinrad(3)=spin(3) !rotation period
 spinrad(4)=spin(4)*deg2rad !initial hour angle of zero meridian [deg] -> [rad]
 
 !determine declination decl and hour angle phi wrt each star
 call planetorient(t,spinrad,r0,r1,r2,dr1,dr2,decl,phi)

 !---------------------------------------------------------------------
 !insolation
 !---------------------------------------------------------------------
 !calculate flux in wave number bands
 do j=1,2 !number of stars
 !find effective temperatures and radii for stars via Zero Age Main Sequence fit
 call zamsmrlt(mass(j),met(j),lum(j),teff(j),rs(j))
 !integrate Planck curve to get band flux
  do i=1,nbands
   call intplanckwn(teff(j)*teffsun,bands(j,i,1:2),bndflux(j,i))
  end do
 end do 
 
 !calculate transit factor (whether one star is hidden behind the other)
 koos(:,3)=0.d0
 koop(3)=0.d0
 koos(1,1:2)=r0
 koos(2,1:2)=r1
 koop(1:2)=r2
 call transitfactor(2,koos,koop,rs,tf)
 
 !calculate the squares of the relative distances between binary stars and the planet
 d2(1)=dot_product(dr1,dr1)
 d2(2)=dot_product(dr2,dr2)

 !diminish flux due to transit and distance 
 do i=1,2
  df(i)=pi*(rs(i)*rsunau)**2/d2(i)              !distance diminishing factor
  bndflux(i,:)=bndflux(i,:)*tf(i)*df(i)   !diminish flux by transit and distance factors
 end do
 
 return
end subroutine

!********************************************************************************************************************************************
subroutine insolation2(mass,teff,rs,semia,e1,man,spin,t,nbands,bands,bndflux,decl,phi) 
!-----------------------------------------------------------------------------------------------------------------------------------------
! Determine the top-of the atmosphere insolation a circumbinary (P-type) planet receives in 
! specified wave-length bands
!
! This routine lets you specify mass, effective temperature and solar radii for each star directly.
! 
!*)planetary and stellar dynamics are solved analytically using relations from Georgakarakos & Eggl 2014
!
!*)spectral fluxes are integrated following Widger & Woodall 1976
!
! written by Siegfried Eggl 20140609
!
! dependencies: pconst, pdist, getspin, planetorient, intplanckwn, transitfactor
!-----------------------------------------------------------------------------------------------------------------------------------------
! Input
!
! real::
! mass(1:3)      		... [Msun] 	m0, m1, m2: 	masses of the primary star, secondary star, and the planet
! teff(1:2)      		... [K] 			effective temperatures of the stars
! rs(1:2)        		... [m] 			stellar (mean) radii 
! semia(1:2)     		... [au]   	a1, a2: 	semimajor axes of the binary (inner) and planetary (outer) orbit. 0<a1<a2  
! e1             		... [] 				eccentricity of binary star's orbit  0 <= e1 < 1 / the planet will always start on a circular orbit      
! man(1:2)       		... [rad]	ma10, ma20:	inital positions of the binary/planet on its orbit (initial mean longitudes), range [0,360[ deg
! spin(1:4)      		... [mixed] 	(1) obl [deg]	obliquity (ref: z-axis = normal to the orbital plane), range [0,180[ deg
! 						(2) pre [deg]	angle of precession (ref: x-axis = along initial pericenter of binary), range [0,360[ deg
! 						(3) pp  [D]     rotation period of the planet
!						(4) eta0 [deg]  initial hour angle of zero-meridian of planet wrt x-axis,   range [0,360[ deg
! t		 		... [D]				current time in (Gaussian) days 
! bands(1:2,1:nbands,1:2)       ... [cm^-1]			wave numbers of insolation band borders (first index: star, second index: band number, third index: lower and upper band limit wave numbers)
! 
! integer::
! nbands			... []				number of spectral bands
!
! 
!
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
! Output
! 
! real::
! bndflux(1:2,1:nbands)		[W m^-2 sr^-1]		...top of the atmosphere flux on the planet for each of the two stars in each band in Watt per square meter per solid angle
! decl(1:2)                    	[rad] 		        ... declination of the planet wrt to each star (angle between equator and direction vector to each star)
! phi (1:2)      		[rad]			... hour angle (zero-meridian vector to direction vector of the respective star)
! 
!----------------------------------------------------------------------------------------------------------------------------------------------------------- 
 implicit none
 !input
 integer,intent(in)::nbands
 real(kind=dp),intent(in)::e1,t
 real(kind=dp),dimension(1:2),intent(in)::teff,rs,semia,man
 real(kind=dp),dimension(1:3),intent(in)::mass
 real(kind=dp),dimension(1:4),intent(in)::spin
 real(kind=dp),dimension(1:2,1:nbands,1:2),intent(in)::bands
 
 !output
 real(kind=dp),dimension(1:2),intent(out)::decl,phi
 real(kind=dp),dimension(1:2,1:nbands),intent(out)::bndflux
 
 !local
 integer::i,j
 real(kind=dp),dimension(1:2)::r0,r1,r2 	!stellar and planetary positions at time t [au]
 real(kind=dp)::cs(1:6),c(1:14) 		!dynamical constants
 real(kind=dp)::e2,w1,w2        		!eccentricity of the outer orbit,argument of pericenter of the inner and outer orbit [rad]
 real(kind=dp)::koos(1:2,1:3),koop(1:3)   	!vectors for transit factor determination
 real(kind=dp),dimension(1:2)::d2,dr1,dr2  	!distances between the planet and the stars [au], relative vectors
 real(kind=dp),dimension(1:4)::spinrad 		!initial spin vector with angles in [rad] 
 real(kind=dp),dimension(1:2)::tf,df 		!transit and distance insolation diminishing factors for each star
 
 !---------------------------------------------------------------------
 !dynamical evolution of the planet and the binary
 !---------------------------------------------------------------------
 call pconst(semia(1),semia(2),e1,man(1)*deg2rad,man(2)*deg2rad,mass(1),mass(2),mass(3),c,cs)
 call pdist(semia(1),semia(2),e1,mass(1),mass(2),mass(3),c,cs,t,.true.,0._dp,r0,r1,r2,e2,w1,w2)

 !get the 3D planetary spin vector in Cartesian coordinates from obliquity, angle of precession and rotation period
 spinrad(1:2)=spin(1:2)*deg2rad !obliquity and angle of precession [deg] -> [rad]
 spinrad(3)=spin(3) !rotation period
 spinrad(4)=spin(4)*deg2rad !initial hour angle of zero meridian [deg] -> [rad]
 
 !determine declination decl and hour angle phi wrt each star
 call planetorient(t,spinrad,r0,r1,r2,dr1,dr2,decl,phi)
 
 !---------------------------------------------------------------------
 !insolation
 !---------------------------------------------------------------------
 !calculate the integral over Planck curve bands
 do j=1,2 !number of stars 
  do i=1,nbands
   call intplanckwn(teff(j),bands(j,i,1:2),bndflux(j,i))
  end do
 end do 
 

!calculate transit factor (whether one star is hidden behind the other)
 koos(:,3)=0.d0
 koop(3)=0.d0
 koos(1,1:2)=r0
 koos(2,1:2)=r1
 koop(1:2)=r2
 call transitfactor(2,koos,koop,rs,tf)
 
  
!calculate the squares of the relative distances between binary stars and the planet
d2(1)=dot_product(dr1,dr1)
d2(2)=dot_product(dr2,dr2)

!diminish flux due to transit and distance 
 do i=1,2
  df(i)=pi*(rs(i)*rsunau)**2/d2(i)              !distance diminishing factor
  bndflux(i,:)=bndflux(i,:)*tf(i)*df(i)   !diminish flux by transit and distance factors
 end do
 
 return
end subroutine


!********************************************************************************************************************************************
subroutine insolation3(mass,teff,rs,semia,e1,man,spin,t,d2,tf,decl,phi,e2) 
!-----------------------------------------------------------------------------------------------------------------------------------------
! Determine the top-of the atmosphere insolation a circumbinary (P-type) planet receives in 
! specified wave-length bands
!
! This routine lets you specify mass, effective temperature and solar radii for each star directly.
! 
!*)planetary and stellar dynamics are solved analytically using relations from Georgakarakos & Eggl 2014
!
!*)spectral fluxes are integrated following Widger & Woodall 1976
!
! written by Siegfried Eggl 20140609
! modified by Eric T Wolf 20190711 to jive with ExoCAM
!
! dependencies: pconst, pdist, getspin, planetorient, intplanckwn, transitfactor
!-----------------------------------------------------------------------------------------------------------------------------------------
! Input
!
! real::
! mass(1:3)      		... [Msun] 	m0, m1, m2: 	masses of the primary star, secondary star, and the planet
! teff(1:2)      		... [K] 			effective temperatures of the stars
! rs(1:2)        		... [m] 			stellar (mean) radii 
! semia(1:2)     		... [au]   	a1, a2: 	semimajor axes of the binary (inner) and planetary (outer) orbit. 0<a1<a2  
! e1             		... [] 				eccentricity of binary star's orbit  0 <= e1 < 1 / the planet will always start on a circular orbit      
! man(1:2)       		... [rad]	ma10, ma20:	inital positions of the binary/planet on its orbit (initial mean longitudes), range [0,360[ deg
! spin(1:4)      		... [mixed] 	(1) obl [deg]	obliquity (ref: z-axis = normal to the orbital plane), range [0,180[ deg
! 						(2) pre [deg]	angle of precession (ref: x-axis = along initial pericenter of binary), range [0,360[ deg
! 						(3) pp  [D]     rotation period of the planet
!						(4) eta0 [deg]  initial hour angle of zero-meridian of planet wrt x-axis,   range [0,360[ deg
! t		 		... [D]				current time in (Gaussian) days 
!
! 
!
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
! Output
! 
! real::
! df(1:2) distance factor
! tf(1:2) transift factor
! decl(1:2)                    	[rad] 		        ... declination of the planet wrt to each star (angle between equator and direction vector to each star)
! phi (1:2)      		[rad]			... hour angle (zero-meridian vector to direction vector of the respective star)
! e2                                                    ... planet eccentricity
!
!----------------------------------------------------------------------------------------------------------------------------------------------------------- 
 implicit none
 !input
 !!! integer,intent(in)::nbands
 real(kind=dp),intent(in)::e1,t
 real(kind=dp),dimension(1:2),intent(in)::teff,rs,semia,man
 real(kind=dp),dimension(1:3),intent(in)::mass
 real(kind=dp),dimension(1:4),intent(in)::spin
 !!! real(kind=dp),dimension(1:2,1:nbands,1:2),intent(in)::bands
 
 !output
 real(kind=dp),dimension(1:2),intent(out)::decl,phi, d2, tf
 real(kind=dp),intent(out)::e2

 !!! real(kind=dp),dimension(1:2,1:nbands),intent(out)::bndflux
 
 !local
 integer::i,j
 real(kind=dp),dimension(1:2)::r0,r1,r2 	!stellar and planetary positions at time t [au]
 real(kind=dp)::cs(1:6),c(1:14) 		!dynamical constants
 real(kind=dp)::w1,w2            		!argument of pericenter of the inner and outer orbit [rad]
 real(kind=dp)::koos(1:2,1:3),koop(1:3)   	!vectors for transit factor determination
 real(kind=dp),dimension(1:2)::dr1,dr2  !d2  	!distances between the planet and the stars [au], relative vectors
 real(kind=dp),dimension(1:4)::spinrad 		!initial spin vector with angles in [rad] 
 real(kind=dp),dimension(1:2)::df  !tf 		!transit and distance insolation diminishing factors for each star
 
 !---------------------------------------------------------------------
 !dynamical evolution of the planet and the binary
 !---------------------------------------------------------------------
 call pconst(semia(1),semia(2),e1,man(1)*deg2rad,man(2)*deg2rad,mass(1),mass(2),mass(3),c,cs)
 call pdist(semia(1),semia(2),e1,mass(1),mass(2),mass(3),c,cs,t,.true.,0._dp,r0,r1,r2,e2,w1,w2)

 !get the 3D planetary spin vector in Cartesian coordinates from obliquity, angle of precession and rotation period
 spinrad(1:2)=spin(1:2)*deg2rad !obliquity and angle of precession [deg] -> [rad]
 spinrad(3)=spin(3) !rotation period
 spinrad(4)=spin(4)*deg2rad !initial hour angle of zero meridian [deg] -> [rad]
 
 !determine declination decl and hour angle phi wrt each star
 call planetorient(t,spinrad,r0,r1,r2,dr1,dr2,decl,phi)
 
 !---------------------------------------------------------------------
 !insolation
 !---------------------------------------------------------------------
 !calculate the integral over Planck curve bands
 !do j=1,2 !number of stars 
 ! do i=1,nbands
 !  call intplanckwn(teff(j),bands(j,i,1:2),bndflux(j,i))
 ! end do
 !end do 
 

!calculate transit factor (whether one star is hidden behind the other)
 koos(:,3)=0.d0
 koop(3)=0.d0
 koos(1,1:2)=r0
 koos(2,1:2)=r1
 koop(1:2)=r2
 call transitfactor(2,koos,koop,rs,tf)
 
  
!calculate the squares of the relative distances between binary stars and the planet
d2(1)=dot_product(dr1,dr1)
d2(2)=dot_product(dr2,dr2)

!diminish flux due to transit and distance 
! do i=1,2
!  df(i)=pi*(rs(i)*rsunau)**2/d2(i)              !distance diminishing factor
!  bndflux(i,:)=bndflux(i,:)*tf(i)*df(i)   !diminish flux by transit and distance factors
! end do
 
 return
end subroutine

!*********************************************************************************************************************************
!						P-TYPE DYNAMICS
!					     (Circumbinary planet)
!*********************************************************************************************************************************
 
subroutine e2av(a1,a2,e1,m0,m1,m2,eav)
!-------------------------------------------------------------------------------------------------------------------------------------
!generates <e^2> estimates for a planet (a2,m2) in a circumbinary orbit 
!GR pericenter advancement for the binary is included - no tides
!-------------------------------------------------------------------------------------------------------------------------------------
!Input:
!real::
!a1,a2               [au]   semimajor axis of the inner and outer orbit (for P-type motion a1:binary, a2:circumbinary planet)
!e1                  []     eccentricity of the inner orbit (P-type motion: eccentricity of the binary)
!m0,m1,m2            [Msun] masses of the primary, secondary, third body (P-type motion: binary star 1, binary star2, planet.)
!            
!-------------------------------------------------------------------------------------------------------------------------------------
!Output:
!real::
!eav                 [] averaged squared eccentricity <e^2> of outer orbit with respect to the binary center of mass
!                   
!
!written by Siegfried Eggl  20140329
!dependencies: none
!-------------------------------------------------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::a1,a2,e1,m0,m1,m2
real(kind=dp)::eav
!local
real(kind=dp)::e12
real(kind=dp)::c(1:11),cs(1:6)


 c(1)=k  !gauss k
 c(2)=clight  !c[au]
 c(3)=c(2)*c(2)        !c^2
 c(4)=pi  !pi

 !mass parameters
 c(5)=m0+m1+m2          !m_tot
 c(6)=m0*m1/(m0+m1)**(4._dp/3._dp)/c(5)**(2._dp/3._dp) !mstar

  
 !mean motions
 c(7)=c(1)*sqrt((m0+m1)/a1**3._dp) !n1
 c(8)=c(1)*sqrt(c(5)/a2**3._dp)    !n2

 !periods
 c(9)=2._dp*c(4)/c(7) !P1
 c(10)=2._dp*c(4)/c(8) !P2
 !period ratio
 c(11)=c(7)/c(8) !X=n1/n2=P2/P1
 
 !e1^2
 e12=e1*e1
  
 !calculate secular amplitudes
 cs(1)=3._dp/8._dp*c(1)*sqrt(c(5))*m0*m1*a1**2*(2._dp+3._dp*e12)/((m0+m1)**2*a2**3.5_dp)
 cs(2)=15._dp/64._dp*c(1)*sqrt(c(5))*m0*m1*(m0-m1)*a1**3*e1*(4._dp+3._dp*e12)/((m0+m1)**3*a2**4.5_dp)
 cs(3)=3._dp/4._dp*c(1)*m2*sqrt(a1**3*(1._dp-e12))/(sqrt(m0+m1)*a2**3) +3._dp*c(1)**3*(m0+m1)**1.5_dp/(c(3)*a1**2.5_dp*(1._dp-e12))
 

 !<e^2>
 eav= c(6)**2/(c(11)**(8._dp/3._dp))*(9._dp/8._dp+27._dp/8._dp*e12+887._dp/64._dp*e12*e12- &
         975._dp/64._dp/c(11)*e12*e12*sqrt(1._dp-e12)+1._dp/c(11)**2*&
         (225._dp/64._dp+6619._dp/64._dp*e12-26309._dp/512._dp*e12*e12-393._dp/64._dp*e12**3))+&
         2._dp*(cs(2)/(cs(1)-cs(3)))**2


 return
 end subroutine
 
!****************************************************************************************************** 
subroutine e2max(a1,a2,e1,m0,m1,m2,emaxsp,emaxsec,emax)
!-----------------------------------------------------------------------------------------------------------------------------------------
!generates maximum eccentricity estimates for a planet (a2,m2) in a circumbinary orbit 
!GR pericenter advancement for the binary is included - no tides
!
! written by Siegfried Eggl  20140329
!
! dependencies: none
!
!-----------------------------------------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!a1,a2               [au]   semimajor axis of the inner and outer orbit (for P-type motion a1:binary, a2:circumbinary planet)
!e1                  []     eccentricity of the inner orbit (P-type motion: eccentricity of the binary)
!m0,m1,m2            [Msun] masses of the primary, secondary, third body (P-type motion: binary star 1, binary star2, planet.)       
!-----------------------------------------------------------------------------------------------------------------------------------------
!Output:
!
!real::
!emaxsp              [] short period contribution to e2 max
!emaxsec             [] secular contribution to e2 max
!emax                [] maximum eccentricity of outer orbit with respect to the binary center of mass
!                   
!----------------------------------------------------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::a1,a2,e1,m0,m1,m2
real(kind=dp),intent(out)::emax

!local
real(kind=dp)::e12,emaxsec,emaxsp
real(kind=dp)::c(1:11),cs(1:6)

 e12=e1*e1 

 c(1)=k  !gauss k
 c(2)=clight  !c[au/D]
 c(3)=c(2)*c(2)        !c^2
 c(4)=pi  !pi

 !mass parameters
 c(5)=m0+m1+m2          !m_tot
 c(6)=m0*m1/(m0+m1)**(4._dp/3._dp)/c(5)**(2._dp/3._dp) !mstar

  
 !mean motions
 c(7)=c(1)*sqrt((m0+m1)/a1**3._dp) !n1
 c(8)=c(1)*sqrt(c(5)/a2**3._dp)    !n2

 !periods
 c(9)=2._dp*c(4)/c(7) !P1
 c(10)=2._dp*c(4)/c(8) !P2
 !period ratio
 c(11)=c(7)/c(8) !X=n1/n2=P2/P1
 
  !calculate secular amplitudes
 cs(1)=3._dp/8._dp*c(1)*sqrt(c(5))*m0*m1*a1**2*(2._dp+3._dp*e12)/((m0+m1)**2*a2**3.5_dp)
 cs(2)=15._dp/64._dp*c(1)*sqrt(c(5))*m0*m1*(m0-m1)*a1**3*e1*(4._dp+3._dp*e12)/((m0+m1)**3*a2**4.5_dp)
 cs(3)=3._dp/4._dp*c(1)*m2*sqrt(a1**3*(1._dp-e12))/(sqrt(m0+m1)*a2**3) +3._dp*c(1)**3*(m0+m1)**1.5_dp/(c(3)*a1**2.5_dp*(1._dp-e12))
 

 !short period contribution to emax
 emaxsp= c(6)*c(11)**(-4._dp/3._dp)*(1.5_dp+17._dp/2._dp*e12+1._dp/c(11)*(3._dp+19._dp*e1+ &
         21._dp/8._dp*e12-1.5_dp*e12*e1))
 
 !secular contribution to emax
 emaxsec=2._dp*cs(2)/(cs(1)-cs(3))

 !emax
 emax=emaxsp+emaxsec

 return
 end subroutine
!***********************************************************************************************

subroutine hw99p(m0,m1,a1,e1,ac)
!------------------------------------------------------------------------------------------------------------
!stability check for P-type restricted 3 body systems following Holman & Wiegert 1999
!
!written by Siegfried Eggl  20140304
!
!dependencies: none
!------------------------------------------------------------------------------------------------------------
!Input: 
!
!real::
!m0,m1                            ...[Msun] masses of the primary and secondary in the binary star system
!a1                               ...[au]   semimajor axes of the inner (binary) 
!e1                               ...[]     eccentricity of the inner orbit
!------------------------------------------------------------------------------------------------------------
!Output
!
!real::
!ac                               ...[au] critical semimajor axis (from the binary's barycenter)
!
!-------------------------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::m0,m1,a1,e1
real(kind=dp),intent(out)::ac
!local
real(kind=dp)::e12,mu,mu2

!dummy variables
e12=e1*e1
mu=min(m0,m1)/(m0+m1)
mu2=mu*mu

!calculate critical semimajor axis
ac=a1*(1.6_dp+5.1_dp*e1-2.22_dp*e12+&
   4.12_dp*mu-4.27_dp*e1*mu-5.09_dp*mu2+4.61_dp*e12*mu2)

return
end subroutine

!**********************************************************************************************************
subroutine pconst(a1,a2,e1,ma10,ma20,m0,m1,m2,c,cs)
!--------------------------------------------------------------------------------------------------------------------------------
!generates constants necessary for describing the outer (planetary) eccentricity evolution
!
! written by Siegfried Eggl  20140304
!
! dependencies: esp
!--------------------------------------------------------------------------------------------------------------------------------
!Input
!
!real::
!a1,a2               [au]   semimajor axes of the inner and outer orbit (for P-type motion a1:binary, a2:circumbinary planet)
!e1                  []     eccentricity of the inner orbit (P-type motion: eccentricity of the binary)
!ma10,ma20           [rad]  initial phases (mean longitude of the inner and mean longitude of the outer orbit)
!m0,m1,m2            [Msun] masses of the primary, secondary, third body (P-type motion: binary star 1, binary star2, planet.)
!--------------------------------------------------------------------------------------------------------------------------------
!Output 
!
!real::
!c(1:14)::           vector containing constants
!c(1)                [sqrt(au^3/Msun/D^2)] Gaussian gravitational constants
!c(2)                [au/D] light speed 
!c(3)                [(au/D)^2] light speed squared
!c(4)                [] PI
!c(5)                [Msun] total mass of the system
!c(6)                [] mass parameter mstar 
!c(7:8)              [rad/D] mean motions of inner and outer orbit
!c(9:10)             [D] periods of inner and outer orbit
!c(11)               [] period ratio
!c(12)               [rad] initial mean longitude of inner orbit: ma10
!c(13)		     [rad] initial mean longitude of the outer orbit ma20
!c(14)		     [rad] initial eccentric anomaly of the inner orbit ea10

!cs(1:6)::           vector containing constants for secular amplitudes and frequencies
!cs(1:3)             K1-K3
!cs(4:6)             C1-C3
!
!---------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::a1,a2,e1,m0,m1,m2,ma10,ma20
real(kind=dp),intent(out)::c(1:14),cs(1:6)
!local
real(kind=dp)::e12,e2spx,e2spy,ea0


 c(1)=k  !gauss k
 c(2)=clight  !c[au]
 c(3)=c(2)*c(2)        !c^2
 c(4)=pi  !pi

 
 !mass parameters
 c(5)=m0+m1+m2          !m_tot
 c(6)=m0*m1/(m0+m1)**(4._dp/3._dp)/c(5)**(2._dp/3._dp) !mstar

  
 !mean motions
 c(7)=c(1)*sqrt((m0+m1)/a1**3._dp) !n1
 c(8)=c(1)*sqrt(c(5)/a2**3._dp)    !n2

 !periods
 c(9)=pix2/c(7) !P1
 c(10)=pix2/c(8) !P2
 !period ratio
 c(11)=c(7)/c(8) !X=n1/n2=P2/P1
 
 !initial eccentric anomaly of inner orbit -> initial mean longitude of inner orbit
 c(12)=ma10
 !initial mean longitude of outer orbit
 c(13)=ma20
 
 !initial eccentric anomaly of the inner orbit (ea0)
 !eccentric anomaly
  call nr(ma10,e1,ea0)
 c(14)=ea0
 
 !calculate secular frequencies and amplitudes
 e12=e1*e1
 
 cs(1)=3._dp/8._dp*c(1)*sqrt(c(5))*m0*m1*a1**2*(2._dp+3._dp*e12)/((m0+m1)**2*a2**3.5_dp)
 cs(2)=15._dp/64._dp*c(1)*sqrt(c(5))*m0*m1*(m0-m1)*a1**3*e1*(4._dp+3._dp*e12)/((m0+m1)**3*a2**4.5_dp)
 cs(3)=3._dp/4._dp*c(1)*m2*sqrt(a1**3*(1._dp-e12))/(sqrt(m0+m1)*a2**3) +3._dp*c(1)**3*(m0+m1)**1.5_dp/(c(3)*a1**2.5_dp*(1._dp-e12))
 
 
 !call short period terms at t=0 to get constants
 call pesp(e1,0._dp,c,e2spx,e2spy)
 
 cs(6)=cs(2)/(cs(1)-cs(3))
 cs(4)=-e2spx-cs(6)
 cs(5)=e2spy 

 return
end subroutine

!******************************************************************************************************

subroutine pdist(a1,a2,e1,m0,m1,m2,c,cs,t,e2evol,e2init,r0,r1,r2,e2,w1,w2)
!-----------------------------------------------------------------------------------------------------------------------------------------
!generates barycentric positions of a binary star and a circumbinary planet using analytical eccentricity vector estimates for the planet 
!and the two body problem + GR pericenter advancement for the binary
!
! written by Siegfried Eggl  20140304
! 
! dependencies: oe2xy,peouter
!------------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!a1,a2               [au]   semimajor axis of the inner and outer orbit (for P-type motion a1:binary, a2:circumbinary planet)
!e1                  []     eccentricity of the inner orbit (P-type motion: eccentricity of the binary)
!m0,m1,m2            [Msun] masses of the primary, secondary, third body (P-type motion: binary star 1, binary star2, planet.)

!c(1:14)::             vector containing constants
!c(1)                [sqrt(au^3/Msun/D^2)] Gaussian gravitational constants
!c(2)                [au/D] light speed 
!c(3)                [(au/D)^2] light speed squared
!c(4)                [] PI
!c(5)                [Msun] total mass of the system
!c(6)                [] mass parameter mstar 
!c(7:8)              [rad/D] mean motions of inner and outer orbit
!c(9:10)             [D] periods of inner and outer orbit
!c(11)               [] period ratio
!c(12)               [rad] initial mean longitude of inner orbit: ma10
!c(13)		     [rad] initial mean longitude of the outer orbit ma20
!c(14)		     [rad] initial eccentric anomaly of the inner orbit ea10

!cs(1:6)::             vector containing constants for secular amplitudes and frequencies
!cs(1:3)             K1-K3
!cs(4:6)             C1-C3

!t                   [D] current time (t-t0)

!e2init              initial eccentricity (only important if e2 evolution is OFF, otherwise e2ini=0!!!)

!logical::         
!e2evol               eccentricity evolution of the planet on/off    
!
!------------------------------------------------------------------------------------------------------------
!Output:
!
!real::
!r0(1:2),r1(1:2)     [au] barycentric vector of star 1 (m0, r0) and star 2 (m1, r1) with respect to the stellar 2body problem
!                    i.e. only the barycenter of the two stars are accounted for                     
!r2(1:2)             [au] vector of the planet with respect to the barycenter of the binary; if m1=0 then it is the heliocentric vector between m0 and m2
!                          
!e2                  [] evolving outer eccentricity
!w1,w2               [rad] argument of pericenter of the inner and outer orbit
!-------------------------------------------------------------------------------------------------------------------------------------
real(kind=dp),intent(in)::a1,a2,e1,m0,m1,m2,t
real(kind=dp),intent(in)::c(1:14),cs(1:6)
real(kind=dp),intent(in)::e2init
logical,intent(in)::e2evol

!local
real(kind=dp)::e12,e2x,e2y,e2,w1,w2
real(kind=dp),dimension(1:2)::rb,r0,r1,r2


!---------------- binary star xy coordinates (unperturbed two body problem) ---------------------------
w1=mod(cs(3)*t,pix2) 			!GR precession of pericenter of the binary
!w1=0._dp
call oe2xy2(a1,e1,c(7),c(12),w1,t,rb)  	!transform orbital elements to heliocentric xy coordinates

r0(:)=-m1/(m0+m1)*rb(:)			!transform to barycentric vectors
r1(:)=m0/(m0+m1)*rb(:)

!---------------- planet xy coordinates -----------------------------------
if(e2evol) then
 call peouter(e1,c,cs,t, e2x,e2y)  	!get Laplace Runge Lenz vector for outer orbit LRL=(e2x,e2y)^T
 w2=atan2(e2y,e2x)
 e2=sqrt(e2x*e2x+e2y*e2y)
 if (w2.lt.0) then
  w2=w2+pix2				 !put w2 from [-PI,PI] to [0,2 PI]
 end if
else
 e2=e2init
 w2=0._dp
end if
 
call oe2xy2(a2,e2,c(8),c(13),w2,t,r2)	!get (binary) barycentric vector of the planet 

return
end subroutine 

!********************************************************************************************************
subroutine peouter(e1,c,cs,t, e2x,e2y)
!-------------------------------------------------------------------------------------------------------
!calculates outer eccentricity vector by linear superposition of short period and secular terms
!
!written by Siegfried Eggl  20140304
!
!dependencies: pesp, pesec
!
!-------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!e1                  []      eccentricity of the inner orbit (P-type motion: eccentricity of the binary)
!c(1:14)             [mixed] vector containing constants
!cs(1:6)             [mixed] containing constants for secular amplitudes and frequencies
!t                   [D]     time
!-------------------------------------------------------------------------------------------------------
!Output:
!
!real::
!e2x                 [] x component of the outer orbit's Laplace Runge Lenz vector
!e2y                 [] y component of the outer orbit's Laplace Runge Lenz vector
!
!-------------------------------------------------------------------------------------------------------
implicit none
!input
real(kind=dp),intent(in)::e1,c(1:14),cs(1:6),t
!output
real(kind=dp),intent(out)::e2x,e2y
!local
real(kind=dp)::e2spx,e2spy,e2secx,e2secy

call pesp(e1,t,c,e2spx,e2spy)
call pesec(cs,t,e2secx,e2secy)

e2x=e2spx+e2secx
e2y=e2spy+e2secy

return
end subroutine

!**************************************************************************************************

subroutine pesec(cs,t,e2secx,e2secy)
!-------------------------------------------------------------------------------------------------
!calculate secular outer eccentricity of a p-type (circumbinary) planet system
!
! written by Siegfried Eggl  20140304
!
! dependencies: none
!----------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!cs(1:6)::           vector containing constants for secular amplitudes and frequencies
!cs(1:3)             K1-K3
!cs(4:6)             C1-C3
!t                   [D]    time
!----------------------------------------------------------------------------------------------------------
!Output:
!
!real::
!e2secx              [] secular x coordinate of outer Laplace Runge Lenz vector
!e2secy              [] secular y coordinate of outer Laplace Runge Lenz vectors
!
!-------------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::cs(1:6),t
real(kind=dp),intent(out)::e2secx,e2secy
!local
real(kind=dp)::cs1t,cs3t,scs1t,ccs1t
 cs1t=cs(1)*t
 cs3t=cs(3)*t
 scs1t=sin(cs1t)
 ccs1t=cos(cs1t)

e2secx=cs(4)*ccs1t+cs(5)*scs1t+cs(6)*cos(cs3t)
e2secy=cs(4)*scs1t-cs(5)*ccs1t+cs(6)*sin(cs3t)

return
end subroutine

!**************************************************************************************************
subroutine pesp(e1,t,c,e2spx,e2spy)
!----------------------------------------------------------------------------------------------------------
!short period outer eccentricity with averaged binary motion (averaged over true anomaly of binary)
!
!
!written by Siegfried Eggl  20140304
!dependencies: none
!----------------------------------------------------------------------------------------------------------
!Input
!
!real::

!e1                  [] binary eccentricity
!t                   [D] time

!c(1:14)::           vector containing constants
!c(1)                [sqrt(au^3/Msun/D^2)] Gaussian gravitational constants
!c(2)                [au/D] light speed 
!c(3)                [(au/D)^2] light speed squared
!c(4)                [] PI
!c(5)                [Msun] total mass of the system
!c(6)                [] mass parameter mstar 
!c(7:8)              [rad/D] mean motions of inner and outer orbit
!c(9:10)             [D] periods of inner and outer orbit
!c(11)               [] period ratio
!c(12)               [rad] initial mean longitude of inner orbit: ma10
!c(13)		     [rad] initial mean longitude of the outer orbit ma20
!c(14)		     [rad] initial eccentric anomaly of the inner orbit ea10
!
!----------------------------------------------------------------------------------------------------------
!Output
!
!real::
!e2spx               [] short period x component of the outer orbits Laplace Runge Lenz vector
!e2spy               [] short period x component of the outer orbits Laplace Runge Lenz vector
!
!----------------------------------------------------------------------------------------------------------
!Local:
!
!real::
!e2aspx               [] averaged short period x component of the outer orbits Laplace Runge Lenz vector
!e2aspy               [] averaged short period y component of the outer orbits Laplace Runge Lenz vector
!
!e2uspx              [] ultra short period terms of x component of the outer orbits Laplace Runge Lenz vector
!e2uspy		     [] ultra short period terms of y component of the outer orbits Laplace Runge Lenz vector
!
!-------------------------------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::e1,c(1:14),t
real(kind=dp),intent(out)::e2spx,e2spy

!local
real(kind=dp)::e12,e13,n1,n2,ea,ea1,l1,l2,b,bm,bp
real(kind=dp)::e2aspx,e2aspy,e2uspx,e2uspy
real(kind=dp),dimension(1:24)::dx,qx,wx,dy,qy,wy
real(kind=dp)::ga,gb
!const

e12=e1*e1
e13=e12*e1
b=sqrt(1._dp-e12)
bm=1._dp-b
bp=1._dp+b
n1=c(7)
n2=c(8)
l2=mod(n2*t+c(13),pix2)
l1=mod(n1*t,pix2)

!global amplitude of P_xy short period terms
ga=c(6)/c(11)**(7._dp/3._dp)
!global amplitude of the other short period terms
gb=c(6)/16._dp/c(11)**(4._dp/3._dp)

!calculate eccentric anomaly for the inner orbit using Newton's method
 call nr(l1,e1,ea) 
  ea1=ea+c(14)

!################ x component ############################  
dx(1)=21._dp/32._dp
dx(2)=3._dp/32._dp
dx(3)=dx(2)
dx(4)=dx(1)
dx(5)=21._dp/96._dp
dx(6)=3._dp/96._dp
dx(7)=dx(2)
dx(8)=dx(6)
dx(9)=105._dp/32._dp
dx(10)=dx(5)
dx(11)=dx(2)
dx(12)=dx(9)
dx(13)=dx(2)
dx(14)=dx(1)
dx(15)=dx(2)
dx(16)=dx(1)
dx(17)=1._dp/64._dp
dx(18)=7._dp/64._dp
dx(19)=3._dp/64._dp
dx(20)=dx(18)
dx(21)=21._dp/64._dp
dx(22)=dx(19)
dx(23)=dx(17)
dx(24)=dx(21)

qx(1)=bm
qx(2)=bm
qx(3)=-bp
qx(4)=-bp
!------------order e1-----------------------
qx(5)=-bm*e1
qx(6)=-bm*e1
qx(7)=(13._dp+5._dp*b)*e1
qx(8)=bp*e1
qx(9)=-bm*e1
qx(10)=bp*e1
qx(11)=-(13._dp-5._dp*b)*e1
qx(12)=bp*e1
!----------- order e2----------------------
qx(13)=(3.5_dp-b)*e12
qx(14)=(0.5_dp-b)*e12
qx(15)=-(3.5_dp+b)*e12
qx(16)=-(0.5_dp+b)*e12
!----------- order e3----------------------
qx(17)=-e13
qx(18)=-e13
qx(19)=e13
qx(20)=e13
qx(21)=e13
qx(22)=-e13
qx(23)=e13
qx(24)=-e13


wx(1)=cos(2._dp*ea1+3._dp*l2)
wx(2)=cos(2._dp*ea1+l2)
wx(3)=cos(2._dp*ea1-l2)
wx(4)=cos(2._dp*ea1-3._dp*l2)
!------------order e1-----------------------
wx(5)=cos(3._dp*ea1+3._dp*l2)
wx(6)=cos(3._dp*ea1+l2)
wx(7)=cos(ea1-l2)
wx(8)=cos(3._dp*ea1-l2)
wx(9)=cos(ea1+3._dp*l2)
wx(10)=cos(3._dp*ea1-3._dp*l2)
wx(11)=cos(ea1+l2)
wx(12)=cos(ea1-3._dp*l2)
!----------- order e2----------------------
wx(13)=wx(2)
wx(14)=wx(1)
wx(15)=wx(3)
wx(16)=wx(4)
!----------- order e3----------------------
wx(17)=wx(6)
wx(18)=wx(10)
wx(19)=wx(7)
wx(20)=wx(5)
wx(21)=wx(9)
wx(22)=wx(11)
wx(23)=wx(8)
wx(24)=wx(12)


e2aspx=gb*(12._dp*cos(l2)+e12*(33._dp* &
        cos(l2)+35._dp*cos(3._dp*(l2))))
          
   
!P_x(t)/X*mstar/X^(4/3)
e2uspx=ga*dot_product(dx(:)*qx(:),wx(:))

!outer exccentricity vector short periodic x-component
e2spx=e2aspx+e2uspx

!################ y component ############################

dy(1)=21._dp/32._dp
dy(2)=3._dp/32._dp
dy(3)=dy(2)
dy(4)=dy(1)
!------------order e1-----------------------
dy(5)=7._dp/32._dp
dy(6)=1._dp/32._dp
dy(7)=3._dp/32._dp
dy(8)=dy(6)
dy(9)=105._dp/32._dp
dy(10)=dy(5)
dy(11)=dy(2)
dy(12)=dy(9)
!----------- order e2----------------------
dy(13)=dy(2)
dy(14)=dy(1)
dy(15)=dy(2)
dy(16)=dy(1)
!----------- order e3----------------------
dy(17)=3._dp/64._dp
dy(18)=7._dp/64._dp
dy(19)=9._dp/64._dp
dy(20)=dy(18)
dy(21)=21._dp/64._dp
dy(22)=dy(19)
dy(23)=dy(17)
dy(24)=dy(21)

qy(1)=bm
qy(2)=-bm
qy(3)=-bp
qy(4)=bp
!------------order e1-----------------------
qy(5)=-bm*e1
qy(6)=bm*e1
qy(7)=-(3._dp-5._dp*b)*e1
qy(8)=bp*e1
qy(9)=-bm*e1
qy(10)=-bp*e1
qy(11)=-(3._dp+5._dp*b)*e1
qy(12)=-bp*e1
!----------- order e2----------------------
qy(13)=(2.5_dp+b)*e12
qy(14)=(0.5_dp-b)*e12
qy(15)=(2.5_dp-b)*e12
qy(16)=(0.5_dp+b)*e12
!----------- order e3----------------------
qy(17)=-e13
qy(18)=e13
qy(19)=-e13
qy(20)=e13
qy(21)=e13
qy(22)=-e13
qy(23)=-e13
qy(24)=e13


wy(1)=sin(2._dp*ea1+3._dp*l2)
wy(2)=sin(2._dp*ea1+l2)
wy(3)=sin(2._dp*ea1-l2)
wy(4)=sin(2._dp*ea1-3._dp*l2)
!------------order e1-----------------------
wy(5)=sin(3._dp*ea1+3._dp*l2)
wy(6)=sin(3._dp*ea1+l2)
wy(7)=sin(ea1-l2)
wy(8)=sin(3._dp*ea1-l2)
wy(9)=sin(ea1+3._dp*l2)
wy(10)=sin(3._dp*ea1-3._dp*l2)
wy(11)=sin(ea1+l2)
wy(12)=sin(ea1-3._dp*l2)
!----------- order e2----------------------
wy(13)=wy(2)
wy(14)=wy(1)
wy(15)=wy(3)
wy(16)=wy(4)
!----------- order e3----------------------
wy(17)=wy(6)
wy(18)=wy(10)
wy(19)=wy(7)
wy(20)=wy(5)
wy(21)=wy(9)
wy(22)=wy(11)
wy(23)=wy(8)
wy(24)=wy(12)


e2aspy=gb*(12._dp*sin(l2)+e12*(3._dp* &
         sin(l2)+35._dp*sin(3._dp*(l2))))
       

!P_y(t)/X*mstar/X^(4/3)
e2uspy=ga*dot_product(dy(:)*qy(:),wy(:))
      
!outer eccentricity vector short periodic y-component
e2spy=e2aspy+e2uspy
 
return
end subroutine

!*************************************************************************************************

subroutine pespav(e1,l0,t,c,e2aspx,e2aspy)
!----------------------------------------------------------------------------------------------------------
! short period outer eccentricity with averaged binary motion (averaged over true anomaly of binary)
!
! written by Siegfried Eggl  20140304
!
! dependencies: none
!----------------------------------------------------------------------------------------------------------
!Input
!
!real::
!
!e1                  [] binary eccentricity
!l0                  [rad] initial outer mean longitude
!t                   [D] time
!
!c(1:14)::           vector containing constants
!c(1)                [sqrt(au^3/Msun/D^2)] Gaussian gravitational constants
!c(2)                [au/D] light speed 
!c(3)                [(au/D)^2] light speed squared
!c(4)                [] PI
!c(5)                [Msun] total mass of the system
!c(6)                [] mass parameter mstar 
!c(7:8)              [rad/D] mean motions of inner and outer orbit
!c(9:10)             [D] periods of inner and outer orbit
!c(11)               [] period ratio
!c(12)               [rad] initial mean longitude of inner orbit: ma10
!c(13)		     [rad] initial mean longitude of the outer orbit ma20
!c(14)		     [rad] initial eccentric anomaly of the inner orbit ea10
!----------------------------------------------------------------------------------------------------------
!Output
!
!real::
!e2aspx               [] short period x component of the outer orbits Laplace Runge Lenz vector
!e2aspy               [] short period x component of the outer orbits Laplace Runge Lenz vector
!
!-------------------------------------------------------------------------------------------------
implicit none
!input
real(kind=dp),intent(in)::e1,c(1:14),l0,t
!output
real(kind=dp),intent(out)::e2aspx,e2aspy
!local
real(kind=dp)::e12,n1,n2,l2
real(kind=dp)::cl2,sl2,a


e12=e1*e1
n1=c(7)
n2=c(8)
l2=n2*t+l0

 cl2=cos(l2)
 sl2=sin(l2)

 !amplitude
 a=c(6)/16._dp/c(11)**(4._dp/3._dp)
 !x component 
 e2aspx=a*(12._dp*cl2+e12*(33._dp*cl2+35._dp*cos(3._dp*l2)))
 !y component       
 e2aspy=a*(12._dp*sl2+e12*(3._dp*sl2+35._dp*sin(3._dp*l2)))
       
 
return
end subroutine

!*********************************************************************************************************************************
subroutine planetorient(t,s,r0,r1,r2,rd1,rd2,delta,phi)
!------------------------------------------------------------------------------------------
! determine declination and hour angle of insolation between the planet and the two stars
!
! written by Siegfried Eggl  20140304
!
! dependencies: none
!-------------------------------------------------------------------------------------------
!
! Input
! real::
! t			...[D] current time in Gaussian days
! s(1:4)      		... [mixed] 	(1) obl [rad]	obliquity (ref: z-axis = normal to the orbital plane),
! 					(2) pre [rad]	angle of precession (ref: x-axis = along initial pericenter of binary) 
! 					(3) pp  [D]     rotation period of the planet
!					(4) eta0 [rad]  initial hour angle of zero-meridian of planet wrt x-axis 
! r0,r1,r2 (1:2)        ...[au] position vectors of the stars (r0,r1) and the planet (r2)
!-------------------------------------------------------------------------------------------
! 
! Output
! real::
! rd1,rd2 (1:2)		...[au] relative distances between the stars and the planet
! delta (1:2)           ...[rad] declination of the planet wrt. the two stars
! phi (1:2)		...[rad] hour angle (angle between zero-meridian and stellar position)
!-------------------------------------------------------------------------------------------

 implicit none
 !input
 real(kind=dp),intent(in)::t,s(1:4)
 real(kind=dp),dimension(1:2),intent(in)::r0,r1,r2

 !output
 real(kind=dp),dimension(1:2),intent(out)::rd1,rd2,delta,phi
 
 !local
 integer::i
 real(kind=dp)::eta ! eta= rotation angle of the planet (zero-meridian to x-axis), sf= spin frequency of the planet [1/D] multiplied by 2PI  

 real(kind=dp),dimension(1:3)::rs1p,rs2p !polar coordinates of relative distance vectors between planet and star
 real(kind=dp),dimension(1:3)::rs1r,rs2r !relative distance vectors between planet and star rotated into the zero meridian planet frame
 real(kind=dp),dimension(1:3,1:3)::rot,irot
 real(kind=dp)::ca,sa,co,so,cs,ss
 real(kind=dp),dimension(1:3)::spinout
 
! To get the declination and hour angles of the stars wrt the planet's equator and zero-meridian, we rotate the relative position vectors into 
! a planetocentric frame, and convert them to polar coordinates.
 
 !obliquity
 co=cos(s(1))
 so=sin(s(1))
 
 !angle of precession
 ca=cos(s(2))
 sa=sin(s(2))
 
 !correct rotation direction of the planet if obliquity is larger than 90° and smaller than 270°
 if(abs(s(1)).gt.pib2.and.abs(s(1)).le.pi3b2) then
  eta=mod(-pix2/s(3)*t+s(4),pix2)
 else
  !planet proper rotation (progression of zero-meridian) 
  eta=mod(pix2/s(3)*t+s(4),pix2)
 end if
 
 cs=cos(eta)
 ss=sin(eta)
 
! Here, irot is the inverse rotation matirx of the transformations that convert barycentric coordinates 
! to planetocentric coordinates using the equator and an arbitrary zero-meridian as a reference.
! The original rotation matrix is constructed as follos: rot = R_z(pre) R_y(obl) R_z(eta). 
!  rot(1,1:3)=(/ca*co*cs - sa*ss, -cs* sa - ca* co* ss, ca* so/)
!  rot(2,1:3)=(/co*cs*sa + ca*ss, ca*cs - co *sa* ss, sa* so/)
!  rot(3,1:3)=(/-cs*so, so*ss, co/)
! where irot is its inverse

 irot(1,1:3)=(/ca*co*cs - sa*ss, co*cs*sa + ca*ss, -cs*so/)
 irot(2,1:3)=(/-cs*sa - ca*co*ss, ca*cs - co*sa*ss, so*ss/)
 irot(3,1:3)=(/ca*so, sa*so, co/)


!relative distance vectors between planet and stars 
 rd1=r0-r2
 rd2=r1-r2

!relative distance vecotrs in the rotated planetocentric frame 
 rs1r=matmul(irot(1:3,1:2),rd1)
 rs2r=matmul(irot(1:3,1:2),rd2)

!conversion to polar coordinates gives declination (delta) and hour angle (phi) 
 call topolar(rs1r,rs1p)
 call topolar(rs2r,rs2p)
 
 delta(1)=rs1p(2)
 delta(2)=rs2p(2)
 
 phi(1)=rs1p(3)
 phi(2)=rs2p(3)

return
end subroutine
 

!************************************************************************************************
!
!					General Purpose
!
!************************************************************************************************

subroutine oe2xy2(a,e,n,m0,w,t,x)
!---------------------------------------------------------------------------
!converts orbital elements to xy coordinates in the planar two body problem 
!
! written by Siegfried Eggl  20140304
!
! dependencies: nr
!----------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!a                  [au] semimajor axis
!e                  [] eccentricity
!n                  [rad/D] mean motion
!m0                 [rad] initial mean longitude
!w                  [rad] argument of pericenter
!t                  [D] time-steps
!----------------------------------------------------------------------------------------------------------
!Output
!
!real::
!x(1:2)             [au] x,y coordinates wrt focus
!

!--------------------------------------------------------------------------
implicit none
real(kind=dp),intent(in)::a,e,n,m0,w,t
real(kind=dp),intent(out)::x(1:2)
!local
real(kind=dp)::rotmat(1:2,1:2),ea,mm,ml,rcosf,rsinf
real(kind=dp),dimension(1:3)::x2,v

!calculate mean longitude as a function of time (Ml=n*t, n=2pi/T)
!since phi is the angle between r1 and r2 rather than r0 and r2 the supplementary angle to phi has to be added to the mean longitude ml
ml=mod(n*t+m0,pix2)

!find mean longitude as a function of time (M=n*t+l0-w, n=2pi/T)
mm=mod(ml-w+pix2,pix2)

!eccentric anomaly
call nr(mm,e,ea)
!true anomaly -> x,y

rcosf=a*(Cos(ea)-e)			!x
rsinf=a*sqrt(1._dp-e*e)*Sin(ea)		!y

!rotate orbit to correct pericenter
rotmat(1,1:2)=(/cos(w),-sin(w)/)
rotmat(2,1:2)=(/sin(w),cos(w)/)

!final xy coordinates
x=matmul(rotmat,(/rcosf,rsinf/))

return
end subroutine

!************************************************************************************************************

subroutine nr(m,ecc,ea)
!-------------------------------------------------------------------------------------------------
!solves Kepler’s equation by applying Newton Raphson Method
!with an accuracy limit of "deat"
!
!----------------------------------------------------------------------------------------------------------
!Input:
!
!m[real].................mean longitude (radians!)
!ecc[real]................Eccentricity (numerical Eccentricity <1!)
!
!----------------------------------------------------------------------------------------------------------
!Output:
!
!ea[real]................Eccentric Anomaly (radians!)
!
!
!written by Siegfried Eggl  20061222
!modified                   20111026
!dependencies: none
!-------------------------------------------------------------------------------------------------
implicit none
        integer(kind=isp)::i,iter
        real(kind=dp)::m,ecc,ea,ea0,dea,deat

ea=0._dp
ea0=1.3421_dp
deat=1.E-8
dea=ABS(ea-ea0)
iter=50


do i=1,iter

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

end do

if (ea.le.1.E-14) then
   ea=0._dp
end if
!if precision is not achieved try with different initial condition
if(dea>deat) then

ea=0._dp
ea0=0.3421_dp
dea=ABS(ea-ea0)
 do i=1,iter

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if


if(dea>deat) then

ea=0._dp
ea0=3.1_dp
dea=ABS(ea-ea0)
 do i=1,iter

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if

if(dea>deat) then
   write(unit=*,fmt=*)'convergence-error in subroutine nr'
   write(unit=*,fmt=*)'target precision:',deat
   write(unit=*,fmt=*)'achieved precision:',dea
end if
return

end subroutine nr
!******************************************************************************
      subroutine topolar(x,p)
!------------------------------------------------------------------------------
!     Changes 3D Cartesian to polar coordinates.
!     Attention: Latitudes {p(2)} are counted from the equator upwards!!!
!     
!     written by Siegfried Eggl 20140611
!
!     dependencies: none
!------------------------------------------------------------------------------
!     Input:
!     real::
!     x(1:3)			...[whatever distance unit = WDU] Cartesian coordinates
!------------------------------------------------------------------------------
!     Output:
!     real::
!     p(1:3)                    ...[mixed] polar coordinate vector,
!     			p(1)	...[WDU] polar distance r
!			p(2)    ...[rad] latitude (declination)
!                       p(3)    ...[rad] longitude
!------------------------------------------------------------------------------
          implicit none
          real(kind=dp),intent(in)::x(1:3) 
          real(kind=dp),intent(out)::p(1:3)
        
          p(1)=sqrt(Dot_Product(x,x)) !r
          p(2)=asin(x(3)/p(1))	      !latitude (declination)
          p(3)=atan2(x(2),x(1))       !longitude
          
          if(p(3).lt.0._dp) then
             p(3)=p(3)+pix2
          end if
     
          return
        end subroutine topolar

!*************************************************************        
        subroutine tocartesian(p,x)
!------------------------------------------------------------------------------
!     Changes 3D polar to Cartesian coordinates.
!     Attention: Latitudes {p(2)} are counted from the equator upwards!!!
!     
!     written by Siegfried Eggl 20140611
!
!     dependencies: none
!------------------------------------------------------------------------------
!     Input:
!     real::
!     p(1:3)                    ...[mixed] polar coordinate vector,
!     			p(1)	...[WDU] polar distance r
!			p(2)    ...[rad] latitude (declination)
!                       p(3)    ...[rad] longitude
!------------------------------------------------------------------------------
!     Output:
!     real::
!     x(1:3)			...[whatever distance unit = WDU] Cartesian coordinates
!------------------------------------------------------------------------------
        implicit none
        !input
        real(kind=dp),intent(in)::p(1:3)
        !output
        real(kind=dp),intent(out)::x(1:3)
       
        x(1) = p(1)*cos(p(2))*cos(p(3))
        x(2) = p(1)*cos(p(2))*sin(p(3))
        x(3) = p(1)*sin(p(2))

        return
        end subroutine
!******************************************************************************


!*********************************************************************************************************************************
!						P-TYPE INSOLATION
!					     (Circumbinary planet)
!*********************************************************************************************************************************

subroutine bxaxb(a,b,c)
!-----------------------------------------------------------------------------------------------
! calculates vector product C = B x (A x B) in 3-space for subroutine transitfactor
! 
!
! written by Siegfried Eggl 20120502
!
! dependencies: none
!
! Input
!
! real::
! a(1:3)                                ... vector a
! b(1:3)				... vector b
!
! Output
!
! real::
! c(1:3)				... vector c
!--------------------------------------------------------------------------------------------------

       implicit none
       !input
       real(kind=dp),intent(in):: a(1:3),b(1:3)
       !output
       real(kind=dp),intent(out)::c(1:3)
       !local
       integer:: i
       real(kind=dp)::d(1:3),projm(1:3,1:3)
       
          
       d(:)=b(:)/Sqrt(Dot_Product(b(:),b(:)))

       projm(1,1:3)=(/1.d0-d(1)*d(1),-d(1)*d(2),-d(1)*d(3)/)      
       projm(2,1:3)=(/-d(2)*d(1),1.d0-d(2)*d(2),-d(2)*d(3)/)
       projm(3,1:3)=(/-d(3)*d(1),-d(3)*d(2),1.d0-d(3)*d(3)/)
 
       c(:)=Matmul(projm(:,:),a(:)) 
 
       return

end subroutine  
!*************************************************************************************************************
 subroutine getspin(spin,spout)
!--------------------------------------------------------------------------------------------------------------
! produces 3D spin vector from obliquity: spin(1), angle of precession: spin(2) and rotation period: spin(3)
!
! written by Siegfried Eggl 20140610
!
! dependencies: none
!
! Input:
! real::
! spin(1:3)          ...initial conditions for the spin vector
!	             ...spin(1) [deg] = obliquity [deg] (wrt. z-axis)
!		     ...spin(2) [deg] = angle of precession [deg] (wrt. x-axis)
!		     ...spin(3) [D]   = period of rotation [D]
!
! Output
! real::
! spout(1:3)         ...3 component spin vector with length = spin frequency = 2xPi/period of rotation
!
!------------------------------------------------------------------------------------------------------------

 implicit none
 !input
 real(kind=dp),dimension(1:3),intent(in)::spin
 !output
 real(kind=dp),dimension(1:3),intent(out)::spout
 !local
 !real(kind=dp),dimension(1:3,1:3)::rot
 real(kind=dp)::cy,cz,sy,sz,s
 
 cy=cos(spin(1)) !cos(obl)
 sy=sin(spin(1)) !sin(obl)
 cz=cos(spin(2)) !cos(pre)
 sz=sin(spin(2)) !sin(pre)
 
 s=pix2/spin(3) !initially the planet's rotation axis points exactly in z direction with frequency 2pi/rotation period
 
 spout=(/s*cz*sy,s*sy*sz,s*cy/) !the final spin vector results from two rotations around the y and z axis spout=[Rot_z(pre).Rot_y(obl)].(0,0,s)^T
 
 return
 end subroutine
 !***************************************************************************************************************
subroutine intplanckwn(t,sigma,flux)
!--------------------------------------------------------------
! integrate Planck curve (wave number) 
! following  Widger & Woodall 1976, Bulletin of the American Meteorological Society, Vol. 57, No. 10, pp. 1217-1219
! 
! written by Siegfried Eggl 20140609
!
! dependencies: none
!----------------------------------------------------------------------------------------------------------
! Input:
!
! real::
! t     ... [K] effective temperature of star
! sigma ... [cm^-1](1:2) lower (1) and upper (2) boundary of wave number  
!----------------------------------------------------------------------------------------------------------
! Output:
!
! real::
! flux   ... [W m^-2 sr^-1] band flux between sigma(1) and sigma(2)
!-----------------------------------------------------------------------------------------------------
implicit none
!input
integer::i,n,iter(1:2)
real(kind=dp),intent(in)::t,sigma(1:2)
!output
real(kind=dp),intent(out)::flux
!local
real(kind=dp)::dn
real(kind=dp),dimension(1:2)::x,x2,x3,summ,iterations

!------------------------------------------------------------------------------------------------------------
!parameters:
! c1=2*h*c^2           ...h: Planck's constant, c: vacuum speed of light, k: Stefan-Boltzmann constant 
! c2=h*c/k
!-----------------------------------------------------------------------------------------------------------
real(kind=dp),parameter::c1=1.1910427584934558E-12 ![W cm^2 sr^-1]
real(kind=dp),parameter::c2=1.4387751601679205  ![K cm]

		  !c1=1.1910427584934558d-16 ![W m^2 sr^-1]
		  !c2=1.4387751601679205d-2  ![K m]


 !integrate Planck curve in region sigma(1), sigma(2) by subtracting sigma(2) to infinity from sigma(1) to infinity   
 !sigma  [1/cm] : wave number = 1/wavelength
 
      summ(:)=0.d0
      
lambd: do i=1,2
      
         x(i)=c2*sigma(i)/t
         x2(i)=x(i)*x(i)
         x3(i)=x2(i)*x(i)
         ! decide how many terms of sum are needed
         iterations(i)=2.d0+20.d0/x(i)
         
         if (iterations(i)>512.d0) then
            iterations(i)=512.d0
         end if

         iter(i)=nint(iterations(i))         
         
         !  add up terms of sum
         do n=1,iter(i)
            dn=1.d0/n
            summ(i)=summ(i)+exp(-real(n)*x(i))*(x3(i)+(3.d0*x2(i)+6.d0*(x(i)+dn)*dn)*dn)*dn
         end do
         
      end do lambd
      
      
  !integral over Planck curve from sigma(1) to sigma(2) [W m^-2 ster^-1]
  flux=c1*(t/c2)**4.d0*abs((summ(2)-summ(1)))*1.d4

return
end subroutine

!******************************************************************************************************** 
 subroutine transitfactor(nrad,koos,koop,rs,tf)
   implicit none
!---------------------------------------------------------------------------------------------------------
! Calculates the transit factor (tf) that reduces the planetary insolation in a multistar system
! due to mutual stellar occultations.
!
! tf:= 1 - (fraction luminous disc area occulted / total disc area of the star)  -> 2D! not surface!
!
! tf=1 means all the flux from this star arrives at the planet
! tf=0 means all the flux is blocked by one of the stars (from point of view of the planet)
!
!
! written by Siegfried Eggl 20120502
!
! dependencies: bxaxb
!----------------------------------------------------------------------------------------------------------
! Input:
!
! integer::
! nrad				...number of radiating bodies in the system
!
! real::
! koos(1:nrad,1:3)		...[au] xyz coordinate vectors of the stars in the system 
! koop(1:3)			...[au] xyz coordinate vectors of the planet
! rs(1:nrad)			...[Rsun] stellar radii
!----------------------------------------------------------------------------------------------------------
! Output:
!
! real::
! tf(1:nrad)			...[] transit factor 
!-----------------------------------------------------------------------------------------------------------

!input
   integer::nrad
   real(kind=dp),intent(in)::kooS(nrad,3),kooP(3),RS(nrad)
!output
   real(kind=dp),intent(out)::TF(nrad)
   
!local   
   integer::i,j,maxd
   real(kind=dp),dimension(1:3)::d,dpp,proj
   real(kind=dp),dimension(nrad)::RSm,RSm2,da
   real(kind=dp)::dproj,covarea,dmax,dproj2

   TF(:)=1.d0 
   covarea=0.d0

!calculate which star is in front and choose the distance vectors correctly
   do i=1,2
    dpp(:)=kooS(i,:)-kooP(:)
    da(i)=sqrt(dot_product(dpp(:),dpp(:)))
   end do

   if(da(1).gt. da(2))then
    d(:)=kooS(1,:)-kooS(2,:)
    dpp(:)=kooS(1,:)-kooP(:)
   elseif (da(1).lt. da(2))then
    d(:)=kooS(2,:)-kooS(1,:)
    dpp(:)=kooS(2,:)-kooP(:)
   end if
     
!calculate d vector's projection onto plane orthogonal to planet's line of sight going through the farthest star 
    call bxaxb(d,dpp,proj)

!length of projected vector in [km]
    dproj=Sqrt(Dot_Product(proj,proj))*au

!radius of stars in [km]
    RSm(:)=RS(:)*rsun

!transit
if (sum(RSm(:)).gt.dproj) then
   
   RSm2(:)=RSm(:)*RSm(:)
   dproj2=dproj*dproj  

  !cases for total disc coverage (passages of smaller star behind larger star, etc)
  if (dproj.le.abs(maxval(RSm(:))-minval(RSm(:)))) then
   if(da(1).gt. da(2)) then
     TF(1)=max(1.d0-(RSm2(2)/RSm2(1)),0.d0)
     TF(2)=1.d0
    elseif(da(1).lt. da(2))then
     TF(2)=max(1.d0-(RSm2(1)/RSm2(2)),0.d0)
     TF(1)=1.d0
   end if
  else 
!calculate covered area following (grazing goat) circle circle intersection. 
    covarea=RSm2(1)*Acos((dproj2+RSm2(1)-RSm2(2))/(2.d0*dproj*RSm(1)))+  &
            RSm2(2)*Acos((dproj2+RSm2(2)-RSm2(1))/(2.d0*dproj*RSm(2)))-  &
            0.5d0*Sqrt((-dproj+RSm(1)+RSm(2))*(dproj+RSm(1)-RSm(2))*(dproj-RSm(1)+RSm(2))*(dproj+RSm(1)+RSm(2)))
  
   if(da(1).gt. da(2))then
    TF(1)=1.d0-(covarea/(PI*RSm2(1)))
   elseif(da(1).lt. da(2))then
    TF(2)=1.d0-(covarea/(PI*RSm2(2)))
   end if
  end if
else 
 TF(:)=1.d0
end if

   return
end subroutine   

!*******************************************************************************************
subroutine zamsmrlt(ms,z,lum,teff,rs)
!----------------------------------------------------------------------------------------------
! Zero Age Main Sequence (ZAMS) relations between mass (ms), metallicity (z), luminosity (lum) and stellar radius (rs)
! this subroutine calculates radius luminosity and effective temperature 
! for a given stellar mass and metallicity following 
! Tout et al. 1996
!
! valid for z : [0.0001, 0.03],  M star : [0.1, 100] 
!
! written by Siegfried Eggl 20140609
!
! dependencies: none
!----------------------------------------------------------------------------------------------------------
! Input:
!
! real:: 
! ms		...[Msun] mass of the star
! z		...[] metallicity (z_sun=0.02)
!----------------------------------------------------------------------------------------------------------
! Output:
!
! real::
! lum         	...[Lsun] stellar luminosity
! teff		...[K] effective temperature
! rs		...[Rsun] stellar radius
!----------------------------------------------------------------------------------------------
implicit none
!input
real(kind=dp),intent(in)::ms,z !mass of the star [Msun], metallicity (Zsun~0.02)
!output
real(kind=dp),intent(out)::lum,teff,rs !luminosity [Lsun], effective temperature [K], radius [Rsun]
!local
integer::i
real(kind=dp),parameter::zsun=0.02d0  		!metallicity of the Sun
real(kind=dp),dimension(1:8,1:5)::lc 		!luminosity fit coefficients
real(kind=dp),dimension(1:10,1:5)::rc 		!radius fit coefficients
real(kind=dp),dimension(1:8)::lcz,pml,ml 	!luminosity fit coefficients (eq. 3), and powers of mass (eq. 1)
real(kind=dp),dimension(1:10)::rcz,pmr,mr 	!luminosity fit coefficients (eq. 4), and powers of mass (eq. 2)
real(kind=dp),dimension(1:5)::logz 		!log(Z/Zsun)

pml=(/5.5d0,11.d0,0.d0,3.d0,5.d0,7.d0,8.d0,9.5d0/)
pmr=(/2.5d0,6.5d0,11.d0,19.d0,19.5d0,0.d0,2.d0,8.5d0,18.5d0,19.5d0/)

do i=1,8
 ml(i)=ms**pml(i)
end do
do i=1,10
 mr(i)=ms**pmr(i)
end do

lc(1,:)=(/0.39704170, -0.32913574,0.34776688,0.37470851,0.09011915/)
lc(2,:)=(/8.52762600,-24.41225973,56.43597107,37.06152575,5.45624060/)
lc(3,:)=(/0.00025546,-0.00123461,-0.00023246,0.00045519,0.00016176/)
lc(4,:)=(/1.d0,0.d0,0.d0,0.d0,0.d0/)
lc(5,:)=(/5.43288900,-8.62157806,13.44202049,14.51584135,3.39793084/)
lc(6,:)=(/5.56357900,-10.32345224,19.44322980,18.97361347,4.16903097/)
lc(7,:)=(/0.78866060,-2.90870942,6.54713531,4.05606657,0.53287322/)
lc(8,:)=(/0.00586685,-0.01704237,0.03872348,0.02570041,0.00383376/)


rc(1,:)=(/1.71535900,0.62246212,-0.92557761,-1.16996966,-0.30631491/)
rc(2,:)=(/6.59778800,-0.42450044,-12.13339427,-10.73509484,-2.51487077/)
rc(3,:)=(/10.08855000,-7.11727086,-31.67119479,-24.24848322,-5.33608972/)
rc(4,:)=(/1.01249500,0.32699690,-0.00923418,-0.03876858,-0.00412750/)
rc(5,:)=(/0.07490166,0.02410413,0.07233664,0.03040467,0.00197741/)
rc(6,:)=(/0.01077422d0,0.d0,0.d0,0.d0,0.d0/)
rc(7,:)=(/3.08223400,0.94472050,-2.15200882,-2.49219496,-0.63848738/)
rc(8,:)=(/17.84778000,-7.45345690,-48.96066856,-40.05386135,-9.09331816/)
rc(9,:)=(/1.d0,0.d0,0.d0,0.d0,0.d0/)
rc(10,:)=(/0.00022582,-0.00186899,0.00388783,0.00142402,-0.00007671/)

logz(1)=1.d0
do i=2,4
logz(i)=logz(i-1)*log10(z/zsun)
end do

do i=1,8
 lcz(i)=dot_product(lc(i,:),logz(:))
end do

do i=1,10
 rcz(i)=dot_product(rc(i,:),logz(:))
end do

lum=dot_product(lcz(1:2),ml(1:2))/dot_product(lcz(3:8),ml(3:8))
rs=dot_product(rcz(1:5),mr(1:5))/dot_product(rcz(6:10),mr(6:10))

teff=(lum/(rs*rs))**0.25d0
return
end subroutine
 

!#############################################################################################
end module
