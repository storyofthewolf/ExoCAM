pro hybrid2pressure,nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi,lev_P,ilev_P

;AUTHOR WOLF E.T.
;7/11/2008
;-------------------------------------------------------------
;PURPOSE:  Convert WACCM hybrid level coordinates to pressure coordinates
;
;NOTES:  for implementation can change lon, lat, lev passing
;parameters to nlon, nlat, nlev to improve speed.  This is meant to be
;used as a function within larger data analysis packages
;
;-------------------------------------------------------------  
;ARGUMENTS:
;nlon                => number of elements, longitude array
;nlat                => number of elements, latitude array
;nlev                => number of elements, hybrid vertical level midpoint coordinates
;PS[lon,lat]         => surface pressure array at each grid point
;P0                  => reference pressre
;hyam                => hybrid A coefficient at layer midpoints
;hybm                => hybrid B coefficient at layer midpoints
;hyai                => hybrid A coefficient at layer interfaces
;hybi                => hybrid B coefficient at layer interfaces
;lev_P[lon,lat,lev]  => returns mid layer pressure coordinate matrix, must be
;                       defined with proper dimensions before calling this procedure
;ilev_P[lon,lat,lev+1] => returns interface pressure coordinate matrix, must be
;                       defined with proper dimensions before calling this procedure
;---------------------------------------------------------------
ilev_P=fltarr(nlon,nlat,nlev+1)
;nlon=n_elements(lon)
;nlat=n_elements(lat)
;nlev=n_elements(lev)

FOR x=0, nlon-1 DO BEGIN
  FOR y=0, nlat-1 DO BEGIN
    FOR z=0,nlev-1 DO BEGIN
      lev_P(x,y,z)=hyam(z)*P0+hybm(z)*PS(x,y)
    ENDFOR
  ENDFOR
ENDFOR

FOR x=0, nlon-1 DO BEGIN
  FOR y=0, nlat-1 DO BEGIN
    FOR z=0,nlev DO BEGIN
      ilev_P(x,y,z)=hyai(z)*P0+hybi(z)*PS(x,y)
    ENDFOR
  ENDFOR
ENDFOR

end
