pro area_weighted_avg_gen,lon,lat,var,weighted_avg 

;-----------------------------------------------------------------
;AUTHOR WOLF E.T.
;7/11/2008
;
;PURPOSE: Calculate area weighted average of geophysical quantities
;from CAM data. 
;
;Error handling: program expects fill_value of -999.0.  Testing to
;ensure error handling capabilities work has not been conducted.
;
;NOTES:  need to subtract duplicate pole
; Kludged to double precision 10/9/17, not sure it makes a difference.
;-----------------------------------------------------------------
;ARGUMENTS:
;lon            =>  longitude array from CAM
;lat            =>  lattitude array from CAM
;var[lon,lat]   =>  atmospheric variable we wish to calculate the area weighted
;                   average of.  
;weighted_avg   =>  returns the area weighted average
;-----------------------------------------------------------------


R_EARTH=1.0                            ; the area_weighted_average is indepedent of planet radius
fill_value=-999.0                      ;set fill-value of variable being used


;Determine number of elements in lat and lon arrays
nlat=n_elements(lat)
nlon=n_elements(lon)
weighted_avg=0.0
missing_area=0.0


area=fltarr(nlon,nlat)
;create staggered lon lat arrays
slon = fltarr(nlon+1)
slat = fltarr(nlat+1)

slon(0:nlon-1) = lon(*)
slon(nlon) = 360.

slat(0)=-90.
slat(nlat)=90
for ya=1,nlat-1 do begin
 slat(ya)=(lat(ya-1)+lat(ya))/2.0
endfor

;calculat area
for xa=0,nlon-1 do begin
  for ya=0,nlat-1 do begin
    area(xa,ya) = double(R_EARTH^2.*!pi/180.*(slon(xa+1)-slon(xa))*(sin(slat(ya+1)*!pi/180)-sin(slat(ya)*!pi/180)))
  endfor
endfor


sumArea=0.0
FOR x=0,nlon-1 DO BEGIN
  FOR y=0, nlat-1 DO BEGIN
    IF(var(x,y) ne fill_value) THEN BEGIN
      weighted_avg=double(weighted_avg+area(x,y)*var(x,y))
      sumArea=double(sumArea+area(x,y))
    ENDIF ELSE BEGIN
      missing_area=double(missing_area+area(x,y))
    ENDELSE
  ENDFOR
ENDFOR

weighted_avg=double(weighted_avg/sumArea)
;print, "weighted_avg", weighted_avg
end
