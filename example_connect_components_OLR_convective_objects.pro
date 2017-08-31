pro example_connect_components_OLR_convective_objects

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;19 July 2013, Bethan White
;14 Feb 2014, Aeron Buchanan - small cluster culling, etc
;04 Mar 2014, Bethan White - edit to process CASCADE OLR data
;
;this procedure reads in an OLR field from WRF output data, reduces it to a binary field according to
;a set threshold value (olrmin), and then finds connected components within that binary field.
;The aim is to identify the degree of organisation of convection. 
; 
; This procedre calls the function bwamb_ccl.pro which performs the fast CCL step
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; read in OLR data
f = '/home/cumulus/cp/Users/white/data/CASCADE/4km_exp/olr_4kmExp_reduced_domain_coarsened_12km.nc'
ncid = NCDF_OPEN(f)
zid = NCDF_VARID(ncid,'olr')
NCDF_VARGET, ncid,  zid, olr
zid = NCDF_VARID(ncid,'x')
NCDF_VARGET, ncid,  zid, x 
zid = NCDF_VARID(ncid,'y')
NCDF_VARGET, ncid,  zid, y
NCDF_CLOSE, ncid

; set threshold OLR for deep convection
olrmin = 150.

; next we need to create a new array where we assign a value of 1 to those elements that are below the threshold and a value 
; of 0 to those elements that are above the threshold    
; copy olr array so we don't overwrite original data
olr_binary = olr
olr_binary[WHERE(olr LE olrmin, /NULL)] = 1
olr_binary[WHERE(olr GT olrmin, /NULL)] = 0

; a secondary threshold as part of "non-max supppression"
olrmin2 = 175.

olr_bin2 = olr
olr_bin2[WHERE(olr LE olrmin2, /NULL)] = 1
olr_bin2[WHERE(olr GT olrmin2, /NULL)] = 0
 
; NEED TO DO THE FOLLOWING FOR UM DATA: 
; rotate lat/lon co-ords from UM  grid to lat and long
; (code taken from John Marsham's plot_umseries_csip_thsw.pro)
x2d = fltarr( (size(x))(1), (size(y))(1) )
y2d = fltarr( (size(x))(1), (size(y))(1) )
for i = 0,(size(y))(1)-1 do x2d(*,i) = x
for i = 0,(size(x))(1)-1 do y2d(i,*) = y

eqtoll, y2d,x2d,  lat2d,long2d, 79.0, 180.0
swaplong = where(long2d gt 180.)

; now convert UM longitudes (0 at Greenwich, 360 going east)
; to standard longitude (0 at Greenwich, + 180 going east, -180 going west)
long2dnew = long2d 
long2dnew(swaplong) = long2d(swaplong) - 360.

; if we're happy, overwrite the original array!
long2d = long2dnew

gridsize = 4. ; horizontal grid length in km; this is for 4 km UM data

; FOR UM DATA:
x_vals = findgen((size(x))(1)) * gridsize
y_vals = findgen((size(y))(1)) * gridsize

times = indgen((size(olr))(3))
printtimes = string(times)

; set plotting symbol to be a filled square (plotsym, sym, symsize):
plotsym, 8, 0.5, /fill

; distance threshold within which objects will be assigned the same label
dmax = 1 ; dmax = 1 means not using the radius merge feature

; plot the connected components labelled regions on a map
window, 0, retain = 2
device, decomposed = 0,retain = 2
loadct, 39
TVLCT,R,G,B,/GET ;Put the table into IDL variables
RR = REVERSE(R)
GG = REVERSE(G)
BB = REVERSE(B)
TVLCT,RR,GG,BB ;Load the reversed colors

; also swap the default foreground/background colors around!
!p.color = 256*256-1
!p.background = 0

; Flags for whether to use certain features of the algorithm:
CULL = 0 ; cull small clusters with sizes below a certain threshold
TWOPASS = 1 ; perform 2-pass CCL
NONMAX = 1 ; perform non-maximum suppression

; if CULL-ing then ignore blobs below this size (area in grid points)
SMALL_THRESHOLD = 9

; label strings for writing out data
LABEL = ''
if NONMAX then LABEL+='NONMAX'
if TWOPASS then LABEL+='TWOPASSCULL'
if CULL then LABEL+='CULL'

if strlen(LABEL) gt 0 then LABEL+='_'

; store the ccl output array
clustered_array = fltarr( (size(olr))(1), (size(olr))(2), (size(olr))(3) )

; loop over time and analyse the convective regions
; we incorporate an option to plot a lat/lon map of the identified convective objects
; in the same time loop, after the object identification step is performed

for t = 0, (size(olr))(3) - 1 do begin
   ; we use label_region to identify connected components using a 4-neighbour search

   input = olr_binary(*,*,t)

	if NONMAX then begin
		; approx non-max suppression
		; spot all areas beyond upper (necessary) threshold
		; then expand areas to include everywhere beyond a lower (sufficient) threshold
		; areas that are wholy beneath upper threshold are ignored
		
		; 1: cluster the lower threshold
		lower_threshold_dmax = dmax
		time_start = systime(1)
		ccl2 = bwamb_ccl(olr_bin2(*,*,t),lower_threshold_dmax)
		time_stop = systime(1) 
		print, 'CCL TIMING (2) (d=',lower_threshold_dmax,'): ', time_stop - time_start, ' (s)'

		; 2. find all the lower-threshold clusters that have the upper-threshold areas within them
		good_labels = input * ccl2
		good_labels = good_labels(sort(good_labels))
		good_labels = good_labels(uniq(good_labels))
		for j = 1, max(ccl2) do begin
			wh = where(good_labels eq j, count)
			if count eq 0 then begin
				; not a good label -> cull
				ccl2(where(ccl2 eq j, /null)) = 0
			endif
		endfor

		input = ccl2
	endif

	cull_count = 0
	if TWOPASS then begin
		pre_dmax = 1 ; adjacent only
		time_start = systime(1)
		ccl = bwamb_ccl(input,pre_dmax) 
		time_stop = systime(1) 
		print, 'CCL TIMING (d=',pre_dmax,'): ', time_stop - time_start, ' (s)'

		; TODO: create a function for culling
		if CULL then begin
			for j = 1, max(ccl) do begin
				wh = where(ccl eq j, count)
				if count le SMALL_THRESHOLD then begin 
					; discard clusters that are too small
					cull_count += 1
					ccl(wh) = 0
					; TODO: remap labels
				endif
			endfor	
		endif
		print, 'cull count = ', cull_count

		input = ccl
	endif
	
   time_start = systime(1)
   ccl = bwamb_ccl(input,dmax)
   time_stop = systime(1) 
   print, 'CCL TIMING (d=',dmax,'): ', time_stop - time_start, ' (s)'
 
   nblobs = max(ccl) ; find the number of blobs

	; CULL step (if used)
	cull_count = 0
	if CULL then begin
		for j = 1, nblobs do begin
			wh = where(ccl eq j, count)
			if count le SMALL_THRESHOLD then begin 
				; discard clusters that are too small
				cull_count += 1
				ccl(wh) = 0
				; TODO: remap labels
			endif
		endfor
		print, 'cull count = ', cull_count
	endif

	print, "nblobs: ", nblobs - cull_count

	;store the clustered blobs
	clustered_array(*,*,t) = ccl

	; do we want to plot the clusters?
	print_plot = 1

	if print_plot EQ 1 then begin

	   ; define start color (as not white!) 
		col = 1.

	    ; sort out string length for times - add leading zeroes so all strings are same length
	    ; values from 0:9
	    if strlen(strtrim(printtimes(t),1)) EQ 1 then printtimes(t) = '00'+strtrim(printtimes(t),1)
	    ; values from 10:99
	    if strlen(strtrim(printtimes(t),1)) EQ 2 then printtimes(t) = '0'+strtrim(printtimes(t),1)
   
		; set plot position and map ranges
	   !p.position=[0.15,0.15,0.95,0.9]

	   latmin = min(lat2d)
	   latmax = max(lat2d)
	   lonmin = min(long2d)
	   lonmax = max(long2d)		

	   map_set,limit = [ latmin, lonmin, latmax, lonmax ]
	   map_continents, /hires
 
	   xr = [lonmin, lonmax]
	   yr = [latmin, latmax]

	   ; loop through the ccl array such that each blob is plotted one by one:
	   for j = 1, nblobs do begin
	      wh = where(ccl eq j, count, complement = c, /null)
			if count eq 0 then continue ; don't try and draw nothing (erased due to culling)
	      wheretomulti, ccl, wh, cols, rows
	        if count gt 1e6 then continue ; don't draw clusters that fill the screen TODO: make sensible
      
		  ; plot the jth blob
		  plot, long2d(wh), lat2d(wh), psym = 8, /noerase, xrange = xr, yrange = yr, xstyle = 1, ystyle = 1, color = col, symsize = 1, xminor = 1, yminor = 1, charsize = 1.	

	      col+=30 ;change color
	   endfor ;looping over blobs

	   ; overplot axes in black (bloody idl)
	   plot, [0], [0], psym = 8, /noerase, xrange = xr, yrange = yr, xstyle = 1, ystyle = 1, color = 255, symsize = 1, xtitle = 'longitude', ytitle = 'latitude', xminor = 1, yminor = 1, charsize = 1, title = 'CASCADE CCL, 4km exp '+strtrim(printtimes(t),1)+' hours, olrmin(1,2) = ('+strtrim(fix(olrmin),1)+','+strtrim(fix(olrmin2),1)+')'

	    saveimage,'/a/cplxfs3/srv/homes/aopp/whiteb/plots/CASCADE/reduced_domains/4km_exp/clustered_olr/CASCADE_4kmexp_coarsened_12km_connected_components_bwamb_OLR_'+LABEL+'d'+strtrim(fix(dmax),1)+'_t'+strtrim(printtimes(t),1)+'.png', PNG=1  

	endif ; only if we want to plot maps of the objects
endfor ;looping over time

; ---------- now write the labelled convective objects to netcdf -------------
print, 'saving data to netcdf'
;create netcdf file
datadir = '/home/cumulus/cp/Users/white/data/CASCADE/4km_exp/'
file_id = ncdf_create(datadir+'CASCADE_4km_exp_coarsened_12km_reduced_domain_olr_connected_components_d_'+strtrim(fix(dmax),1)+LABEL+'_olrmin_'+strtrim(fix(olrmin),1)+'_olrmin2_'+strtrim(fix(olrmin2),1)+'.nc', /clobber)

;define dimensions
dim1 = ncdf_dimdef(file_id,'xdim',(size(olr))(1))
dim2 = ncdf_dimdef(file_id,'ydim',(size(olr))(2))
dim3 = ncdf_dimdef(file_id,'times',(size(olr))(3))
;define the array
vid = ncdf_vardef(file_id, 'clusters',[dim1,dim2,dim3],/float)
; Exit 'define mode', entry 'data mode'
ncdf_control, file_id, /endef
; Write the variable to the file
ncdf_varput, file_id, vid, clustered_array
; Close the file
ncdf_close, file_id
;stop

print, 'ALL DONE!'
stop
end
