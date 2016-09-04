;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;*NAME: BWAMB_CLUSTERINDEX
;*PURPOSE: determine a 'clustering index' for data that has been clustered
;	(e.g. with bwamb_ccl.pro). 
;
;*CALLING SEQUENCE:
;   RESULT = BWAMB_FRAGINDEX(ARRAY [, MAX_MEANINGFUL_DISTANCE, [SCALE_INVARIANCE_PARAMETER]])
;
;*PARAMETERS:
; INPUTS:
; ARRAY =	2D data field. Elements that are non-zero must already have been labelled in increasing 
;			numerical order by a clustering algorithm (such as bwamb_ccl.pro)
; MAX_MEANINGFUL_DISTANCE =	The distance metric (in grid points) by which the inverse 
;							interaction potential will be scaled, in order to avoid scale 
;							invariance of the clustering metric. 
;							DEFAULT VALUE = 1, but note that this leads to scale invariance.
; SCALE_INVARIANCE_PARAMETER =	the parameter k by which the scale invariance treatement is weighted.
;								k = 0: the clustering metric is completely scale invariant
;								k = 1/2, 2 etc: takes into account the 'max meaningful interaction distance'
;								k = 1000: the metric is almost binary: everything will either interact
;										  with everything else, or will never interact.									
;								DEFAULT VALUE = 0 (scale invariant)
;
; CALCULATIONS:
;	1.	INVERSE INTERACTION POTENTIAL, DIIP(I,J), scaled by MAX_MEANINGFUL_DISTANCE
;		and SCALE_INVARIANCE_PARAMETER to avoid scale invariance (unless default values
;		are chosen, in which case DIIP is scale invariant). The interaction potential, IP 
;		is the inverse of DIIP. 
;	2.  CLUSTERING INDEX, CLUSTERI 
;
; OUTPUTS:
;   1. CLUSTERI: a clustering index value for the time frame passed in.
;
; REQUIRES: WHERETOMULTI.PRO
;
;*VERSION: 1.0
;
;*AUTHOR: BETHAN WHITE, March 2014, creator
;         AERON BUCHANAN, March 2014, development discussion
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function bwamb_clusterindex, _input_array, _input_max_distance, _input_scale_invariance_param

print, 'BWAMB CLUSTERINDEX V. 1.0'

; CHECKS:
; first check input parameters, and use defaults where required
if n_params() eq 0 then begin
    print, 'CALLING SEQUENCE: RESULT = BWAMB_FRAGINDEX(ARRAY [, MAX_MEANINGFUL_DISTANCE [, SCALE_INVARIANCE_PARAMETER]])'
    retall
endif

input_array = _input_array

; defaults
input_max_distance = 1
input_scale_invariance_param = 0

if n_params() ge 2 then input_max_distance = _input_max_distance
if n_params() ge 3 then input_scale_invariance_param = _input_scale_invariance_param


; and check that the input array is 2D:
if (size(input_array))(0) ne 2 $
    || ~ isa(input_array, /number, /array) $
then begin
    print, 'INPUT ARRAY MUST BE A 2D NUMERIC ARRAY'
    retall
endif

;+++++++++++++++++++++ GET CENTRAL POINTS OF EACH CLUSTER +++++++++++++++++++++++++++++++++++++++++++++++++++++
; get total number of clusters
nclusters = max(input_array)

; if we only have one cluster, set cluster index to 0 and exit!
IF nclusters LE 1. then begin

	clusteri = 0.
	return, clusteri

ENDIF ELSE BEGIN ;run full code as long as there is more than one cluster!

; find central point of each cluster
cluster_populations = fltarr(nclusters)
cluster_centres_x = fltarr(nclusters) ;x, y co-ords of cluster centre
cluster_centres_y = fltarr(nclusters)

for n = 1, nclusters do begin ; ignore large empty cluster!
    ; get cluster population
    cluster_pop = where(input_array eq n, /null, count)
    cluster_populations(n-1) = count

    ; find central point of cluster by summing the coordinates of every point in the cluster
    ; and dividing by the total cluster population
    cluster_coords = where(input_array eq n, /null, count)
    wheretomulti, input_array, cluster_coords, col, row

    xcentre = total(col) / cluster_populations(n-1)
    ycentre = total(row) / cluster_populations(n-1)

    cluster_centres_x(n-1) = xcentre
    cluster_centres_y(n-1) = ycentre
endfor

;+++++++++++++++++++++++ INVERSE INTERACTION POTENTIAL CALCULATION +++++++++++++++++++++++++++++++++++++++++++++
; Calculate the INVERSE INTERACTION POTENTIAL (diip) using the disstances from each point to every other point
; and accounting for the parameters used to avoid scale invariance
; naive scale-invariant inverse interaction potential between two clusters is:
;
; diip(i,j) = d(i,j)* sqrt(!pi) / sqrt(ci) + sqrt(cj) SCALE-INVARIANT VERSION
;
; where ci, cj are populations of each cluster and d(i,j) is the distance between the clusters
; to account for variations of the interaction potential with scale, we introduce two parameters:
; the maximum interaction distance (in grid points), D, and the scale invariance parameter, k
; these are introduced as follows: 
;
; diip (i,j) = ((d(i,j)* sqrt(!pi) / sqrt(ci) + sqrt(cj) ) * (d(i,j) / D)^k )^1/(k+1)  

; total number of connections
nc = (1./2.)*nclusters*(nclusters-1) ; total number of connections is (N(N-1))/2

distances_set = fltarr(nc) ; total number of connections is (N(N-1))/2
diip_set = fltarr(nc); same for diip
diip_naive_set = fltarr(nc) ;TODO: remove this when testing of non-scale-invariant version is complete
ip_set = fltarr(nc) ; interaction potential ( IP = 1/DIIP)

count = 0

d = input_max_distance
k = input_scale_invariance_param

;loop over total number of clusters but don't need to use the final one as all connections have already been counted
for i = 1, nclusters - 1  do begin 
    for j = i + 1, nclusters  do begin ; look to every other cluster that's not the current one

            xi = cluster_centres_x(i-1)
            xj = cluster_centres_x(j-1)
            yi = cluster_centres_y(i-1)
            yj = cluster_centres_y(j-1)

            ; get distance between clusters
            distance_ij = sqrt( (xi - xj)^2 + (yi - yj)^2 )
            distances_set(count) = distance_ij

			; get inverse interaction potential between clusters
            ci = cluster_populations(i-1)
            cj = cluster_populations(j-1)

            ; get naive scale-invariant diip
            diip_ij_naive = (distance_ij*sqrt(!PI) ) / ( sqrt(ci) + sqrt(cj) )
            diip_naive_set(count) = diip_ij_naive

            ; get diip accounting for scale invariance of first formulation
            diip = ( (distance_ij* sqrt(!PI) / ( sqrt(ci) + sqrt(cj) )) * (distance_ij / d)^k )
            diip_set(count) = diip

			ip = 1./diip
			ip_set(count) = ip

            count+=1
    endfor
endfor

;+++++++++++++++++ CLUSTERING INDEX +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; clustering index, clusteri, is calculated from the sum of the interaction potentials divided by the 
; total number of possible connection. For values of maximum interaction distance greater than 1 and scale
; invariance parameter greater than 0, the clustering index should be non-scale-invariant

clustering_index = total(ip_set) / ( (1./2.)*nclusters*(nclusters - 1))

clusteri = clustering_index

; pass out the clustering index
return, clusteri

ENDELSE ; more than one cluster
end ; function
