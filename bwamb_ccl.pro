;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;*NAME: BWAMB_CCL
;
;*PURPOSE: connected-component or cluster labelling using a within radius test.
;	It has been written with execution speed in mind so should generally run 
;	quite fast. See examples for more information. The distance threshold test
;	is based on the L2 norm (pythagorean distance) by default, but any norm, 
;	or indeed any arbitrary inclusion mask can be used. NB: the generalizations
;	regarding distance metric and inclusion mask 'come for free' with the 
;	fast implementation, i.e. no speed compromize for these features.
;
;*CALLING SEQUENCE:
;	RESULT = BWAMB_CCL(ARRAY [, DISTANCE_THRESHOLD [, DISTANCE_METRIC])
;
;*PARAMETERS:
; INPUTS:
;	ARRAY = 2D data field. Elements that are non-zero will be labelled
;	DISTANCE_THRESHOLD = Two non-zero elements within this radius will
;	                     be given the same label. DEFAULT VALUE = 1
;	                     NB: If an array is given then it is used as the
;	                     inclusion mask.
;	DISTANCE_METRIC = How the distance between two elements is calculated.
;	                  'l1' = L1 norm (manhatten distance, diamond)
;	                  'l2' = DEFAULT L2 norm (pythagorean distance, circle)
;	                  'linf' = L-infinity norm (max coord distance, square)
;	                  scalar = arbitrary norm (e.g. 1 for l1, 2 for l2, 1.5, 10)
;	                  This argument is ignored if an array is given above
;	
; OUTPUTS:
;	1. CCL: an array with the connected components labelled
;
; EXAMPLES:
;	BWAMB_CCL(ARRAY, 1) gives the 4-neighbour connected component labelling
;	BWAMB_CCL(ARRAY, 1.5) gives the 8-neighbour connected component labelling
;	BWAMB_CCL(ARRAY, 50) clusters together elements such that there is a path
;		between any two elements with the same label, using only elements with
;		that label, where all jumps are less than 50 units. 1 unit is the distance
;		between adjacent 4-neighbours. For any two elements with different labels
;		no such path exists.
;
;*VERSION: 3.32
;
;*AUTHOR: BETHAN WHITE, July 2013, creator
;         AERON BUCHANAN, December 2013, added distance threshold capability
;         AERON BUCHANAN, December 2013, significant speed improvements
;         AERON BUCHANAN, December 2013, generalizations
; 
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function bwamb_ccl, _input_array, _input_threshold, _input_norm 

print, 'BWAMB CCL V. 3.32'

;t1 = systime(1) ; DEBUG

; CHECKS

; first check input parameters, and use defaults where required
; plus, we must protect the caller from pass-by-reference

if n_params() eq 0 then begin
	print, 'CALLING SEQUENCE: RESULT = BWAMB_CCL(ARRAY [, DISTANCE_THRESHOLD [, NORM]])'
	retall
endif

input_array = _input_array

; defaults
input_threshold = 1
input_norm = 'l2'

if n_params() ge 2 then input_threshold = _input_threshold
if n_params() ge 3 then input_norm = _input_norm
;if n_params() ge 4 then input_param = _input_param

; and check that the input array is 2D:
if (size(input_array))(0) ne 2 $
	|| ~ isa(input_array, /number, /array) $ 
then begin
	print, 'INPUT ARRAY MUST BE A 2D NUMERIC ARRAY'
	retall 
endif

; save some typing
N = (size(input_array))(1)
M = (size(input_array))(2)

; set up storage array for labels
connected_labellings = lonarr(N, M) ; i, j

; find all the elements which need processing and label them
indices = where( input_array ne 0, count )
; protect ourselves from over-runs 
input_array(indices) = 1

; check whether there is anything to do
if count eq 0 then begin
	; there is nothing to cluster
	return, connected_labellings
endif

; initialize all the labels
num_labels = n_elements(indices)
; note we reserve the zero label for background, so start label numbers at 1
initial_labels = indgen(num_labels, /long) + 1L
; assign
connected_labellings(indices) = initial_labels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CREATE INCLUSION MASK
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;t2 = systime(1) ; DEBUG

; note that beause distance is commutative d(a,b) = d(b,a) only two quadrants
; of the inclusion mask are actually needed to have comparisons between all
; pairs of elements. However, using the full inclusion mask means more
; equivalences are found more quickly, speeding up the first pass.

; check whether the threshold argument is actually a user-supplied mask

; if (size(input_threshold))(0) eq 2 then begin
if isa(input_threshold, /number, /array) then begin
	; user supplied mask 

	i_diameter = (size(input_threshold))(1)
	j_diameter = (size(input_threshold))(2)
	
	need_to_resize = 0
	if long(i_diameter / 2L) * 2L - i_diameter eq 0 then begin ; modulo even
		need_to_resize = 1
		i_diameter += 1L
	endif
	if long(j_diameter / 2L) * 2L - j_diameter eq 0 then begin ; modulo even
		need_to_resize = 1
		j_diameter += 1L
	endif
	
	if need_to_resize then begin
		temp = dblarr(i_diameter, j_diameter)
		temp[0, 0] = input_threshold
		input_threshold = temp
	endif
	
	inclusion_mask = lonarr(i_diameter, j_diameter)	
	inclusion_mask( where( input_threshold , /null ) ) = 1L
	
	i_radius = ( i_diameter - 1L ) / 2L
	j_radius = ( j_diameter - 1L ) / 2L
	
endif else begin
	; convert input_threshold to be sure of sufficient precision
	input_threshold = double(input_threshold)
	
	base_radius = floor(input_threshold)
	base_diameter = 2L * base_radius + 1L

	i_radius = base_radius
	j_radius = base_radius
	
	i_diameter = base_diameter
	j_diameter = base_diameter

	; create an array holding the distance of each element from the center of the mask
	distances = lonarr(i_diameter, j_diameter)

	; create matrices of horizontal and vertical distance from center
	dis = ( indgen(i_diameter, /long) # make_array(1, j_diameter, /long, value=1) ) - i_radius
	djs = ( make_array(i_diameter, 1, /long, value=1) # indgen(j_diameter, /long) ) - j_radius

	threshold = 1L

	; if isa(input_norm, /scalar) then begin ; FAIL?
	if (size(input_norm))(0) eq 0 and ( (size(input_norm))(1) eq 2 or (size(input_norm))(1) eq 3 or (size(input_norm))(1) eq 4 or (size(input_norm))(1) eq 5 ) then begin
		p = double(input_norm)
		input_norm = 'lp'
	endif

	case input_norm of
		'l1': begin ; l1 norm, i.e. diamond
				distances = abs(dis) + abs(djs)
				threshold = input_threshold
			end
		'l2': begin ; l2 norm, i.e. circle
				distances = dis * dis + djs * djs
				threshold = input_threshold * input_threshold
			end
		'lp': begin ; lp norm, e.g. rounded square
				distances = abs(dis) ^ p + abs(djs) ^ p
				threshold = input_threshold ^ p
			end
		'linf': begin ; l-infinity norm, i.e. square
				; actually, there is nothing to be done here, but below is technically the correct calculation, i.e. max
				; indices = where( abs(djs) gt abs(dis) , /null )
				; distances = dis & distances( indices ) = djs( indices ) ; max
				; threshold = input_threshold
			end
		else: begin
				print, "UNKNOWN DISTANCE METRIC: ", input_norm
				retall
			end
	endcase
	
	; create a mask indicating which elements are within radius
	inclusion_mask = lonarr(i_diameter, j_diameter)
	inclusion_mask( where( distances le threshold , /null ) ) = 1L

endelse

; ensure the central element is set
inclusion_mask(i_radius, j_radius) = 1 

; DEBUG
;print, "CHECK MASK!"
;print, inclusion_mask
;PAUSE

; double check whether there is anything to do
if count gt n_elements(input_array) - i_radius * N + 1L and $
	count gt n_elements(input_array) - j_radius * M + 1L $
then begin
	; not enough empty elements to separate groups by enough rows/columns
	; therefore everything is connected
	connected_labellings(indices) = 1
	return, connected_labellings
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PASS 1 - determine equivalence of labels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;t3 = systime(1) ; DEBUG
;debug_count = 0L ; DEBUG

; also, to help updates to the labels in the connected_labellings array
; (it helps efficiencies if there are as few different label numbers
; being processed as possible) create coordinate reference arrays
i_coords = indgen(N, /long) # make_array(1, M, /long, value=1)
j_coords = make_array(N, 1, /long, value=1) # indgen(M, /long)
; double check because wtf idl dimension consistency
;if (size(i_coords))(1) ne N or (size(i_coords))(2) ne M or i_coords(123, 234) ne 123 then print, "FAILURE (case I)"
;if (size(j_coords))(1) ne N or (size(j_coords))(2) ne M or j_coords(123, 234) ne 234 then print, "FAILURE (case J)"

; to determine the equivalences efficiently, we're going to use an array
; where the n-th element is a number representing canonical label in the
; equivalence group of the label n
equiv_refs = indgen(num_labels + 1, /long) ; "+1" because labelling starts at '1' but array indices start at '0'
; start by setting the equivalence group of each label to be itself

; array to store values of neighbouring labels
equiv_labels = lonarr(num_labels) 

; keep track of the number of labels found during each search
label_count = 0L

; IDL stores matrices with (i,j) and (i+1,j) next to each other in memory
for j = 0L, M - 1 do begin
	
	for i = 0L, N - 1 do begin
	
		; if the pixel is not in the background:
		if input_array(i, j) ne 0 then begin
		
			; determine the index bounds of the region we are looking at
			min_i_raw = i - i_radius
			min_j_raw = j - j_radius
			min_i = min_i_raw
			min_j = min_j_raw
			if min_i_raw lt 0 then min_i = 0
			if min_j_raw lt 0 then min_j = 0
			max_i_raw = i + i_radius
			max_j_raw = j + j_radius
			max_i = max_i_raw
			max_j = max_j_raw
			if max_i_raw ge N then max_i = N - 1
			if max_j_raw ge M then max_j = M - 1
			; the 'if's are slightly faster than max/min function calls
			
			; get the correct section of the distance mask
			mmin_i = 0L
			if min_i_raw lt 0 then mmin_i = mmin_i - min_i_raw
			mmin_j = 0L
			if min_j_raw lt 0 then mmin_j = mmin_j - min_j_raw
			mmax_i = 2L * i_radius
			if max_i_raw ge N then mmax_i = mmax_i - max_i_raw + N - 1
			mmax_j = 2L * j_radius
			if max_j_raw ge M then mmax_j = mmax_j - max_j_raw + M - 1
			
			; find coords of all within-radius non-background elements
			ws = where( input_array(min_i:max_i, min_j:max_j) * inclusion_mask(mmin_i:mmax_i, mmin_j:mmax_j) , /null )
			
			; gdl work-around, although a minor speed up too perhaps
			; ( gdl 0.9.2 shift([1], 1) core dumps )
			if n_elements(ws) gt 1 then begin
				is = (i_coords(min_i:max_i, min_j:max_j))(ws)
				js = (j_coords(min_i:max_i, min_j:max_j))(ws)
				
				; grab the section of the matrix we are working with
				labels_vector = connected_labellings(is, js)
				
				; disard duplicates (it is faster to perform this unique search
				; than it is to deal with possibly large numbers of duplicates
				; in the while loop below)
				; the uniq function is slower than it should be, so 
				; find all the unique labels by first sorting them
				labels_vector = labels_vector(sort(labels_vector))
				; then picking out the ends of runs
				unique_indices = where( labels_vector - shift(labels_vector, 1), label_count)
				; and saving those values
				if label_count eq 0 then begin
					labels_vector = [labels_vector(0)]
					label_count = 1
				endif else begin
					labels_vector = labels_vector(unique_indices) ; at least 2 as there can't be only one boundary
				endelse
			endif else begin
				; must be 1 because input_array says so!
				label_count = 1
				; no need to set labels_vector
			endelse
						
			; if only one, then there is nothing to do: skip to next
			if label_count le 1 then continue
			
			; DEBUG
			;debug_count += 1
			
			; DEBUG
			;if n_elements(labels_vector) ne label_count then print, "FAILURE (case Z)"
			
			; lookup the label references of all the labels in this set and
			; find the canonical label for this group, i.e. the smallest label
			min_label = min(equiv_refs(labels_vector))
			
			; DEBUG
			;if min_label eq 0 then print, "FAILURE (case A)"
			
			; save in our "extendable" array
			equiv_labels[0] = labels_vector
			
			; go through the current set and make all labels reference the label
			; that represents this group, i.e. min_label
			; NOTE: in order to propagate this group's canonical label, all
			; referenced labels must also be updated. To achieve this they are
			; added to this set during the loop. To account for this, we need the while
			index = 0L ; index variable
			while index lt label_count do begin ; while there are still set elements to see (length might change)
				this_label = equiv_labels(index) ; for easy referencing
				ref_label = equiv_refs(this_label)
				; is the ref label in our to-do list?
				x = where( ref_label eq equiv_labels, count )
				if count eq 0 then begin
					; this label's reference is another label that we will need to update
					; so add it to the current set
					equiv_labels(label_count) = ref_label
					; update length variable
					label_count += 1
					; check whether this should be the canonical label
					if ref_label lt min_label then begin
						min_label = ref_label
					endif
				endif
				
				; DEBUG
				;if this_label lt min_label then print, "FAILURE (case B)"
				
				; look at next label in set
				index += 1
			endwhile
			
			; DEBUG
			;if min_label eq 0 then print, "FAILURE (case C)"
			
			; update the label references
			equiv_refs(equiv_labels(0:label_count-1)) = min_label
			
			; DEBUG
			;if equiv_refs(0) ne 0 then print, "FAILURE (case E1) ", i, ",", j
			
			; update
			connected_labellings(is, js) = min_label
			
		endif
	endfor 
endfor 

; DEBUG
;print, "equiv_refs"
;print, equiv_refs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; INTERMEDIATE STEP: record label equivalences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;t4 = systime(1) ; DEBUG

; now for some intermediate processing
; the equiv_refs hold the canonical labels, which now need to be replaced
; with a set of contiguous label numbers starting at 1
contig_refs = lonarr(num_labels + 1L)
new_label_count = 0L
for label = 1L, num_labels do begin
	; every label reference in equiv_refs is that label's group's canonical label
	; i.e. every entry in equiv_refs must either point to itself or 
	; recursively point to a lower entry that eventually points to itself
	
	if equiv_refs(label) eq label then begin
		; points to itself, therefore canonical
		new_label_count += 1
		new_label = new_label_count
	endif else begin
		; we need to do some searching 
		refs = lonarr(label) ; keep track of where we have been 
		refs(0) = equiv_refs(label) ; start at the beginning
		recurse_count = 0L ; count how many times we have stepped down the chain
		while equiv_refs(refs(recurse_count)) ne refs(recurse_count) do begin
			; this relies on the condition explained above
			recurse_count += 1
			; follow the references down one more step
			refs(recurse_count) = equiv_refs(refs(recurse_count - 1))
		endwhile
		; we have reached the reference that refers to itself, i.e. the canonical label
		min_label = refs(recurse_count)
		; go back and make a note of the final reference to make future searches faster
		if recurse_count gt 0 then equiv_refs(refs(0L:recurse_count - 1)) = min_label
		
		; DEBUG
		;if equiv_refs(0) ne 0 then print, "FAILURE (case E2) ", label
		
		; is this a label that we have seen?
		if contig_refs(min_label) eq 0 then begin
			; no! increase the new label count
			new_label_count += 1
			new_label = new_label_count
		endif else begin
			; yes! we know what to do
			new_label = contig_refs(min_label)
		end
		
		; DEBUG
		;if new_label eq 0 then print, "FAILURE (case D)"
		; a condition was contravened! ...or maybe something much worse :-(
	endelse
	
	; update the reference list
	contig_refs(label) = new_label
	
endfor 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PASS 2 - merge labels and use the smallest label value for any given region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;t5 = systime(1) ; DEBUG

; turn the connected_labellings array of labels into a 1D array
; then use that as look-up indices into the equiv_refs
; translation vector and finally reshape back into the matrix
connected_labellings = reform( contig_refs(reform(connected_labellings, n_elements(connected_labellings), 1, /overwrite)), N, M, /overwrite )

;t6 = systime(1) ; DEBUG

; DEBUG
;print, "-----"
;print, "checks: ", t2 - t1, " sec"
;print, "mask: ", t3 - t2, " sec"
;print, "initial labels: ", t4 - t3, " sec (debug count = ", debug_count, ")" 
;print, "minimization: ", t5 - t4, " sec"
;print, "final labels: ", t6 - t5, " sec"

; pass out the ccl array
return, connected_labellings  
end ; function 
