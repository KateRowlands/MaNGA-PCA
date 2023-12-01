;*** MASKSKY
; masks sky emission lines [OI], N2+
; vw: updated Sep 14th to widen 5578 line
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function masksky, wave,flux,pixel = pixel,silent=silent

;5578.5486       4279.2039       6301.7423       6365.7595
;emlines = $
;[[5574,5582],$
;[4276,4282],$
;[6297,6305],$
;[6364,6368]]

;; update 14/09
emlines = $
[[5574,5590],$
[4276,4282],$
[6297,6305],$
[6364,6368]]



nlines = n_elements(emlines)/2

n = 0
index2 = lonarr(1000)
mask_arr = fltarr(1000)
for i=0, nlines-1 do begin

    index = where(wave gt emlines[0,i] and wave lt emlines[1,i],count)

    if count gt 0 then begin
       ind_mask = where((wave gt emlines[0,i]-75 and wave lt emlines[0,i]) or (wave gt emlines[1,i] and wave lt emlines[1,i]+75))


        mask_arr[n:n+count-1] = median(flux[ind_mask])


        index2[n:n+count-1] = index
        n = n+count
    endif

endfor



if n ne 0 then begin
    index2 = index2[0:n-1]

    pixel = uintarr(n_elements(wave))
    pixel[index2]=1

    mask_arr = mask_arr[0:n-1]
    newflux=flux
    newflux[index2] = mask_arr
endif else begin
    if n_elements(silent) eq 0 then print,'masksky.pro:no em lines masked'
    newflux = flux
    pixel = uintarr(n_elements(wave))
endelse

return,newflux

end
