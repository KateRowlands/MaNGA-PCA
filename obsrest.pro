;; OBSREST

;; convert observed frame wavelength array into rest frame
;; Optional OUTPUTS: shift flux and error to per-restframe-AA

;; vw: 2008-01-24 updated to more sensible method

FUNCTION OBSREST, wave, z,flux=flux,error=error

wave_rest = wave - (round(alog10(1.+z)*10000.))/10000.

if n_elements(flux) ne 0 then flux = flux*10^wave/10^wave_rest

if n_elements(error) ne 0 then error = error*10^wave/10^wave_rest

return, wave_rest

END
