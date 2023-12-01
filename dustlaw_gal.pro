;***DUSTLAW_GAL
; to correct given flux, with wavelengths in AA and measured g-band
; extinction, for Galactic redenning. Uses Pei 1992 Galactic redenning
; curve.
;
; UPDATE 20/04/05 now uses Cardelli dust to be consistent with SDSS
; UPDATE 22/01/2014 to accept ebv instead of ext_g and output mags if required

FUNCTION DUSTLAW_GAL,flux,wave,ext_g,corr=corr,dust=dust,lam_dust=lam_dust,ebv=e_bv,mag=mag

if NOT(KEYWORD_SET(dust)) then begin
    DIR2 = '/data/kate/REDDEN/'
    infile = 'cardelli_gal.ext'
    tmp = fltarr(4,4751)
    openr,1,DIR2+infile
    readf,1,tmp
    close,1
    lam_dust = reform(tmp[0,*])
    dust = reform(tmp[1,*])
endif

if min(wave) lt min(lam_dust) then message, 'wavelength wrong'
if max(wave) gt max(lam_dust) then message, 'wavelength wrong'

f_dust = interpol(dust,lam_dust,wave)

if n_elements(e_bv) eq 0 then E_BV=ext_g/3.793               ;E(B-V)=A_lambda/R_lambda


R = 3.1                         ;R_V
A = E_BV*(f_dust + R)           ;extinction in mags

if keyword_set(mag) then return,A

;now we know what the extinction is as function of lambda, so divide
;this out of the flux

flux_corr = flux/10^(-A/2.5)    ;corrected flux
corr = 10^(-A/2.5)              ;to correct error array if needed


return, flux_corr

END
