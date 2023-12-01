FUNCTION PCA_PROJECTMaNGA,info,runno,plotspec=plotspec,nrecon=nrecon,snr=snr,rms=rms_snr,$
                         espec=espec,normval=norm,pcerr=pcs_err,cov=cov, rootdir=rootdir, roothome=roothome

c = 299792458.0 / 1e3 ;km/s

dataDIR = rootdir
figDIR = roothome+'FIGS/'
outDIR = roothome+'PCA_results_DR17_IDL/'
MaNGA_HYB10_dir = rootdir+'analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/'

;;-- restore eigenvectors
restore, roothome+'VO/PCARUN/ESPEC/pcavo_espec_25.sav'
wave_espec = wave
npix = n_elements(wave_espec)

;;-- dust correction
tmp = fltarr(4,4751)
openr,1, roothome+'REDDEN/cardelli_gal.ext'
readf,1,tmp
close,1
lam_dust = reform(tmp[0,*])
dust = reform(tmp[1,*])

restwave = [3727.092, 3934.777, 3969.588, 4102.89, $
              4341.68, 4364.436, 4862.68, 4960.295,$
              5008.240, 5176.7, 5895.6, 6302.046, $
             6549.86, 6564.61, 6585.27, 6718.29, $
             6732.672]

ngal = n_elements(info)

pcs_store = replicate( {plateifu:'a', z:0., pcs:fltarr(nrecon,16500), $
  pcserr:fltarr(nrecon,16500), norm_store:fltarr(16500), $
  snr_4000A:fltarr(16500)}, ngal )
pcs_store.plateifu = info.plateifu
pcs_store.z = info.z
;(pcerr_store = fltarr(nrecon,ngal))

; ps1, figDIR+'Plotpsbs_MaNGA_test.ps'
; device,ysize=200,xsize=15
; !p.multi=[0,3,50,0,0]

nn = 0L
lost = 0L

print, 'Projecting spectra....',ngal
time = systime(1)

for i = 0, ngal-1 do begin

  print, i, ' ', info[i].plateifu

  filename = MaNGA_HYB10_dir+strcompress(info[i].plate, /rem)+'/'+strcompress(info[i].ifudsgn, /rem)+'/manga-'+strcompress(info[i].plate, /rem)+'-'+strcompress(info[i].ifudsgn, /rem)+'-LOGCUBE-HYB10-MILESHC-MASTARSSP.fits.gz'

  file_exist = FILE_TEST(filename)

  if (file_exist eq 1) then begin

    primaryhdu = mrdfits(filename,0,hdu, /silent)
    fluxcube = mrdfits(filename,1,/silent)
    ivar = mrdfits(filename,2,/silent) ;inverse variance
    errorcube = 1./sqrt(ivar)
    ;maskcube = mrdfits(filename,3,/silent)
    wave = mrdfits(filename,5,/silent)

    ;Read velocity map
    ;print, 'Reading stellar velocity map'
    stellar_vel = mrdfits(MaNGA_HYB10_dir+strcompress(info[i].plate, /rem)+'/'+strcompress(info[i].ifudsgn, /rem)+'/manga-'+strcompress(info[i].plate, /rem)+'-'+strcompress(info[i].ifudsgn, /rem)+'-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz', 15, hdu, /silent)
    stellar_vel_ivar = mrdfits(MaNGA_HYB10_dir+strcompress(info[i].plate, /rem)+'/'+strcompress(info[i].ifudsgn, /rem)+'/manga-'+strcompress(info[i].plate, /rem)+'-'+strcompress(info[i].ifudsgn, /rem)+'-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz', 16, hdu, /silent)

    arrsize = size(stellar_vel, /dim)
    stellar_vel_central = stellar_vel[(arrsize[0]-1)/2, (arrsize[0]-1)/2]

    ;binid = mrdfits('/Users/katerowlands/mount/sdss4/manga/mpl5/HYB10-GAU-MILESHC/8155/12701/manga-8155-12701-MAPS-HYB10-GAU-MILESHC.fits.gz',6,hdu,/silent)

    ;Project spaxel position into a 1D array. values of xout are wrong (pixscale wrong input) but don't need to use these.
    REFORM_CUBE, fluxcube, xout, yout, imaout, imaout_spec=flux_spec
    REFORM_CUBE, errorcube, xouterr, youterr, errout, imaout_spec=error_spec
    ;REFORM_CUBE, binid, xoutbinid, youtbinid, binidout
    REFORM_CUBE, stellar_vel-stellar_vel_central, xoutstellar_vel, youtstellar_vel, stellar_velout ;correct for non-zero central velocity
    REFORM_CUBE, stellar_vel_ivar, xoutstellar_vel_ivar, youtstellar_vel_ivar, stellar_velout_ivar

    stellar_vel_err = 1./sqrt(stellar_vel_ivar)

    ;1 PCA results file per cube.
    pcs = (pcerr = fltarr(nrecon,n_elements(xout)))
    snr_4000A = fltarr(n_elements(xout))
    norm_store = fltarr(n_elements(xout))

    ;Run on Voronoi binned LOGCUBE (modelcube in Marvin language), each bin S/N~10 in r-band.
    for k = 0, n_elements(xout)-1 do begin

      flux = flux_spec[*,k] ;Shape of spectrum length, 74*74 i.e. unravelled cube position
      error = error_spec[*,k]

      ;Need to take out galaxy rotation
      z_map = (1. + info[i].z) * (1. + (stellar_velout[k] / c)) - 1. ;Need to subtract central velocity first

      wave_rest = obsrest(alog10(wave), z_map, flux=flux, error=error)

      if ( (total(flux) ne 0.) or (stellar_vel_err[k] le 500.) ) then begin
      flux = dustlaw_gal(flux, wave, ebv=info[i].ebv, corr=corr, dust=dust, lam_dust=lam_dust)
      error = error/corr

      ;Check this
      flux = masksky(wave,flux,pixel=mask)
      tmp = where(mask eq 1)
      if tmp[0] ne -1 then error[tmp]=0.

      masksize = 5.0
      ind = where(10^wave_rest gt 4102.9-masksize and 10^wave_rest lt 4102.9+masksize)
      if ind[0] ne -1 then error[ind]=0.

      ind = where(10^wave_rest gt 3971.195-masksize and 10^wave_rest lt 3971.195+masksize)
      if ind[0] ne -1 then error[ind]=0.

      ind = where(10^wave_rest gt 3890.151-masksize and 10^wave_rest lt 3890.151+masksize)
      if ind[0] ne -1 then error[ind]=0.

      masksize = 2.5

      ind = where(10^wave_rest gt 3836.472-masksize and 10^wave_rest lt 3836.472+masksize)
      if ind[0] ne -1 then error[ind]=0.

      ind = where(10^wave_rest gt 3798.976-masksize and 10^wave_rest lt 3798.976+masksize)
      if ind[0] ne -1 then error[ind]=0.

      ;plot, 10d^wave_rest, flux, yr=[-0.1, max(flux)], /xs, /ys, title=info[i].plateifu
      ;oplot, 10d^wave_rest, error,color=cgcolor('green')

      ;interpolate onto SDSS eigenbasis
      linterp, 10d^wave_rest, flux, wave_espec, newflux
      linterp, 10d^wave_rest, error^2, wave_espec, newerror_sq

      ;; bad pixel masks get screwed up by interpolation
      ind = where((finite(newflux) eq 0) or (newflux eq 0.))
      if ind[0] ne -1 then begin
         newerror_sq[ind] = 0.
         newflux[ind] = 0.
      endif
      newerror = sqrt(newerror_sq)

      ;; grow bad pixel regions:
      ind = where(newerror eq 0,nn)
      dpix = 3
      if nn gt 0 then for j=0, nn-1 do begin
         if ind[j] gt dpix and ind[j] lt npix-dpix-1 then newerror[ind[j]-dpix:ind[j]+dpix]=0.0
         if ind[j] le dpix then newerror[0:ind[j]+dpix]=0.0
         if ind[j] ge npix-dpix-1 then newerror[ind[j]:npix-1]=0.0
      endfor

      ind = where(wave_rest gt alog10(min(wave_espec)) and $
                  wave_rest lt alog10(max(wave_espec)) and error ne 0, count)

      if count gt 1 then snr_4000A[k] = median(flux[ind]/error[ind]) ;SNR before interpolation

      ;if more than 20% of the pixels are bad
      if n_elements(where(newerror eq 0.))/float(npix) gt 0.2 then begin
          ;print,'excess bad pixels, continue ', filename, info[i].z
          lost = lost+1
          continue
      endif

      ;;-- calculate PCS

           pcs[*,k] = vwpca_normgappy(newflux, newerror, espec[*,0:nrecon-1], meanarr, norm=norm, cov=cov)
           for p=0,2 do pcerr[p,k]=sqrt(cov[p,p])
           norm_store[k] = norm

            ; if ( (pcerr[1,k] lt 0.15) ) then begin
            ;   plot, 10d^wave_rest, flux, xr=[3750,4150], /xs, title=string(k)+string(pcs[*,k]*[1,-1,1],form='(3(F0.2,1X))')
            ;   oplot, 10d^wave_rest, error,color=cgcolor('green')
            ;   oplot, wave_espec, (vwpca_reconstruct(pcs[*,k],espec[*,0:nrecon-1])+meanarr)*norm[0],color=cgcolor('red')
            ;
            ;   for p = 0, n_elements(restwave)-1 do oplot, [(restwave[p]),(restwave[p])], $
            ;     [MIN(!Y.crange),MAX(!Y.crange)], linestyle=1, color=cgcolor('BLUE')
            ;   ;xyouts, restwave, 4.1, name, /data, alignment = 0.7, charthick=3, charsize=0.8
            ; endif

         endif ;Only do PCA on non-zero arrays
      endfor ;Loop over spaxels

      name = info[i].plateifu
      outfile = outDIR+'pca_results_PSB_MaNGA_'+strcompress(name, /rem)+'.sav'

      save, name, pcs, pcerr, snr_4000A, norm_store, file=outfile ;Save one file per galaxy

  pcs_store[i].pcs[*,0:k-1] = pcs
  pcs_store[i].pcserr[*,0:k-1] = pcerr
  pcs_store[i].snr_4000A[0:k-1] = snr_4000A
  pcs_store[i].norm_store[0:k-1] = norm_store

endif ;file test

endfor ;Loop over ngal

;
;------------------------------------------------------------------
print, 'time for projection',systime(1)-time
print, 'number lost ',lost,'out of ',ngal, lost/float(ngal)

return, pcs_store

END
