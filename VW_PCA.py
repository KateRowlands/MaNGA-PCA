#!/usr/bin/env python
# encoding: utf-8
#
# VW_PCA.py
#

import numpy as np
import math
from scipy.io.idl import readsav
from scipy import interpolate
from astropy.io import ascii
import matplotlib.patches as patches
from matplotlib.path import Path
from marvin.utils.dap.bpt import get_masked
from marvin.utils.dap.bpt import get_snr

dir = '/data/kate/' #'/home/ksearle2/' #'/data/kate/'

class VW_PCA():

    def SFR(self, maps, distance):
        #calculate log10(SFR) from the Halpha flux, following Kennicutt et al. (1998).
        # Assumes a Salpeter IMF.

        snr = 3.

        hb = get_masked(maps, 'hb_4862', snr=get_snr(snr, 'hb'))
        ha = get_masked(maps, 'ha_6564', snr=get_snr(snr, 'ha'))

        # Av derivation from A_Ha following Jorge from CALIFA paper (Catalán-Torrecilla et al. (2015)
        K_Ha = 2.53
        K_Hb = 3.61

        balmerdec = np.log10((ha/hb)/2.87)
        A_Ha = balmerdec * K_Ha/(-0.4*(K_Ha-K_Hb))
        Ha_cor = ha*(10.**(0.4*A_Ha))

        Lha = Ha_cor*1e-17*4.*math.pi*((3.0857e18*1E6*distance)**2)
        K98 = 7.9*1e-42
        SFR = np.log10(K98*Lha)

        return SFR

    def dustlaw_gal(self, wave, flux, ebvgal):
        #Galactic dust correction. Based on dustlaw_gal.pro by Vivienne Wild.

        data = ascii.read(dir+'REDDEN/cardelli_gal.ext')

        lam_dust = data['col1']
        dust = data['col2']

        if np.amin(wave) < np.amin(lam_dust):
            print('wavelength wrong')
        if np.amax(wave) > np.amax(lam_dust):
            print('wavelength wrong')

        f_dust = np.interp(wave, lam_dust, dust) #interpolate

        R = 3.1                         #R_v
        A = ebvgal*(f_dust + R)           #extinction in mags

        #now we know what the extinction is as function of lambda, so divide
        #this out of the flux
        flux_corr = flux/10**(-A/2.5)    #corrected flux
        corr = 10**(-A/2.5)              #to correct error array if needed

        return corr

    def dustlaw_interp(self):
        #Emission line dust correction. Based on dust_balmerdec.pro by Vivienne Wild.
        #Balmerdec is linear (not log10).

        wave_ha = 6564.610
        wave_hb = 4862.683

        R_v = 3.1  #R_v is not specified in this fit, have to state it

        data = ascii.read(dir+'REDDEN/cardelli_gal.ext')

        lam_dust = data['col1']
        dust = data['col2']

        ll = lam_dust/1.e4         #microns
        x = 1./ll                 #wave number
        y = x-1.82

        ax = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 \
            -1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8

        bx = 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 \
            +11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8

        dust = np.zeros(len(x))
        ind = (x > 1.1) & (x < 3.3)
        dust[ind] = R_v * ( (ax[ind] + bx[ind]/R_v) -1. )

        ind = (x > 0.3) & (x < 1.1)
        dust[ind] = R_v * ( ((0.574*(x[ind])**1.61)+(-0.527*(x[ind])**1.61)/R_v)-1. )

        f = interpolate.interp1d(lam_dust, dust)
        dust_ha = f(wave_ha)
        dust_hb = f(wave_hb)

        return dust_ha, dust_hb, lam_dust, dust


    def dust_balmerdec(self, dust_ha, dust_hb, lam_dust, dust, balmerdec, a_lambda):
        #Emission line dust correction. Based on dust_balmerdec.pro by Vivienne Wild.
        #Balmerdec is linear (not log10).

        ratio_uncorr = balmerdec * 2.87
        ratio_corr = 2.87

        R_v = 3.1

        E_BV = -2.5*np.log10(ratio_uncorr.flatten()/ratio_corr)/(dust_ha-dust_hb)

        a_flux_tmp = np.empty(shape=[len(a_lambda), len(E_BV)] )

        for p in range(len(E_BV)):
            arr = 1./(10**(-E_BV[p]*(dust+R_v)/2.5))

            f2 = interpolate.interp1d(lam_dust, arr)
            a_flux_tmp[:,p] = f2(a_lambda)

        a_flux = a_flux_tmp.reshape((len(a_lambda),ratio_uncorr.shape[0],ratio_uncorr.shape[1]))

        return a_flux


    def masksky(self, wave, flux):
        # masks sky emission lines [OI], N2+. Based on masksky.pro witten by Vivienne Wild

        emlines = [[5574.,5590.],[4276.,4282.],[6297.,6305.],[6364.,6368.]]
        nlines = np.shape(emlines)[0]

        n = 0
        index2 = [0] * 1000
        mask_arr = np.zeros(1000)

        for i in range(nlines):

            index = np.where((wave > emlines[i][0]) & (wave < emlines[i][1]))
            count = len(wave[index])

            if count > 0:
                ind_mask = np.where(((wave > emlines[i][0]-75) & (wave < emlines[i][0])) | \
                                   ((wave > emlines[i][1]) & (wave < emlines[i][1]+75)))
                mask_arr[n:n+count] = np.median(flux[ind_mask])
                index2[n:n+count] = np.squeeze(index)
                n = n + count

        if n != 0:
            index2_new = np.array(index2[0:n])
            index2_new.astype(int)

            pixel = np.zeros(len(wave))
            pixel.astype(int)
            pixel[index2_new] = 1

            newflux=flux
            newflux[index2_new] = mask_arr[0:n]

        if n == 0:
            if n_elements(silent) == 0:
                print('masksky: no em lines masked')
                newflux = flux
                pixel = np.zeros(len(wave))

        return newflux, pixel

    def pca_prepro(self, waverest, z, flux, error, ewave, ebvgal, masksky=False, mask_emlines=False,):
        #Does preprocessing of MaNGA spectra needed before doing PCA, e.g.
        #Dust correction, skyline masking and emission line masking.

        waveobs = waverest*(1.+z)

        #Dust correct
        if ebvgal>0.:
            corr = self.dustlaw_gal(waveobs, flux, ebvgal)
        elif ebvgal<=0.:
            corr = 1.
            print('no dust correction')

        fluxdustcorr = flux/corr
        errordustcorr = error/corr

        if masksky==True:
            fluxdustcorr_new, pixel = self.masksky(waveobs, fluxdustcorr)
            skymask = np.where(pixel == 1)
            errordustcorr[skymask] = 0.
        elif masksky==False:
            fluxdustcorr_new = fluxdustcorr
            print('no sky line masking')

        if mask_emlines==True:
            #Masking emission lines
            print('Masking emission lines')
            masksize = 5.0
            indmask = np.where((waverest > 4102.92-masksize) & (waverest < 4102.92+masksize)) #Hdelta
            errordustcorr[indmask]=0.
            indmask = np.where((waverest > 3971.195-masksize) & (waverest < 3971.195+masksize)) #Hepsilon
            errordustcorr[indmask]=0.
            indmask = np.where((waverest > 3890.151-masksize) & (waverest < 3890.151+masksize)) #H8
            errordustcorr[indmask]=0.

            masksize = 2.5
            indmask = np.where((waverest > 3836.472-masksize) & (waverest < 3836.472+masksize)) #H9
            errordustcorr[indmask]=0.
            indmask = np.where((waverest > 3798.976-masksize) & (waverest < 3798.976+masksize)) #H10
            errordustcorr[indmask]=0.

            ind = np.where((waverest > 3869-5.) & (waverest < 3869+5.)) #[NeIII] - common in AGN
            errordustcorr[ind]=0.

        #Not sure this is needed
        #elif mask_emlines==False:

        #    indbad = np.where(np.isnan(fluxdustcorr_new) == True)
        #    fluxdustcorr_new[indbad] = 0.
        #    errordustcorr[indbad] = 0.

        #Interpolate onto eigenbasis
        newflux = np.interp(ewave, waverest, fluxdustcorr_new)
        newerror_sq = np.interp(ewave, waverest, errordustcorr**2)

        #bad pixel masks get screwed up by interpolation
        indmask = np.where((np.isnan(newflux) == True) | (newflux == 0.))
        newerror_sq[indmask] = 0.
        newflux[indmask] = 0.
        newerror = np.sqrt(newerror_sq)

        ind = np.where(error==0)[0]
        dpix = 3
        npix = len(ewave)
        if (len(ind) > 0):
            for j in range(len(ind)):
                if ((ind[j] > dpix) & (ind[j] < npix-dpix-1)):
                    error[ind[j]-dpix:ind[j]+dpix] = 0.0
                if (ind[j] <= dpix):
                    error[0:ind[j]+dpix] = 0.0
                if (ind[j] >= npix-dpix-1):
                    error[ind[j]:npix-1] = 0.0

        return newflux, newerror

    def PCA_classify(self, points, vertices):
        #Classify spaxels into different PCA classes.

        # SF
        sf_verts = [
            (vertices['sf_vert1'], vertices['junk_y_lower']), # left, bottom
            (vertices['sf_vert2'], vertices['sb_cut']), # left, top
            (vertices['psb_peak_x'], vertices['psb_cut']), # center, top
            (vertices['sf_vert3'], vertices['psb_cut3']), # top, right
            (vertices['sf_vert4'], vertices['junk_y_lower']) # right, bottom
            ]

        sf_path = Path(sf_verts)
        SF = sf_path.contains_points(points)

        # PSB
        psb_verts = [
            (vertices['psb_peak_x'], vertices['psb_cut']), # center, top
            (vertices['sf_vert2'], vertices['sb_cut']), # right, top
            (vertices['left_cut'], vertices['sb_cut']), # left, bottom
            (vertices['left_cut'], vertices['junk_y_upper']), # left, top
            (vertices['green_vert1'], vertices['junk_y_upper']), # right, top
            (vertices['green_vert1'], vertices['psb_cut2']) # right, bottom
            ]

        psb_path = Path(psb_verts)
        PSB = psb_path.contains_points(points)

        # SB
        sb_verts = [
            (vertices['left_cut'], vertices['sb_vert1']), # left, bottom
            (vertices['left_cut'], vertices['sb_cut']), # left, top
            (vertices['sf_vert2'], vertices['sb_cut']), # right, top
            (vertices['sf_vert1'], vertices['junk_y_lower']), # right, bottom
            (vertices['sf_vert2'], vertices['junk_y_lower']), # right, bottom
            (vertices['sf_vert2'], vertices['sb_vert1']) # right, bottom
            ]

        sb_path = Path(sb_verts)
        SB = sb_path.contains_points(points)

        # Green valley
        green_verts = [
            (vertices['sf_vert4'], vertices['junk_y_lower']), # left, bottom
            (vertices['sf_vert3'], vertices['psb_cut3']), # left, top
            (vertices['green_vert1'], vertices['psb_cut2']), # right, top
            (vertices['green_vert2'], vertices['junk_y_lower']) # right, bottom
            ]

        green_path = Path(green_verts)
        Green = green_path.contains_points(points)

        # Red
        red_verts = [
            (vertices['green_vert2'], vertices['junk_y_lower']), # left, bottom
            (vertices['green_vert1'], vertices['psb_cut2']), # left, mid
            (vertices['green_vert1'], vertices['psb_cut2']+0.5), # left, top
            (vertices['right_cut'], vertices['junk_y_upper']), # right, top
            (vertices['right_cut'], vertices['junk_y_lower']) # right, bottom
            ]

        red_path = Path(red_verts)
        Red = red_path.contains_points(points)

        # Junk
        junk_verts = [
            (vertices['sf_vert2'], vertices['junk_y_lower2']), # left, bottom
            (vertices['sf_vert2'], vertices['junk_y_lower']), # left, top
            (vertices['right_cut'], vertices['junk_y_lower']), # right, top
            (vertices['right_cut'], vertices['junk_y_lower2']) # right, bottom
            ]

        junk_path = Path(junk_verts)
        Junk = junk_path.contains_points(points)

        class_map = np.empty(len(points[:,0]))
        class_map.fill(-99)
        class_map[Red] = 1
        class_map[SF] = 2
        class_map[SB] = 3
        class_map[Green] = 4
        class_map[PSB] = 5
        class_map[Junk] = 0
        class_map[np.where(np.add(points[:,0], points[:,1]) == 0.)] = -99

        return class_map
