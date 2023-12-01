#!/usr/bin/env python
# coding: utf-8

# In[1]:


#This code is the pipeline for turning .sav files of the PCA results into .fits files with the appropriate information.
#Also used for generating the central PSB, SF and SB samples.


# In[114]:


get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import marvin
import re
from scipy import interpolate
#from marvin.tools.image import getImagesByList
from marvin.tools.image import Image
from marvin.utils.general.images import showImage
from marvin.tools.maps import Maps
from marvin.tools.quantities.map import Map
from marvin.tools.cube import Cube
from marvin.tools.modelcube import ModelCube
from marvin.utils.dap.bpt import get_masked
from marvin.utils.dap.bpt import get_snr
import marvin.utils.plot.map as mapplot
from scipy.io.idl import readsav
from scipy import stats
import pdb
from operator import add
from matplotlib import colors
import pandas as pd
import math
from MaNGA_Utils import spx_map
from VW_PCA import VW_PCA
import time
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants as c
import pickle
from astropy.convolution import convolve, Box1DKernel
import os
from matplotlib.colors import LogNorm
from astropy.table import Table
import matplotlib.patches as patches
from matplotlib.path import Path
from pathlib import Path as filepath
from subprocess import call
from random import randint

vwPca = VW_PCA()

#from sdss_access import RsyncAccess
#rsync = RsyncAccess()
#rsync.remote()

from marvin import config

release = 'DR17'

# by default the mode is set to 'auto', but let's set it explicitly to local/remote.
#config.mode = 'auto'
#config.download = True
##config.access = 'collab' #'public'
config.setRelease(release)

# login to receive a token
#config.login(refresh=True)

# see token
#config.token
#print(config.access, config.release)

if (release=='DR15'):
    release='MPL-7'
print(release)

dir = '/home/idies/workspace/Storage/KateRowlands/persistent/'
fig_dir = dir+'FIGS/'+release+'/'
#PCA_results_DR17_IDL

label_size = 16
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 14

# In[115]:


#Open cat
if (release=='MPL-7'):
    drpallfile = 'drpall-v2_4_3_gals.fits'
if (release=='MPL-8'):
    drpallfile = 'drpall-v2_5_3.fits'
if (release=='DR17'):
    drpallfile = '/home/idies/workspace/sdss_sas/dr17/manga/spectro/redux/v3_1_1/drpall-v3_1_1.fits'

cat = fits.open(drpallfile)
tbdata = cat[1].data
cat.close()

#Remove bad IFUs with critical flags, potential bad things. Missing DAp maps will be skipped in the loop.
ind_good = np.where((tbdata.field('DRP3QUAL') < 20000) & (tbdata.field('nsa_z') > 0.) & \
                    (tbdata.field('nsa_elpetro_th50_r') > 0.) )

#DAPQUAL also check this https://trac.sdss.org/wiki/MANGA/TRM/TRM_ActiveDev/DAPMetaData#Maskbits

mangaid = np.array(tbdata.field('mangaid')[ind_good])
nsa_id = np.array(tbdata.field('NSA_NSAID')[ind_good])
plateifu = tbdata.field('plateifu')[ind_good]
plate = tbdata.field('plate')[ind_good]
ifu = tbdata.field('ifudsgn')[ind_good]
z = tbdata.field('nsa_z')[ind_good]
ra = tbdata.field('objra')[ind_good]
dec = tbdata.field('objdec')[ind_good]
ra_plate = tbdata.field('ifura')[ind_good]
dec_plate = tbdata.field('ifudec')[ind_good]
petro_r_50 = tbdata.field('nsa_elpetro_th50_r')[ind_good]
nsa_elpetro_ba = tbdata.field('nsa_elpetro_ba')[ind_good]
nsa_elpetro_phi = tbdata.field('nsa_elpetro_phi')[ind_good]
nsa_elpetro_mass = tbdata.field('nsa_elpetro_mass')[ind_good]
nsa_elpetro_mass_log = np.log10(nsa_elpetro_mass)
ebvgal = tbdata.field('ebvgal')[ind_good]
nsa_sersic_n = tbdata.field('nsa_sersic_n')[ind_good]
nsa_elpetro_absmag = tbdata.field('nsa_elpetro_absmag')[ind_good]
nsa_elpetro_absmag_r = nsa_elpetro_absmag[4]

cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
distance = cosmo.luminosity_distance(z) #Mpc
angular_diameter_distance = cosmo.arcsec_per_kpc_comoving(z) #arcsec / kpc
arc2pc = 1000./angular_diameter_distance.value #area of each spaxel in pc

spaxel_size = 0.5 # "/pix pixel_scale
spaxel_area = (spaxel_size / np.array(angular_diameter_distance))**2 #kpc^2
ngal = len(plateifu)
print(ngal)

MaNGA_HYB10_dir = '/home/idies/workspace/sdss_sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/'


# In[116]:


petro_frac_out = np.array([2.])

ellip_ap_radius_as_out_major = petro_r_50 / petro_frac_out  #arcsec

ellip_ap_radius_pix_out_major = ellip_ap_radius_as_out_major / spaxel_size #arcsec -> pix
ellip_ap_radius_pix_out_minor = nsa_elpetro_ba * ellip_ap_radius_pix_out_major #Axis ratio b/a

theta_rad = np.radians(nsa_elpetro_phi-90) #radians - use for getAperture in pixel coords
theta = nsa_elpetro_phi #deg - use for getAperture in sky coords
#The rotation angle in radians of the semimajor axis from the
#positive x axis. The rotation angle increases counterclockwise.
#nsa_elpetro_phi is Angle (E of N) - needs -90 to make the right way up.


# In[117]:


vertices_lowmass = {
"psb_cut":0.24, \
"psb_cut2":-0.4, \
"psb_cut3":-0.1, \
"psb_peak_x":-4.7, \
"sb_cut":0.15, \
"sb_vert1":-2.3, \
"sf_vert1":-6.3, \
"sf_vert2":-5.5, \
"sf_vert3":-2.6, \
"sf_vert4":-3.2, \
"green_vert1": -1.0, \
"green_vert2": -2.1, \
"junk_y_lower": -1.4, \
"junk_y_lower2": -3., \
"junk_y_upper": 2., \
"left_cut": -7.1, \
"right_cut": 2., \
}

vertices_highmass = {
"psb_cut":0.24, \
"psb_cut2":0.18, \
"psb_cut3":0.2, \
"psb_peak_x":-4.7, \
"sb_cut":0.0, \
"sb_vert1":-1.9, \
"sf_vert1":-5.5, \
"sf_vert2":-5.2, \
"sf_vert3":-2.2-0.07, \
"sf_vert4":-3.2-0.07, \
"green_vert1": -0.1, \
"green_vert2": -2.0-0.07, \
"junk_y_lower": -1.2, \
"junk_y_lower2": -3., \
"junk_y_upper": 2., \
"left_cut": -7.1-0.07, \
"right_cut": 2., \
}

snr_min = 3.

tmp = []

total_spaxels_ell = np.zeros((ngal))
red_spaxels_ell = np.zeros((ngal))
blue_spaxels_ell = np.zeros((ngal))
sb_spaxels_ell = np.zeros((ngal))
green_spaxels_ell = np.zeros((ngal))
psb_spaxels_ell = np.zeros((ngal))
psb_sb_spaxels_ell = np.zeros((ngal))
masked_spaxels_ell = np.zeros((ngal))
total_spaxels_ell_all = np.zeros((ngal))
#junk_spaxels_ell = np.zeros((ngal))

red_frac_ell = np.zeros((ngal))
blue_frac_ell = np.zeros((ngal))
sb_frac_ell = np.zeros((ngal))
green_frac_ell = np.zeros((ngal))
psb_frac_ell = np.zeros((ngal))
psb_sb_frac_ell = np.zeros((ngal))
masked_frac_ell = np.zeros((ngal))
#junk_frac_ell = np.zeros((ngal))

area_annulus_ell = np.zeros((ngal)) #In spaxels

for i in range(ngal):

    name = MaNGA_HYB10_dir+str(plate[i]).ljust(4,'0')+'/'+str(ifu[i])+'/manga-'+plateifu[i]+\
    '-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz'

#    try:
    print(i, plateifu[i])
    if (os.path.isfile(name) == True):
        maps = Maps(filename=name, bintype='HYB10')
    elif (os.path.isfile(name) == False):
        print(name, ' missing')
        continue

    #Read in PCA results
    s = readsav(dir+'PCA_results_'+release+'_IDL/pca_results_PSB_MaNGA_'+str(plateifu[i])+'.sav')

    pc1 = s.pcs[:,0]
    pc2 = s.pcs[:,1]
    pc3 = s.pcs[:,2]

    pc1err = s.pcerr[:,0].copy()
    pc2err = s.pcerr[:,1].copy()
    pc3err = s.pcerr[:,2].copy()

    sidelen = int(np.sqrt(np.shape(pc2)))

    #SNR array
    snr_4000A = np.reshape(s.snr_4000a, (sidelen,sidelen))

    #norm array
    norm = (np.reshape(s.norm_store, (sidelen,sidelen)))

    binid = maps.binid_binned_spectra.value
    spx_bin_mask = spx_map(binid)

    pc1_masked = pc1.copy()
    pc2_masked = pc2.copy()
    pc3_masked = pc3.copy()

    pc1_masked[np.isnan(pc1)] = -99.
    pc2_masked[np.isnan(pc2)] = -99.
    pc3_masked[np.isnan(pc3)] = -99.

    #Flag zeros where PCA failed
    pc1_masked[np.add(pc1, pc2) == 0.] = -99.
    pc2_masked[np.add(pc1, pc2) == 0.] = -99.
    pc3_masked[np.add(pc1, pc2) == 0.] = -99.

    #Get things in the right shape
    pc1_map_reshaped = (np.reshape(pc1_masked, (sidelen,sidelen)))
    pc2_map_reshaped = (np.reshape(pc2_masked, (sidelen,sidelen)))
    pc3_map_reshaped = (np.reshape(pc3_masked, (sidelen,sidelen)))

    #pcerr
    pc1err[np.isnan(pc1)] = -99.
    pc2err[np.isnan(pc2)] = -99.
    pc3err[np.isnan(pc3)] = -99.

    #Flag zeros
    pc1err[np.add(pc1, pc2) == 0.] = -99.
    pc2err[np.add(pc1, pc2) == 0.] = -99.
    pc3err[np.add(pc1, pc2) == 0.] = -99.

    #Get things in the right shape
    pc1err_map_reshaped = (np.reshape(pc1err, (sidelen,sidelen)))
    pc2err_map_reshaped = (np.reshape(pc2err, (sidelen,sidelen)))
    pc3err_map_reshaped = (np.reshape(pc3err, (sidelen,sidelen)))

    #PCA classification of spaxels. points needs to be a (ngal, 2) array
    points = np.stack((pc1_masked, -pc2_masked), axis=1) #Flip PC2

    #Classification boundaries depend on stellar mass
    if nsa_elpetro_mass[i] < 1E10:
        class_map = vwPca.PCA_classify(points, vertices_lowmass)
    elif nsa_elpetro_mass[i] >= 1E10:
        class_map = vwPca.PCA_classify(points, vertices_highmass)

    class_map_reshaped = np.reshape(class_map, (sidelen,sidelen))

    hb = get_masked(maps, 'hb_4862', snr=get_snr(snr_min, 'hb'))
    ha = get_masked(maps, 'ha_6564', snr=get_snr(snr_min, 'ha'))

    balmerdec_tmp = np.ma.log10((ha/hb)/2.86)
    balmerdec = balmerdec_tmp.filled([0.]) #If masked due to low SNR, fill with 0. so 10**0 = 1 (i.e. 2.86 for case B recombination).
    #Otherwise class map is masked where low SNR in ha or hb which might bias the results.

    #Quality mask - see https://www.sdss.org/dr15/algorithms/bitmasks/#MANGA_DRP3QUAL
    ha_qual = maps['emline_gflux_ha_6564']
    nocov = ha_qual.pixmask.get_mask('NOCOV')
    lowcov = ha_qual.pixmask.get_mask('LOWCOV')
    donotuse = ha_qual.pixmask.get_mask('DONOTUSE')
    deadfiber = ha_qual.pixmask.get_mask('DEADFIBER')
    forestar = ha_qual.pixmask.get_mask('FORESTAR')


    qualmask = np.ma.MaskedArray(pc1_map_reshaped.T, mask=(snr_4000A.T < 4.) | (pc1_map_reshaped.T < -10.) | \
                                (nocov) | (lowcov) | (donotuse) | (deadfiber) | (forestar) )
    #                             | (pc2err_map_reshaped > 0.20) | (10**balmerdec > 1.8))
    qualmask_new = np.zeros((sidelen, sidelen))
    qualmask_new[qualmask.mask.T] = 1

    mstar = np.zeros((sidelen, sidelen))
    mstar.fill(np.log10(nsa_elpetro_mass[i]))

    master = pd.DataFrame({ 'pca_store': class_map_reshaped.flatten(), \
                          'pc1_store': pc1_map_reshaped.flatten(), \
                          'pc2_store': pc2_map_reshaped.flatten(), \
                          'pc1err_store': pc1err_map_reshaped.flatten(), \
                          'pc2err_store': pc2err_map_reshaped.flatten(), \
                          'balmerdec_store': balmerdec.flatten(), \
                          'spxmap_store': spx_bin_mask.flatten(), \
                          'mstar_store': mstar.flatten(), \
                          'snr_4000A_store': snr_4000A.flatten(), \
                          'qualmask_store': qualmask_new.flatten() \
                          })
    tmp.append(master)


    primhdr = maps.header

    hdr_new = fits.HDUList([fits.PrimaryHDU(),
                   fits.ImageHDU(data=pc1_map_reshaped.T, name='PC1'),
                   fits.ImageHDU(data=-pc2_map_reshaped.T, name='PC2'), #Flip PC2
                   fits.ImageHDU(data=pc3_map_reshaped.T, name='PC3'),
                   fits.ImageHDU(data=pc1err_map_reshaped.T, name='PC1ERR'),
                   fits.ImageHDU(data=pc2err_map_reshaped.T, name='PC2ERR'),
                   fits.ImageHDU(data=pc3err_map_reshaped.T, name='PC3ERR'),
                   fits.ImageHDU(data=qualmask_new.T, name='qualmask'),
                   fits.ImageHDU(data=snr_4000A.T, name='snr4000A'),
                   fits.ImageHDU(data=norm.T, name='norm'),
                   fits.ImageHDU(data=class_map_reshaped.T, name='classmap'),
                   fits.ImageHDU(data=spx_bin_mask.T, name='spx_bin_mask')])

    hdr_new.writeto(dir+'PCA_results_'+release+'_fits/manga-'+plateifu[i]+'_PCA.fits.gz', overwrite=True)
    hdr_new.close()


    ######### Elliptical aperture measurements ############

    #Apply SNR and pixel-spaxel mask so binned spaxels don't get double counted.
    #class_map_masked = np.ma.MaskedArray(class_map_reshaped.T, mask=spx_bin_mask.T < 1)
    #I think that we don't want to apply the spaxel mask, we want to calculate the surface area fraction of spaxels
    class_map_masked = np.ma.MaskedArray(class_map_reshaped.T)

    #Quality mask.
    class_map_snr_masked_balmerdec = np.ma.MaskedArray(class_map_masked, mask=( (qualmask_new.T==1) ))

    aperture = maps.getAperture([(ra[i], dec[i])], (ellip_ap_radius_as_out_major[i], \
                        nsa_elpetro_ba[i]*ellip_ap_radius_as_out_major[i], theta[i]), \
                        coord_type='sky', aperture_type='elliptical')

    ap = aperture.to_pixel(maps.wcs)
    mask_ellipical = ap.to_mask()[0]

     #Because mask.cutout is a square and we need to account for the circle edges by multiplying by mask.data
    data_cutout_balmerdec_tmp = mask_ellipical.cutout(class_map_snr_masked_balmerdec)
    data_cutout_balmerdec = np.ma.MaskedArray(data_cutout_balmerdec_tmp, mask=data_cutout_balmerdec_tmp <= 0.)

    masked_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec.mask==True)] )
    total_spaxels_ell_all[i] = np.ma.sum( mask_ellipical.data ) #All spaxels in aperture, masked or unmasked

    if (data_cutout_balmerdec.mask.all() == False):
        total_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec > 0)] )
        red_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec == 1)] )
        blue_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec == 2)] )
        sb_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec == 3)] )
        green_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec == 4)] )
        psb_spaxels_ell[i] = np.ma.sum( mask_ellipical.data[np.ma.where(data_cutout_balmerdec == 5)] )

        #Checks
        sum1 = total_spaxels_ell[i]+masked_spaxels_ell[i]
        sum2 = red_spaxels_ell[i]+blue_spaxels_ell[i]+sb_spaxels_ell[i]+green_spaxels_ell[i]+psb_spaxels_ell[i]

        sumtest1 = np.isclose(sum1, total_spaxels_ell_all[i], rtol=0.05)
        sumtest2 = np.isclose(sum2, total_spaxels_ell[i], rtol=0.05)

        assert (sumtest1==True), "spaxels don't add up 1"
        assert (sumtest2==True), "spaxels don't add up 2"

    #Check if meet min spaxel number, otherwise array value will stay empty.
    if (float(total_spaxels_ell[i]) > 3):
        if (float(red_spaxels_ell[i]) > 0):
            red_frac_ell[i] = float(red_spaxels_ell[i]) / float(total_spaxels_ell[i])

        if (float(blue_spaxels_ell[i]) > 0):
            blue_frac_ell[i] = float(blue_spaxels_ell[i]) / float(total_spaxels_ell[i])

        if (float(sb_spaxels_ell[i]) > 0):
            sb_frac_ell[i] = float(sb_spaxels_ell[i]) / float(total_spaxels_ell[i])

        if (float(green_spaxels_ell[i]) > 0):
            green_frac_ell[i] = float(green_spaxels_ell[i]) / float(total_spaxels_ell[i])

        if (float(psb_spaxels_ell[i]) > 0):
            psb_frac_ell[i] = float(psb_spaxels_ell[i]) / float(total_spaxels_ell[i])

        if (float(psb_sb_spaxels_ell[i]) > 0):
            psb_sb_frac_ell[i] = float(psb_sb_spaxels_ell[i]) / float(total_spaxels_ell[i])

    if (float(masked_spaxels_ell[i]) > 0):
        masked_frac_ell[i] = float(masked_spaxels_ell[i]) / float(total_spaxels_ell_all[i])

    elif (float(total_spaxels_ell[i]) <= 3):
        total_spaxels_ell[i] = 0.

    elif (data_cutout_balmerdec.mask.all() == True):
        total_spaxels_ell[i] = 0.
        masked_frac_ell[i] = 1.0

spaxel_properties_master = pd.concat(tmp)

f = open(dir+'OUTPUT/spaxel_properties_master_'+release+'.p', 'wb')
pickle.dump(spaxel_properties_master, f)          # dump data to f
f.close()

#Output
f2 = open(dir+'OUTPUT/elliptical_radii_params_'+release+'.p', 'wb')
pickle.dump([total_spaxels_ell, red_spaxels_ell, blue_spaxels_ell, sb_spaxels_ell, green_spaxels_ell, \
            psb_spaxels_ell, red_frac_ell, blue_frac_ell, \
            sb_frac_ell, green_frac_ell, psb_frac_ell, total_spaxels_ell_all, masked_spaxels_ell, masked_frac_ell], \
             f2)  # dump data to f
f2.close()
