# MaNGA-PCA
Code to work with MaNGA PCA maps

This repository contains code related to the project to classify all MaNGA spaxels into five spectral classes. Following the method outlined in Wild et al. (2007), we define two spectral indices which are based on a Principal Component Analysis (PCA) of the 3750-4150A region of the spectra. PC1 is related to the strength of the 4000A break (equivalent to the D_{n}4000 index), and PC2 is the excess Balmer absorption (of all Balmer lines simultaneously) over that expected for the 4000A break strength. The eigenbasis that defines the principal components is taken from Wild et al. (2007), and was built using Bruzual and Charlot (2003) model spectra.

The classification into the five classes of quiescent, starforming, starburst, green valley and post-starburst are based on the PC1, PC2 values. The classification boundaries are slightly different depending on whether a galaxy has M*<10^10 or M*>10^10 Msun.

    PCA plotting example.ipynb: Notebook showing how to read the maps of PCA amplitudes, classifications and quality masks, make maps of the PCA classes, reclassify the spaxels using different boundaries in PC1, PC2, and count how many spaxels there are in each class in an aperture.
    MaNGA_Utils.py: Class containing the code which makes the spx_bin_mask to account for the double counting of spaxels which have the same values because they are in the same Voronoi bin.
    VW_PCA.py: Class containing the classification code, and other useful functions such as deriving SFR from the Halpha flux, Galaxtic extinction correction, sky line masking and pre-processesing of spectra required before performing the PCA.
