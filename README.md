# MaNGA-PCA
This respository contains code for running the PCA on MaNGA data. The IDL files contains the PCA code, and the Python directory contains post-processing and analysis codes, that runs on the output of the IDL code. For details of the IDL code see https://github.com/SEDMORPH/VWPCA.

MaNGA spaxels are classified into five spectral classes. Following the method outlined in Wild et al. (2007), we define two spectral indices which are based on a Principal Component Analysis (PCA) of the 3750-4150A region of the spectra. PC1 is related to the strength of the 4000A break (equivalent to the D_{n}4000 index), and PC2 is the excess Balmer absorption (of all Balmer lines simultaneously) over that expected for the 4000A break strength. The eigenbasis that defines the principal components is taken from Wild et al. (2007), and was built using Bruzual and Charlot (2003) model spectra. More details are given in Rowlands et al. (2018) and Boardman et al. (2023).

The plotting routines depend on the Marvin code: https://github.com/sdss/marvin

run_make_PCA_maps.py is the code used to turn the IDL output files into fits files on Sciserver (essentially a copy of what is in make_PCA_maps.py)	

analyse_PCA_maps.ipynb does some useful plotting of the PCA maps, and PCA space.

The classification into the five classes of quiescent, starforming, starburst, green valley and post-starburst are based on the PC1, PC2 values. The classification boundaries are slightly different depending on whether a galaxy has M*<10^10 or M*>10^10 Msun.

MaNGA_Utils.py: Class containing the code which makes the spx_bin_mask to account for the double counting of spaxels which have the same values because they are in the same Voronoi bin.

VW_PCA.py: Class containing the classification code, and other useful functions such as deriving SFR from the Halpha flux, Galaxtic extinction correction, sky line masking and pre-processesing of spectra required before performing the PCA.

Kate Rowlands
krowlands [at] stsci.edu
