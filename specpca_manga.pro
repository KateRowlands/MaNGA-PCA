pro specpca_MaNGA

;roothome = '/uufs/chpc.utah.edu/common/home/sdss42/mangawork/users/u6048365/'
roothome = '/uufs/chpc.utah.edu/common/home/u6048365/'
rootdir = '/uufs/chpc.utah.edu/common/home/sdss42/mangawork/manga/spectro/'
dataDIR = rootdir
figDIR = roothome+'FIGS/'
outDIR = roothome+'PCA_results_DR17_IDL/'

MaNGA_cat = mrdfits(dataDIR+'redux/v3_1_1/drpall-v3_1_1.fits', 1, hdu, /silent)

;Flag critical IFUs
ind = where( (MaNGA_cat.nsa_z gt 0.), ngal ) ;and (MaNGA_cat.DRP3QUAL lt 300) )

MaNGA_cat = MaNGA_cat[ind]

ngal = n_elements(MaNGA_cat)

info = replicate({plateifu:'a',ifudsgn:0,plate:0,z:0.,EBV:0.},ngal)
info.plateifu = MaNGA_cat.plateifu
info.ifudsgn = MaNGA_cat.ifudsgn
info.plate = MaNGA_cat.plate
info.z = MaNGA_cat.nsa_z
info.EBV = MaNGA_cat.ebvgal

;Run PCA on MaNGA spaxels
pcs_manga_store = pca_projectmanga(info, nrecon=3, pcerr=pcerr, rootdir=rootdir, roothome=roothome)

outfile = outDIR+'pca_results_DR17.sav'
save, pcs_manga_store, file=outfile

END
