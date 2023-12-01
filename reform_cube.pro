PRO reform_cube, ima, xout, yout, imaout, imaout_spec=imaout_spec, xpixscale=xpixscale, ypixscale=ypixscale, xcenter=xcenter, ycenter=ycenter

; programme that takes a 2D/3D velocity field or image/cube and reformats it into long
; arrays, that can be read by Voronoi binning, pPXF codes, kinemetry, etc.

; written by Anne-Marie Weijmans, Anstruther, 1 June 2015

; ---------------------------------

; determine pixel scales, if given
xpixscale=0.2
ypixscale=0.2

xcenter=0.
ycenter=0.

; determine dimensions of cube
sz = SIZE(ima, /dim)

nx = sz[0]*1.
ny = sz[1]*1.

if (n_elements(sz) gt 2) then begin
  speclength = sz[2]
  imaout_spec = fltarr(speclength, nx*ny)
endif

; construct output arrays
xout = fltarr(nx*ny)
yout = fltarr(nx*ny)
imaout = fltarr(nx*ny)

FOR i = 0., nx-1. DO BEGIN
    FOR j = 0., ny-1. DO BEGIN
        xout[i*ny+j] = (i - xcenter)*xpixscale
        yout[i*ny+j] = (j - ycenter)*ypixscale
        imaout[i*ny+j] = ima[i,j]
        if (n_elements(sz) gt 2) then imaout_spec[*, i*ny+j] = ima[i, j, *]
    ENDFOR
ENDFOR

END
