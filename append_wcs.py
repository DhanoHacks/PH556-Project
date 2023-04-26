from astropy.io import fits
from astropy.wcs import WCS
filename = "20230415202558-023-RA.science.fits"
hdulist = fits.open("wcs.fits")
w = WCS(hdulist[0].header)
hdulist2 = fits.open(filename)
header2=hdulist2[0].header
for i in w.to_header():
    header2[i]=w.to_header()[i]
hdulist2.writeto(filename.split(".")[0]+".wcs.proc.fits",overwrite=True)