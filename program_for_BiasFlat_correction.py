#credits - Vishwajeet Swain IITB
from astropy.io import fits
import numpy as np
import glob
from photutils.aperture import CircularAperture


def load_fits(image, extension=0):
    with fits.open(image) as hdu:
        header = hdu[extension].header
        data = hdu[extension].data
    return data, header


images = glob.glob("*-RA.fits")
for image in images:
    date = image[:8]
    sci_image = image
    mBias = f"masterbias-{date}.fits"
    mFlat = f"masterflat-{date}.fits"

    data, header = load_fits(sci_image)
    data = np.float32(data)
    aperture = CircularAperture((data.shape[0]//2,data.shape[1]//2),r=data.shape[0]//2)
    mask = aperture.to_mask(method='exact')
    data[mask.to_image(data.shape) == 0] = np.NaN

    mBais_data, _ = load_fits(mBias)
    mFlat_data, _ = load_fits(mFlat)
    corr_data = (data-mBais_data)/mFlat_data
    header.add_history('Bias corrected and flat-fielded')
    fits.writeto(sci_image.split(".")[0]+".science.fits", corr_data, header, overwrite=True)
