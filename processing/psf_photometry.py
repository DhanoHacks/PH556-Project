import numpy as np
import numpy.ma as ma
import os, re
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import subprocess
import warnings
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture
import glob

warnings.filterwarnings("ignore")

def get_table_from_ldac(filename, frame=1):
    from astropy.table import Table
    if frame>0:
        frame = frame*2
    tbl = Table.read(filename, hdu=frame)
    return tbl

def psf_photometry(idx_image, imageName):
    print(imageName)
    # Move to the data directory for our analysis
    f = fits.open(imageName)
    data = f[0].data  #This is the image array
    header = f[0].header
    band = header["filter"]
    tar_ra = header["TARRA"]
    tar_dec = header["TARDEC"]

    mask_corners = data[0][0]
    data[data == mask_corners] = np.NaN

    #Compute some image statistics for scaling the image plot
    mean, median, sigma = sigma_clipped_stats(data)

    #strong the image WCS into an object
    w = WCS(header)
    
    [center_brightsource_x, center_brightsource_y] = w.all_world2pix(293.64,19.77,1)
    aperture = CircularAperture((center_brightsource_x,center_brightsource_y),r=800)
    mask = aperture.to_mask(method='exact')
    data[mask.to_image(data.shape) == 1] = np.NaN
    if ".wcs.fits" in imageName:
        mask_img_name = re.sub(".wcs.fits",".masked.wcs.fits",imageName)
    else:
        mask_img_name = re.sub(".wcs.proc.fits",".masked.wcs.proc.fits",imageName)
    f.writeto(mask_img_name, overwrite=True)

    #Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0]/2, data.shape[1]/2, 1)

    #Set the box size to search for catalog stars
    boxsize = 30 # arcminutes

    #Magnitude cut-offs of sources to be cross-matched against
    maxmag = 20

    catNum = 'II/349'#This is the catalog number of PS1 in Vizier
    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a radius of %.4f arcmin in %s band'%(catNum, raImage, decImage, boxsize, band))


    #You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
    v = Vizier(columns=['*'], column_filters={f"{band}mag":"<%.2f"%maxmag, "Nd":">6", f"e_{band}mag":"< 0.005"}, row_limit=-1, timeout = 120)
    Q = v.query_region(SkyCoord(ra = raImage, dec = decImage, unit = (u.deg, u.deg)), radius = str(boxsize)+'m', catalog=catNum, cache=False)
    #query vizier around (ra, dec) with a radius of boxsize

    #Convert the world coordinates of these stars to image coordinates
    ps1_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)

    #Another round of filtering where we reject sources close to the edges
    good_cat_stars = Q[0][np.where((ps1_imCoords[0] > 500) & (ps1_imCoords[0] < 3500) & (ps1_imCoords[1] > 500) & (ps1_imCoords[1] < 3500))]
    ps1CatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
    print(len(ps1CatCoords))
    ps1_imCoords = w.all_world2pix(good_cat_stars['RAJ2000'],good_cat_stars['DEJ2000'], 1)

    #Sextractor commands
    configFile = 'photomCat.sex'
    catalogName = mask_img_name+'.cat'
    paramName = 'photomCat.param'
    command = 'sex -c %s %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (configFile, mask_img_name, catalogName, paramName)
    print('Executing command: %s' % command)
    rval = subprocess.run(command.split(), check=True)


    #psfex commands
    psfConfigFile = 'psfex_conf.psfex'
    command = 'psfex -c %s %s' % (psfConfigFile, catalogName)
    print('Executing command: %s' % command)
    rval = subprocess.run(command.split(), check=True)

    psfModelHDU = fits.open('moffat_'+mask_img_name+'.fits')[0]
    psfModelData = psfModelHDU.data
    mean, median, std = sigma_clipped_stats(psfModelData)

    #computing psf values for sextractor catalog
    psfName = mask_img_name + '.psf'
    psfcatalogName = mask_img_name+'.psf.cat'
    psfparamName = 'photomPSF.param' #This is a new set of parameters to be obtained from SExtractor, including PSF-fit magnitudes
    #We are supplying SExtactor with the PSF model with the PSF_NAME option
    command = 'sex -c %s %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (configFile, mask_img_name, psfcatalogName, psfName, psfparamName)
    print("Executing command: %s" % command)
    rval = subprocess.run(command.split(), check=True)

    psfsourceTable = get_table_from_ldac(psfcatalogName)

    #Selecting the clean sources away from image edges as before 
    cleanPSFSources = psfsourceTable[(psfsourceTable['FLAGS']==0) & (psfsourceTable['FLAGS_MODEL']==0)  & (psfsourceTable['FWHM_WORLD'] < 2) & (psfsourceTable['XMODEL_IMAGE']<3500) & (psfsourceTable['XMODEL_IMAGE']>500) &(psfsourceTable['YMODEL_IMAGE']<3500) &(psfsourceTable['YMODEL_IMAGE']>500)]
    # print(cleanPSFSources.keys())
    # fig = plt.figure(figsize=(10,10))
    # ax = fig.gca()
    # plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma)
    # #plotting circles on top of all detected sources
    # circles = [plt.Circle((source['XWIN_IMAGE'], source['YWIN_IMAGE']), radius = 5, edgecolor='r', facecolor='None') for source in cleanPSFSources]
    # for c in circles:
    #     ax.add_artist(c)


    # plt.show()

    cleanPSFSourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHAWIN_J2000'], dec=cleanPSFSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
    PSFSourceCatCoords = SkyCoord(ra=psfsourceTable['ALPHAWIN_J2000'], dec=psfsourceTable['DELTAWIN_J2000'], frame='icrs', unit='degree')
    print(len(cleanPSFSourceCatCoords))
    #Now cross match sources
    #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
    photoDistThresh = 0.6
    idx_psfimage, idx_psfps1, d2d, d3d = ps1CatCoords.search_around_sky(cleanPSFSourceCatCoords, photoDistThresh*u.arcsec)

    print('Found %d good cross-matches'%len(idx_psfimage))


    # plt.figure(figsize=(8,8))
    # #Plotting instrumental magnitude for aperture sizes of 5.0, 6.0 and 7.0 pixels
    # plt.plot(good_cat_stars['rmag'][idx_psfps1], cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage] , 'r.', markersize=14)
    # #plt.ylim(-16, -7.5)
    # plt.xlabel('PS1 magnitude', fontsize=15)
    # plt.ylabel('Instrumental PSF-fit magnitude', fontsize=15)
    # plt.show()

    #ZP Derivation 
    psfoffsets = ma.array(good_cat_stars[f'{band}mag'][idx_psfps1] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
    #Compute sigma clipped statistics
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

    # fig = plt.figure(figsize=(10,10))
    # ax = fig.gca()
    # mean, median, sigma = sigma_clipped_stats(data)
    # plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma)
    # #plotting circles on top of all detected sources
    # circles = [plt.Circle((source['XWIN_IMAGE'], source['YWIN_IMAGE']), radius = 5, edgecolor='r', facecolor='None') for source in psfsourceTable]
    # circles += [plt.Circle((w.world_to_pixel(SkyCoord(ra=tar_ra,dec=tar_dec,frame="icrs",unit="degree"))), radius = 50, edgecolor='b', facecolor='None'),]
    # for c in circles:
    #     ax.add_artist(c)
    # plt.savefig(f"../figures/sextracted_targets_{imageName}.png")
    # plt.show()

    # Calculating PSF Mag of target
    tar_coords = SkyCoord(ra=[tar_ra,], dec=[tar_dec,], frame='icrs', unit='degree')
    print(tar_coords)
    idx_tar, idx_psf_tar, d2d, d3d = PSFSourceCatCoords.search_around_sky(tar_coords, 2*u.arcsec)
    try:
        print('Found the source at index %d'%idx_psf_tar[0])

        tar_psfinstmag = psfsourceTable[idx_psf_tar]['MAG_POINTSOURCE'][0]
        tar_psfinstmagerr = psfsourceTable[idx_psf_tar]['MAGERR_POINTSOURCE'][0]

        tar_psfmag = zero_psfmed + tar_psfinstmag
        tar_psfmagerr = np.sqrt(tar_psfinstmagerr**2 + zero_psfstd**2)

        print('PSF-fit magnitude of tar is %.2f +/- %.2f'%(tar_psfmag, tar_psfmagerr))
        os.system(f'echo "{header["JD"]},{tar_psfmag:.2f},{tar_psfmagerr:.2f}" >> ../results.txt')
    except IndexError:
        print(f"Could not find source in image {imageName}")
    return
os.system(f'echo "JD,tar_psfmag,tar_psfmagerr" > results.txt')
os.chdir("data")
# images = glob.glob("../../*.wcs.fits")
# images = ["20230410204702-229-RA.wcs.proc.fits","20230411204011-067-RA.wcs.proc.fits","20230412203312-198-RA.wcs.proc.fits","20230414203330-318-RA.wcs.proc.fits","20230415202558-023-RA.wcs.proc.fits","20230421200129-113-RA.wcs.proc.fits","20230422195357-965-RA.wcs.proc.fits"]
# images = ["20230410204702-229-RA.wcs.proc.fits","20230411204011-067-RA.wcs.proc.fits","20230414203330-318-RA.wcs.proc.fits","20230415202558-023-RA.wcs.proc.fits","20230422195357-965-RA.wcs.proc.fits"]
images = ["20230317225123-039-RA.wcs.proc.fits","20230321211314-542-RA.wcs.proc.fits","20230322231059-990-RA.wcs.proc.fits","20230326214525-048-RA.wcs.proc.fits","20230404211143-368-RA.wcs.proc.fits","20230409205209-789-RA.wcs.proc.fits","20230413203519-882-RA.wcs.proc.fits","20230416222652-985-RA.wcs.proc.fits","20230417204141-924-RA.wcs.proc.fits"]
for i in range(len(images)):
    psf_photometry(i, os.path.basename(images[i]))