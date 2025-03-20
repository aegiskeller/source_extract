# here we develop a function that obtains an image from DSS and saves it in a file
# the function is called getDSSImage
# the function takes as input the coordinates of the object and the size of the image
# the function returns the name of the file where the image was saved
# the function uses the astroquery module to access the DSS database
# the function uses the astropy module to manipulate the image
# the function uses the matplotlib module to display the image
# the function uses the os module to save the image in a file

# import the necessary modules
from astroquery.skyview import SkyView
from astropy.io import fits
from astroquery.vizier import Vizier
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.profiles import CurveOfGrowth
from photutils.aperture import ApertureStats
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import string
import os
import shutil
import time
import logging

logger = None
def set_logger(external_logger):
    global logger
    logger = external_logger

def create_logger():
    # create a logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler("coreSource.log")
    handler.setLevel(logging.DEBUG)
    # create a logging format
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    # add the handler to the logger
    logger.addHandler(handler)  # add the handler to the logger
    # also handle the logging to the console
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(formatter)
    logger.addHandler(consoleHandler)
    return logger

# define the function getDSSImage
def getDSSImage(ra, dec, size):
    """
    This function gets an image from the DSS database
    and saves it in a file.
    ra: the right ascension of the object in degrees
    dec: the declination of the object in degrees
    size: the size of the image in pixels
    """
    # get the image from the DSS database
    images = SkyView.get_images(
        position=f"{ra} {dec}", survey="DSS", coordinates="J2000", pixels=size
    )
    # get the data from the image
    data = images[0][0].data
    # save the image in a file
    filename = "DSSImage.fits"
    fits.writeto(filename, data, overwrite=True)
    # return the name of the file
    return filename


# define a function to display the image
def displayImage(filename):
    # open the image file
    hdul = fits.open(filename)
    # get the data from the image
    data = hdul[0].data
    # display the image
    plt.imshow(data, origin="lower", cmap="gray")
    plt.show()
    # save the image in a png file
    # append .png to the filename
    pngfile = filename + ".png"
    plt.imsave(pngfile, data, cmap="gray")
    # close the image file
    hdul.close()


# a function to print the image header for an image
def printImageHeader(filename):
    # open the image file
    hdul = fits.open(filename)
    # print the header
    logger.debug(hdul[0].header)
    # close the image file
    hdul.close()


# a function to run astrometry.net on an image
def runAstrometry(filename):
    # create a uniquely named directory by selecting a random string
    random_string = "".join(random.choices(string.ascii_uppercase + string.digits, k=6))
    try:
        # create the directory
        os.mkdir(f"{random_string}")
        # run astrometry.net on the image
        start_time = time.time()
        #os.system(f"solve-field --guess-scale --no-plots -D {random_string} {filename}")
        os.system(f"solve-field -D {random_string} {filename}")
        duration = time.time() - start_time
        logger.info(f"Astrometry.net ran in {duration} seconds")
        # create the filename of the new solved image
        # this is <base name> +'.new'
        # split the filename on '.'
        newfile = filename.split(".")[0]+".new"
        # move the file up a directory
        shutil.move(f"{random_string}/{newfile}", f"./{newfile.replace('.new', '_solved.fits')}")
        # remove the directory
        shutil.rmtree(f"{random_string}")
    except Exception as e:
        logger.error(f"Error: {e}")
        # remove the directory and contents
        shutil.rmtree(f"{random_string}")
        return "failed"
    return f"success-{duration}"


# a function that executes wcsinfo on an image
# and extracts ra_center and dec_center, the pixscale, orientation, fieldw, fieldh, fieldunits
def wcsinfo(filename):
    # run wcsinfo on the image
    try:
        wcsinfo_output = os.popen(f"wcsinfo {filename}").read()
    except Exception as e:
        logger.error(f"Error: {e}")
        return "failed"
    # extract the values from the output
    wcsinfo_output = wcsinfo_output.split("\n")
    ra_center = wcsinfo_output[18].split()[1]
    dec_center = wcsinfo_output[19].split()[1]
    pixscale = wcsinfo_output[16].split()[1]
    orientation = wcsinfo_output[17].split()[1]
    fieldw = wcsinfo_output[33].split()[1]
    fieldh = wcsinfo_output[34].split()[1]
    fieldunits = wcsinfo_output[35].split()[1]
    # return the values
    return ra_center, dec_center, pixscale, orientation, fieldw, fieldh, fieldunits


def getUCAC4Objects(filename):
    """
    This function takes a image with a defined wcs and
    returns the objects in the UCAC4 catalog that are
    within the field of view of the
    image. The function returns a pandas dataframe
    with the objects.
    """
    # get the wcs information from the image
    try:
        ra_center, dec_center, pixscale, orientation, fieldw, fieldh, fieldunits = (
            wcsinfo(filename)
        )
    except Exception as e:
        logger.error(f"Error: {e}")
        return "failed - no wcs info"

    # Create a SkyCoord object for the center of the field
    center_coord = SkyCoord(
        ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame="icrs"
    )

    # Define the search radius
    search_radius = max(float(fieldw), float(fieldh)) / 2 * u.arcmin

    # Query the UCAC4 catalog
    v = Vizier(columns=["*"], catalog="I/322A")
    result = v.query_region(center_coord, radius=search_radius)

    if result:
        ucac4_table = result[0]
        ucac4_df = ucac4_table.to_pandas()
        # drop rows with missing Vmag - too restrictive?
        # ucac4_df = ucac4_df.dropna(subset=["Vmag"])
        # drop rows with bad object quality
        ucac4_df = ucac4_df[ucac4_df["of"] == 0]
        # open the image file to get the WCS information
        try:
            hdul = fits.open(filename)
            wcs = WCS(hdul[0].header)
        except Exception as e:
            logger.error(f"Error: {e}")
            hdul.close()
            return "failed - no WCS info"
        # Create a SkyCoord object for the objects in the catalog
        ucac4_coords = SkyCoord(
            ra=ucac4_df["RAJ2000"],
            dec=ucac4_df["DEJ2000"],
            unit=(u.deg, u.deg),
            frame="icrs",
        )
        # Project the coordinates to the image
        x, y = wcs.world_to_pixel(ucac4_coords)
        ucac4_df["x"] = x
        ucac4_df["y"] = y
        return ucac4_df
    else:
        return "No objects found"


def removeBackground(filename):
    """
    a function that fits a background to the image
    and removes it. This is required prior to aperture photometry
    """
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    try:
        hdul = fits.open(filename)
        data = hdul[0].data
        bkg = Background2D(
            data,
            (50, 50),
            filter_size=(3, 3),
            sigma_clip=sigma_clip,
            bkg_estimator=bkg_estimator,
        )
        data_nobkg = data - bkg.background
        hdul[0].data = data_nobkg
        hdul.writeto(filename.replace(".fits", "_nobkg.fits"), overwrite=True)
        hdul.close()
    except Exception as e:
        logger.error(f"Error: {e}")
        hdul.close()
        return "failed"
    return "success"


# create a function to perform curve of growth analysis
def curveOfGrowth(ucac4_df, filename):
    """
    a function to project the ra,dec of ucac4_df to the image
    coordinates (x,y) and perform curve of growth analysis
    """
    # generate an aperture for each object in ucac4_df[x,y]
    x = ucac4_df["x"]
    y = ucac4_df["y"]
    hdul = fits.open(filename)
    # check that hdul[0].data is not None
    if hdul[0].data is None:
        hdul.close()
        return "failed - no image data"
    positions = list(zip(x, y))
    # screen out objects with x,y outside the image
    # either greater than the image extent in x,y or
    # less than zero
    x_max, y_max = hdul[0].data.shape
    positions = [pos for pos in positions if pos[0] < x_max and pos[1] < y_max]
    positions = [pos for pos in positions if pos[0] > 0 and pos[1] > 0]
    # this could be quite fragile, need to check the radii
    # as these vary depending on the image and pixel scale
    # obtain the header keyword for pixscale
    pixscale = float(wcsinfo(filename)[2])
    if pixscale > 1.0:
        radii = np.arange(1.5, 8, 0.5)
    else:
        radii = np.arange(3, 15, 1)
    # create a curve of growth object
    cog = CurveOfGrowth(hdul[0].data, positions, radii=radii)
    # normalise the curve of growth
    cog.normalize(method="max")
    # option to create a plot of the curve of growth using matplotlib
    # and save to a png file
    cog.plot()
    plt.savefig("curve_of_growth.png")
    plt.close()
    # now determine the best apertures
    ee_vals = [0.9, 0.99]
    photAp, skyApInner = cog.calc_radius_at_ee(ee_vals)
    skyApOuter = skyApInner + 5
    hdul.close()
    return photAp, skyApInner, skyApOuter


def doPhotometry(ucac4_df, filename):
    """
    a function to project the ra,dec of ucac4_df to the image
    coordinates (x,y) and perform aperture photometry
    """
    # open the image file to get the WCS information
    try:
        hdul = fits.open(filename)
        wcs = WCS(hdul[0].header)
        data = hdul[0].data
    except Exception as e:
        print(f"Error: {e}")
        hdul.close()
        return "failed - no WCS info"

    # use the curveofgrowth function to determine the best apertures
    photAp, skyApInner, skyApOuter = curveOfGrowth(ucac4_df, filename)

    # Project the coordinates to the image
    x, y = ucac4_df["x"], ucac4_df["y"]
    # generate an aperture for each object
    positions = list(zip(x, y))
    # screen out objects with x,y outside the image
    # either greater than the image extent in x,y or
    # less than zero
    x_max, y_max = hdul[0].data.shape
    positions = [pos for pos in positions if pos[0] < x_max and pos[1] < y_max]
    positions = [pos for pos in positions if pos[0] > 0 and pos[1] > 0]
    # get the apertures from curveofgrowth
    photAp, skyApInner, skyApOuter = curveOfGrowth(ucac4_df, filename)

    apertures = CircularAperture(positions, r=photAp)
    # use a circular annulus for the background
    annulus_apertures = CircularAnnulus(positions, r_in=skyApInner, r_out=skyApOuter)
    aperstats = ApertureStats(data, annulus_apertures)
    bkg_mean = aperstats.mean

    # perform the photometry in the circular aperture
    phot_table = aperture_photometry(data, apertures)
    # the area of the aperture
    aperture_area = np.pi * photAp**2
    # subtract the background
    phot_bkgsub = phot_table["aperture_sum"] - bkg_mean * aperture_area
    # convert the aperture sum to magnitudes
    phot_bkgsub = -2.5 * np.log10(phot_bkgsub)
    #drop the column with the aperture sum
    phot_table.remove_column("aperture_sum")
    # add the column with the aperture sum with background
    phot_table["inst_mag"] = phot_bkgsub
    # convert the astropy table to a pandas dataframe
    phot_table = phot_table.to_pandas()
    return phot_table

def examine_zeropoint(ucac4, phot, photband='f.mag', iplots=False):
    """
    a function to examine the zeropoint of the photometry

    use the iplots flag to generate plots
    """
    # join the ucac4 and phot dataframes on the ucac4.x and phot.xcenter
    # columns
    merged = pd.merge(ucac4, phot, left_on="x", right_on="xcenter")
    # find the number of rows with valid f.mag
    # and valid inst_mag
    nvalidfmag = len(merged.dropna(subset=["f.mag", "inst_mag"]))
    # and for Vmag
    nvalidVmag = len(merged.dropna(subset=["Vmag", "inst_mag"]))
    # and for Jmag  
    nvalidJmag = len(merged.dropna(subset=["Jmag", "inst_mag"]))
    logger.info(f"Number of rows with valid f.mag: {nvalidfmag}")
    logger.info(f"Number of rows with valid Vmag: {nvalidVmag}")
    logger.info(f"Number of rows with valid Jmag: {nvalidJmag}")
    logger.info(f"The requested photband is: {photband}")
    # calculate the zeropoint
    merged['zeropoint'] = merged[photband] - merged["inst_mag"]
    if iplots:
        # use matplotlib to create a png plot of merged["Vmag"] vs zeropoint
        plt.scatter(merged[photband], merged["inst_mag"])
        plt.xlabel(f"UCAC4 {photband}")
        plt.ylabel("Instrumental Magnitude")
        plt.savefig("comparePhot.png")
        plt.close()
        # and also generate a plot of zeropoint vs Vmag
        plt.scatter(merged[photband], merged["zeropoint"])
        plt.xlabel(f"UCAC4 {photband}")
        plt.ylabel("Zeropoint")
        plt.savefig("zeropoint.png")
        plt.close()
    # return the zeropoint
    # keep only columns x,y, xcenter, ycenter
    # merged = merged[["x", "y", "xcenter", "ycenter", "f.mag", "inst_mag", "zeropoint"]]   
    return merged