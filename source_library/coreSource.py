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
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import random
import string
import os
import shutil
import time

# define the function getDSSImage
def getDSSImage(ra, dec, size):
    # get the image from the DSS database
    images = SkyView.get_images(position=f'{ra} {dec}', survey='DSS', coordinates='J2000', pixels=size)
    # get the data from the image
    data = images[0][0].data
    # save the image in a file
    filename = 'DSSImage.fits'
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
    plt.imshow(data, origin='lower', cmap='gray')
    plt.show()
    # save the image in a png file
    # append .png to the filename
    pngfile = filename + '.png'
    plt.imsave(pngfile, data, cmap='gray')
    # close the image file
    hdul.close()

#a function to print the image header for an image
def printImageHeader(filename):
    #open the image file
    hdul = fits.open(filename)
    #print the header
    print(hdul[0].header)
    #close the image file
    hdul.close()

# a function to run astrometry.net on an image
def runAstrometry(filename):
    # create a uniquely named directory by selecting a random string
    random_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
    try: 
        # create the directory  
        os.mkdir(f'{random_string}') 
        # run astrometry.net on the image
        start_time = time.time()
        os.system(f'solve-field --guess-scale --no-plots -D {random_string} {filename}')
        duration = time.time() - start_time
        print(f'Astrometry.net ran in {duration} seconds')
        # get the name of the new image
        newfile = filename.replace('.fits', '.new')
        # return the name of the new image
        os.rename(f'{random_string}/{newfile}', newfile.replace('.new', '_solved.fits'))
        # remove the directory
        shutil.rmtree(f'{random_string}')
    except Exception as e:
        print(f'Error: {e}')
        return('failed')
    return (f'success-{duration}')

# a function that executes wcsinfo on an image
# and extracts ra_center and dec_center, the pixscale, orientation, fieldw, fieldh, fieldunits
def wcsinfo(filename):
    # run wcsinfo on the image
    try:
        wcsinfo_output = os.popen(f'wcsinfo {filename}').read()
    except Exception as e:
        print(f'Error: {e}')
        return('failed')
    # extract the values from the output
    wcsinfo_output = wcsinfo_output.split('\n')
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
        ra_center, dec_center, pixscale, orientation, fieldw, fieldh, fieldunits = wcsinfo(filename)
    except Exception as e:
        print(f'Error: {e}')
        return 'failed - no wcs info'
    
    # Create a SkyCoord object for the center of the field
    center_coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame='icrs')
    
    # Define the search radius
    search_radius = max(float(fieldw), float(fieldh)) / 2 * u.arcmin
    
    # Query the UCAC4 catalog
    v = Vizier(columns=['*'], catalog='I/322A')
    result = v.query_region(center_coord, radius=search_radius)
    
    if result:
        ucac4_table = result[0]
        ucac4_df = ucac4_table.to_pandas()
        # drop rows with missing Vmag
        ucac4_df = ucac4_df.dropna(subset=['Vmag'])
        # drop rows with bad object quality
        ucac4_df = ucac4_df[ucac4_df['of'] == 0]
        return ucac4_df
    else:
        return 'No objects found'
