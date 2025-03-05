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
import matplotlib.pyplot as plt
import os

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

#a function to proint the image header for an image
def printImageHeader(filename):
    #open the image file
    hdul = fits.open(filename)
    #print the header
    print(hdul[0].header)
    #close the image file
    hdul.close()

# a function to run astrometry.net on an image
def runAstrometry(filename):
    # run astrometry.net on the image
    os.system(f'solve-field --guess-scale -D solver-output --no-plots --overwrite {filename}')
    # get the name of the new image
    newfile = filename.replace('.fits', '-wcs.fits')
    # return the name of the new image
    return newfile