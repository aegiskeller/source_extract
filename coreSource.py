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
    # create a new figure
    plt.figure()
    # display the image
    plt.imshow(data, cmap='gray')
    # save the image in a file
    filename = 'DSSImage.fits'
    fits.writeto(filename, data, overwrite=True)
    # return the name of the file
    return filename

# test the function by obtaining an image of NGC 330
ra = 10.721
dec = -72.080
size = 512
filename = getDSSImage(ra, dec, size)
print(f'The image was saved in the file {filename}')
# save the resulting plot as a png file
plt.savefig('DSSImage.png')

#a function to proint the image header for an image
def printImageHeader(filename):
    #open the image file
    hdul = fits.open(filename)
    #print the header
    print(hdul[0].header)
    #close the image file
    hdul.close()