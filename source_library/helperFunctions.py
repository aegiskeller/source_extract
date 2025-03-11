import astropy.io.fits as fits

# create a function that prints the contents of a fits table
def printTable(filename):
     # open the table file
    hdul = fits.open(filename)
    # print the header
    print(hdul[0].header)
    # print the data
    print(hdul[1].data)
    # close the table file
    hdul.close()