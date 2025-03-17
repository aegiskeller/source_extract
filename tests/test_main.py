from source_library.coreSource import getDSSImage
from source_library.helperFunctions import resolve_coordinates
from source_library.helperFunctions import printTable
import os
import subprocess
def test_hello_world():
    assert 1 + 1 == 2

# def test_printTable():
#     filename = "tests/data/2mass-atlas-990502s-j1430240.fits"
#     assert printTable(filename) == None

def test_resolve_coordinates():
    target = "NGC330"
    assert resolve_coordinates(target) == ('00 56 18.55', '-72 27 45.1')

def test_getDSSImage():
    ra = '00 56 18.55'
    dec = '-72 27 45.1'
    size = 50
    assert getDSSImage(ra, dec, size) == 'DSSImage.fits'

def test_runAstrometry():
    # wget the test image
    # "https://g5555.neocities.org/input-image.jpg"
    url = "https://g5555.neocities.org/input-image.jpg"
    filename = "input-image.jpg"
    subprocess.run(["wget", url, "-O", filename], check=True)
    #check that the return has 'success' contained in it
    assert runAstrometry(filename).find('success') != -1
     

    