from source_library.coreSource import getDSSImage, runAstrometry, wcsinfo, getUCAC4Objects, removeBackground
from source_library.coreSource import curveOfGrowth, doPhotometry
from source_library.helperFunctions import resolve_coordinates
from source_library.helperFunctions import printTable
import os

def test_resolve_coordinates():
    target = "NGC330"
    # check that the return is correct to 5 decimal places
    result = resolve_coordinates(target)
    assert round(result[0], 5) == 14.07729
    assert round(result[1], 5) == -72.46253

def test_getDSSImage():
    ra = '00 56 18.55'
    dec = '-72 27 45.1'
    size = 500
    assert getDSSImage(ra, dec, size) == 'DSSImage.fits'

def test_runAstrometry():
    filename = "DSSImage.fits"
    assert runAstrometry(filename).find('success') != -1
     
def test_wcsinfo():
    # check that the return is correct ish
    filename = "DSSImage_solved.fits"
    result = wcsinfo(filename)
    # test if the wcsinfo function returns the correct values
    # to within plus / minus 10% of the expected values
    print(result)
    assert abs((float(result[0]) / 0.93938) - 1) < 0.1
    assert abs((float(result[1]) / -72.462944) - 1) < 0.1
    assert abs((float(result[4]) / 14.17) - 1) < 0.1
    assert abs((float(result[5]) / 14.17) - 1) < 0.1

def test_getUCAC4Objects():
    filename = "DSSImage_solved.fits"
    # check that a dataframe is returned
    assert getUCAC4Objects(filename) is not None
    # check that the dataframe has the correct number of rows
    assert len(getUCAC4Objects(filename)) > 10
    # check that the dataframe has the correct number of columns
    assert len(getUCAC4Objects(filename).columns) == 24
    # check that the dataframe has the correct column names
    assert getUCAC4Objects(filename).columns[0] == 'UCAC4'
    assert getUCAC4Objects(filename).columns[1] == 'RAJ2000'

def test_removeBackground():
    filename = "DSSImage_solved.fits"
    removeBackground(filename)
    assert os.path.isfile("DSSImage_solved_nobkg.fits")

def test_curveOfGrowth():
    filename = "DSSImage_solved_nobkg.fits"
    df = getUCAC4Objects(filename)
    photAp, skyApInner, skyApOuter = curveOfGrowth(df, filename)
    assert abs(photAp/3.85 -1) <0.1
    assert abs(skyApInner/6.14 -1) <0.1
    assert abs(skyApOuter/11.15 -1) <0.1
    assert photAp < skyApInner
    assert skyApInner < skyApOuter

def test_doPhotometry():
    filename = "DSSImage_solved_nobkg.fits"
    df = getUCAC4Objects(filename)
    phot_df = doPhotometry(df, filename)
    assert phot_df is not None
    assert len(phot_df) == len(df)
    assert 'aperture_sum_bkgsub' in phot_df.columns

