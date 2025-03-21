#!/usr/bin/env python
from source_library.coreSource import create_logger, getDSSImage
from source_library.helperFunctions import resolve_coordinates
from source_library.coreSource import runAstrometry, wcsinfo, getUCAC4Objects,\
    removeBackground, curveOfGrowth, doPhotometry, examine_zeropoint
from source_library.coreSource import set_logger, getSMphotometry

logger = create_logger()
set_logger(logger)

if __name__ == "__main__":
    # target = "NGC330"
    # size = 500

    # # Resolve the coordinates of the target
    # ra, dec = resolve_coordinates(target)
    # print(f"RA: {ra}, DEC: {dec}")

    # # Get the DSS image
    # filename = getDSSImage(ra, dec, size)
    # print(f"Image saved as: {filename}")


# The photometric properties of the DSS photographic plates are 
# highly non-linear: basically, a comparison of the magnitudes of refcat to
# those derived form the DSS show a general trend of brighter objects having 
# more counts but the relationship is not linear.

# So I will now use an image from astrometry.net

    filename = "erty.fits"

    # #run astrometry.net on the image
    # rtn = runAstrometry(filename)
    # #has the image been solved? if it has not raise error
    # if rtn.find('success') == -1:
    #     raise ValueError(f"Image {filename} could not be solved")

    # # append _solved to the filename
    filename = filename.replace('.fits', '_solved.fits')
    rtn = wcsinfo(filename)
    logger.info(f"wcsinfo: {rtn}")

    #refcat = getUCAC4Objects(filename)
    refcat = getSMphotometry(filename)
    logger.info(f"Reference Catalog: {refcat}")

    removeBackground(filename)
    #append _nobkg to the filename
    filename = filename.replace('.fits', '_nobkg.fits')

    curveOfGrowth(refcat, filename)
    

    phot_df = doPhotometry(refcat, filename)  
    logger.info(f"Photometry: {phot_df}")
    merged = examine_zeropoint(refcat, phot_df, photband='gPSF', iplots=True)
    logger.info(f"Merged: {merged}")
    logger.info("Done")