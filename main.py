#!/usr/bin/env python
from source_library.coreSource import getDSSImage
from source_library.helperFunctions import resolve_coordinates
from source_library.coreSource import runAstrometry

if __name__ == "__main__":
    target = "NGC330"
    size = 500

    # Resolve the coordinates of the target
    ra, dec = resolve_coordinates(target)
    print(f"RA: {ra}, DEC: {dec}")

    # Get the DSS image
    filename = getDSSImage(ra, dec, size)
    print(f"Image saved as: {filename}")

    #run astrometry.net on the image
    runAstrometry(filename)