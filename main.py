from source_library.coreSource import getDSSImage

# Example usage of the functions
ra = '10h00m00s'
dec = '+02d00m00s'
size = 500

# Get the DSS image
filename = getDSSImage(ra, dec, size)
print(f'Image saved as: {filename}')
