#!/bin/bash
# Install astrometry.net 
sudo apt-get update
sudo apt-get install -y astrometry.net 
# this will install the index files for the astrometry.net
# in the directory /usr/share/astrometry
# although tycho2 files suffice for some images
# it could be the case we need to install more index files
# sudo apt-get install astrometry-data-tycho2
# the 4000 series index files
# check if ths directory /usr/share/astrometry exists
# if not create it
if [ ! -d /usr/share/astrometry ]; then
    sudo mkdir /usr/share/astrometry
fi
cd /usr/share/astrometry
for ((i=7; i<20; i++)); do
    I=$(printf %02i $i)
    sudo wget http://data.astrometry.net/4100/index-41$I.fits
done
# index-5000-*
for ((i=0; i<48; i++)); do
    I=$(printf %02i $i)
    sudo wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5200-$I.fits
done

for ((i=0; i<48; i++)); do
    I=$(printf %02i $i)
    sudo wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5201-$I.fits
done

for ((i=0; i<48; i++)); do
    I=$(printf %02i $i)
    sudo wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5202-$I.fits
done