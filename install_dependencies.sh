#!/bin/bash
# Install astrometry.net 
sudo apt-get update
sudo apt-get install -y astrometry.net 
# this will install the index files for the astrometry.net
# in the directory /usr/share/astrometry
# although tycho2 files suffice for some images
# it could be the case we need to install more index files
sudo apt-get install astrometry-data-tycho2
