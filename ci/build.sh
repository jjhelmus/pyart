#!/bin/bash

cp -r /home/jhelmus/python/pyart/* .
rm pyart/version.py
ls -la
#export SRC_DIR=/home/jhelmus/python/pyart/
export RSL_PATH=$PREFIX
$PYTHON setup.py install

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
