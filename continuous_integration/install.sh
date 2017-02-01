#!/bin/bash

set -e
# use next line to debug this script
set -x

# Install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda3/bin:$PATH
conda config --set always_yes yes
conda config --set show_channel_urls true
conda update -q conda

## Create a testenv with the correct Python version
conda env create -f continuous_integration/environment.yml
source activate testenv
conda install -c conda-forge pyproj   # KLUDGE to replace jjhelmus::pyproj 

# Install Py-ART dependencies
#conda install -q -c conda-forge numpy scipy matplotlib netcdf4 nose hdf4
#conda install -q -c conda-forge basemap
#conda install -q -c conda-forge trmm_rsl

#if [[ $PYTHON_VERSION == '2.7' ]]; then
    #conda install -q -c http://conda.anaconda.org/jjhelmus cbc cylp
    #conda install -q -c http://conda.anaconda.org/jjhelmus glpk pyglpk
    #conda install -q -c http://conda.anaconda.org/jjhelmus cvxopt_glpk

    # wradlib and dependencies
    #conda install -q h5py
    # KLUDGE libgdal does not report its version dependency on geos which
    # causes either gdal or basemap to break, force the exact libgdal version
    # see: https://github.com/ContinuumIO/anaconda-issues/issues/584
    #conda install -q gdal basemap libgdal=2.0.0=0 krb5 proj4
    #conda install --no-deps -q -c conda-forge wradlib
#fi

# install coverage modules
pip install nose-cov
if [[ "$COVERALLS" == "true" ]]; then
    pip install python-coveralls
fi

# install Py-ART
export RSL_PATH=~/miniconda3/envs/testenv

if [[ "$FROM_RECIPE" == "true" ]]; then
    source deactivate
    conda install -q conda-build
    conda install -q jinja2 setuptools
    conda config --add channels http://conda.anaconda.org/jjhelmus
    source activate testenv
    conda build --no-test -q conda_recipe/
   
    export CONDA_PACKAGE=`conda build --output conda_recipe/ | grep bz2`
    conda install -q $CONDA_PACKAGE
    mkdir foo   # required so source directory not picked up during tests
    cd foo
else
    python setup.py build_ext --inplace
fi
