#!/bin/bash
# This script is adapted from the install.sh script from the scikit-learn
# project: https://github.com/scikit-learn/scikit-learn

# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

set -e
set -x

sudo apt-get update -qq

# Use Miniconda to provide a Python environment.  This allows us to perform
# a conda based install of the SciPy stack on multiple versions of Python
# as well as use conda and binstar to install additional modules which are not
# in the default repository.
wget http://repo.continuum.io/miniconda/Miniconda-3.5.2-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda/bin:$PATH
conda update --yes conda
conda update --yes conda

# Create a testenv with the correct Python version
conda create -n testenv --yes pip python=$PYTHON_VERSION
source activate testenv

# Install Py-ART dependencies
conda install --yes numpy scipy matplotlib netcdf4 nose
conda install --yes -c http://conda.binstar.org/jjhelmus trmm_rsl
sudo apt-get install -qq libglpk-dev
conda install --yes -c http://conda.binstar.org/jjhelmus pyglpk

if [[ $PYTHON_VERSION == '2.7' ]]; then
    conda install --yes basemap 
    conda install --yes -c http://conda.binstar.org/jjhelmus cbc cylp
    conda install --yes -c http://conda.binstar.org/jjhelmus cvxopt_glpk
fi

# install coverage modules
pip install nose-cov
if [[ "$COVERALLS" == "true" ]]; then
    pip install python-coveralls
fi

# install Py-ART
export RSL_PATH=~/miniconda/envs/testenv

if [[ "$FROM_RECIPE" == "true" ]]; then
    source deactivate
    conda install --yes conda-build
    conda config --add channels http://conda.binstar.org/jjhelmus
    source activate testenv
    export CONDA_PACKAGE=`conda build --output conda_recipe/`
    conda build -b -q conda_recipe/
    echo $CONDA_PACKAGE
    conda install $CONDA_PACKAGE
    mkdir foo   # required so source directory not picked up during tests
    cd foo
else
    python setup.py build_ext --inplace
fi
