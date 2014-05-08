# remove and rebuild the _fourdd_interface.so module
rm  _fourdd_interface.so
#cython _fourdd_interface.pyx
#export RSL_PATH=/Users/jhelmus/anaconda/
python setup.py build_ext -i
#install_name_tool -change librsl.1.dylib /Users/jhelmus/anaconda/lib/librsl.1.dylib _fourdd_interface.so
