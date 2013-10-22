# remove and rebuild _echo_classification module
rm _echo_classification.so
rm -r build
python setup.py build_ext -i
