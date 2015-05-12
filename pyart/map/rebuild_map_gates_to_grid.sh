# remove and rebuild the _load_nn_field_data module
rm  _map_gates_to_grid.so
cython -a _map_gates_to_grid.pyx
python setup.py build_ext -i
