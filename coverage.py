# Configuration for coverage.py

[run]
branch = True
source = pyart
include = */pyart/*
omit =
    */setup.py

[report]
exclude_lines =
    if self.debug:
    if debug:
    raise NotImplementedError
