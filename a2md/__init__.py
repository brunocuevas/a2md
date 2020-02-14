import json
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))



def open_file(name):
    with open(name) as f:
        try:
            FILE = json.load(f)
        except json.JSONDecodeError:
            raise IOError("could not read amd params file")
        return FILE

LIBRARY_PATH = os.path.dirname(os.path.abspath(__file__))
PARAMETERS_PATH = LIBRARY_PATH + '/parameters/'
AMD_PARAMS_18_FILE= PARAMETERS_PATH + 'amd_params_18'
AMD_PARAMS_19_FILE= PARAMETERS_PATH + 'amd_params_19'
RADIUS_FILE = PARAMETERS_PATH + 'radius.json'
BOND_FILE = PARAMETERS_PATH + 'bond_types.json'


AMD_PARAMS_18 = open_file(AMD_PARAMS_18_FILE)
AMD_PARAMS_19 = open_file(AMD_PARAMS_19_FILE)
RADIUS = open_file(RADIUS_FILE)
BOND_TYPES = open_file(BOND_FILE)
with open(LIBRARY_PATH + '/parameters/symmetry_params.json') as f:
    SYMMETRY_PARAMETERS = json.load(f)