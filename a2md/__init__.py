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

TOPO_RESTRICTED_FILE = PARAMETERS_PATH + "a2md_topo_bonded_model.json"
TOPO_RESTRICTED_PARAMS = open_file(TOPO_RESTRICTED_FILE)
with open(LIBRARY_PATH + '/parameters/symmetry_params.json') as f:
    SYMMETRY_PARAMETERS = json.load(f)

LEBEDEV_DESIGN = dict(
    tight=PARAMETERS_PATH + "lebedev_101.txt",
    medium=PARAMETERS_PATH + "lebedev_053.txt",
    coarse=PARAMETERS_PATH + "lebedev_027.txt"
)