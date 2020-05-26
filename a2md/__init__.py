import json
import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def open_file(name):
    with open(name) as f:
        try:
            FILE = json.load(f)
        except json.JSONDecodeError:
            raise IOError("could not read amd params file")
        return FILE

# LIBRARY_PATH = os.path.dirname(os.path.abspath(__file__))
LIBRARY_PATH = Path(__file__).parent
PARAMETERS_PATH = LIBRARY_PATH / 'parameters/'

TOPO_RESTRICTED_FILE = PARAMETERS_PATH / "a2md_topo_bonded_model.json"
TOPO_RESTRICTED_PARAMS = open_file(TOPO_RESTRICTED_FILE)
HARMONIC_TOPO_RESTRICTED_FILE = PARAMETERS_PATH / "a2md_htopo_bonded_model_extended.json"
HARMONIC_TOPO_RESTRICTED_PARAMS = open_file(HARMONIC_TOPO_RESTRICTED_FILE)
EXTENDED_TOPO_RESTRICTED_FILE = PARAMETERS_PATH / "a2md_topo_bonded_model_extended.json"
EXTENDED_TOPO_RESTRICTED_PARAMS = open_file(EXTENDED_TOPO_RESTRICTED_FILE)
SPHERICAL_FILE = PARAMETERS_PATH / "a2md_spherical_model.json"
SPHERICAL_PARAMS = open_file(SPHERICAL_FILE)
# SYMMETRY_PARAMETERS = open_file(LIBRARY_PATH / '/parameters/symmetry_params.json')

LEBEDEV_DESIGN = dict(
    tight=PARAMETERS_PATH / "lebedev_101.txt",
    medium=PARAMETERS_PATH / "lebedev_053.txt",
    coarse=PARAMETERS_PATH / "lebedev_027.txt"
)