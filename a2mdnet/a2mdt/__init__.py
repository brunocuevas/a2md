import json
import os
LIBRARY_PATH = os.path.dirname(os.path.abspath(__file__))
with open(LIBRARY_PATH + r'/a2md_topo_bonded_model.old.json') as f:
    OLD_A2MD_MODEL=json.load(f)

with open(LIBRARY_PATH + r'/a2md_topo_bonded_model.json') as f:
    A2MD_MODEL=json.load(f)