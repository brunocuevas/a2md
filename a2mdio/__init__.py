import os
import time
import json
from pathlib import Path
from mendeleev import get_table

path2files = Path(__file__).parent


class A2MDlib:
    def __init__(
            self,
            name='aAMDlib',
            verbose=False,
    ):
        """

        :param name:
        :param verbose:
        """
        self.__name = name
        self.__v = verbose

    def log(self, message):
        """

        :param message:
        :return:
        """
        if self.__v:
            ltm = time.localtime()
            mssg_head = "[{:02d}:{:02d}:{:02d}] ".format(ltm.tm_hour, ltm.tm_min, ltm.tm_sec)
            print(mssg_head + message)

with open(path2files / 'parameters/standard_topology.json') as f:
    PDB_PROTEIN_TOPOLOGY = json.load(f)
with open(path2files / 'parameters/standard_labels.json') as f:
    PDB_PROTEIN_TYPES = json.load(f)
with open(path2files / 'parameters/standard_charges.json') as f:
    PDB_PROTEIN_CHARGES= json.load(f)
with open(path2files / 'parameters/standard_charges_by_atom.json') as f:
    PDB_PROTEIN_TYPE_CHARGES = json.load(f)
with open(path2files / 'parameters/babel2standard.json') as f:
    BABEL2STANDARD = json.load(f)

PERIODIC_TABLE = get_table('elements')
MAP_AN2SYMBOL = PERIODIC_TABLE[['atomic_number', 'symbol']].set_index('atomic_number')
MAP_SYMBOL2AN = PERIODIC_TABLE[['atomic_number', 'symbol']].set_index('symbol')

with open(path2files / 'parameters/wfn_symmetry_index.json') as f:
    WFN_SYMMETRY_INDEX = json.load(f)


def get_atomic_number(symbol):
    an = MAP_SYMBOL2AN.loc[symbol]['atomic_number']
    return an


def get_symbol(atomic_number):
    atomic_number = int(atomic_number)
    s = MAP_AN2SYMBOL.loc[atomic_number]['symbol']
    return s
