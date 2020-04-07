import os
import time
import json
from pathlib import Path

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
