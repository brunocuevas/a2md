from a2mdnet import LIBRARY_PATH
from pathlib import Path
import json
import os

class TestDataset:
    def __init__(self, dataset_dict):
        try:
            self.environ_path = Path(os.environ[dataset_dict['path']])
        except KeyError:
            print(("you need to speficy the environment variable {:s} pointing towards the folder "
                                    "where data is stored".format(dataset_dict['path'])))
        self.molecular_path = self.environ_path / dataset_dict['molecular_path']
        self.parameters_path = self.environ_path / dataset_dict['parameters_path']
        self.density_path = self.environ_path / dataset_dict['density_path']
        with open(self.environ_path / dataset_dict['idx']) as f:
            self.idx = json.load(f)
        self.energy = self.environ_path / dataset_dict['energy']
        self.other_idx = dict(
            (key, value) for key, value in dataset_dict.items() if key[-3:] == 'idx'
        )

aniset_dict = dict(
    path="ANISETDIR",
    molecular_path=r'mol2/',
    parameters_path=r'ppp/',
    density_path=r'npy/',
    idx=r'aniset_all.json',
    validation_idx=r'aniset_validation.json',
    training_idx=r'aniset_training.json',
    test_idx=r'aniset_test.json',
    energy=r'aniset_all_energies.json'
)
fdaset_dict = dict(
    path="FDASETDIR",
    molecular_path=r'mol2/',
    parameters_path=r'ppp/',
    density_path=r'npy/',
    idx=r'fdaset.json',
    energy=r'fdaset_all_energies.json'
)

aniset = TestDataset(aniset_dict)
fdaset = TestDataset(fdaset_dict)
