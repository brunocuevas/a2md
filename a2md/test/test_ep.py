import json
import numpy as np
from a2md.preprocessor import preprocessor
from a2md.amd import A2MD
from a2md.utils import element2an
from a2md.preprocessor import symmetry_features
from a2md import SYMMETRY_PARAMETERS
from a2mdtests import gdb_test

MOL2_FILE = gdb_test['gdb_000214']['path'] / gdb_test['gdb_000214']['mol2']
TRAINING_FILE = gdb_test['gdb_000214']['path'] / list(gdb_test['gdb_000214']['surfaces'])[3]
DX_FILE = r'C:\Users\Bruno\ownCloud\main\a2md\a2md\test\gdb_000214\gdb_000214_ep.dx'
OPT_PARAMS = r'C:\Users\Bruno\ownCloud\main\a2md\a2md\test\gdb_000214\test_ep.ppp'

if __name__ == '__main__':

    bps_featurizer = symmetry_features(file=MOL2_FILE)
    radial_feats = bps_featurizer.get_radial_features(SYMMETRY_PARAMETERS['radial'])

    pp = preprocessor(
        file=MOL2_FILE,
        input_format='mol2'
    )
    pp.read_parameters(r'C:\Users\Bruno\ownCloud\main\a2md\a2md\parameters\amd_params_19')


    wx, wc, top, pars = pp.parametrize(kind='frozen_core_gaussian')
    an = np.array([element2an(i) for i in pp.get_labels()])

    density_model = A2MD(
        coordinates=wx,
        charge=wc,
        topology=top,
        parameters=pars,
        atomic_numbers=an,
    )

    training = np.loadtxt(TRAINING_FILE, delimiter=',', skiprows=1)
    density_model.optimize(
        target_coordinates=training[:, :3],
        target_density=training[:, 3],
        weigths=np.ones(training.shape[0]) / training.shape[0],
        method='restricted',
        atom_features=radial_feats

    )
    optimized_parameters = density_model.get_parametrization()
    with open(OPT_PARAMS, 'w') as f:
        json.dump(optimized_parameters, f, indent=4, sort_keys=True)

    dx_density = density_model.eval_volume(extend=3.0, resolution=0.15, field='ep')
    dx_density.write(DX_FILE)

    density_model.save_model(file='gdb_000214.ca2md')
