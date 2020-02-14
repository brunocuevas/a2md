from a2md.preprocessor import preprocessor, symmetry_features
from a2md import SYMMETRY_PARAMETERS
from a2md.amd import A2MD
from a2mdtests import gdb_test
import numpy as np
from a2md.utils import element2an, create_all2all_topology
import json

MOL2_FILE= gdb_test['gdb_000214']['path'] / gdb_test['gdb_000214']['mol2']
TRAINING_FILE=gdb_test['gdb_000214']['path'] / list(gdb_test['gdb_000214']['surfaces'].keys())[3]
MODEL_FILE=r'C:\Users\Bruno\ownCloud\main\a2md\a2md\parameters\a2md_topo_bonded_model.json'
OPT_PARAMS=r'C:\Users\Bruno\ownCloud\main\a2md\a2md\test\gdb_000214\model_file.ppp'
DX_FILE=r'C:\Users\Bruno\ownCloud\main\a2md\a2md\test\gdb_000214\model_file_ep.dx'

if __name__ == '__main__':

    ppp = preprocessor(file=MOL2_FILE, input_format='mol2')
    _, _, _, fc_pars = ppp.parametrize(kind='frozen_core_gaussian')

    ppp.set_topology(create_all2all_topology(12))
    out_coords, out_charge, out_topo, out_pars = ppp.parametrize_from_modelfile(MODEL_FILE)

    an = np.array([element2an(i) for i in ppp.get_labels()])
    bps_featurizer = symmetry_features(file=MOL2_FILE)
    radial_feats = bps_featurizer.get_radial_features(SYMMETRY_PARAMETERS['radial'])
    print("Test Preprocessor, Completed")
    density_model = A2MD(
        coordinates=out_coords,
        charge=out_charge,
        topology=out_topo,
        parameters=out_pars,
        atomic_numbers=an
    )

    print("Test A2MD, Completed")

    training = np.loadtxt(TRAINING_FILE, delimiter=',', skiprows=1)

    density_model.optimize(
        target_coordinates=training[:, :3],
        target_density=training[:, 3],
        weigths=np.ones(training.shape[0]) / training.shape[0],
        method='restricted',

    )

    with open(OPT_PARAMS, 'w') as f:
        json.dump(density_model.get_parametrization(), f, indent=4, sort_keys=True)

    dx = density_model.eval_volume(extend=3.0, resolution=0.25, field='ep')

    dx.write(DX_FILE)