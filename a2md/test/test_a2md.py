from a2md.preprocessor import preprocessor
from a2md.preprocessor import symmetry_features
from a2md.amd import A2MD
from a2mdtests import gdb_test
from a2md import SYMMETRY_PARAMETERS
import numpy as np
import json

MOL2 = gdb_test['gdb_000214']['path'] / gdb_test['gdb_000214']['mol2']
TRAIN = gdb_test['gdb_000214']['path'] / list(gdb_test['gdb_000214']['surfaces'].keys())[0]
OUT = './gdb_000214/a2md_test.ppp'
DX = './gdb_000214/a2md_test.dx'

if __name__ == '__main__':
    pp = preprocessor(
        file=MOL2,
        input_format='mol2'
    )
    bps_featurizer = symmetry_features(file=MOL2)
    radial_feats = bps_featurizer.get_radial_features(SYMMETRY_PARAMETERS['radial'])
    elements = bps_featurizer.get_atomic_nums()
    wx, wc, top, pars = pp.parametrize(kind='frozen_core_gaussian')

    ibu_a2md = A2MD(
        coordinates=wx,
        charge=wc,
        topology=top,
        parameters=pars,
        atomic_numbers=elements
    )

    training = np.loadtxt(TRAIN, delimiter=',', skiprows=1)
    ibu_a2md.optimize(
        target_coordinates=training[:,:3],
        target_density=training[:,3],
        weigths=np.ones(training.shape[0])/training.shape[0],
        method='restricted',
    )

    optimized_parameters = ibu_a2md.get_parametrization()
    with open(OUT, 'w') as f:
        json.dump(optimized_parameters, f, indent=4, sort_keys=True)

    ibu_dx = ibu_a2md.eval_volume(extend=2.0, resolution=0.25)
    ibu_dx.write('./ibuprofen/ibuprofen_ns.dx')

    print(ibu_a2md.eval_potential_ne())
    print("HOLA!")