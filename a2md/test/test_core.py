import numpy as np

from a2md.preprocessor import preprocessor, symmetry_features
from a2md.amd import A2MD
from a2md import SYMMETRY_PARAMETERS


if __name__ == '__main__':

    ppp = preprocessor()
    ppp.read_molecule(file='zinc897291/ZINC000000897291.mol2', input_format='mol2')
    x, q, t, pars = ppp.parametrize(kind='frozen_core_gaussian')
    an = ppp.get_atomic_numbers()

    training_data = np.loadtxt(fname='zinc897291/ZINC000000897291.csv', delimiter=',', skiprows=1)

    instance = A2MD(
        coordinates=x, charge=q, topology=t,
        parameters=pars, verbose=False,
        atomic_numbers=an
    )

    bpsf = symmetry_features(file='zinc897291/ZINC000000897291.mol2')
    sym_feats = bpsf.get_radial_features(SYMMETRY_PARAMETERS['radial'])

    opt_indicators = instance.optimize(
        target_coordinates=training_data[:,:3],
        target_density=training_data[:, 3],
        weigths=np.ones(training_data.shape[0]) /training_data.shape[0],
        method='restricted', atom_features=sym_feats
    )

    print(instance.eval_core(instance.get_coordinates()))
