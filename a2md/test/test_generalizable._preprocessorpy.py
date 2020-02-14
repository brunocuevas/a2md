from a2md.amd import A2MD
from a2md.preprocessor import generalizable_preprocessor, preprocessor
from a2md.preprocessor import a2md_general_parameters
import pandas as pd
import numpy as np

test_params = a2md_general_parameters(version=1.0,verbose=True, filename="test_params.json")
test_params.read_parameters()

gp = generalizable_preprocessor(verbose=True)
gp.read_custom_parameters(
    filename="test_params.json"
)
# --------------------------------------------------- #
gp.read_molecule(file='../input/water_gas.aamd')
wx, wc, wt, pars = gp.parametrize(kind='trigonometric')
# --------------------------------------------------- #
targets = pd.read_csv('../output/water_gas_s2.csv').values
# --------------------------------------------------- #
water_aamd = A2MD(coordinates=wx, charge=wc, topology=wt, parameters=pars)
water_aamd.read(pars)
water_aamd.optimize(
    target_coordinates=targets[:,:3],
    target_density=targets[:,3],
    weigths=np.ones(targets.shape[0])
)
