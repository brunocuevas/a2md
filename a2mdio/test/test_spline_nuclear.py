from a2mdio.molecules import Mol2, UNITS_TABLE
from a2mdio.utils import splined_nuclear_charge, eval_volume, eval_charge
from a2mdtest.a2mdtests import benzene
from a2md.models import a2md_from_mol
import numpy as np
import matplotlib.pyplot as plt
if __name__ == '__main__':

    mm = Mol2(benzene.mol2)
    # r = mm.get_coordinates()
    # q = mm.get_absolute_charges()
    # dx = eval_volume(
    #     lambda x : splined_nuclear_charge(x, r, q, 0.13 * UNITS_TABLE['angstrom']['au'] ), size=60, dims=0.13, shift=[0,0,0]
    # )
    # dx.write('splined.dx')

    dm = a2md_from_mol(mm)

    dm.parametrize()

    chgm = lambda x : eval_charge(x, dm, 0.13*UNITS_TABLE['angstrom']['au'])
    dx = eval_volume(chgm, size=60, resolution=0.13, shift=[0, 0, 0])
    dx.write('charge.dx')