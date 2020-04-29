from a2md.integrate import polar_integral, split_space
from a2mdtest.a2mdtests import water
from a2mdio.molecules import Mol2
from a2md.models import a2md_from_mol
import numpy as np
import time
if __name__ == '__main__':
    START = time.time()
    fun = lambda x : np.exp(-2 * np.linalg.norm(x, axis=1))/np.pi
    integral = polar_integral(fun, r_max=10.0, radial_res=200, angular_res=20)
    print(integral)
    print("done!")
    print("TE = {:8.4f}".format(time.time() - START))

    START = time.time()
    water_molecule = Mol2(water.mol2)
    water_dm = a2md_from_mol(water_molecule)
    water_dm.parametrize()
    surf = np.loadtxt(water.surfaces[1], skiprows=1, delimiter=',')
    water_dm.optimize(surf[:, :3], surf[:, 3], optimization_mode='unrestricted')

    fun = lambda x: water_dm.eval(x)
    integral = 0

    for fx in split_space(water_molecule, fun):
        integral += polar_integral(
            fun=fx,
            r_max=15.0,
            radial_res=100,
            angular_res=20
        )

    print(integral)
    print("done")
    print("TE = {:8.4f}".format(time.time() - START))
