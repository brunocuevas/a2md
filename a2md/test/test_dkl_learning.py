import numpy as np
from a2md.integrate import integrate_density_functional_gradient, dkl_gradient_functional, kullback_leibler_functional
from a2md.integrate import integrate_density_functional
from a2md.models import a2md_from_mol
from a2md.utils import RBFSymmetryCluster
from a2mdio.molecules import Mol2
from a2mdio.qm import WaveFunction
from a2mdtest.a2mdtests import water
from a2md.utils import project
from scipy.optimize import minimize


if __name__ == '__main__':

    water_mol2 = Mol2(water.mol2)
    water_a2md = a2md_from_mol(water_mol2)
    water_wfn = WaveFunction(file=water.wfn)
    water_density_sample = np.loadtxt(water.surfaces[1], skiprows=1, delimiter=',')
    cm = RBFSymmetryCluster(verbose=False)
    water_a2md.parametrize()
    water_a2md.clustering(cm.cluster)
    water_a2md.optimize(
        training_coordinates=water_density_sample[:, :3],
        training_density=water_density_sample[:, 3],
        optimization_mode='restricted'
    )


    n = water_a2md.get_number_optimizable_functions()
    x = water_a2md.get_unfrozen_integrals()
    q = np.array(water_a2md.atom_charges).sum() - water_a2md.get_frozen_integrals().sum()




    for i in range(10):

        dkl_fun = kullback_leibler_functional(water_wfn.eval, water_a2md.eval)
        gdkl_fun = dkl_gradient_functional(water_wfn.eval, water_a2md)
        dkl = integrate_density_functional(dkl_fun, water_mol2)
        gdkl = integrate_density_functional_gradient(gdkl_fun, water_mol2, nfuns=n, res=50)
        gdkl = project(gdkl, x)
        print("{:12.6e} {:12.6e}".format(
            dkl, np.linalg.norm(gdkl)
        ))

        c = water_a2md.get_opt_coefficients()
        cp = c - 0.1 * gdkl
        check = (cp * x).sum()
        water_a2md.set_opt_coefficients(cp)


    water_a2md.inflate()

    dx = water_a2md.eval_volume(spacing=3.0, resolution=0.2)
    dx.write("dkl.dx")
    print("DONE")