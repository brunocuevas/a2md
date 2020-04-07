from a2md.utils import RBFSymmetryCluster
from a2mdio.molecules import Mol2
from a2mdtest.a2mdtests import benzene
from a2md.models import a2md_from_mol2
import numpy as np

if __name__ == "__main__":

    ptl = Mol2(benzene.mol2)
    coords = ptl.get_coordinates(units="angstrom")
    labels = ptl.get_symbols()
    topo = ptl.get_bonds()
    cluster_machine = RBFSymmetryCluster()

    model = a2md_from_mol2(ptl)
    model.parametrize()
    model.clustering(cluster_machine.cluster)

    training = np.loadtxt(benzene.surfaces[0], delimiter=',', skiprows=1)
    training_coords = training[:, :3]
    training_density = training[:, 3]

    model.optimize(
        training_density=training_density,
        training_coordinates=training_coords
    )

    model.inflate()

    inflated_params = model.get_parametrization()
    model.read(inflated_params)

    dx = model.eval_volume(spacing=2.0, resolution=0.25)
    dx.write("benzene.dx")

    print("molecule loaded!")