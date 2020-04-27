from a2md.models import a2md_from_mol
from a2md.utils import RBFSymmetryCluster
from a2mdio.molecules import Mol2
from a2mdtest.a2mdtests import benzene
import json
import numpy as np
if __name__ == "__main__":

    print("a2md/molecules")
    print("---")
    m = Mol2(file=benzene.mol2)

    dm = a2md_from_mol(m)


    with open(r"C:\Users\Bruno\ownCloud\projects\a2md\a2md\parameters\a2md_topo_bonded_model_extended.json") as f:
        dm.parametrize(json.load(f))

    training_coords  = np.loadtxt(benzene.surfaces[2], skiprows=1, delimiter=',')
    training_density = training_coords[:,  3]
    training_coords  = training_coords[:, :3]
    cluster_machine = RBFSymmetryCluster(verbose=False)
    dm.clustering(cluster_machine.cluster)

    c = dm.optimize(
        training_coordinates=training_coords, training_density=training_density,
        optimization_mode='restricted'
    )

    dm.inflate()

    q1 = dm.integrate()
    with open("benzene.ppp", "w") as f:
        json.dump(dm.get_parametrization(), f, indent=4)
    with open("benzene.ppp") as f:
        dm.read(json.load(f))
    q2 = dm.integrate()
    print("{:4.3f}, {:4.3f}".format(q1, q2))
    print("DONE!")
