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


    with open(r"C:\Users\Bruno\ownCloud\projects\a2md\a2md\parameters\a2md_htopo_bonded_model_extended.json") as f:
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

    x = np.random.rand(100, 3)

    q1 = dm.integrate()
    r1 = dm.eval(x)
    ppp = dm.get_parametrization()
    with open('benzene_harmonic.ppp', 'w') as f:
        json.dump(ppp, f, indent=4)

    with open('benzene_harmonic.ppp') as f:
        dm.read(json.load(f))
    r2 = dm.eval(x)
    ppy = dm.get_parametrization()

    test1 = np.equal(r1, r2)
    # with open("benzene.ppp", "w") as f:
    #     json.dump(dm.get_parametrization(), f, indent=4)
    # with open("benzene.ppp") as f:
    #     dm.read(json.load(f))
    # q2 = dm.integrate()
    # print("{:4.3f}, {:4.3f}".format(q1, q2))
    # print("DONE!")

    # dx = dm.eval_volume(spacing=6.0, resolution=0.25, kind='ep')
    # dx.write("benzene_harmonic_ep.dx")
    dx = dm.eval_volume(spacing=4.0, resolution=0.25, kind='density')
    dx.write("benzene_harmonic.dx")