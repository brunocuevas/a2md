import torch
import h5py
from a2mdnet.utils import IntegrationGrid
from a2mdnet.density_models import GenAMD
from a2mdnet.modules import QMDensityFun
from a2mdnet.data import Coordinates
from a2mdio.params import AMDParameters
from a2mdio.units import bohr
from a2mdtest.a2mdtests import benzene
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import math

TEST_TESELLATION = False
TEST_GRID = False
TEST_AMD_INT = False
TEST_QM_INT = True

if __name__ == '__main__':

    ig = IntegrationGrid(
        device=torch.device('cuda:0'), dtype=torch.float, grid='coarse', radial_resolution=25
    )

    if TEST_TESELLATION:

        print('-- tesellation test for a two centers body')
        r = torch.tensor([[[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]]], dtype=torch.float, device=torch.device('cuda:0'))

        x = torch.stack(
            [
                torch.zeros(1000, dtype=torch.float),
                torch.zeros(1000, dtype=torch.float),
                ((torch.arange(1000, dtype=torch.float) / 1000) - 0.5)*2.0
            ], dim=1
        ).to(torch.device('cuda:0'))
        x = x.unsqueeze(0)
        i = torch.zeros(1, 1000, dtype=torch.long, device=torch.device('cuda:0'))

        r = Coordinates(values=r, units=bohr)
        r.distance_matrix_()
        x = Coordinates(values=x, units=bohr)

        w = ig.becke_tesellation(x, r.get_cartessian(unit=bohr), r.get_distance_matrix(unit=bohr), i)

        plt.plot(w[0, :].data.cpu().numpy())
        plt.show()
        print('-- done')

    if TEST_GRID:

        print('-- integration concentric spheres')

        r = Coordinates([[0.0, 0.0, 0.0]], dtype=torch.float, device=torch.device, units=bohr)
        r.distance_matrix_()
        _, z, w = ig.spheric_integration_grid()

        amdp = AMDParameters.from_file('../a2mdt/test/spherical.json')
        gamd = GenAMD(amdp, device=torch.device('cuda:0'), dtype=torch.float)

        # HYDROGEN

        print('-- \t integrating a hydrogen AMD')
        labels = torch.zeros(1, 1, dtype=torch.long, device=torch.device('cuda:0'))
        p = gamd.protodensity(z, labels, r)
        integral = (p * w).sum() - gamd.integrate(labels).sum()
        print('-- \t\t integral error : {:f}'.format(integral.item()))

        # CARBON

        print('-- \t integrating a carbon AMD')
        labels = torch.tensor([[2]], dtype=torch.long, device=torch.device('cuda:0'))
        p = gamd.protodensity(z, labels, r)
        integral = (p * w).sum() - gamd.integrate(labels).sum()
        print('-- \t\t integral error : {:f}'.format(integral.item()))

        del z, w

    if TEST_AMD_INT:

        print('-- integration grid for a 32 center body')
        r = Coordinates(torch.rand(1, 32, 3), dtype=torch.float, device=torch.device('cuda:0'), units=bohr)
        l = torch.zeros(1, 32, dtype=torch.long, device=torch.device('cuda:0'))
        start = time.time()
        z, w = ig.integration_grid(coords=r)
        p = gamd.protodensity(z, l, r)
        integral = (p * w).sum() - 32.0
        end = time.time()
        print('-- \t time ellapsed : {:12.6f} s'.format(end - start))
        print('-- \t integration error : {:12.6f}'.format(integral))
        print('-- done!')

    if TEST_QM_INT:

        with h5py.File(benzene.path / 'benzene.wfn.hdf5', 'r') as f:
            benzene_qm = QMDensityFun(f['benzene.wfn'], dtype=torch.float, device=torch.device('cuda:0'))
        r = benzene_qm.coords
        r = Coordinates(r.unsqueeze(0), units=bohr)
        r.distance_matrix_()
        start = time.time()
        z, w = ig.integration_grid(coords=r)
        p = benzene_qm(
            z.get_cartessian(unit=bohr).squeeze(0)
        )
        integral = (p * w).sum() - 42.0
        end = time.time()
        print('-- \t time ellapsed : {:12.6f} s'.format(end - start))
        print('-- \t integration error : {:12.6f}'.format(integral))
        print('-- done!')
        p = benzene_qm(z.get_cartessian(unit=bohr).squeeze(0)).pow(4.0 / 3.0)
        integral = (p * w).sum()
        print('-- \t dirac exchange : {:12.6f}'.format(integral))
