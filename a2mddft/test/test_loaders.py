from a2mddft.modules import WaveFunctionReader, DensityFunction, MLDensityFunctional
from a2md.integrate import integrate_density_functional
from a2mdio.molecules import Mol2
from a2mdio.qm import WaveFunction
import torch
import time
import numpy as np
example_wfn =[ "F:/aniset/wfn/rep_01/" + i for i in [
    "ani_000001_0000.wfn",
    "ani_000001_0001.wfn",
    "ani_000001_0002.wfn",
    "ani_000001_0003.wfn",
    "ani_000001_0004.wfn",
    "ani_000002_0000.wfn",
    "ani_000002_0001.wfn",
    "ani_000002_0002.wfn",
    "ani_000002_0003.wfn",
    "ani_000002_0004.wfn"
]]

if __name__ == '__main__':
    mm = Mol2(file="F:/aniset/mol2/ani_000001_0000.mol2")
    wfn = WaveFunction(file=example_wfn[0])

    x = integrate_density_functional(
        functional=lambda r: wfn.eval(r),
        mol=mm, grid='coarse', res=50
    )
    print(x)

    # wfnr = WaveFunctionReader(example_wfn)
    # wfnd = DensityFunction(device=torch.device('cuda:0'))
    # mldf = MLDensityFunctional(architecture=[2, 25, 25, 1], device=torch.device('cuda:0'))
    # r = np.zeros((2, 3), dtype='float64')
    # z = integrate_density_functional(
    #     functional= lambda x: np.ones(x.shape[0])
    # )
    #
    # for r, c, xp, s, dm in wfnr.generate_wfns():
    #
    #
    #     wfnd.parametrize(r, c, xp, s, dm)
    #     x, w = mldf.generate_lebdenev_grid(r)
    #     p, _ = wfnd.forward(x)
    #     z = (p * w).sum()
    #     print(z)
    #     break
    #     # p, _ = wfnd.forward(x)
    #     # z = (u*p).sum()
    #     # end = time.time()
    #     # print("{:8.4f}".format(end - start_))
    # print("DONE!")
