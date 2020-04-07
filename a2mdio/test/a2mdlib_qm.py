from a2mdio.qm import WaveFunction, GaussianLog, WaveFunctionGPU
from a2mdtest.a2mdtests import methane
import torch
import numpy as np
import time
if __name__ == '__main__':
    print("testing to read a wfn")

    random_coords = np.random.rand(10000, 3)

    wfn = WaveFunction(
        verbose=True,
        file=methane.wfn,
        batch_size=10000000
    )

    wfn_gpu = WaveFunctionGPU(
        file=methane.wfn,
        device='cuda:0', dtype=torch.double
    )

    sampled_gpu = torch.tensor(random_coords, device=torch.device('cuda:0'), dtype=torch.double)

    T1 = time.time()
    ro_numpy = torch.from_numpy(wfn.eval(sampled_gpu.data.cpu().numpy()))
    T2 = time.time()

    T3 = time.time()
    ro_torch = wfn_gpu.eval(sampled_gpu)
    T4 = time.time()

    print("numpy : {:8.4f}".format(T2 - T1))
    print("torch : {:8.4f}".format(T4 - T3))

    ro_numpy = torch.tensor(ro_numpy, dtype=torch.double, device=torch.device('cuda:0'))
    print("RMSE {:12.4e}".format((ro_torch - ro_numpy).pow(2.0).sum().div(ro_torch.size(0)).sqrt().item()))