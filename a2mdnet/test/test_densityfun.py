import torch
import time
import os
from a2mdnet.modules import QMDensityFun
from a2mdio.qm import WaveFunction, WaveFunctionHDF5
from a2mdtest.a2mdtests import methane, benzene


if __name__ == '__main__':

    print("---")
    print("-- qm density fun test")

    device = torch.device('cuda:0')
    start = time.time()
    benzene_wfn = WaveFunction.from_file(benzene.wfn, program='g09')
    methane_wfn = WaveFunction.from_file(methane.wfn, program='g09')
    wfnh5 = WaveFunctionHDF5('.wfn.h5py', mode='w')
    wfnh5.add(key='benzene', wfn=benzene_wfn)
    wfnh5.add(key='methane', wfn=methane_wfn)
    wfnh5.close()
    end = time.time()
    print("-- {:12.4f}".format(end - start))
    print("-- done!")

    load_qm = lambda cwfn : QMDensityFun(group=cwfn, dtype=torch.float, device=device)
    wfnh5 = WaveFunctionHDF5('.wfn.h5py', mode='r', wfn_init=load_qm)

    for _, torch_wfn in wfnh5.iterall():
        for i in range(10):
            start = time.time()
            x = torch.rand(100000, 3, dtype=torch.float, device=device)
            torch_res = torch_wfn.forward(x)
            end = time.time()
            print("TE : {:12.6f}".format(end - start))

    wfnh5.close()
    os.remove('.wfn.h5py')
    print("DONE!")

    print("---")
    print("-- advanced test")
    wfnh5 = WaveFunctionHDF5('C:/scratch/.aaa.wfn.hdf5py', mode='r', wfn_init=load_qm)
    for _, torch_wfn in wfnh5.iterall():
        for i in range(1):
            start = time.time()
            x = torch.rand(10000, 3, dtype=torch.float, device=device)
            torch_res = torch_wfn.forward(x)
            end = time.time()
            print("TE : {:12.6f}".format(end - start))
