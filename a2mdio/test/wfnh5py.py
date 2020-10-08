from a2mdio.qm import WaveFunction, WaveFunctionHDF5
from a2mdtest.a2mdtests import methane, benzene
import numpy as np
import time
import os
if __name__ == '__main__':
    # STORING

    print("-- storing")
    start = time.time()
    benzene_wfn = WaveFunction.from_file(benzene.wfn, program='g09')
    methane_wfn = WaveFunction.from_file(methane.wfn, program='g09')
    wfnh5 = WaveFunctionHDF5('.wfn.h5py', mode='w')
    wfnh5.add(key='benzene', wfn=benzene_wfn)
    wfnh5.add(key='methane', wfn=methane_wfn)
    wfnh5.close()
    end = time.time()
    print("-- done")
    print("-- TE : {:12.6f}".format(end - start), end="\n\n")
    random_sample = np.random.rand(1000, 3)

    # GET ITEM

    print("-- loading")
    start = time.time()
    wfnh5 = WaveFunctionHDF5('.wfn.h5py', mode='r')
    _, benzene_nwfn = wfnh5['benzene']
    wfnh5.close()
    p = benzene_wfn(random_sample)
    pn = benzene_nwfn(random_sample)
    l2 = np.power(p - pn, 2.0).sum()
    print("l2 : {:12.4f}".format(l2))
    print("-- done")
    end = time.time()
    print("-- TE : {:12.6f}".format(end - start), end="\n\n")

    # ITERATE

    print("-- loading")
    start = time.time()
    wfnh5 = WaveFunctionHDF5('.wfn.h5py', mode='r')
    for key, nwfn in wfnh5.iterall():
        nwfn(random_sample)
        print("--- {:s}".format(key))
    wfnh5.close()
    print("-- done")
    end = time.time()
    print("-- TE : {:12.6f}".format(end - start), end="\n\n")

    os.remove('.wfn.h5py')
