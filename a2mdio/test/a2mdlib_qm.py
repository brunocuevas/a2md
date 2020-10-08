from a2mdio.qm import GaussianLog
from a2mdio.parsers import forces
from a2mdtest.a2mdtests import aca
from a2mdio.qm import WaveFunction
from a2mdtest.a2mdtests import ammonium
from a2mdio.wfx import WaveFunctionX
from a2mdio.utils import eval_volume
import numpy as np
if __name__ == '__main__':
    sample = np.loadtxt(ammonium.surfaces[0], delimiter=',', skiprows=1)
    r = sample[:, :3]
    p = sample[:, 3]
    wfn = WaveFunction.from_file(ammonium.wfn, program='g09', prefetch_dm=True)
    pr = wfn.eval(r)
    l2 = ((p - pr) ** 2).sum()
    print("-- l2 = {:6.4f}".format(l2))
    print("done!")

    # gl = GaussianLog(
    #     file='urea_c_000968.g09.output', method='dft-B3LYP',
    #     charges='MK', ep=False
    # )
    # c, f = gl.seek(forces)
    # for (x, y, z) in f:
    #     print("{:8.4f} {:8.4f} {:8.4f}".format(x, y, z))
    # print("DONE!")
    # out_dict = gl.read()
    #
    # gl = GaussianLog(
    #     file=aca.out, method='MP2',
    #     charges='NPA', ep=False
    # )
    # out_dict = gl.read()
    #
    # print('DONE!')
    # wfx = WaveFunctionX(file='crambin.dft.orca.molden.wfx')