from a2mdio.qm import GaussianLog
from a2mdio.parsers import forces
from a2mdtest.a2mdtests import aca
from a2mdio.wfx import WaveFunctionX
from a2mdio.utils import eval_volume
import numpy as np
if __name__ == '__main__':

    gl = GaussianLog(
        file='urea_c_000968.g09.output', method='dft-B3LYP',
        charges='MK', ep=False
    )
    c, f = gl.seek(forces)
    for (x, y, z) in f:
        print("{:8.4f} {:8.4f} {:8.4f}".format(x, y, z))
    print("DONE!")
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