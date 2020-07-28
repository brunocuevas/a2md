from a2mdio.molecules import Mol2, QmSetUp
from a2mdtest.a2mdtests import aca

# In this test we expect that:
# - Mol2 files can be read
# - Mol2 files can be written
# - Mol2 files can be read back
# - G09 can be written by using Mol2 coordinates
# - AAMD can be written by using Mol2 coordinates after using a change of units
if __name__ == '__main__':
    # print("test a2md library")
    #
    # # Testing on ANIset
    #
    # qmstp = QmSetUp(
    #     basis='STO-3G', method='B3LYP', calculation_type='single', nprocs=1, output='aca.wfn',
    #     additional_commands=['pop=npa', 'output=wfn', 'density=current']
    # )
    #
    # aca = Mol2(
    #     file=aca.mol2
    # )
    #
    # qmstp.write_g09("aca.g09", aca)
    #
    # # aca.write("aca.mol2")
    #
    # print("DONE!")
    mm = Mol2('F:/solv/a2md/solv_000051.mol2')
    print("DONE!")