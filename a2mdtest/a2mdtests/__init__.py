import json
from pathlib import Path

path2files = Path(__file__).parent.parent

with open(Path(__file__).parent.parent / 'a2mdtests/gdb_test_files.json') as f:
    gdb_test = json.load(f)

for k, item in gdb_test.items():
    item['path'] = path2files / item['path']


class A2mdTest:
    def __init__(self, name, reference):
        self.name = name
        self.path = path2files / reference['path']
        self.mol2 = self.path / reference['mol2']
        self.out = self.path / reference['out']
        self.pdb = self.path / reference['pdb']
        self.surfaces = [self.path / i for i in reference['surfaces'].keys()]
        self.wfn = self.path / reference['wfn']
    def __str__(self):
        return "a2md test, {:s}".format(self.name)

benzene = A2mdTest("benzene", gdb_test['gdb_000214'])
methane = A2mdTest("methane", gdb_test['gdb_000001'])
ammonium = A2mdTest("ammonium", gdb_test['gdb_000002'])
water = A2mdTest("water", gdb_test['gdb_000003'])
methanethiol = A2mdTest("methanethiol", gdb_test['methanethiol'])
all20peptide = A2mdTest("all20peptide", gdb_test['all20peptide'])
aca = A2mdTest("aca", gdb_test['aca'])
prup3 = A2mdTest("prup3", gdb_test['prup3'])
# ethane = A2mdTest("ethane", gdb_test['gdb_000007'])
#propanotriol = A2mdTest("propanotriol", gdb_test['gdb_000397'])