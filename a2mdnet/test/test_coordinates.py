from a2mdio.units import angstrom, bohr, Unit
from a2mdio.molecules import Mol2
from a2mdnet.data import Coordinates
import torch

if __name__ == '__main__':
    m = Unit(fullname='metre', shortname='m', reference='metre', value=1.0)
    x = torch.tensor([[1.0, 1.0, 1.0]], device=torch.device('cuda:0'), dtype=torch.float)
    c = Coordinates(x, units=bohr)
    u = c.get_cartessian(unit=angstrom)
    print(u)

