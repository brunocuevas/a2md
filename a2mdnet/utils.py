from a2mdio.molecules import Mol2
from a2mdnet.data import convert_label2tensor
from a2md.utils import integrate_from_dict
from a2mdnet.data import match_fun_names
import warnings
import torch


def get_charges(l, t, i_iso, i_aniso, iso, aniso, device):
    """

    :param l:
    :param t:
    :param i_iso:
    :param i_aniso:
    :param iso:
    :param aniso:
    :return:
    """

    # Defining problems
    n_batch = l.size(0)
    n_atoms = l.size(1)
    n_bond = t.size(1)
    n_bond_funs = int(i_aniso.size(2)/2)

    # Calculating function charges
    charge_iso = (i_iso * iso)
    charge_aniso = (i_aniso * aniso)
    q_iso = charge_iso.sum(2)
    # Calculating charges per atom
    # To do so, the topology tensor is expanded to match the shape of the anisotric charge tensor
    # Then a range tensor is added to allow to map index in flat vectors
    # Finally, two operations are employed to map back the anisotropic charge to each atom
    #       1. scatter add. Sums all the positions related to the same index and stores them in that given index
    #       2. masked_scatter. Gets all the non zero values and stores them in the output

    t_r = t.reshape(n_batch, n_bond_funs * n_bond, 1).expand(n_batch, n_bond_funs * n_bond, n_bond_funs).reshape(
        n_batch, n_bond, 2 * n_bond_funs)
    r = (torch.arange(0, n_batch, device=device) * n_atoms).unsqueeze(1).unsqueeze(2).expand(n_batch, n_bond,
                                                                                             2 * n_bond_funs)
    t_r.masked_scatter_(t_r != -1, t_r.masked_select(t_r != -1) + r.masked_select(t_r != -1))
    t_r[t_r == -1] = 0
    q_buffer = torch.zeros_like(charge_aniso).flatten()
    t_rf = t_r.flatten()
    ca_f = charge_aniso.flatten()
    q_buffer = q_buffer.scatter_add(0, t_rf, ca_f)
    q_aniso = torch.zeros_like(q_iso)
    q_aniso = q_aniso.masked_scatter(l != - 1, q_buffer.masked_select(q_buffer != 0.0))

    q_per_atom = q_aniso + q_iso
    return q_per_atom

class Parametrizer:

    def __init__(self, model, device):
        self.model = model
        self.device = device

    def parametrize(self, mol, params):
        """

        :param mol:
        :param params
        :return:
        """
        if not isinstance(mol, Mol2):
            try:
                mol = Mol2(file=mol)
            except ValueError:
                raise IOError("unknown format for mol. Please, use Mol2 or str")
        coords = torch.tensor(mol.get_coordinates(), device=self.device, dtype=torch.float).unsqueeze(0)
        labels = convert_label2tensor(mol.get_atomic_numbers(), device=self.device).unsqueeze(0)
        connectivity = torch.tensor(mol.get_bonds(), dtype=torch.long, device=self.device).unsqueeze(0)
        charge = torch.tensor(mol.charges + mol.atomic_numbers, dtype=torch.float, device=self.device).unsqueeze(0)
        natoms = mol.get_number_atoms()
        nbonds = mol.get_number_bonds()
        connectivity -= 1

        int_iso = torch.zeros(1, natoms, 2)
        int_aniso = torch.zeros(1, nbonds, 4)

        for fun in params:
            center = fun['center']
            funtype, pos = match_fun_names(fun)
            if funtype == 'core':
                charge[0, center] -= integrate_from_dict(fun)
            elif funtype in 'iso':
                int_iso[0, center, pos] = integrate_from_dict(fun)
            elif funtype in 'aniso':
                idx, col = self.match_bond(connectivity, fun)
                int_aniso[0, idx, col + pos] = integrate_from_dict(fun)

        int_iso = int_iso.to(self.device)
        int_aniso = int_aniso.to(self.device)

        _, _, iso_out, aniso_out = self.model.forward_coefficients(labels, connectivity, coords, charge, int_iso, int_aniso)

        if iso_out.is_cuda:
            iso_out = iso_out.squeeze(0).cpu().data.numpy()
            aniso_out = aniso_out.squeeze(0).cpu().data.numpy()
        else:
            iso_out = iso_out.squeeze(0).data.numpy()
            aniso_out = aniso_out.squeeze(0).data.numpy()

        for fun in params:
            center = fun['center']
            funtype, pos = match_fun_names(fun)
            if funtype == 'core':
                continue
            if funtype == 'iso':
                fun['coefficient'] = iso_out[center, pos].item()
            elif funtype == 'aniso':
                idx, col = self.match_bond(connectivity, fun)
                fun['coefficient'] = aniso_out[idx, col + pos].item()

        return params

    @staticmethod
    def match_bond(topology, fun):
        topology = topology.squeeze(0)
        for j in range(topology.shape[0]):
            c1 = (topology[j, 0] == fun['center']) and (
                    topology[j, 1] == fun['bond']
            )
            c2 = (topology[j, 1] == fun['center']) and (
                    topology[j, 0] == fun['bond']
            )
            if c1:

                return j, 0
            if c2:

                return j, 0 + 2

        warnings.warn("bond was not matched")
