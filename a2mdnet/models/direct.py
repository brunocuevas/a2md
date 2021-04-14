from a2mdnet import modules
# from a2mdnet.functions import APEV
from a2mdnet.a2mdt.modules import A2MDt
from a2mdnet.a2mdt import A2MD_MODEL
import json
import torch
import torch.nn as nn
import warnings
warnings.filterwarnings('ignore')
DEVICE = torch.device('cuda:0')
DEFAULT_ARCHITECTURE = dict(
    common_net=[384, 128, 112, 96, 36, 6],
    atom_net=[6, 2],
    bond_net=[22, 2],
    subnet=0
)
PAIR_FEATURES = dict(
    EtaR = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
    ShfR = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5],
    Rc = 6.0
)
LR = 1e-5
EPOCHS = 100
BETAS = (0.9, 0.999)
WEIGHT_DECAY = 1e-2


class DensityCoupled(nn.Module):

    def __init__(self, architecture, device):
        from a2mdnet.a2mdt.modules import A2MDt
        from a2mdnet.a2mdt import A2MD_MODEL
        super(DensityCoupled, self).__init__()

        self.device = device
        try:

            common_net_architecture = architecture['common_net']
            atom_net_architecture = architecture['atom_net']
            bond_net_architecture = architecture['bond_net']
            feature_extraction_layer = architecture['fe_layer']
            feature_extraction_net = architecture['fe_net']
            coefficients_distribution = architecture['coefficients_distribution']

        except KeyError:
            raise IOError("architecture dict did not contain the right fields")

        self.common_net = modules.TorchElementSpecificA2MDNN(
            nodes=common_net_architecture, elements=[1, 6, 7, 8], device=self.device
        )
        self.atom_net = modules.TorchElementSpecificA2MDNN(
            nodes=atom_net_architecture, elements=[1, 6, 7, 8], device=self.device,
            distribution=coefficients_distribution['isotropic']
        )
        self.bond_net = modules.TorchPairSpecificA2MDNN(
            nodes=bond_net_architecture, elements=[1, 6, 7, 8], device=self.device,
            distribution=coefficients_distribution['anisotropic']
        )

        if feature_extraction_layer is None:
            self.feature_extractor = modules.SymFeats()
        else:
            self.feature_extractor = modules.TorchaniFeats(
                net=feature_extraction_net, feats_layer=feature_extraction_layer,
                device=self.device
            )

        for param in self.feature_extractor.parameters():
            param.requires_grad = False


        self.normalizer = modules.QPChargeNormalization(device=self.device)

        self.density_model = A2MDt(
            params=A2MD_MODEL,
            device=self.device
        )

    @staticmethod
    def reverse_connectivity(conn):
        """

        :param conn:
        :return:
        """
        return torch.flip(conn, dims=(2, ))

    def forward_common_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, coordinates = x

        labels, ta_feats = self.feature_extractor(labels, coordinates)
        labels, common_layer = self.common_net(labels, ta_feats)

        return labels, common_layer

    def forward_atom_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, common_layer = x

        labels, atom_targets = self.atom_net(labels, common_layer)

        return labels, atom_targets


    def forward_bond_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, connectivity, common_layer = x

        reverse_con = self.reverse_connectivity(connectivity)

        labels, connectivity, bond_targets_forward = self.bond_net(labels, connectivity, common_layer)
        labels, connectivity, bond_targets_reverse = self.bond_net(
            labels, reverse_con, common_layer
        )

        bond_targets = torch.cat([bond_targets_forward, bond_targets_reverse], dim=2)

        return labels, connectivity, bond_targets

    def forward_coefficients(self, *x):

        labels, connectivity, coordinates, charges, int_iso, int_aniso = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)
        atom_targets, bond_targets = self.normalizer.forward(
            atom_targets, bond_targets, int_iso, int_aniso, charges
        )
        return labels, connectivity, atom_targets, bond_targets

    def forward(self, *x):
        """

        :param x:
        :return:
        """

        labels, connectivity, coordinates, charges, int_iso, int_aniso, sampling_coords = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)
        atom_targets, bond_targets = self.normalizer.forward(
            atom_targets,
            bond_targets,
            int_iso,
            int_aniso,
            charges
        )
        assert type(self.density_model) is A2MDt
        density_prediction = self.density_model.forward(
            labels, connectivity, coordinates, atom_targets, bond_targets, sampling_coords
        )

        return labels, connectivity, density_prediction, atom_targets, bond_targets

    def get_atom_feats(self, *x):
        """

        :param x:
        :return:
        """
        labels, coordinates = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        return labels, atom_targets

    def get_bond_feats(self, *x):
        """

        :param x:
        :return:
        """
        labels, connectivity, coordinates = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, connectivity, atom_targets = self.forward_bond_net(labels, connectivity, common_layer)
        return labels, connectivity, atom_targets

    def to(self, new_device):
        # Following https://stackoverflow.com/questions/54706146/moving-member-tensors-with-module-to-in-pytorch
        new_self = super(DensityCoupled, self).to(new_device)
        new_self.device = new_device
        new_self.normalizer.device = new_device
        new_self.feature_extractor.device = new_device
        new_self.atom_net.device = new_device
        new_self.bond_net.device = new_device
        new_self.common_net.device = new_device
        new_self.density_model.device = new_device
        return new_self


class MultTaskDensityCoupled(nn.Module):

    def __init__(self, architecture, device):
        from a2mdnet.a2mdt.modules import A2MDt
        super(MultTaskDensityCoupled, self).__init__()

        self.device = device
        try:

            common_net_architecture = architecture['common_net']
            atom_net_architecture = architecture['atom_net']
            mol_net_architecture = architecture['mol_net']
            bond_net_architecture = architecture['bond_net']
            feature_extraction_layer = architecture['fe_layer']
            feature_extraction_net = architecture['fe_net']
            density_function_iso_order = architecture['iso_order']
            density_function_aniso_order = architecture['aniso_order']
            coefficients_distribution = architecture['coefficients_distribution']

        except KeyError:
            raise IOError("architecture dict did not contain the right fields")

        self.common_net = modules.TorchElementSpecificA2MDNN(
            nodes=common_net_architecture, elements=[1, 6, 7, 8], device=self.device
        )

        self.mol_net = modules.TorchElementSpecificA2MDNN(
            nodes=mol_net_architecture, elements=[1, 6, 7, 8], device=self.device,
        )

        self.atom_net = modules.TorchElementSpecificA2MDNN(
            nodes=atom_net_architecture, elements=[1, 6, 7, 8], device=self.device,
            distribution=coefficients_distribution['isotropic']
        )

        self.bond_net = modules.TorchPairSpecificA2MDNN(
            nodes=bond_net_architecture, elements=[1, 6, 7, 8], device=self.device,
            distribution=coefficients_distribution['anisotropic']
        )

        if feature_extraction_layer is None:

            self.feature_extractor = modules.SymFeats()

        else:

            self.feature_extractor = modules.TorchaniFeats(
                net=feature_extraction_net,
                feats_layer=feature_extraction_layer,
                device=self.device
            )

        for param in self.feature_extractor.parameters():
            param.requires_grad = False

        self.normalizer = modules.QPChargeNormalization(device=self.device)
        self.density_model = A2MDt(
            params=A2MD_MODEL,
            device=self.device
        )

    @staticmethod
    def reverse_connectivity(conn):
        """

        :param conn:
        :return:
        """
        return torch.flip(conn, dims=(2, ))

    def forward_common_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, coordinates = x

        labels, ta_feats = self.feature_extractor(labels, coordinates)
        labels, common_layer = self.common_net(labels, ta_feats)

        return labels, common_layer

    def forward_atom_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, common_layer = x

        labels, atom_targets = self.atom_net(labels, common_layer)

        return labels, atom_targets

    def forward_bond_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, connectivity, common_layer = x

        reverse_con = self.reverse_connectivity(connectivity)

        labels, connectivity, bond_targets_forward = self.bond_net(labels, connectivity, common_layer)
        labels, connectivity, bond_targets_reverse = self.bond_net(
            labels, reverse_con, common_layer
        )

        bond_targets = torch.cat([bond_targets_forward, bond_targets_reverse], dim=2)

        return labels, connectivity, bond_targets

    def forward_coefficients(self, *x):
        """

        :param x:
        :return:
        """

        labels, connectivity, coordinates, charges, int_iso, int_aniso = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)
        atom_targets, bond_targets = self.normalizer.forward(
            atom_targets,
            bond_targets,
            int_iso,
            int_aniso,
            charges
        )
        return labels, connectivity, atom_targets, bond_targets

    def forward_mol_property(self, *x):
        """

        :param x:
        :return:
        """
        labels, connectivity, coordinates, charges, int_iso, int_aniso = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, mol_targets = self.forward_mol_net(labels, common_layer)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)
        atom_targets, bond_targets = self.normalizer.forward(
            atom_targets,
            bond_targets,
            int_iso,
            int_aniso,
            charges
        )
        return labels, connectivity, atom_targets, bond_targets, mol_targets

    def forward_mol_net(self, *x):
        """

        :param x:
        :return:
        """
        labels, common_layer = x

        labels, mol_prop = self.mol_net(labels, common_layer)

        return labels, mol_prop.sum(2).sum(1)


    def forward(self, *x):
        """

        :param x:
        :return:
        """

        labels, connectivity, coordinates, charges, int_iso, int_aniso, sampling_coords = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, mol_targets = self.forward_mol_net(labels, common_layer)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)
        atom_targets, bond_targets = self.normalizer.forward(
            atom_targets,
            bond_targets,
            int_iso,
            int_aniso,
            charges
        )
        assert type(self.density_model) is A2MDt
        density_prediction = self.density_model.forward(
            labels, connectivity, coordinates, atom_targets, bond_targets, sampling_coords
        )

        return labels, connectivity, density_prediction, atom_targets, bond_targets, mol_targets

    def get_atom_feats(self, *x):
        """

        :param x:
        :return:
        """
        labels, coordinates = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        return labels, atom_targets

    def get_bond_feats(self, *x):
        """

        :param x:
        :return:
        """
        labels, connectivity, coordinates = x
        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, connectivity, atom_targets = self.forward_bond_net(labels, connectivity, common_layer)
        return labels, connectivity, atom_targets

    def to(self, new_device):
        # Following https://stackoverflow.com/questions/54706146/moving-member-tensors-with-module-to-in-pytorch
        new_self = super(MultTaskDensityCoupled, self).to(new_device)
        new_self.device = new_device
        new_self.normalizer.device = new_device
        new_self.feature_extractor.device = new_device
        new_self.atom_net.device = new_device
        new_self.bond_net.device = new_device
        new_self.common_net.device = new_device
        new_self.density_model.device = new_device
        return new_self