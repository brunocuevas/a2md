from a2mdnet import modules
import torch
import torch.nn as nn

DEVICE = torch.device('cpu')
DEFAULT_ARCHITECTURE = dict(
    common_net=[96, 36, 6],
    atom_net=[6, 2],
    bond_net=[12, 2],
    subnet=0
)
LR = 1e-5
EPOCHS = 100
BETAS = (0.9, 0.999)
WEIGHT_DECAY = 1e-2


class A2MDnetDensity(nn.Module):

    def __init__(self, architecture, normalize=False):

        super(A2MDnetDensity, self).__init__()

        # self.normalize = normalize

        try:

            common_net_architecture = architecture['common_net']
            atom_net_architecture = architecture['atom_net']
            bond_net_architecture = architecture['bond_net']
            subnet = architecture['subnet']

        except KeyError:
            raise IOError("architecture dict did not contain the right fields")

        self.common_net = modules.TorchElementSpecificA2MDNN(
            nodes=common_net_architecture, elements=[1, 6, 7, 8], device=DEVICE
        )
        self.atom_net = modules.TorchElementSpecificA2MDNN(
            nodes=atom_net_architecture, elements=[1, 6, 7, 8], device=DEVICE
        )
        self.bond_net = modules.TorchPairSpecificA2MDNN(
            nodes=bond_net_architecture, elements=[1, 6, 7, 8]
        )
        self.feature_extractor = modules.TorchaniFeats(net=subnet, device=DEVICE)

        for param in self.feature_extractor.parameters():
            param.requires_grad = False

        # self.normalizer = modules.ChargeNormalization()

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

    def forward(self, *x):
        """

        :param x:
        :return:
        """

        # charges = None
        # integrals_iso = None
        # integrals_aniso = None

        # if self.normalize:
        #     labels, connectivity, coordinates, charges, integrals_iso, integrals_aniso = x
        # else:
        #     labels, connectivity, coordinates = x

        labels, connectivity, coordinates = x

        labels, common_layer = self.forward_common_net(labels, coordinates)
        labels, atom_targets = self.forward_atom_net(labels, common_layer)
        labels, connectivity, bond_targets = self.forward_bond_net(labels, connectivity, common_layer)

        # if self.normalize:
        #
        #     factor, atom_targets, bond_targets = self.normalizer(
        #         charges,
        #         atom_targets, bond_targets,
        #         integrals_iso, integrals_aniso
        #     )
        #     return labels, connectivity, atom_targets, bond_targets, factor
        # else:

        return labels, connectivity, atom_targets, bond_targets

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


if __name__ == "__main__":

    from a2mdnet.data import CompleteSetDensityParams
    from torch.utils import data
    import torch

    device = torch.device('cpu')
    params = dict(
        batch_size=64,
        shuffle=False,
        num_workers=0
    )

    csdp_t = CompleteSetDensityParams(
        device=device, dtype=torch.float, normalize=True, kind='training', number=10, charges=True, integrals=True
    )

    csdp_t_dl = data.DataLoader(csdp_t, **params)
    model = A2MDnetDensity(architecture=DEFAULT_ARCHITECTURE, normalize=True)

    # sgd_opt = torch.optim.SGD(lr=1e-5, params=model.parameters())
    lr = 1e-4
    gamma = 1e-2
    for l_, c_, coords, charges, int_iso, int_aniso, iso, aniso in csdp_t_dl:

        _, _, iso_targets, aniso_targets, factor = model.forward(l_, c_, coords, charges, int_iso, int_aniso)
        l2_factor = torch.pow(factor - 1, 2).sum() * gamma
        l2 = torch.pow(iso_targets - iso, 2.0).sum() \
             + torch.pow(aniso_targets - aniso, 2.0).sum() \
             + torch.pow(factor - 1, 2).sum() * gamma
        l2.backward()

        # sgd_opt.step()

        with torch.no_grad():
            for p in model.parameters():
                if p.grad is None:
                    pass
                else:
                    p -= lr * p.grad

        print(l2.item())
        model.zero_grad()
