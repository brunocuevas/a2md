import torch
import torchani
import torch.nn as nn
from a2mdnet import ELEMENT2NN


class TorchElementA2MDNN(nn.Module):

    def __init__(self, nodes, means=None, stds=None, untrained=False):
        """

        Basic NN. Uses CELU activations as TorchANI. Just requires the number of nodes in each layer

        :param nodes: list of layer terms
        :type nodes: list

        Example:

            TorchElementA2MDNN([384, 128, 96, 2])

        """
        super(TorchElementA2MDNN, self).__init__()

        if means:
            self.means = nn.Parameter(torch.tensor(means), requires_grad=False)
        else:
            self.means = None
        if stds:
            self.stds = nn.Parameter(torch.tensor(stds), requires_grad=False)
        else:
            self.stds = None

        self.layers = nn.ModuleList()

        for i in range(1, len(nodes)-1):

            self.layers.append(nn.Linear(nodes[i-1], nodes[i], bias=True))
            self.layers.append(nn.CELU())

        self.layers.append(nn.Linear(nodes[-2], nodes[-1], bias=True))
        self.n_layers = len(self.layers)
        self.untrained = untrained

    def forward(self, x):
        """

        :param x: input tensor
        :type x: torch.Tensor
        :return: nn output
        :rtype torch.Tensor

        """

        if self.untrained:
            print("WARNING. Possibly untrained nn is being used to make a prediction. Results may be inconsistent")

        x_m = x

        for i in range(self.n_layers):

            x_m = self.layers[i](x_m)

        if self.stds is not None:
            x_m = x_m * self.stds

        if self.means is not None:
            x_m = x_m + self.means

        return x_m


class TorchElementSpecificA2MDNN(nn.Module):

    def __init__(self, nodes, elements, device, distribution=None):
        """

        This module allows to pass the features corresponding to each element to a different and
        specific network.


        :param nodes: list of the number of nodes at each layer
        :param elements: allowed atoms
        :param device: torch device


        """
        from a2mdnet import ELEMENT2SYMBOL
        super(TorchElementSpecificA2MDNN, self).__init__()
        self.allowed_elements = elements
        self.subnets = nn.ModuleList()
        if distribution is None:
            for _ in range(len(elements)):
                self.subnets.append(
                    TorchElementA2MDNN(nodes=nodes).to(device)
                )
        else:
            for i, element in enumerate(sorted(elements)):

                symbol = ELEMENT2SYMBOL[element]
                self.subnets.append(
                    TorchElementA2MDNN(
                        nodes=nodes,
                        means=distribution[symbol]['mean'],
                        stds=distribution[symbol]['std'],
                    ).to(device)
                )

        self.output_size = nodes[-1]
        self.device = device

    def forward(self, *u):
        """

        Input features must be (n, M, f) where.
            - n is the batch size.
            - M is the number of atoms of each batch instance.
            - f is the number of features/attributes
        Labels must be (n,M)


        :param u: tensor containing elements labels, tensor containing coordinates/features
        :return:

        """

        elements, features = u

        assert isinstance(features, torch.Tensor)
        assert isinstance(elements, torch.Tensor)

        if len(features.size()) != 3:
            raise IOError("check the input shape to the features. Should be n,M,f")
        if len(elements.size()) != 2:
            raise IOError("check the input shape of the labels. Should be n, M")

        # Flattening
        elements_ = elements.flatten()
        present_elements = elements
        present_elements = present_elements.unique()
        if present_elements[0].item() == -1:
            present_elements = present_elements[1:]
        present_elements, _ = present_elements.sort()
        features = features.flatten(0, 1)

        # Declaring output tensor
        y = torch.zeros(elements_.size()[0], self.output_size, dtype=torch.float, device=self.device)

        # Iteration through elements
        for i, n in enumerate(present_elements):

            mask_input = (elements_ == n).to(self.device)
            x = features.index_select(0, mask_input.nonzero().squeeze())
            mask_output = mask_input.unsqueeze(1).expand(mask_input.size(0), self.output_size)
            y.masked_scatter_(mask_output, self.subnets[n.item()](x))


        y = y.reshape(elements.size()[0], elements.size()[1], self.output_size)
        return elements, y


class TorchPairSpecificA2MDNN(nn.Module):
    def __init__(self, nodes, elements, device, distribution=None):
        """

        This module organizes input by its connectivity, concatenating the input.

        :param nodes:
        :param elements:
        """
        from a2mdnet import ELEMENT2SYMBOL
        super(TorchPairSpecificA2MDNN, self).__init__()
        self.allowed_elements = elements
        sub = nn.ModuleList()

        if distribution is None:
            for _ in elements:
                for _ in elements:
                    sub.append(TorchElementA2MDNN(nodes=nodes).to(device))
        else:
            for i, element1 in enumerate(sorted(elements)):
                for j, element2 in enumerate(sorted(elements)):
                    symbol1 = ELEMENT2SYMBOL[element1]
                    symbol2 = ELEMENT2SYMBOL[element2]
                    joint_name = '{:s}{:s}'.format(symbol1, symbol2)
                    if joint_name not in distribution.keys():
                        sub.append(TorchElementA2MDNN(nodes=nodes, untrained=True).to(device))
                    else:
                        sub.append(
                            TorchElementA2MDNN(
                                nodes=nodes,
                                means=distribution[joint_name]['mean'],
                                stds=distribution[joint_name]['std'],
                            ).to(device)
                        )

        self.subnets = sub
        self.nodes = nodes
        self.output_size = nodes[-1]
        self.device = device

    def forward(self, *u):
        """

        :param u:
        :return:
        """

        elements, connectivity, features = u

        assert isinstance(elements, torch.Tensor)
        assert isinstance(connectivity, torch.Tensor)
        assert isinstance(features, torch.Tensor)

        n_connectivity = connectivity.size()[1]
        n_batch = connectivity.size()[0]
        n_atoms = features.size(1)
        n_elements = len(self.allowed_elements)
        n_batch_arange = torch.arange(n_batch, dtype=torch.long, device=self.device) * n_atoms
        connectivity_f = connectivity.flatten(0, 1)
        features_f = features.flatten(0, 1)
        elements_f = elements.flatten(0, 1).unsqueeze(1)

        # Updating values on C to match the new indexing
        con_mask = connectivity_f[:, 0] != -1
        con_mask_features = con_mask.unsqueeze(1).expand(connectivity_f.size(0), 2*features_f.size(1))
        con_mask_elements = con_mask.unsqueeze(1).expand(connectivity_f.size(0), 2)
        connectivity_f = connectivity_f + n_batch_arange.reshape(-1, 1).repeat(1, n_connectivity).reshape(-1, 1).repeat(1, 2)

        # Creating combinations of the features and the elements

        features_c = torch.zeros(connectivity_f.size(0), 2 * features_f.size(1), dtype=torch.float, device=self.device)
        elements_c = torch.ones(connectivity_f.size(0), 2, dtype=torch.long, device=self.device) * -1

        features_c.masked_scatter_(
            con_mask_features,
            torch.cat(
                (
                    features_f.index_select(dim=0, index=connectivity_f[con_mask, 0]),
                    features_f.index_select(dim=0, index=connectivity_f[con_mask, 1]),
                ), dim=1
            )
        )

        elements_c.masked_scatter_(
            con_mask_elements,
            torch.cat(
                (
                    elements_f.index_select(dim=0, index=connectivity_f[con_mask, 0]),
                    elements_f.index_select(dim=0, index=connectivity_f[con_mask, 1]),
                ), dim=1
            )
        )

        # finding the present bonds, as done by torchani people

        nn_index = elements_c[:, 0] * n_elements
        nn_index += elements_c[:, 1]

        present_nn_index = nn_index.unique(sorted=True)

        if present_nn_index[0].item() == -5:
            present_nn_index = present_nn_index[1:]

        y = torch.zeros(features_c.size()[0], self.output_size, dtype=torch.float, device=self.device)

        for i in present_nn_index:

            mask_input = (nn_index == i)
            x = features_c.index_select(0, mask_input.nonzero().squeeze())
            mask_output = mask_input.unsqueeze(1).expand(mask_input.size(0), self.output_size)
            y.masked_scatter_(mask_output, self.subnets[i](x).squeeze())

        y = y.reshape(n_batch, n_connectivity, self.output_size)

        return elements, connectivity, y

class TorchPairDistances(nn.Module):
    def __init__(self, nodes, elements, device):
        """

        This module organizes input by its connectivity, concatenating the input.

        :param nodes:
        :param elements:
        """
        super(TorchPairDistances, self).__init__()
        self.allowed_elements = elements
        sub = nn.ModuleList()

        for _ in elements:
            for _ in elements:
                sub.append(TorchElementA2MDNN(nodes=nodes).to(device))

        self.subnets = sub
        self.nodes = nodes
        self.output_size = nodes[-1]
        self.device = device

    def forward(self, *u):
        """

        :param u:
        :return:
        """

        elements, connectivity, sym_features, pair_features  = u

        assert isinstance(elements, torch.Tensor)
        assert isinstance(connectivity, torch.Tensor)
        assert isinstance(sym_features, torch.Tensor)
        assert isinstance(pair_features, torch.Tensor)

        n_connectivity = connectivity.size()[1]
        n_batch = connectivity.size()[0]
        n_atoms = sym_features.size(1)
        n_elements = len(self.allowed_elements)
        n_batch_arange = torch.arange(n_batch, dtype=torch.long, device=self.device) * n_atoms

        connectivity_f = connectivity.flatten(0, 1)
        features_f = sym_features.flatten(0, 1)
        elements_f = elements.flatten(0, 1).unsqueeze(1)
        pair_f = pair_features.flatten(0, 1)

        # Updating values on C to match the new indexing
        con_mask = connectivity_f[:, 0] != -1
        con_mask_features = con_mask.unsqueeze(1).expand(connectivity_f.size(0), 2*features_f.size(1))
        con_mask_elements = con_mask.unsqueeze(1).expand(connectivity_f.size(0), 2)
        connectivity_f = connectivity_f + n_batch_arange.reshape(-1, 1).repeat(
            1, n_connectivity
        ).reshape(-1, 1).repeat(1, 2)

        # Creating combinations of the features and the elements

        features_c = torch.zeros(connectivity_f.size(0), 2 * features_f.size(1), dtype=torch.float, device=self.device)
        elements_c = torch.ones(connectivity_f.size(0), 2, dtype=torch.long, device=self.device) * -1

        features_c.masked_scatter_(
            con_mask_features,
            torch.cat(
                (
                    features_f.index_select(dim=0, index=connectivity_f[con_mask, 0]),
                    features_f.index_select(dim=0, index=connectivity_f[con_mask, 1]),
                ), dim=1
            )
        )

        elements_c.masked_scatter_(
            con_mask_elements,
            torch.cat(
                (
                    elements_f.index_select(dim=0, index=connectivity_f[con_mask, 0]),
                    elements_f.index_select(dim=0, index=connectivity_f[con_mask, 1]),
                ), dim=1
            )
        )

        # including pairs

        features_cp = torch.cat((features_c, pair_f), dim=1)

        # finding the present bonds, as done by torchani people

        nn_index = elements_c[:, 0] * n_elements
        nn_index += elements_c[:, 1]

        present_nn_index = nn_index.unique(sorted=True)

        if present_nn_index[0].item() == -5:
            present_nn_index = present_nn_index[1:]

        y = torch.zeros(features_c.size()[0], self.output_size, dtype=torch.float, device=self.device)

        for i in present_nn_index:

            mask_input = (nn_index == i)
            x = features_cp.index_select(0, mask_input.nonzero().squeeze())
            mask_output = mask_input.unsqueeze(1).expand(connectivity_f.size(0), self.output_size)
            y.masked_scatter_(mask_output, self.subnets[i](x).squeeze())

        y = y.reshape(n_batch, n_connectivity, self.output_size)

        return elements, connectivity, y

class TorchaniFeats(nn.Module):
    def __init__(self, net=0, feats_layer=6, device=torch.device('cpu')):
        """

        Torch ANI Feats

        The main purpose of this class is to extract the internal layers from TorchANI
        network to use them as features in a feature-extraction paradigm.

        :param net: TorchANIx is an ensamble of 8 networks. Specify which.
        :param feats_layer: Uses that layer as features.
        """
        if net > 7:
            raise IOError("the TorchANI ensamble has only 8 networks")
        super(TorchaniFeats, self).__init__()
        torchani_main = torchani.models.ANI1x()
        self.torchani_model_sf = torchani_main[net][0]
        self.torchani_model_nn = torchani_main[net][1]
        self.feats_layer = feats_layer
        if feats_layer not in [5, 6]:
            raise NotImplementedError("feature extraction at lower layers would lead to different input sizes")
        self.device = device

    def forward(self, *x):
        labels_tensor = x[0]
        coords_tensor = x[1]
        # print(self.torchani_model_sf.sizes)
        species, aevs = self.torchani_model_sf((labels_tensor, coords_tensor))

        n_batch = species.size(0)
        n_atoms = species.size(1)

        try:
            n_feats = self.torchani_model_nn[0][self.feats_layer - 1].out_features
        except AttributeError:
            n_feats = self.torchani_model_nn[0][self.feats_layer - 2].out_features

        final_layer = torch.zeros([n_batch, n_atoms, n_feats], dtype=torch.float, device=self.device)
        final_layer_ = final_layer.flatten(0, 1)

        species_ = species.flatten()
        present_species = torchani.utils.present_species(species)
        aevs = aevs.flatten(0, 1)

        for i in present_species:

            mask = (species_ == i)
            input_ = aevs.index_select(0, mask.nonzero().squeeze())
            x = input_

            for j, layer in enumerate(self.torchani_model_nn[i]):

                if j == self.feats_layer:

                    final_layer_[mask, :] = x
                    break

                x = layer(x)

        final_layer = final_layer_.view_as(final_layer)
        return labels_tensor, final_layer

class SymFeats(nn.Module):
    def __init__(self):
        from torchani.aev import AEVComputer
        from torchani.neurochem import Constants
        from a2mdnet import AEV_PARAMETERS

        super(SymFeats, self).__init__()
        self.aev = AEVComputer(
            **Constants(filename=AEV_PARAMETERS)
        )

    def forward(self, *x):
        labels_tensor = x[0]
        coords_tensor = x[1]
        species, aevs = self.aev((labels_tensor, coords_tensor))
        return species, aevs


class ChargeNormalization(nn.Module):
    def __init__(self):
        super(ChargeNormalization, self).__init__()

    def forward(self, *x):
        charge = x[0]
        coefficients_iso = x[1]
        coefficients_aniso = x[2]
        integrals_iso = x[3]
        integrals_aniso = x[4]

        assert isinstance(charge, torch.Tensor)
        assert isinstance(coefficients_iso, torch.Tensor)
        assert isinstance(coefficients_aniso, torch.Tensor)
        assert isinstance(integrals_iso, torch.Tensor)
        assert isinstance(integrals_aniso, torch.Tensor)

        charge = charge.sum(dim=1)

        ciso = coefficients_iso.flatten(1)
        caniso = coefficients_aniso.flatten(1)

        integrals_iso = integrals_iso.flatten(1)
        integrals_aniso = integrals_aniso.flatten(1)

        coefficients = torch.cat((ciso, caniso), dim=1)
        function_integrals = torch.cat((integrals_iso, integrals_aniso), dim=1)

        estimated_charge = coefficients * function_integrals
        estimated_charge = torch.sum(estimated_charge, dim=1)

        factor = charge / estimated_charge
        factor = factor.unsqueeze(dim=1)
        factor = factor.unsqueeze(dim=1)
        return factor, coefficients_iso * factor, coefficients_aniso * factor


class QPChargeNormalization(nn.Module):
    def __init__(self, device):
        super(QPChargeNormalization, self).__init__()
        self.device = device

    def forward(self, *x):
        """

        Least squares with restraints solving.
        This module allows to normalize the input coefficients so it reproduces the integral of the charge.
        To do so:

            1.  Merges the coefficients into a single tensor
            2.  Defines a linear system that represents the derivative of the lagrangian:

                    L(x, u) = (x-c)^T (x-c) - u (qx - Q)

            3. Solves this system and reshapes the output into the input shape

        :param x:
        :return:
        """
        c_iso = x[0]
        c_aniso = x[1]
        int_iso = x[2]
        int_aniso = x[3]
        charge = x[4]

        assert isinstance(c_iso, torch.Tensor)
        assert isinstance(c_aniso, torch.Tensor)
        assert isinstance(int_iso, torch.Tensor)
        assert isinstance(int_aniso, torch.Tensor)
        assert isinstance(charge, torch.Tensor)

        nbatch = c_iso.size()[0]
        natoms = c_iso.size()[1]
        nbonds = c_aniso.size()[1]
        niso_functions = c_iso.size()[2]
        naniso_functions = c_aniso.size()[2]

        c_iso_f = c_iso.reshape(nbatch, -1)
        c_aniso_f = c_aniso.reshape(nbatch, -1)
        int_iso_f = int_iso.reshape(nbatch, -1)
        int_aniso_f = int_aniso.reshape(nbatch, -1)

        coefficients = torch.cat([c_iso_f, c_aniso_f], dim=1)
        integrals = torch.cat([int_iso_f, int_aniso_f], dim=1)

        nfunctions = coefficients.size()[1]

        lqoperator = (torch.eye(nfunctions, device=self.device) * 2).repeat(nbatch, 1, 1)
        integrals_flat = integrals.reshape(-1, nfunctions, 1)
        lqoperator = torch.cat([lqoperator, -integrals_flat], dim=2)

        integrals_flat = integrals.reshape(-1, 1, nfunctions)
        integrals_flat = torch.cat(
            [
                integrals_flat, torch.zeros(nbatch, 1, 1, dtype=torch.float, device=self.device)
            ], dim=2
        )
        lqoperator = torch.cat([lqoperator, integrals_flat], dim=1)

        coefficients_flat = coefficients.reshape(nbatch, nfunctions, 1)
        charge_flat = charge.sum(1).reshape(nbatch, 1, 1)
        operator_problem = torch.cat(
            [
                2 * coefficients_flat,
                charge_flat
            ], dim=1
        )

        operator_solution, lu = torch.solve(operator_problem, lqoperator)
        solutions_iso, solutions_aniso, lagmult = torch.split(
            operator_solution, [natoms * niso_functions, nbonds * naniso_functions, 1], dim=1
        )

        assert isinstance(solutions_iso, torch.Tensor)
        assert isinstance(solutions_aniso, torch.Tensor)

        solutions_iso = solutions_iso.view_as(c_iso)
        solutions_aniso = solutions_aniso.view_as(c_aniso)

        return solutions_iso, solutions_aniso


# -------------------------------------------------------------------------------------------------------------------- #
# This part here is still experimental
# -------------------------------------------------------------------------------------------------------------------- #


class NearestNeighbourModule(nn.Module):
    def __init__(self, alpha, optimize_alpha=True):
        """

        :param alpha:
        :param optimize_alpha:
        """
        super(NearestNeighbourModule, self).__init__()
        self.alpha = nn.Parameter(
            torch.tensor([alpha], requires_grad=optimize_alpha, dtype=torch.float),
            requires_grad=True
        )

    @staticmethod
    def distance_matrix(*x):
        x0 = x[0]
        x1 = x[1]
        n = x0.size(0)  # this is the number of references
        m = x1.size(0)  # this is the number of training points
        d = x0.size(1)  # this should be n dim
        x0p = x0.unsqueeze(1).expand(n, m, d)
        x1p = x1.unsqueeze(0).expand(n, m, d)
        return torch.sqrt(torch.pow(x0p - x1p, 2).sum(2))

    def exp_kernel(self, *x):
        """

        :param x:
        :return:
        """
        input_feats = x[0]
        n = input_feats.size(0)
        m = input_feats.size(1)
        p = torch.exp(-input_feats * self.alpha)
        z = p.sum(dim=1).reshape(-1, 1)
        z = z.expand(n, m)
        w = p / z
        return w

    def forward(self, *x):
        """

        :param x:
            x[0] - training feats
            x[1] - reference feats
            x[2] - reference targets
        :return:
        """
        a = self.distance_matrix(x[0], x[1])
        w = self.exp_kernel(a)
        return w.mm(x[2].reshape(-1, 1)).reshape(-1)


class TorchElementN3(nn.Module):

    def __init__(self, nodes):
        """

        :param nodes:
        """
        super(TorchElementN3, self).__init__()
        self.layers = nn.ModuleList()
        for i in range(1, len(nodes)-1):
            self.layers.append(nn.Linear(nodes[i-1], nodes[i], bias=True))
            self.layers.append(nn.CELU())
        self.layers.append(nn.Linear(nodes[-2], nodes[-1], bias=True))
        self.layers.append(NearestNeighbourModule(alpha=1e2, optimize_alpha=True))
        self.n_layers = len(self.layers)

    def forward(self, *method_set):
        """

        :param method_set:
        :return:
        """

        x_ref = method_set[0][0]
        x_train = method_set[1]
        y_ref = method_set[0][1]

        for i in range(self.n_layers - 1):
            x_ref = self.layers[i](x_ref)
            x_train = self.layers[i](x_train)

        return self.layers[-1](x_train, x_ref, y_ref)


class TorchElementSpecificN3(nn.Module):
    def __init__(self, nodes, elements):
        super(TorchElementSpecificN3, self).__init__()
        self.allowed_elements = elements
        self.subnets = nn.ModuleList()
        for _ in range(len(elements)):
            self.subnets.append(
                TorchElementN3(nodes=nodes)
            )
        self.output_size = nodes[-1]

    def forward(self, *method_sets):
        x_ref = method_sets[0][0]
        y_ref = method_sets[0][1]
        l_ref = method_sets[0][2]

        x_train = method_sets[1][0]
        l_train = method_sets[1][1]

        assert isinstance(x_ref, torch.Tensor)
        assert isinstance(l_ref, torch.Tensor)
        assert isinstance(y_ref, torch.Tensor)
        assert isinstance(x_train, torch.Tensor)
        assert isinstance(l_train, torch.Tensor)

        l_ref_ = l_ref.flatten()
        l_train_ = l_train.flatten()

        present_l_ref = l_ref.unique()
        present_l_train = l_train.unique()

        if present_l_ref[0].item() == 0:
            present_l_ref = present_l_ref[1:]

        if present_l_train[0].item() == 0:
            present_l_train = present_l_train[1:]

        present_l_ref, _ = present_l_ref.sort()
        present_l_train, _ = present_l_train.sort()

        y = torch.zeros(l_train_.size()[0], self.output_size, dtype=torch.float)
        x_ref_ = x_ref.flatten(0, 1)
        x_train_ = x_train.flatten(0, 1)
        # y_ref_ = y_ref.flatten(0, 1)

        for i, n in enumerate(present_l_train):
            mask_input_train = (l_train_ == n)
            mask_input_ref = (l_ref_ == n)

            if not torch.any(mask_input_ref):
                raise RuntimeError("the reference has not the element requested")

            m_ref = mask_input_ref.nonzero().squeeze()
            m_train = mask_input_train.nonzero().squeeze()

            u_ref = x_ref_.index_select(0, m_ref)
            u_train = x_train_.index_select(0, m_train)
            v_ref = y_ref.index_select(0, m_ref)
            y.masked_scatter_(mask_input_train, self.subnets[ELEMENT2NN[n.item()]]((u_ref, v_ref), u_train))
        y.reshape(l_train.size()[0], l_train.size()[1], self.output_size)
        return l_train, y
