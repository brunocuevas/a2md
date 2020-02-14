
from a2mdnet.modules import *
from a2mdnet import FIELDS


class TorchA2MDNetConfiguration:

    def __init__(self, file=None, *kws):
        if file is None and len(kws) == 0:
            raise IOError("missing information")
        elif file is not None:
            print("reading file {:s}".format(file))
            self.conf = self.__read_from_file(file)
        elif file is None and len(kws) > 0:
            print("reading from list")

    @staticmethod
    def __read_from_file(file):

        import json
        with open(file) as f:
            conf = json.load(f)
        for i, item in enumerate(conf):
            assert isinstance(item, dict)
            check_fields = all([k in FIELDS for k in item.keys()])
            if not check_fields:
                raise IOError("missing field in instance number {:d}".format(i))
        return conf

    def get_configuration(self):

        return self.conf


class TorchRadialSupport(nn.Module):

    def __init__(self, nodes, elements, device):
        """

        :param nodes:
        :param elements:
        """
        super(TorchRadialSupport, self).__init__()
        self.net = TorchElementSpecificA2MDNN(
            nodes=nodes,
            elements=elements,
            device=device
        )

    def forward(self, u):
        """

        :param u:
        :return:
        """
        return self.net(u)


class TorchAngularSupport(nn.Module):

    def __init__(self, nodes, elements):
        """

        :param nodes:
        :param elements:
        """
        super(TorchAngularSupport, self).__init__()
        self.net = TorchPairSpecificA2MDNN(
            nodes=nodes, elements=elements
        )

    def forward(self, u):
        """

        :param u:
        :return:
        """
        return self.net(u)


AVAILABLE_NETS = dict(
    radial=TorchRadialSupport,
    angular=TorchAngularSupport
)


class TorchA2MDNet(nn.Module):
    def __init__(self, configuration, nodes, elements):
        """

        :param configuration:
        :param nodes:
        :param elements:
        """
        super(TorchA2MDNet, self).__init__()
        if isinstance(configuration, TorchA2MDNetConfiguration):
            self.conf = configuration.get_configuration()
        self.nodes = nodes
        self.elements = elements
        self.net = self.__build()

    def __build(self):
        """

        :return:
        """
        nets = []
        for i, fun in enumerate(self.conf):
            print("including function : {:s}".format(fun['NAME_']))
            if fun['TYPE_'] in AVAILABLE_NETS.keys():
                # Just to go faster
                if fun['TYPE_'] == 'radial':
                    net = TorchRadialSupport(
                        nodes=fun['NODES_'],
                        elements=fun['ELEMENTS'],
                        device=fun['DEVICE']
                    )
                elif fun['TYPE_'] == 'angular':
                    net = TorchAngularSupport(
                        nodes=fun['NODES_'],
                        elements=fun['ELEMENTS']
                    )
                else:
                    raise IOError
                nets.append(net)
            else:
                raise KeyError("could not found such a {:s} function".format(fun['TYPE_']))
        return nets

    def forward(self, u):
        """

        :param u:
        :return:
        """
        res = []
        for i, net in enumerate(self.net):
            res.append(net(u))
        return res


class A2MDN3(nn.Module):
    def __init__(self, architecture, elements, feats_ref, values_ref, labels_ref):
        """

        :param architecture:
        :param elements:
        """
        super(A2MDN3, self).__init__()
        self.architecture = architecture
        self.elements = elements
        self.net = self.__build()

        self.reference_coords = feats_ref
        self.feats = None
        self.values = values_ref

        self.labels = labels_ref

    def __build(self):

        net = nn.ModuleList([])
        net.append(
            TorchaniFeats(net=0, feats_layer=6)
        )
        net.append(
            TorchElementSpecificN3(
                nodes=self.architecture, elements=self.elements
            )
        )
        return net

    def forward_reference(self):

        max_size = 0
        buffer_feats = []
        buffer_labels = []
        for i in range(len(self.reference_coords)):
            ta_feats = self.net[0](
                self.labels[i].unsqueeze(0), self.reference_coords[i].unsqueeze(0)
            )
            if ta_feats.size(0) > max_size:
                max_size = ta_feats.size(0)
            buffer_feats.append(ta_feats)

        for i in range(len(self.reference_coords)):
            tmp = self.pad(buffer_feats[i], dims=(max_size, 96))
            buffer_feats[i] = tmp
            buffer_labels.append(
                self.pad(
                    self.labels[i].unsqueeze(1),
                    dims=(max_size, 1),
                    dtype=torch.long)
            )

        assert isinstance(self.values, list)
        self.feats = torch.stack(buffer_feats)
        self.labels = torch.stack(buffer_labels)
        self.labels = self.labels.squeeze(2)
        self.values = torch.stack(self.values)

    def forward(self, *x):

        feats = x[1]
        labels = x[0]
        output = []

        for i in range(len(feats)):

            feats_instance = feats[i].unsqueeze(0)
            labels_instance = labels[i].unsqueeze(0)

            # ta_feats_ref = self.net[0](labels_ref, feats_ref)
            ta_feats_train = self.net[0](labels_instance, feats_instance)

            y = self.net[1](
                (self.feats, self.values, self.labels),
                (ta_feats_train, labels_instance),
            )

            output.append(y)

        return torch.tensor(output, dtype=torch.float, device=torch.device('cpu')).reshape(-1, 1)

    @staticmethod
    def pad(tnsr, dims, dtype=torch.float):
        u = torch.zeros(dims, dtype=dtype, device=torch.device('cpu'))
        s = tnsr.size()
        u[:s[0], :] = tnsr
        return u


class A2MDNearestNeighbors:

    def __init__(self, reference_set, reference_values, n_neighbors=1):
        """

        The purpose of this class is making an efficient nearest neighbours using Sklearn libraries
        over PyTorch Autoencoders.

        :param reference_set:
        :param reference_values:
        """
        from sklearn.neighbors import NearestNeighbors
        self.reference_set = reference_set
        self.reference_values = reference_values
        self.estimator = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto')
        self.estimator = self.estimator.fit(self.reference_set)

    def eval(self, u):
        """

        :param u:
        :return:
        """
        distance, indices = self.estimator.kneighbors(u)
        predicted_values = self.reference_values[indices, :]
        return predicted_values.mean(axis=1)


if __name__ == '__main__':

    pass
