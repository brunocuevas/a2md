import torch
import numpy as np
from torch import nn
from a2mdnet.data import convert_label2tensor
from a2mdnet.a2mdt.functions import distance, distance_vectors, angle
from a2mdnet.a2mdt.functions import gaussian_kernel, exponential_kernel, xexponential_kernel
from a2mdnet.a2mdt.functions import select_distances, select_labels
from a2mdnet.a2mdt.functions import expand_parameter
from a2md.utils import element2an
from a2mdio.molecules import UNITS_TABLE


def convert_params(pars):
    new_pars = dict()
    for key, item in pars.items():
        nkey = key.lower().replace('3', '')
        new_pars[nkey] = item
    return new_pars


class A2MDtFun:

    def __init__(self, elements, params, device, frozen=False):
        """

        Base clase for the density functions

        :param elements:
        :param params:
        :param device:
        :param frozen:
        """
        self.device = device
        self.elements = convert_label2tensor(
            [element2an(i) for i in elements], device=device
        )
        self.raw_params = params
        self.params = dict()
        for elem, item in params.items():
            idx = elements.index(elem)
            for par, value in item.items():
                try:
                    self.params[par][idx] = value
                except KeyError:
                    self.params[par] = torch.zeros(len(elements), dtype=torch.float, device=device)
                    self.params[par][idx] = value
        if type(frozen) is bool:
            self.frozen = frozen
        else:
            raise IOError("frozen must be bolean")

    def mask_input(self, labels):
        missing = [i for i in labels.unique() if not i in self.elements]
        mask = labels.clone()
        for i in missing:
            mask[labels == i] = -1
        for i, j in enumerate(self.elements):
            mask[labels == j] = i
        return mask

    @staticmethod
    def filter_args(fun, args):

        args_list = fun.__code__.co_varnames
        tmp_args = dict()
        for i in args_list:
            try:
                tmp_args[i] = args[i]
            except KeyError:
                pass
        #tmp_args = dict((i, args[i]) for i in args_list)
        return tmp_args


class A2MDtIso(A2MDtFun):

    def __init__(self, elements, params, device, frozen):
        """

        Isotropic density function

        :param elements: use the torch tensor format
        :param params: dictionary of parameters
        :param device:
        :param frozen: use True if no coefficient must be applied
        """
        A2MDtFun.__init__(self, elements, params, device, frozen)


    def forward(self, l, d):
        """

        Returns density

        :param l: labels
        :param d: distances
        :return:
        """


        mask = self.mask_input(l)

        expanded_params = dict(
            (i, expand_parameter(mask, self.params[i])) for i in self.params.keys()
        )

        return exponential_kernel(
            d, **self.filter_args(exponential_kernel, expanded_params)
        )


class A2MDtAniso(A2MDtFun):

    def __init__(self, elements, params, device, frozen):
        """

        Anisotropic density function

        :param elements: use the torch tensor format
        :param params: dictionary of parameters
        :param device:
        :param frozen: use True if no coefficient must be applied
        """
        A2MDtFun.__init__(self, elements, params, device, frozen)


    def forward(self, l, d, z):
        """

        Returns density

        :param l: labels
        :param d: distance
        :param z: angles
        :return:
        """

        mask = self.mask_input(l)

        expanded_params = dict(
            (i, expand_parameter(mask, self.params[i])) for i in self.params.keys()
        )

        return xexponential_kernel(
            d, **self.filter_args(xexponential_kernel, expanded_params)
        ) * gaussian_kernel(
            z, **self.filter_args(gaussian_kernel, expanded_params)
        )


OLD_FUNCTION_DICT = dict(
    AG01=A2MDtAniso,
    AG02=A2MDtAniso,
    ORC=A2MDtIso,
    ORV=A2MDtIso,
    ORCV=A2MDtIso,
    OR=A2MDtIso
)

FUNCTION_DICT = dict(
    B01=A2MDtAniso,
    B02=A2MDtAniso,
    CR=A2MDtIso,
    CVR=A2MDtIso,
    VR=A2MDtIso,
    hVR=A2MDtIso
)

class A2MDt:

    def __init__(self, params, device):
        model_info = params['_FEATURES']
        model = params['_MODEL']
        functions_dict = dict()
        functions_pos = dict()

        self.functions = []
        self.positions = []
        self.device = device
        frozen_dict = dict()
        for elem, elem_list in model.items():

            for fun in elem_list:
                name = fun['_NAME']

                frozen_dict[name] = fun['_FROZEN']
                if name in functions_dict.keys():
                    functions_dict[name][elem] = convert_params(fun['_PARAMS'])
                    functions_pos[name].append(fun['_TENSOR_POS'])
                else:
                    functions_dict[name] = dict()
                    functions_dict[name][elem] = convert_params(fun['_PARAMS'])
                    functions_pos[name] = [fun['_TENSOR_POS']]

        for fun_key, pos in functions_pos.items():

            fun_mapped = FUNCTION_DICT[fun_key]
            fun = functions_dict[fun_key]
            try:
                pos = int(np.unique(pos))
            except TypeError:
                pos = None
            self.functions.append(
                fun_mapped(
                    elements=list(fun.keys()),
                    params=fun, device=device,
                    frozen=frozen_dict[fun_key]
                )
            )
            self.positions.append(pos)

    def forward_core(self, l, r, x):

        dv = distance_vectors(x, r, l, self.device)
        d = distance(dv)
        p = torch.zeros(
            (x.size(0), x.size(1)), dtype=torch.float, device=self.device
        )
        for fun in self.functions:

            if fun.frozen:
                if type(fun) is A2MDtIso:
                    p += fun.forward(l, d).sum(2)

        return p

    def forward(self, l, t, r, c_i, c_a, x):
        """

        :param l: labels
        :param t: topology
        :param r: molecular coordinates (angstroms)
        :param c_i: isotropic coefficients
        :param c_a: anisotropic coefficients
        :param x: sample coordinates (atomic units)
        :return:
        """

        r = r * UNITS_TABLE['angstrom']['au']
        c_i = torch.split(c_i, 1, 2)
        c_a = torch.chunk(c_a, 2, 2)
        c_a_f = torch.split(c_a[0], 1, 2)
        c_a_b = torch.split(c_a[1], 1, 2)

        dv = distance_vectors(x, r, l, self.device)
        d = distance(dv)
        z1 = angle(dv, d, x, r, t, self.device)
        t_flipped = t.flip(2)
        z2 = angle(dv, d, x, r, t_flipped, self.device)

        p = torch.zeros(
            (x.size(0), x.size(1)), dtype=torch.float, device=self.device
        )


        for fun, pos in zip(self.functions, self.positions):

            if fun.frozen:
                if type(fun) is A2MDtIso:
                    p += fun.forward(l, d).sum(2)
                    continue
                else:
                    raise NotImplementedError("only frozen isotropic functions are considered")

            if type(fun) is A2MDtIso:

                c_i_ = c_i[pos].transpose(1, 2)
                p += (c_i_ * fun.forward(l, d)).sum(2)


            elif type(fun) is A2MDtAniso:

                c_a_forward = c_a_f[pos].transpose(1, 2)
                c_a_reverse = c_a_b[pos].transpose(1, 2)

                d_selected_forward, d_selected_reverse = select_distances(
                    d, t, device=self.device
                )
                l_selected_forward, l_selected_backwards = select_labels(
                    l, t, device=self.device
                )

                p += (c_a_forward * fun.forward(
                    l_selected_forward, d_selected_forward, z1)
                      ).sum(2)
                p += (c_a_reverse * fun.forward(
                    l_selected_backwards, d_selected_reverse, z2)
                      ).sum(2)


        return p


# UNTESTED
    def eval_volume(self, *mol, resolution=0.2, extend=3.0):
        """

        :param resolution:
        :param extend:
        :return:
        """
        from a2mdio.qm import ElectronDensity
        l, t, r, c_i, c_a = mol
        extend = extend * UNITS_TABLE['au']['angstrom']
        resolution = resolution * UNITS_TABLE['au']['angstrom']
        n_batch = l.size(0)
        # r = r * UNITS_TABLE['angstrom']['au']
        for i in range(n_batch):
            l_slice = l[i, :]
            t_slice = t[i, :, :]
            r_slice = r[i, :, :]
            c_i_slice = c_i[i, :, :]
            c_a_slice = c_a[i, :, :]

            coords_min = [
                r_slice[:, 0].min() - extend,
                r_slice[:, 1].min() - extend,
                r_slice[:, 2].min() - extend
            ]

            coords_max = [
                r_slice[:, 0].max() + extend,
                r_slice[:, 1].max() + extend,
                r_slice[:, 2].max() + extend
            ]

            xx = torch.arange(
                coords_min[0], coords_max[0], resolution,
                dtype=torch.float, device=self.device
            )
            yy = torch.arange(
                coords_min[1], coords_max[1], resolution,
                dtype=torch.float, device=self.device
            )
            zz = torch.arange(
                coords_min[2], coords_max[2], resolution,
                dtype=torch.float, device=self.device
            )

            nx = xx.size(0)
            ny = yy.size(0)
            nz = zz.size(0)

            xx = xx.unsqueeze(1).unsqueeze(2).expand(nx, ny, nz)
            yy = yy.unsqueeze(0).unsqueeze(2).expand(nx, ny, nz)
            zz = zz.unsqueeze(0).unsqueeze(1).expand(nx, ny, nz)
            dxr = torch.stack([xx, yy, zz], dim=3).reshape(-1, 3)
            dxr *= UNITS_TABLE['angstrom']['au']
            dx = self.forward(
                l_slice.unsqueeze(0),
                t_slice.unsqueeze(0),
                r_slice.unsqueeze(0),
                c_i_slice.unsqueeze(0),
                c_a_slice.unsqueeze(0),
                dxr.unsqueeze(0)
            ).reshape(nx, ny, nz)
            if dx.type == 'cpu':
                dx = dx.data.numpy()
            else:
                dx = dx.data.cpu().numpy()

            vol_density = ElectronDensity(verbose=False)
            vol_density.set_r0(np.array(coords_min))
            vol_density.set_basis(np.identity(3) * resolution)
            vol_density.set_volume(dx)
            yield vol_density