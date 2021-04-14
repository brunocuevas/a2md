import torch


def distance_vectors(sample_coords, mol_coords, labels, device):
    """

    Takes the coordinates of the molecule and the sampling coordinates,
    and calculates the distance vectors between them. If n is the batch size,
    m the number of density samples, Ma the number of molecule coordinates;
    the result is a tensor of range 4, dimensions: n, m, Ma, 3

    :param sample_coords:
    :type sample_coords: torch.Tensor
    :param mol_coords:
    :type mol_coords: torch.Tensor
    :param labels:
    :type labels: torch.Tensor
    :param device:
    :type device: torch.device
    :return: Tensor (n, m, Ma, 3, dtype=torch.float)
    """
    # dimensions
    n = sample_coords.size(0)
    m_sample = sample_coords.size(1)
    m_mol = mol_coords.size(1)
    # output tensor
    dv = torch.zeros(
        n, m_sample, m_mol, 3,
        device=device, dtype=torch.float
    )
    # distance operation takes place by each element of the batch
    # the operation is performed by expanding the molecule and the
    # sampling
    for i in range(n):
        mask = (labels[i, :] != -1).unsqueeze(1).expand(m_mol, 3)
        sliced_mol = mol_coords[i, :, :].reshape(m_mol, 3)
        masked_mol = sliced_mol.masked_select(mask).reshape(-1, 3)
        n_ = masked_mol.size(0)
        masked_mol = masked_mol.unsqueeze(0)
        masked_mol = masked_mol.expand(m_sample, n_ , 3)
        sliced_sample = sample_coords[i, :, :].unsqueeze(1)
        sliced_sample = sliced_sample.expand(m_sample, n_, 3)
        dv[i, :,:n_, :] = sliced_sample - masked_mol
    return dv

def distance(dv):
    """

    Takes a (n, m, Ma, 3) tensor, and returns the p2 norm along the
    third dimensions of the tensor. The resulting tensor is
    (n, m, Ma).

    :param dv:
    :type dv: torch.Tensor
    :return:
    """
    return dv.norm(dim=3, p=2)

def select_labels(l, t, device):
    """

    Takes a (n, Ma) tensor containing the labels of each atom of a molecule, and
    a (n, Mb, 2) tensor containing the topology of a molecule. Returns two tensors,
    both of dimensions (n, Mb, 1). The first tensor is a selection of the elements
    participating in a bond by one side, while the second tensor is a selection of
    the elements participating by the other

    :param l:
    :type l: torch.Tensor
    :param t:
    :type t: torch.Tensor
    :param device:
    :type device: torch.device
    :return:
    """
    # Dimensions
    n_batch = l.size(0)
    n_connectivity = t.size(1)
    n_atoms = l.size(1)
    # Flatten tensors
    t_f = t.flatten(0, 1)
    l_f = l.flatten(0, 1)
    # The way this operation is performed is:
    #   1- flatten everything
    #   2- update the topology to match the right index in the flat version of the labels
    #   3- select those index and stores them in the two output tensors
    n_batch_arange = torch.arange(n_batch, device=device, dtype=torch.long) * n_atoms
    l_forward = torch.zeros(n_batch * n_connectivity, dtype=torch.long, device=device)
    l_backward = torch.zeros(n_batch * n_connectivity, dtype=torch.long, device=device)
    mask = t_f[:, 0] != -1
    t_f = t_f + n_batch_arange.reshape(-1, 1).repeat(1, n_connectivity).reshape(-1, 1)
    l_forward.masked_scatter_(
        mask, l_f.index_select(dim=0, index=t_f[mask, 0])
    )
    l_backward.masked_scatter_(
        mask, l_f.index_select(dim=0, index=t_f[mask, 1])
    )
    return l_forward.reshape(n_batch, n_connectivity), l_backward.reshape(n_batch, n_connectivity)


def select_distances(d, t, device):
    """

    Takes a (n, Ma) tensor containing the distances between the sample points and
    each atom of a molecule, and a (n, Mb, 2) tensor containing the topology of
    a molecule. Returns two tensors, both of dimensions (n, Mb, 1). The first
    tensor is a selection of the distances between the samples and the atom
    participating in a bond by one side, while the second tensor is
    a selection of the distances between the samples and the atom
    by the other side.

    :param d:
    :type d: torch.Tensor
    :param t:
    :type t: torch.Tensor
    :param device:
    :type device: torch.device
    :return:
    """
    n_batch = d.size(0)
    n_connectivity = t.size(1)
    n_atoms = d.size(2)
    n_batch_arange = torch.arange(n_batch, device=device, dtype=torch.long) * n_atoms
    n_samples = d.size(1)
    t_f = t.flatten(0, 1)
    mask = t_f[:, 0] != -1
    t_f = t_f + n_batch_arange.reshape(-1, 1).repeat(1, n_connectivity).reshape(-1, 1)
    d_f = d.transpose(1,2).flatten(0, 1)
    d_forward = torch.zeros(
        n_batch * n_connectivity, n_samples,
        device=device, dtype=torch.float
    )
    d_reverse = torch.zeros(
        n_batch * n_connectivity, n_samples,
        device=device, dtype=torch.float
    )
    mask_expanded = mask.unsqueeze(1).expand(mask.size(0), n_samples)
    d_forward.masked_scatter_(
        mask_expanded,
        d_f.index_select(dim=0, index=t_f[mask, 0])
    )
    d_reverse.masked_scatter_(
        mask_expanded,
        d_f.index_select(dim=0, index=t_f[mask, 1])
    )

    d_forward = d_forward.reshape(n_batch, n_connectivity, n_samples)
    d_reverse = d_reverse.reshape(n_batch, n_connectivity, n_samples)
    d_forward = d_forward.transpose(1, 2)
    d_reverse = d_reverse.transpose(1, 2)
    return d_forward, d_reverse


def angle(dv, v, sample_coords, mol_coords, connectivity, device):
    """

    Calculates the angle between a bond-vector and a vector joining an
    atomic centre and the sampling point. Takes the (n, m, Ma, 3) distance
    tensor, the (n, Ma, m) distance norm tensor, the (n, Ma, 3) molecule
    coordinates, the (n, Mb, 2) molecule topology. Returns a (n, m, Mb)
    tensor. Angles are in radians, and should be between 0 and 3.14

    :param dv: [n, m, Ma, 3]
    :type dv: torch.Tensor
    :param v: [n, m, Ma]
    :type v: torch.Tensor
    :param sample_coords: [n, m, 3]
    :type sample_coords: torch.Tensor
    :param mol_coords: [n, Ma, 3]
    :type mol_coords: torch.Tensor
    :param connectivity: [n, Mb, 2]
    :param device:
    :return:
    """
    n = sample_coords.size(0)
    m_sample = sample_coords.size(1)
    m_con = connectivity.size(1)

    z = torch.zeros(
        n, m_sample, m_con,
        device=device, dtype=torch.float
    )

    for i in range(n):
        # Selection of the coordinates of the atoms at both
        # sides of each bond
        mask = connectivity[i, :, 0] != -1
        connectivity_masked_i1 = connectivity[i, mask.nonzero(), 0].flatten()
        connectivity_masked_i2 = connectivity[i, mask.nonzero(), 1].flatten()
        term1 = mol_coords[i, :, :].index_select(0, connectivity_masked_i1)
        term2 = mol_coords[i, :, :].index_select(0, connectivity_masked_i2)
        # distance vectors
        du = term2 - term1
        u = du.pow(2.0).sum(1).add(1e-18).sqrt()
        du = du.unsqueeze(0)
        u = u.unsqueeze(0)
        # selecting the pertinent distance vectors
        sdv = dv[i, :, :, :].index_select(1, connectivity_masked_i1)
        sv = v[i, :, :].index_select(1, connectivity_masked_i1)
        # use of the scalar product to obtain the angle
        dusdv = (du * sdv).sum(2)
        usv = (u * sv)
        expanded_mask = mask.unsqueeze(0).expand(m_sample, m_con)
        ratio = dusdv/usv
        ratio[ratio< -1.0] = -1.0
        ratio[ratio > 1.0] = 1.0
        z[i, :, :].masked_scatter_(expanded_mask, torch.acos(ratio))

    return z

def expand_parameter(labels, param):
    """

    Takes a tensor containing parameters, and expands it to match the
    molecule atoms

    :param labels:
    :param param:
    :return:
    """
    labels_copy = labels.clone()
    labels_copy[labels == -1] = 0
    output = torch.zeros_like(labels, dtype=torch.float)

    output.masked_scatter_(
        labels != -1,
        param.index_select(0, labels_copy.flatten()).view_as(labels)[labels != - 1]
    )
    return output

def exponential_kernel(d, a, b):
    """

    Applies a*exp(-B*d)

    :param d:
    :param a:
    :param b:
    :return:
    """
    a_usq = a.unsqueeze(1)
    b_usq = b.unsqueeze(1)
    buffer = torch.exp(-d * b_usq)
    buffer = a_usq * buffer
    return buffer

def xexponential_kernel(d, a, b):
    """

    Applies u*d*exp(-g*d)

    :param d:
    :param a:
    :param b:
    :return:
    """

    buffer =  torch.exp(-d * b.unsqueeze(1))
    buffer = a.unsqueeze(1) * buffer * d
    return buffer

def gaussian_kernel(z, alpha):
    """

    Applies exp(-alpha * (z**2))

    :param z:
    :param alpha:
    :return:
    """
    buffer =  torch.exp(-alpha.unsqueeze(1) * z.pow(2.0))
    return buffer