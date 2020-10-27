import torch
import math


def kl_divergence(p: torch.Tensor, q: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
    """
    Kullback-Leibler divergence

    Implemented as dkl = int p * log(p/q) dr = H(P,Q) - H(P)

    ¡It is important to notice the assymetry!
    ¡Also important to scale input densities!

    :param p: optimizable function
    :param q: reference function
    :param w: weights
    :return: kullback-leibler divergence. Units: nats
    """
    dkl = (p * (p/q).log() * w).sum()
    return dkl


def mean_squared_error(p: torch.Tensor, q: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
    """
    Mean Squared Error functional

    Implemented as MSE = int q * (p - q)^2 dr = E_{q} [ (p-q)^2 ]

    :param p: optimizable function
    :param q: reference function
    :param w: weights
    :return: mean squared error. Units: (electrons/(bohr^3))^2
    """
    mse = (q * (p - q).pow(2.0) * w).sum()
    return mse


def dirac_exchange(p: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
    """

    :param p:
    :param w:
    :return:
    """
    cx = (3.0/4.0) * math.pow((3.0/math.pi), 1.0 / 3.0)
    kd = - cx * (p.pow(4.0/3.0) * w).sum()
    return kd


def thomas_fermi_kinetic(p: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
    """

    :param p:
    :param w:
    :return:
    """
    cf = 3 * math.pow(3 * math.pi * math.pi, 2.0/3.0) / 10.0
    ttf = cf * (p.pow(5.0/3.0) * w).sum()
    return ttf
