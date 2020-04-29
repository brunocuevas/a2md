import numpy as np
from typing import Callable
from a2mdio.molecules import Mol2

def voronoi(x : np.ndarray, r : np.ndarray):
    """
    voronoi
    ---
    tessellates x points by calculating their distances to the different
    points from r
    :param x:
    :param r:
    :return:
    """
    d = np.zeros((x.shape[0], r.shape[0]))
    for i, ir in enumerate(r):
        d[:, i] = np.linalg.norm(x - r[i, :], axis=1)
    t = np.argmin(d, axis=1)
    return t

def split_space(mm : Mol2, fun : Callable):
    """
    split space
    ---
    divides fun in n functions centered around atoms
    :param mm:
    :param fun:
    :return:
    """
    r = mm.get_coordinates(units='au')
    n = mm.get_number_atoms()
    for i in range(n):
        r0 = r[i, :]
        yield lambda x : fun(x + r0) * (voronoi(x + r0, r) == i)

def polar_integral(fun : Callable, r_max=10.0, radial_res=100, angular_res=10):
    """
    polar integral
    ---
    performs a middle point integral using polar coordinates
    :param fun:
    :param r_max:
    :param radial_res:
    :param angular_res:
    :return:
    """
    u = np.log(r_max + 1) / radial_res
    dphi = 2 * np.pi / angular_res
    drho = np.pi / angular_res
    phi_grid = np.arange(dphi, 2 * np.pi + dphi, dphi)
    rho_grid = np.arange(drho, np.pi + drho, drho)
    r_grid = np.exp(np.arange(1, radial_res + 1) * u) - 1
    integral = 0.0
    for i, r in enumerate(r_grid):
        for j, h in enumerate(rho_grid):

            dr = r - (np.exp(u * i) - 1)

            dv = ((r ** 3) / 3) - (((r - dr) ** 3) / 3)
            dv *= (np.cos(h - drho) - np.cos(h))
            dv *= dphi

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = (r - dr) * np.sin(h - drho) * np.cos(phi_grid - dphi)
            coords[:, 1] = (r - dr) * np.sin(h - drho) * np.sin(phi_grid - dphi)
            coords[:, 2] = (r - dr) * np.cos(h - drho)
            f000 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = (r - dr) * np.sin(h - drho) * np.cos(phi_grid)
            coords[:, 1] = (r - dr) * np.sin(h - drho) * np.sin(phi_grid)
            coords[:, 2] = (r - dr) * np.cos(h - drho)
            f010 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = (r - dr) * np.sin(h) * np.cos(phi_grid - dphi)
            coords[:, 1] = (r - dr) * np.sin(h) * np.sin(phi_grid - dphi)
            coords[:, 2] = (r - dr) * np.cos(h)
            f001 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = (r - dr) * np.sin(h) * np.cos(phi_grid)
            coords[:, 1] = (r - dr) * np.sin(h) * np.sin(phi_grid)
            coords[:, 2] = (r - dr) * np.cos(h)
            f011 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = r * np.sin(h - drho) * np.cos(phi_grid - dphi)
            coords[:, 1] = r * np.sin(h - drho) * np.sin(phi_grid - dphi)
            coords[:, 2] = r * np.cos(h)
            f100 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = r * np.sin(h - drho) * np.cos(phi_grid)
            coords[:, 1] = r * np.sin(h - drho) * np.sin(phi_grid)
            coords[:, 2] = r * np.cos(h)
            f110 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = r * np.sin(h) * np.cos(phi_grid - dphi)
            coords[:, 1] = r * np.sin(h) * np.sin(phi_grid - dphi)
            coords[:, 2] = r * np.cos(h)
            f101 = fun(coords)

            coords = np.zeros((angular_res, 3), dtype='float64')
            coords[:, 0] = r * np.sin(h) * np.cos(phi_grid)
            coords[:, 1] = r * np.sin(h) * np.sin(phi_grid)
            coords[:, 2] = r * np.cos(h)
            f111 = fun(coords)

            integral += (dv * (f000 + f001 + f010 + f011 + f100 + f101 + f110 + f111)/8.0).sum()

    return integral