import numpy as np
from typing import Callable
from a2mdio.molecules import Mol2
from a2md.models import Molecule

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


def pi_lebedev(fun : Callable, r_max=10.0, radial_res=100, grid='medium'):
    """
    polar integral
    ---
    performs a middle point integral using polar coordinates
    :param fun:
    :param r_max:
    :param radial_res:
    :param grid: either medium, coarse, tight
    :return:
    """
    from a2md import LEBEDEV_DESIGN
    lebedevdesgin = np.loadtxt(LEBEDEV_DESIGN[grid])

    lebedev = np.zeros((lebedevdesgin.shape[0], 3), dtype='float64')
    lebedev[:, 0] = np.cos(np.deg2rad(lebedevdesgin[:, 0])) * np.sin(np.deg2rad(lebedevdesgin[:, 1]))
    lebedev[:, 1] = np.sin(np.deg2rad(lebedevdesgin[:, 0])) * np.sin(np.deg2rad(lebedevdesgin[:, 1]))
    lebedev[:, 2] = np.cos(np.deg2rad(lebedevdesgin[:, 1]))

    w = lebedevdesgin[:, 2]
    u = np.log(r_max + 1) / radial_res
    r_grid = np.exp(np.arange(1, radial_res + 1) * u) - 1
    integral = 0.0
    for i, r in enumerate(r_grid):
        dr = r - (np.exp(u * i) - 1)
        r2 = r
        r1 = r - 0.5 * dr
        r0 = r - dr

        dv = ((r ** 3) / 3) - (((r - dr) ** 3) / 3)

        coords = r0 * lebedev
        f00 = fun(coords)
        coords = r1 * lebedev
        f01 = fun(coords)
        coords = r2 * lebedev
        f02 = fun(coords)

        f = (f00 + 4 * f01 + f02) / 6
        f = f * w * dv * 4 * np.pi
        integral +=  f.sum()

    return integral

def pi_lebedev_m(fun : Callable, n_out : int, r_max=10.0, radial_res=100, grid='coarse'):
    from a2md import LEBEDEV_DESIGN
    lebedevdesgin = np.loadtxt(LEBEDEV_DESIGN[grid])

    lebedev = np.zeros((lebedevdesgin.shape[0], 3), dtype='float64')
    lebedev[:, 0] = np.cos(np.deg2rad(lebedevdesgin[:, 0])) * np.sin(np.deg2rad(lebedevdesgin[:, 1]))
    lebedev[:, 1] = np.sin(np.deg2rad(lebedevdesgin[:, 0])) * np.sin(np.deg2rad(lebedevdesgin[:, 1]))
    lebedev[:, 2] = np.cos(np.deg2rad(lebedevdesgin[:, 1]))

    w = lebedevdesgin[:, 2].reshape(-1, 1)
    u = np.log(r_max + 1) / radial_res
    r_grid = np.exp(np.arange(1, radial_res + 1) * u) - 1
    integral = np.zeros(n_out, dtype='float64')
    for i, r in enumerate(r_grid):
        dr = r - (np.exp(u * i) - 1)
        r2 = r
        r1 = r - 0.5 * dr
        r0 = r - dr

        dv = ((r ** 3) / 3) - (((r - dr) ** 3) / 3)

        coords = r0 * lebedev
        f00 = fun(coords)
        coords = r1 * lebedev
        f01 = fun(coords)
        coords = r2 * lebedev
        f02 = fun(coords)

        f = (f00 + 4 * f01 + f02) / 6
        f = f * w * dv * 4 * np.pi
        integral += f.sum(0)

    return integral


def integrate_density_functional(functional:Callable, mol:Mol2, grid='coarse', res=100):
    functional_value = 0.0
    for fx in split_space(mm=mol, fun=functional):
        functional_value += pi_lebedev(fx, radial_res=res, grid=grid)
    return functional_value

def integrate_density_functional_gradient(functional:Callable, mol:Mol2, nfuns:int, grid='coarse', res=100):
    r = mol.get_coordinates(units='au')
    n = mol.get_number_atoms()
    functional_value = np.zeros(nfuns, dtype='float64')
    for i in range(n):
        r0 = r[i, :]
        fx = lambda x : functional(x + r0) * (voronoi(x + r0, r) == i).reshape(-1, 1)
        functional_value += pi_lebedev_m(fx, nfuns, radial_res=res, grid=grid)
    return functional_value

def kinetic_energy_functional(fun:Callable):
    cf = (3/10) * ((3 * np.pi**2)**(2/3))
    ke = lambda x : cf * np.power(fun(x),5/3)
    return ke

def exchange_energy_functional(fun:Callable):
    cx = -(3/4)*((3/np.pi)**(1/3))
    xe = lambda x: cx * np.power(fun(x), 4 / 3)
    return xe

def mse_functional(ref : Callable, fun : Callable):
    def mse(x):
        rxf = ref(x)
        fxf = fun(x)
        return rxf * np.power(fxf - rxf, 2.0)
    return mse

def mlse_functional(ref : Callable, fun : Callable):
    def mlse(x):
        rxf = ref(x)
        fxf = fun(x)
        return rxf * np.log(np.power(fxf - rxf, 2.0))
    return mlse

def dkl_functional(ref : Callable, fun : Callable):
    def dkl(x):
        rxf = ref(x)
        fxf = fun(x)
        return fxf * np.log(fxf/rxf)
    return dkl

def vdwvolume_functional(ref: Callable, eps=1e-3):
    def vdwvol(x):
        rxf = ref(x)
        return (rxf > eps).astype(float)
    return vdwvol

def dkl_gradient_functional(ref : Callable, model : Molecule):
    nopt = model.get_number_optimizable_functions()
    nfuns = model.get_number_functions()
    def dkl_gradient(x):
        u = np.zeros((x.shape[0], nopt), dtype='float64')
        pref = ref(x)
        cand = model.eval(x)
        lp = np.log(cand / pref) + 1
        indx = [i for i in range(nfuns) if not model.map_frozenfunctions[i]]
        for idx, fun_idx in enumerate(indx):
            u[:, idx] = model.functions[fun_idx].eval(x) * lp
        return u
    return dkl_gradient