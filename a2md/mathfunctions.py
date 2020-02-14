import numpy as np
from scipy.special import erfi, erf, expi

INVERSE_DIST_DUMPING = 1e-3

SINEINTEGRAL = [
    np.pi/2,
    0,
    0,
    0,
]
COSINEINTEGRAL = [
    0,
    -2/3,
    0,
    -2/15,
]

TRIGOEXPANSION = [
    lambda u,v : 2*u + ((4 + np.pi)*v/2),
    lambda u,v : 4*u/3 + 2*v,
    lambda u,v : 2*(u+v),
    lambda u,v : (28/15)*u + 2*v
]
NUM2ELE = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
ELE2NUM = dict(H = 0, He = 1, Li = 2, Be = 3, B = 4 , C = 5, N = 6, O = 7, F = 8, Ne = 9 )

def filter_min(d):
    if isinstance(d, float) or isinstance(d, int):
        return np.max([float(d), INVERSE_DIST_DUMPING])
    elif isinstance(d, np.ndarray):
        d = d.copy()
        d[d < INVERSE_DIST_DUMPING] = INVERSE_DIST_DUMPING
        return d
    else:
        raise IOError

# Angle calculation

def get_angle(x, center, ref_frame = None):
    """

    :param x:
    :param center:
    :param ref_frame:
    :return:
    """
    x = x - center
    if ref_frame is None:
        d = np.linalg.norm(x, axis=1)
        return np.zeros(x.shape[0], dtype='float64'), d
    else:
        x = ref_frame.dot(x.T).T
        d = np.linalg.norm(x, axis=1)
        module = np.copy(d)
        module[np.abs(module) < 1e-12] = 1  # This is a random value. It does not matter
        # since in those cases r = 0
        z = x[:, 2] / module
        z[np.isnan(z)] = 0.0
        # For the sake of numerical integrity ----------------------------
        z[z > 1.0] = 1.0
        z[z < -1.0] = -1.0
        # ----------------------------------------------------------------
        z = np.arccos(z)
        return z, d

# DENSITY FUNCTIONS

#   Exponential functions
# ------------------------------------
def expfun(A,B, d):
    """
    exponential function A*exp(-B*x)
    :param A: Coefficient
    :param B: Exponent
    :param d: Distance
    :type A: float
    :type B: float
    :type d: np.ndarray
    :return: exponential function
    :rtype: np.ndarray
    """
    return A*np.exp(-B*d)

def xexpfun(A, B, d):
    """
    xexponential function x*F*exp(-G*x)
    :param A: Coefficient
    :param B: Exponent
    :param d: Distance
    :type A: float
    :type B: float
    :type d: np.ndarray
    :return: xexponential function
    :rtype: np.ndarray
    """
    return d * A * np.exp(-B * d)

def x2expfun(A, B, d):
    """

    x2exponential funcion (x^2)*F*exp(-G*x)
    :param A:
    :param B:
    :param d:
    :return:
    """
    return A * (d * d) * np.exp(-B * d)

#   Gaussian functions
# ------------------------------------

def gaussian(alpha, h, d):
    """

    :param alpha:
    :param h:
    :param d:
    :return:
    """
    return h * np.exp(-alpha * (d**2))

#   Trigonometric functions
# ------------------------------------

def nonefun(z):
    return 1.0

def sinfun(psi, z):
    """

    :param psi:
    :param z:
    :return:
    """
    return np.sin(psi*z)

def cosfun(psi, z):
    """

    :param psi:
    :param z:
    :return:
    """
    return np.cos(psi*z)


def cosssin(u,v,s,z):
    """
    sine + cosine function v*sin(s*z)+u*(cos(s*z))
    :param u: cosine coefficient
    :param v: sine coefficient
    :param s: period
    :param z: angle
    :type u: float
    :type v: float
    :type s: int
    :type z: np.ndarray
    :return: sin + cos function
    :rtype: np.ndarray
    """
    s = float(s)
    return (u*(1.0 + np.cos(s*z))) + (v*(1.0 + np.sin(s*z)))

# INTEGRAL FUNCTIONS

# Integral of exponential functions
# ------------------------------------------

def expfun_integral(A,B):
    """
    integral of A*exp(-B*x) in R3 from 0 to inf
    :param A: coefficient
    :param B: exponent
    :type A: float
    :type B: float
    :return: the integral
    :rtype: float
    """
    return 4*np.pi * A / (B**3)

def xexpfun_integral(A,B):
    """
    integral of A*exp(-B*x) in R3 from 0 to inf
    :param A: coefficient
    :param B: exponent
    :type A: float
    :type B: float
    :return: the integral
    :rtype: float
    """
    return 12.0 * np.pi * A / (B**4)


def x2expfun_integral(A,B):
    """

    :param A:
    :param B:
    :return:
    """
    return 48.0 * np.pi * A/(B**5)

def x2expfun_integral_trigo(fun, G, u, k):
    """

    :param fun:
    :param G:
    :param u:
    :param k:
    :return:
    """
    if fun == np.sin :
        return np.pi * 48.0 * u * (G ** -5) * SINEINTEGRAL[k]
    elif fun == np.cos :
        return np.pi * 48.0 * u * (G ** -5) * COSINEINTEGRAL[k]
    else :
        raise NotImplementedError

def xexpfun_integral_trigo(fun, G,u,k):
    """

    :param fun:
    :param G:
    :param u:
    :param k:
    :return:
    """
    if fun == np.sin :
        return np.pi * 12 * u * (G ** -4) * SINEINTEGRAL[k]
    elif fun == np.cos :
        return np.pi * 12 * u * (G ** -4) * COSINEINTEGRAL[k]
    else :
        raise NotImplementedError

def expfun_integral_trigo(fun, G,u,k):
    """

    :param fun:
    :param G:
    :param u:
    :param k:
    :return:
    """
    if fun == np.sin :
        return np.pi * 8 * u * (G ** -3) * SINEINTEGRAL[k] # need to check it
    elif fun == np.cos :
        return np.pi * 8 * u * (G ** -3) * COSINEINTEGRAL[k] # need to check it
    else :
        raise NotImplementedError

def xexpfun_integral_vk(F,G,u,v,k):
    """
    integral of F*x*exp(-G*x)*sum_psi^Psi(u*(1+cos(psi*z))+v*(1+sin(psi*z)) in R3
    :param F: Radial coefficient
    :param G: Radial exponent
    :param u: Array of sine coefficients
    :param v: Array of cosine coefficients
    :param k: Size of the expansion
    :type F: float
    :type G: float
    :type u: np.ndarray
    :type v: np.ndarray
    :type k: int
    :return: the integral
    :rtype: float
    """
    U = np.zeros(k, dtype='float64')
    for psi in range(k):
        U[psi] = TRIGOEXPANSION[psi](u[psi], v[psi])
    return np.pi*12*F*(G**-4)*U.sum()

# Integral of gaussian functions
# ------------------------------------------

def angular_gaussian_integral(G, u, alpha):
    """

    :param G:
    :param u:
    :param alpha:
    :return:
    """

    factor1 = (np.exp(-1 / (4 * alpha)) * np.sqrt(np.pi)) / (4 * np.sqrt(alpha))
    factor2 = 0
    factor2 += 2 * erfi(1 / (2 * np.sqrt(alpha)))

    inside_term_1 = 1 - (np.pi * alpha * (0 + 2j))
    inside_term_2 = 1 + (np.pi * alpha * (0 + 2j))
    factor2 -= erfi(inside_term_1 / (2 * np.sqrt(alpha)))
    factor2 -= erfi(inside_term_2 / (2 * np.sqrt(alpha)))

    factor3 = 12 * np.pi/(G**4) * u

    return np.real(factor1*factor2*factor3)

def nonfun_integral():
    return 1.0

def sinfun_integral(psi):
    """

    :param psi:
    :return:
    """
    return SINEINTEGRAL[int(psi) - 1]

def cosfun_integral(psi):
    """

    :param psi:
    :return:
    """
    return COSINEINTEGRAL[int(psi) - 1]

# First moments integral

def dipole_gaussian(G, u, alpha):
    """

    :param G:
    :param u:
    :param alpha:
    :return:
    """

    factor1 = (np.exp(-1 / alpha) * np.sqrt(np.pi)) / (8 * np.sqrt(alpha))
    factor2 = 0
    factor2 += 2 * erfi(1 / np.sqrt(alpha))

    inside_term_1 = 1 - (np.pi * alpha * (0+1j))
    inside_term_2 = 1 + (np.pi * alpha * (0+1j))
    factor2 -= erfi(inside_term_1 / (np.sqrt(alpha)))
    factor2 -= erfi(inside_term_2 / (np.sqrt(alpha)))

    factor3 = 12 * np.pi / (G ** 4) * u

    return np.real(factor1 * factor2 * factor3)

# Electrostatic potential expressions

def electrostatic_potential_exp(A, B, d):

    d = filter_min(d)
    term1 = A * 4 * np.pi /(np.power(B, 3) * d)
    term2 = (-2 * np.exp(-B*d)) + 2 - (B * d * np.exp(-B * d))

    return term1 * term2


def ep_xg_radial0(G, d):
    """

    :param G:
    :param d:
    :return:
    """
    # SHORT
    d = filter_min(d)
    y = G * d
    poly = (-6 - np.power(y, 3) - 3 * np.power(y, 2) - 6 * y)
    short = (6 + (np.exp(-y) * poly)) / (np.power(G, 4.0) * d)
    # LONG
    poly = 2 + (2 * G * d) + ((G * d) ** 2)
    long = np.exp(-G * d) * poly * (G ** -3)
    return short + long

def ep_xg_radial1(G, d):
    """

    :param G:
    :param d:
    :return:
    """
    # SHORT
    d = filter_min(d)
    y = G * d
    poly = (-24.0 - np.power(y, 4.0) - 4.0 * np.power(y, 3.0) - 12.0 * np.power(y, 2.0) - 24.0 * y)
    short = ((24 + (np.exp(-y) * poly)) / (np.power(G, 5.0) * np.power(d, 2.0)))
    long = np.exp(-y)*d*(1+y)/(G**2)
    return short + long

def ep_xg_radial2(G, d):
    """

    :param G:
    :param d:
    :return:
    """
    d = filter_min(d)
    y = G * d
    factor = 1 / (np.power(G, 6.0) * np.power(d, 3.0))
    poly = -120.0 - 120 * y - 60.0 * np.power(y, 2.0) \
           - 20.0 * np.power(y, 3.0) - 5.0 * np.power(y, 4.0) - np.power(y, 5.0)
    short = factor * (120 + np.exp(-y) * poly)
    long = np.exp(-y)*(d**2)/G
    return short + long

def ep_xg_radial3(G, d):
    """

    :param G:
    :param d:
    :return:
    """
    d = filter_min(d)
    y = d * G
    # SHORT
    p1 = -np.exp(-y)
    y2 = y * y
    y3 = y2 * y
    y4 = y3 * y
    y5 = y4 * y
    y6 = y5 * y
    # p2 = (-720.0 - y * (720.0 + y * (360.0 + y * (120.0 + y * (30.0 + y * (6 + y))))))
    p2 = 720.0 + 720.0 * y + 360.0 * y2 + 120 * y3 + 30.0 * y4 + 6.0 * y5 + y6
    p3 = (p1 * p2) + 720
    den = ((G**7) * (d**4))
    p3 /= den
    # LONG
    u1 = d ** 3
    u2 = -expi(-y)
    u3 = u1 * u2
    return p3 + u3


def ep_xg_angular0(alpha, z):
    """

    :param alpha:
    :param z:
    :return:
    """
    factor1 = (np.exp(-1 / (4 * alpha)) * np.sqrt(np.pi)) / (4 * np.sqrt(alpha))
    factor2 = 0
    factor2 += 2 * erfi(1 / (2 * np.sqrt(alpha)))

    inside_term_1 = 1 - (np.pi * alpha * (0 + 2j))
    inside_term_2 = 1 + (np.pi * alpha * (0 + 2j))
    factor2 -= erfi(inside_term_1 / (2 * np.sqrt(alpha)))
    factor2 -= erfi(inside_term_2 / (2 * np.sqrt(alpha)))

    return factor1 * factor2 * 0.5


def ep_xg_angular1(alpha, z):
    factor1 = np.exp(-1 / alpha) * np.sqrt(np.pi) / (8 * np.sqrt(alpha))
    it1 = 1 / np.sqrt(alpha)
    it2 = (1 - ((0 + 1j) * alpha * np.pi)) / np.sqrt(alpha)
    it3 = (1 + ((0 + 1j) * alpha * np.pi)) / np.sqrt(alpha)
    factor2 = 2 * erfi(it1) - erfi(it2) - erfi(it3)
    return (3 / 2) * np.cos(z) * factor1 * factor2


def ep_xg_angular2(alpha, z):
    factor = 5 * np.sqrt(np.pi) / (256.0 * np.sqrt(alpha) * np.exp(9 / (4*alpha)))
    factor *= (1.0 + 3.0 * np.cos(2 * z))

    sqa = 2.0 * np.sqrt(alpha)

    it1 = (3.0j + (2.0 * alpha * np.pi)) / sqa
    it2 = (3.0 + (2.0 * alpha * np.pi * 1j)) / sqa
    it3 = (1.0 - (np.pi * 2j * alpha)) / sqa
    it4 = (1.0 + (np.pi * 2j * alpha)) / sqa

    sum1 = 6 * erfi(3.0 / sqa) + (3 * 1j * erf(it1)) - (3 * erfi(it2))
    sum2 = -2 * erfi(1.0 / sqa) + erfi(it3) + erfi(it4)
    return factor * (sum1 + np.exp(2/alpha) * sum2)

def ep_xg_angular3(alpha, z):

    factor = 7 * np.exp(-4.0 / alpha) * np.sqrt(np.pi) * np.cos(z) / (256.0 * np.sqrt(alpha))
    factor *= (-3 + 5*np.power(np.cos(z), 2.0))

    it1 = (2.0j + alpha * np.pi) / np.sqrt(alpha)
    it2 = 1/np.sqrt(alpha)
    it3 = 2/np.sqrt(alpha)
    it4 = (1 - 1.0j*alpha*np.pi)/np.sqrt(alpha)
    it5 = (1 + 1.0j * alpha * np.pi) / np.sqrt(alpha)
    it6 = 2.0 + (1.0j * alpha * np.pi)

    sum1 = 5.0j * erf(it1)
    sum1 -= (4.0 * np.exp(3/alpha))*erfi(it2) - 10.0 * erfi(it3)
    sum1 += 2*np.exp(3/alpha) * (erfi(it4)+erfi(it5))
    sum1 -= 5.0 * erfi(it6)
    return factor*sum1


def electrostatic_potential_xexp_gaussian(G, alpha, d, z):
    """

    :param d:
    :param z:
    :param G:
    :param alpha:
    :return:
    """
    # Solving the electrostatic potential requires calculating by one side
    # the integral of the radial part, and by the other, the integral of the gaussian part

    radial_terms = [ep_xg_radial0, ep_xg_radial1, ep_xg_radial2, ep_xg_radial3]
    angular_terms = [ep_xg_angular0, ep_xg_angular1, ep_xg_angular2, ep_xg_angular3]
    electrostatic_potential_buffer = np.zeros(d.shape[0], dtype='float64')
    P = lambda order : (4 * np.pi)/((2 * order) + 1)
    for l, (rad, ang) in enumerate(zip(radial_terms, angular_terms)):
        u_rad = rad(G,d)
        u_ang = np.real(ang(alpha, z))
        u_p = P(l)
        electrostatic_potential_buffer = electrostatic_potential_buffer + (u_rad * u_ang * u_p)
    return electrostatic_potential_buffer

def trapezoidal_3d_integral(fun, domain, resolution, sign_operator=False):
    """

    :param fun:
    :param domain:
    :param resolution:
    :param sign_operator:
    :return:
    """
    h = resolution
    xx = np.arange(domain[0], domain[1] + h, h)
    yy = np.arange(domain[2], domain[3] + h, h)
    zz = np.arange(domain[4], domain[5] + h, h)

    zz_m = np.array([domain[4], domain[4] + h], dtype='float64')
    I = 0.0
    N = 0.0
    h3 = h * h * h / 8
    #    X, Y, Z = np.meshgrid(xx, yy, zz_m)
    P = np.zeros((xx.size, yy.size, zz_m.size))
    XYZ = np.zeros((xx.size, 3), dtype='float64')
    for iy in range(yy.size):
        XYZ[:, 0] = xx
        XYZ[:, 1] = yy[iy]
        XYZ[:, 2] = zz[0]
        P[:, iy, 0] = fun(XYZ)
        XYZ[:, 2] = zz[1]
        P[:, iy, 1] = fun(XYZ)

    for iz in range(1, zz.size):

        P_av = np.zeros((P.shape[0] - 1, P.shape[1] - 1), dtype='float64')
        P_av[:, :] += P[:-1, :-1, -1]  # lll
        P_av[:, :] += P[:-1, :-1, 1]  # llu
        P_av[:, :] += P[:-1, 1:, -1]  # lul
        P_av[:, :] += P[:-1, 1:, 1]  # luu
        P_av[:, :] += P[1:, :-1, -1]  # ull
        P_av[:, :] += P[1:, :-1, 1]  # ulu
        P_av[:, :] += P[1:, 1:, -1]  # uul
        P_av[:, :] += P[1:, 1:, 1]  # uuu

        if sign_operator:
            N += np.sum(P_av[P_av < 0.0]) * h3
        I += np.sum(P_av) * h3

        for iy in range(yy.size):
            XYZ[:, 0] = xx
            XYZ[:, 1] = yy[iy]
            XYZ[:, 2] = zz[iz]
            if iz % 2 == 0:
                P[:, iy, 0] = fun(XYZ)
            else:
                P[:, iy, 1] = fun(XYZ)

    if sign_operator:
        return I, N
    else:
        return N

RADIAL_PARAMS = {
    expfun: ['A', 'B'],
    xexpfun: ['A', 'B'],
    x2expfun: ['A', 'B']
}

ANGULAR_PARAMS = {
    cosssin: ['u', 'v', 's'],
    sinfun: ['psi'],
    cosfun: ['psi'],
    nonefun: [],
    gaussian : ['h', 'alpha']
}

RADIAL_INTEGRAL_FUNCTIONS = {
    expfun: expfun_integral,
    xexpfun: xexpfun_integral,
    x2expfun: x2expfun_integral
}

ANGULAR_INTEGRAL_FUNCTIONS = {
    sinfun: sinfun_integral,
    cosfun: cosfun_integral,
    nonefun: nonfun_integral
}

