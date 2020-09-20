import numpy as np
from scipy.special import erfi, erf, expi, gamma, gammaincc, exp1

INVERSE_DIST_DUMPING = 1e-3
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

def get_polar_rep(x, center, ref_frame = None):
    """
    calculates the angle of the vector x-center against
    the z axis and module of such vector.

    the reference frame allows to rotate the system to
    an arbitrary z vector (for instance, a bonding vector).

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


def generalized_exponential(A, B, d, P=0):
    """
    generalized exp
    Formula:
    A Exp[-B*d] * (d^p)

    :param A: Coefficient
    :param B: Exponent
    :param d: distance
    :param P: polynomial degree
    :return:
    """
    return A * np.power(d, P) * np.exp(-B * d)

def generalized_exponential_integral(A, B, P=0):
    """
    generalized exponential integral
    Formula:
    int[A d^{P+2} Exp[-B d] dd] = A(P+2!)/(B^{P+1})
    """
    return (A * np.power(B, -3-P)) * gamma(3 + P)


def gaussian(alpha, h, d):
    """
    gaussian function
    Formula:
    h Exp[-alpha * (d^2)]
    """
    return h * np.exp(-alpha * (d**2))

# Trigonometric functions

def nonefun(z):
    return 1.0

def angular_gaussian_integral(alpha):
    """
    angular gaussian integral
    Formula:
    int[G[alpha, theta]*Sin(theta) dtheta dphi, {dtheta, 0, Pi}, {dphi, 0, 2Pi}]

    """

    factor1 = (np.exp(-1 / (4 * alpha)) * np.sqrt(np.pi)) / (4 * np.sqrt(alpha))
    factor2 = 0
    factor2 += 2 * erfi(1 / (2 * np.sqrt(alpha)))

    inside_term_1 = 1 - (np.pi * alpha * (0 + 2j))
    inside_term_2 = 1 + (np.pi * alpha * (0 + 2j))
    factor2 -= erfi(inside_term_1 / (2 * np.sqrt(alpha)))
    factor2 -= erfi(inside_term_2 / (2 * np.sqrt(alpha)))

    return np.real(factor1*factor2)

def nonfun_integral():
    return 1.0

# First moments integral

def dipole_gaussian(G, u, alpha):
    """
    Obtains the first moment of function:

    rho[d] = (u d Exp[-G d - alpha * (theta^2)])

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


# Electrostatic potential of harmonic functions
# harmonics
def yl1m0(t):
    """
    harmonic l=1, m=0
    """
    return np.cos(t)

def yl2m0(t):
    """
    harmonic l=2, m=0
    """
    x = np.cos(t)
    return 0.5 * (3 * (x * x) - 1)

def yl3m0(t):
    """
    harmonic l=3, m=0
    """
    x = np.cos(t)
    return 0.5 * (5 * (x ** 3) - 3 * x)

def inc_gamma(a, x):
    """
    incomplete gamma definition
    """
    return exp1(x) if a == 0 else gamma(a)*gammaincc(a, x)

def spherical_harmonic(t, l):
    if l == 0:
        return np.ones(t.size, dtype='float64')
    elif l == 1:
        return yl1m0(t)
    elif l == 2:
        return yl2m0(t)
    elif l == 3:
        return yl3m0(t)
    else:
        raise NotImplementedError("only implemented up to l = 3")

# short and long terms
def short_generalized(r, B, l, P):
    term1 = r ** (2 + P)
    term2 = (B * r) ** (-3 - l - P)
    term31 = gamma(3 + l + P)
    term32 = gammaincc(3 + l + P, B * r) * gamma(3 + l + P)
    return term1 * term2 * (term31 - term32)

def long_generalized(r, B, l, P):
    if 2 - l + P > 0:
        factor = (r ** l) / (B ** 2)
        term1 = (B ** (l - P)) * gamma(2 - l + P)
        term21 = (r ** (-l + P)) * ((B * r)**(l - P))
        term22 = -gamma(2 - l + P) + gammaincc(2-l+P, B*r) * gamma(2-l+P)
        return (term1 + term21 * term22) * factor
    else:
        if 2 - l + P == 0:
            return (r ** l) * inc_gamma(0, B * r)
        if 2 - l + P == -1:
            return (r ** l) * ((np.exp(-B*r)/r) - B * inc_gamma(0, B*r))

def pe_harmonic(r, t, l, P, B):
    factor = (4*np.pi)/((2 * l) + 1)
    radial = short_generalized(r, B, l, P) + long_generalized(r, B, l, P)
    u = np.cos(t)
    if l == 1 :
        return u * factor * radial
    elif l == 2 :
        return factor * radial * ((3.0 * u**2) - 1.0) / 2.0
    elif l == 3 :
        return factor * radial * ((5.0 * u**3) - (3.0 * u)) / 2.0
    else:
        raise NotImplementedError("only implemented for 0 < l < 4")