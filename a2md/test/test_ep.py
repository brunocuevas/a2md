import json
import numpy as np
from a2md.mathfunctions import short_generalized, long_generalized, pe_harmonic
from a2md.mathfunctions import spherical_harmonic
import matplotlib.pyplot as plt

if __name__ == '__main__':
    t = np.arange(0.0, 2.1, 0.1) * np.pi
    sph0 = spherical_harmonic(t, l=0)
    sph1 = spherical_harmonic(t, l=1)
    sph2 = spherical_harmonic(t, l=2)
    sph3 = spherical_harmonic(t, l=3)

    r = np.arange(1.0, 5.0, 0.1)
    B = 1.0
    l = 1.0
    P = 0.0
    sg1 = short_generalized(r, B, l, P)
    sg2 = long_generalized(r, B, l, P)

    plt.plot(r, sg1, marker='.')
    plt.plot(r, sg2, marker='x')
    plt.show()