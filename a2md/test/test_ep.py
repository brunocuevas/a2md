import json
import numpy as np
from a2md.mathfunctions import short_generalized, long_generalized, pe_harmonic
from a2md.mathfunctions import spherical_harmonic
from a2md.mathfunctions import electrostatic_potential_exp
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

    # CASE 1
    r = 2.25
    P = 1.0
    l = 1
    B = 2.24
    t = np.pi/3

    V = pe_harmonic(r,t, l, P, B)
    print(
        "r={:3.2f} t={:3.2f} P=1, l=1, B={:3.2f} -> V={:18.12f}".format(
            r, t, B, V
        )
    )

    r = 2.25
    P = 0.0
    l = 2
    B = 2.24
    t = np.pi/3

    V = pe_harmonic(r,t, l, P, B)
    print(
        "r={:3.2f} t={:3.2f} P=0, l=2, B={:3.2f} -> V={:18.12f}".format(
            r, t, B, V
        )
    )


    r = 2.25
    P = 1
    l = 3
    B = 2.24
    t = np.pi/3

    V = pe_harmonic(r,t, l, P, B)
    print(
        "r={:3.2f} t={:3.2f} P=1, l=3, B={:3.2f} -> V={:18.12f}".format(
            r, t, B, V
        )
    )

    V = electrostatic_potential_exp(A=1.0, B=B, d=r)
    print(
        "r={:3.2f} t={:3.2f} P=0, l=0, B={:3.2f} -> V={:18.12f}".format(
            r, t, B, V
        )
    )