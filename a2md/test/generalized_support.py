import numpy as np
from a2md.support import support_agnostic
from a2md.mathfunctions import RADIAL_INTEGRAL_FUNCTIONS
from a2md.mathfunctions import RADIAL_PARAMS
from a2md.mathfunctions import expfun, xexpfun
from a2md.mathfunctions import cosfun, sinfun
from a2md.mathfunctions import nonefun



def expsum(d,A1,B1,A2,B2):
    return A1 * np.exp(-B1*d) + A2 * np.exp(-B2*d)

def expsum_integral(A1,B1,A2,B2):
    return 8 * np.pi * ((A1/(B1**3))+(A2/(B2**3)))

RADIAL_PARAMS[expsum] = ['A1', 'A2', 'B1', 'B2']
RADIAL_INTEGRAL_FUNCTIONS[expsum] = expsum_integral


E1 = lambda args : support_agnostic(name='E1', angular=nonefun, radial=expsum, **args)
E2 = lambda args : support_agnostic(name='E2', angular=nonefun, radial=expfun, **args)
F1C = lambda args : support_agnostic(name='F1C', angular=cosfun, radial=xexpfun, **args)
F1S = lambda args : support_agnostic(name='F1S', angular=sinfun, radial=xexpfun, **args)
F2C = lambda args : support_agnostic(name='F2C', angular=cosfun, radial=xexpfun, **args)
F2S = lambda args : support_agnostic(name='F2S', angular=sinfun, radial=xexpfun, **args)