from a2md.support import SupportAngular, SupportRadial
import numpy as np


carbon_coords = np.array([2.64363, 0.0, 0.0])
hydrogen_coords = np.array([4.68918, 0.0, 0.0])

carbon_core = SupportRadial(
    coordinates=carbon_coords,
    A3 = 128.049,
    B3 = 11.85
)

carbon_core_valence = SupportRadial(
    coordinates=carbon_coords,
    A3 = -2.535,
    B3 = 3.508
)

carbon_valence = SupportRadial(
    coordinates=carbon_coords,
    A3 = 2.041,
    B3 = 2.099
)

hydrogen_core = SupportRadial(
    coordinates=hydrogen_coords,
    A3=0.62,
    B3=2.5
)

carbon_carbon1 = SupportAngular(
    coordinates=carbon_coords,
    Alpha=3.5,
    U=1.29,
    G=2.1,
    Psi=0
)
carbon_carbon2 = SupportAngular(
    coordinates=carbon_coords,
    Alpha=1.75,
    U=1.29,
    G=2.1,
    Psi=1
)

carbon_hydrogen1 = SupportAngular(
    coordinates=carbon_coords,
    Alpha=3.5,
    U=1.29,
    G=2.1,
    Psi=0
)

carbon_hydrogen2 = SupportAngular(
    coordinates=carbon_coords,
    Alpha=1.75,
    U=1.29,
    G=2.1,
    Psi=1
)

hydrogen_carbon1 = SupportAngular(
    coordinates=hydrogen_coords,
    Alpha=0.6,
    U=0.92,
    G=2.95,
    Psi=1
)

hydrogen_carbon2 = SupportAngular(
    coordinates=hydrogen_coords,
    Alpha=0.3,
    U=0.92,
    G=2.95,
    Psi=1
)


carbon_hydrogen1.set_reference_frame(np.array([1.0, 0.0, 0.0]))
carbon_hydrogen2.set_reference_frame(np.array([1.0, 0.0, 0.0]))
hydrogen_carbon1.set_reference_frame(np.array([-1.0, 0.0, 0.0]))
hydrogen_carbon2.set_reference_frame(np.array([-1.0, 0.0, 0.0]))
carbon_carbon1.set_reference_frame(-np.array([np.cos(np.pi/3), np.sin(np.pi/3), 0.0]))
carbon_carbon2.set_reference_frame(-np.array([np.cos(np.pi/3), np.sin(np.pi/3), 0.0]))

probe = np.zeros((1,3), dtype='float64')


ep_functions = [
    carbon_core, carbon_core_valence, carbon_valence,
    hydrogen_core,
    carbon_hydrogen1, carbon_hydrogen2,
    carbon_carbon1, carbon_carbon2,
    hydrogen_carbon1, hydrogen_carbon2
]

# Ccore, Ccv, Cval, Hcore, CH1, CH2, CC1, CC2, HC1, HC2
coeffs = np.array(
    [1.0, #  Ccore
    0.3463157062754157, #  Ccv
    0.5859931990427911,  #  Cval
    0.6594062329968, # Hcore
    0.33566877777753296,  # CH1
    0.6124181453216823,  # CH2
    0.12611667604194682 * 2,  # CC1
    0.711623485801416 * 2,   # CC2
    0.31466244934340654,  # HC1
    0.11741622036555388], dtype='float64' # HC2
)

names =  'Ccore Ccv Cval Hcore CH1 CH2 CC1 CC2 HC1 HC2'.split()
eps = np.array([i.eval_ep(probe) for i in ep_functions]).flatten()
pos_eps = (1.0/np.linalg.norm(hydrogen_coords - probe)) + (6.0/np.linalg.norm(carbon_coords - probe))
mult = eps*coeffs
print("negative ep is {:18.12f}".format(mult.sum()))
print("positive ep is {:18.12f}".format(pos_eps))
print("total ep is {:18.12f}".format(pos_eps - mult.sum()))

for name, value in zip(names, mult):
    print("{:<6s} {:18.6f}".format(name, value))


