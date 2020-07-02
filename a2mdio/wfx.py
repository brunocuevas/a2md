from a2mdio.qm import WaveFunction, symetry_index
import re
import numpy as np


def parse_scalar_text(contents):
    return contents[0].strip()

def parse_scalar_int_number(contents):
    return int(contents[0].strip())

def parse_vector_int_number(contents):
    values = []
    for line in contents:
        values += line.strip().split()
    values = [int(i) for i in values]
    return values

def parse_vector_float_number(contents):
    values = []
    for line in contents:
        values += line.strip().split()
    values = [float(i) for i in values]
    return values

def parse_vector_text(contents):
    values = []
    for line in contents:
        values += line.strip().split()
    return values

def parse_orbital_coefficients(contents):
    coefficients = []
    coefficients_split = '\n'.join(contents).split('<MO Number>')[1:]
    for i, block in  enumerate(coefficients_split):
        block = block.split('\n')[3:]
        x = parse_vector_float_number(block)
        coefficients.append(x)
    return coefficients

def parse_coordinates(contents):
    coordinates = []
    for line in contents:
        coordinates.append([float(i) for i in line.split()])
    return coordinates

fields = dict(
    title = dict(label='Title', fun=parse_scalar_text),
    ncenters = dict(label='Number of Nuclei', fun=parse_vector_int_number),
    norb = dict(label='Number of Occupied Molecular Orbitals', fun=parse_scalar_int_number),
    nelec = dict(label='Number of Electrons', fun=parse_scalar_int_number),
    naelec = dict(label='Number of Alpha Electrons', fun=parse_scalar_int_number),
    nbelec = dict(label='Number of Beta Electrons', fun=parse_scalar_int_number),
    mult = dict(label='Electron Spin Multiplicity', fun=parse_scalar_int_number),
    ncelec = dict(label='Number of Core Electrons', fun=parse_scalar_int_number),
    types = dict(label='Nuclear Names', fun=parse_vector_text),
    atomn = dict(label='Nuclear Charges', fun=parse_vector_float_number),
    coordinates = dict(label='Nuclear Cartesian Coordinates', fun=parse_coordinates),
    nprims = dict(label='Number of Primitives', fun=parse_scalar_int_number),
    centers = dict(label='Primitive Centers', fun=parse_vector_int_number),
    sym = dict(label='Primitive Types', fun=parse_vector_int_number),
    exp = dict(label='Primitive Exponents', fun=parse_vector_float_number),
    occupation = dict(label='Molecular Orbital Occupation Numbers', fun=parse_vector_float_number),
    energies = dict(label='Molecualr Orbital Energies', fun=parse_vector_float_number),
    spins = dict(label='Molecular Orbital Spin Types', fun=parse_vector_text),
    coeff = dict(label='Molecular Orbital Primitive Coefficients', fun=parse_orbital_coefficients)
)

def split_by_labels(contents):
    flag = False
    block = None
    label = None
    for i, line in enumerate(contents):
        line = line.strip()
        if line[0] == '<' and line[1] != '/' and flag == False:
            flag = True
            block = []
            label = line.strip().replace('<', '').replace('>', '')
        elif line == '</{:s}>'.format(label) and flag == True:
            flag = False
            yield label, block
        elif flag:
            block.append(line)

def match_block(label, block, fields_dict):
    for key, item in fields_dict.items():
        if item['label'] == label:
            return key, item['fun'](block)
    return None, None

class WaveFunctionX(WaveFunction):

    def __init__(self, file, verbose=False, prefetch_dm=False, batch_size=10000):
        WaveFunction.__init__(self, file=file, verbose=verbose, prefetch_dm=True, batch_size=batch_size)
        if prefetch_dm:
            self.calculate_density_matrix()

    def read(self):
        with open(self._WaveFunction__file) as f:
            contents = f.readlines()

        contents = [i.strip() for i in contents]
        contents = [i for i in contents if i[0] != '#']
        wavefunction = dict()
        for label, block in split_by_labels(contents):
            l, u = match_block(label, block, fields)
            if u is not None:
                wavefunction[l] = u

        coords = np.array(wavefunction['coordinates'], dtype='float64')
        exp = np.array(wavefunction['exp'], dtype='float64')
        cent = np.array(wavefunction['centers'], dtype='int64') - 1
        coeff = np.array(wavefunction['coeff'], dtype='float64')
        types = wavefunction['types']
        nprims = wavefunction['nprims']
        occ = wavefunction['occupation']
        norbs = wavefunction['norb']
        charges = wavefunction['atomn']
        sym = np.array(wavefunction['sym'], dtype='int64') - 1
        ncenters = wavefunction['ncenters']
        return coeff, occ, exp, sym, cent, coords, types, charges, nprims, norbs, ncenters