import re

def std_coordinates(lines):
    coordinates_flag = False
    coordinates_barrier = 0
    coordinates = []
    atom_numbers = []
    for line in lines:
        line = line.strip()
        if re.match('Standard\sorientation', line):

            coordinates_flag = True
            continue
        if coordinates_flag:
            if re.match('-{5}', line):
                if coordinates_barrier == 2:
                    return atom_numbers, coordinates
                else:
                    coordinates_barrier += 1
            elif re.match('Center\s*Atomic\s*Atomic\s*Coordinates\s\(Angstroms\)', line):
                pass
            elif re.match('Number\s*Number\s*Type\s*X\s*Y\s*Z', line):
                pass
            else:
                cnter, atmnmbr, atmtype, x, y, z = line.split()
                atmnmbr = int(atmnmbr)
                coordinates.append([
                    float(x),
                    float(y),
                    float(z)
                ])
                atom_numbers.append(atmnmbr)

def mk_charges(lines):
    flagged_line = None
    charge = []
    for i, line in enumerate(lines):
        line = line.strip()
        if re.search(r"Fitting", line):
            flagged_line = i
            break
    if flagged_line is None:
        raise RuntimeError("MK charges not found")
    flagged_line += 4

    for line in lines[flagged_line:]:
        line = line.strip()
        if line[0] == "-":
            break
        charge.append(float(line.split()[2]))
    return charge

def npa_charges(lines):

    ncharges = []
    flag = False
    for line in lines:
        line = line.strip()
        if re.match('\s*Atom\s*No\s*Charge\s*Core\s*Valence\s*Rydberg\s*Total', line):
            flag=True
        elif re.match('-{5}', line) and flag:
            pass
        elif re.match('={5}', line) and flag:
            return ncharges
        elif flag:
            atmsymbl, no, charge, core, valence, rydberg, total = line.split()
            ncharges.append(float(charge))
    raise RuntimeError("could not find npa charges")

def hf_energy(lines):
    for line in lines:
        if re.match('SCF\sDone:\s*E\(RHF\)\s=(.*)\s*A\.U\.\safter\s*\d*\scycles', line):
            m = re.match('SCF\sDone:\s*E\(RHF\)\s=(.*)\s*A\.U\.\safter\s*\d*\scycles', line)
            return float(m.group(1).replace('D', 'E'))
    raise RuntimeError("could not find HF energy")

def dft_energy(lines, functional):
    x = 'SCF\sDone:\s*E\(R{:s}\)\s=(.*)\s*A\.U\.\safter\s*\d*\scycles'.format(functional)
    print(x)
    for line in lines:
        line = line.strip()
        if re.match(x, line):
            m = re.match(x, line)
            return float(m.group(1).replace('D', 'E'))
    raise RuntimeError("could not find DFT energy")

def mp2_energy(lines):
    for line in lines:
        line = line.strip()
        if re.match('E2\s=\s*(.*)\sEUMP2\s=\s*(.*)', line):
            m = re.match('E2\s=\s*(.*)\sEUMP2\s=\s*(.*)', line)
            return  float(m.group(2).replace('D', 'E'))
    raise RuntimeError("could not find MP2 energy")

def energy_decomposition(lines):
    for i, line in enumerate(lines):
        line = line.strip()
        if re.match(r'N-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*E-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*KE=\s*(-?\d*.\d*D[+,-]\d{2})\s*',
                    line):
            m = re.match(r'N-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*E-N=\s*(-?\d*.\d*D[+,-]\d{2})\s*KE=\s*(-?\d*.\d*D[+,-]\d{2})\s*',
                         line)

            nn = float(m.group(1).replace('D', 'E'))
            ne = float(m.group(2).replace('D', 'E'))
            kin = float(m.group(3).replace('D', 'E'))

            energy_terms = dict(
                nuclei_nuclei_potential=nn,
                nuclei_electron_potential=ne,
                kinetic=kin
            )

            return energy_terms

    raise RuntimeError("could not find the decomposition of energy")

def dipole(lines):
    for i, line in enumerate(lines):
        line = line.strip()
        if re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line):
            m = re.match(r'^\s*X=\s*(-?\d*.\d*)\s*Y=\s*(-?\d*.\d*)\s*Z=\s*(-?\d*.\d*)', line)
            dx = float(m.group(1))
            dy = float(m.group(2))
            dz = float(m.group(3))
            dipole = [dx, dy, dz]
            return dipole
    raise RuntimeError("could not find dipole")

def electrostatic_potential(lines):
    flag = False
    ep = []
    start = None
    for i, line in enumerate(lines):
        line = line.strip()
        if re.match('Electrostatic Properties \(Atomic Units\)', line):
            start = i
            flag = True
            break
    if start is None:
        raise RuntimeError("could not find ep")

    if flag:
        for line in lines[start + 6:]:
            line = line.strip()
            if re.search(r'Atom', line):
                continue
            elif re.search(r'------', line):
                break
            else:
                ep.append(float(line.split()[1]))

    return ep

def forces(lines, symmetry=True):
    flag = False
    force_ = []
    center_ = []
    start = None
    for i, line in enumerate(lines):
        line = line.strip()
        if re.match('Calling FoFJK', line):
            if symmetry:
                start = i + 6
            else:
                start = i + 5
            flag = True
            break
    if flag:
        for line in lines[start:]:
            if re.search(r'---', line):
                break
            else:
                cnt, _, x, y, z = line.split()
                x = float(x)
                y = float(y)
                z = float(z)
                force_.append([x, y, z])
                center_.append(int(cnt))

    return center_, force_