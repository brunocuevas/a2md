from a2mdio.units import Unit, bohr, angstrom, nanometer


if __name__ == '__main__':
    m = Unit(fullname='metre', shortname='m', reference='metre', value=1.0)
    km = Unit(fullname='kilometre', shortname='km', reference='metre', value=1000.0)
    print('1 km is {:8.4f} m'.format(km / m))
    print('1 m is {:8.4f} km'.format(m / km))

    print('1 angstrom is : {:8.4e} bohrs'.format(angstrom / bohr))
    print('1 bohr is : {:8.4e} bohrs'.format(bohr / bohr))
    print('1 bohr is : {:8.4e} angstroms'.format(bohr / angstrom))
    print('1 {:s} is : {:8.4e} angstroms'.format(nanometer, nanometer / angstrom))