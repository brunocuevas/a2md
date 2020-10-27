from typing import TypeVar, List
T = TypeVar('T')


class Unit:
    def __init__(self, fullname: str, shortname: str, reference: T, value: float):
        """
        A2MDio Unit. Its purpose is to ensure unit converssion by type checking.
        :param fullname: name of the unit
        :param reference: reference of the unit
        :param value: value of the unit in reference units

        For instance:
        meter = Unit(fullname='metre', reference='meter', value=1.0)
        km = Unit(fullname='kilometre', reference='meter', value=1000)

        """

        self.fullname = fullname
        self.shortname = shortname
        self.reference_unit = reference
        self.value = value

    def __str__(self):

        return self.fullname

    def __format__(self, format_spec):
        return format(str(self), format_spec)

    def __truediv__(self, other) -> float:
        if self.reference_unit == other.reference_unit:
            ratio = self.value / other.value
            return ratio

    def __mul__(self, other) -> float:
        if self.reference_unit == other.reference_unit:
            ratio = self.value * other.value
            return ratio


class UnitCollection:
    def __init__(self, name, units: List[Unit]):
        self.units = units
        self.name = name
        self.keys = [str(i) for i in self.units]
    def unit_from_name(self, name):
        try:
            return self.units[self.keys.index(name)]
        except ValueError:
            KeyError("unknown unit in collection")


bohr = Unit(fullname='bohr', shortname='bohr', reference='bohr', value=1.0)
angstrom = Unit(fullname='angstrom', shortname='A', reference='bohr', value=1.889725989)
nanometer = Unit(fullname='nanometer', shortname='nm', reference='bohr', value=18.89725989)


atomic_length = UnitCollection(
    name='atomic lenght units', units=[bohr, angstrom, nanometer]
)