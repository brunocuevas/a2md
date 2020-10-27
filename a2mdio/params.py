import json
from typing import Dict, Union
from pathlib import Path
from a2mdio.units import atomic_length


class AMDParameters:

    primary_template = {
        "_NAME": str,
        "_VERSION": str,
        "_DESCRIPTION": str,
        "_MODEL": dict,
        "_MAXFUNS": int,
        "_NELEMENTS": int,
        "_UNITS": str
    }

    function_template = {
        "_NAME": str,
        "_PARAMS": dict,
        "_FROZEN": bool,
        "_CONNECT": str,
        "_TYPE": str
    }

    def __init__(self, contents):
        AMDParameters.check_fields(contents)
        self.contents = contents

    def get_maxfunctions(self):
        return self.contents['_MAXFUNS']

    def get_included_elements(self):
        return list(self.contents['_MODEL'].keys())

    def get_nelements(self):
        return int(self.contents['_NELEMENTS'])

    def get_element(self, element: str):
        try:
            return self.contents['_MODEL'][element]
        except KeyError:
            raise IOError("element {:s} is not included".format(element))

    def get_units(self):
        return atomic_length.unit_from_name(self.contents['_UNITS'])

    def iter_element(self, element: str):
        try:
            funs = self.contents['_MODEL'][element]
        except KeyError:

            raise IOError("element {:s} is not included".format(element))
        for f in funs:
            yield f['_PARAMS']

    def keep_frozen(self):
        keep = dict()
        max_elements = 0
        for element, fun_list in self.contents['_MODEL'].items():
            keep[element] = []
            n_elements = 0
            for f in fun_list:
                if f['_FROZEN']:
                    keep[element].append(f)
                    n_elements += 1
            if max_elements < n_elements:
                max_elements = n_elements

        contents = self.contents.copy()
        contents['_MODEL'] = keep
        contents['_MAXFUNS'] = max_elements
        return AMDParameters(contents)

    def remove_frozen(self):
        keep = dict()
        max_elements = 0
        for element, fun_list in self.contents['_MODEL'].items():
            keep[element] = []
            n_elements = 0
            for f in fun_list:
                if not f['_FROZEN']:
                    keep[element].append(f)
                    n_elements += 1
            if max_elements < n_elements:
                max_elements = n_elements

        contents = self.contents.copy()
        contents['_MODEL'] = keep
        contents['_MAXFUNS'] = max_elements
        return AMDParameters(contents)

    @staticmethod
    def from_file(filename: Union[str, Path]):

        with open(filename) as f:

            contents = json.load(f)

        return AMDParameters(contents)

    @staticmethod
    def check_fields(contents: Dict):

        for key, item_type in AMDParameters.primary_template.items():

            try:
                item = contents[key]
            except KeyError:
                raise IOError("missing field {:s} in primary keys".format(key))

            if type(item) is not item_type:
                raise TypeError("item {:s} should be type {:s}".format(key, str(item_type)))

        model = contents['_MODEL']
        for element_symbol, element_model in model.items():
            for fun in element_model:
                for key, item_type in AMDParameters.function_template.items():

                    try:
                        item = fun[key]
                    except KeyError:
                        raise IOError("missing field {:s} within models".format(key))

                    if type(item) is not item_type:
                        raise TypeError("item {:s} should be type {:s}".format(key, str(item_type)))
        return True
