from a2md.baseclass import A2MDBaseClass
import numpy as np
from a2md.mathfunctions import get_polar_rep
from a2md.mathfunctions import generalized_exponential, generalized_exponential_integral
from a2md.mathfunctions import gaussian, angular_gaussian_integral
from a2md.mathfunctions import electrostatic_potential_exp, electrostatic_potential_xexp_gaussian
from a2md.mathfunctions import spherical_harmonic, pe_harmonic


class Support(A2MDBaseClass) :
    def __init__(self, **kwargs):
        r"""
        Support functions. It requires:
        - coordinates : where function is centered
        - radial function type : 0 - sum of exp, 1 - exp, 2 - xexp
        - angular function type : 0 - isotropic, 1 - expansion of size 1, 2- expansion of size 2...
        :param \**kwargs['coordinates']:
            See below
        :keyword Arguments:
        * *coordinates* *(``np.ndarray``)
        * *radial_function_type* *(``int``)
        * *angular_function_type* *(``int``)

        """
        A2MDBaseClass.__init__(self, name ='support function ', verbose=False)
        self.coordinates = kwargs['coordinates']
        self.coordinate_system = None
        self.eval_method = None
        self.eval_ep_method = None
        self.integral_method = None
        self.support_type = None
        self.params_kw = None
        self.anisotropic = False

    def eval(self, x):
        """
        provides evaluation at coordinates x

        :param x: cartesian coordinates
        :type x: ndarray
        :return: density values
        :rtype: ndarray
        """
        if self.eval_method is None:
            raise NotImplementedError
        else:
            return self.eval_method(x)

    def eval_ep(self, x):
        """
        provides evaluation at coordinates x

        :param x: cartesian coordinates
        :type x: ndarray
        :return: density values
        :rtype: ndarray
        """
        if self.eval_ep_method is None:
            raise NotImplementedError
        else:
            return self.eval_ep_method(x)

    def get_params(self):
        """

        Returns a dictionary with the parameters of the model and their names.
        Note: This names must be stored in params_kw. If for some reason the variables
        that are needed are different than those specified (as for instance, is coordinates
        of bonding atom are needed, it may be possible just to change the arguments of this list

        :return: a dictionary of paramters
        :rtype: dict
        """
        if self.params_kw is None:
            raise NotImplementedError
        else:
            params = dict(
                (i, getattr(self, '_' + self.__class__.__name__ + '__' + i)) for i in self.params_kw
            )
            return params

    def integral(self):
        """
        provides integral
        :return: density integral
        :rtype: np.ndarray
        """
        if self.integral_method is None:
            raise NotImplementedError
        else:
            return self.integral_method()

    def is_anisotropic(self):
        return self.anisotropic

    def set_reference_frame(self, w):
        """
        creates a new coordinates reference system to ease calculations
        :param w: vector to align coordinates
        :type w: np.ndarray
        :return: True
        """
        w /= np.linalg.norm(w)
        B = np.zeros((3, 3), dtype='float64')
        new_x = np.array([-w[1]-w[2], w[0], w[0]], dtype='float64')
        new_x /= np.linalg.norm(new_x)
        new_y = np.cross(w, new_x)
        B[:, 0] = new_x
        B[:, 1] = new_y
        B[:, 2] = w
        try:
            self.coordinate_system = np.linalg.inv(B)
        except np.linalg.LinAlgError:
            new_x = np.cross(w, np.array([1, 0, 0], dtype='float64'))
            new_y = np.cross(w, new_x)
            B[:, 0] = new_x
            B[:, 1] = new_y
            B[:, 2] = w
            self.coordinate_system = np.linalg.inv(B)
        return True


class SupportRadial(Support):
    """
    support - outer radial function
        \ro(r) = A3 * e^-(B3*r)
    """
    def __init__(self, **kwargs):
        r"""
        Support functions. It requires:
        - coordinates : where function is centered
        - A : coefficient
        - B : exponent
        :param \**kwargs['coordinates']:
            See below
        :keyword Arguments:
        - A: coefficients
        - B: exponents
        """
        Support.__init__(self, **kwargs)
        self.set_name(
            'radial function at %04.2f:%04.2f:%04.2f' % (
                kwargs['coordinates'][0], kwargs['coordinates'][1], kwargs['coordinates'][2]
            )
        )
        try:
            self.__A = kwargs['A']
            self.__B = kwargs['B']
        except KeyError:
            raise IOError("missing parameters A3|B3|P")
        self.eval_method = self.__eval_outer
        self.integral_method = self.__integral_outer
        self.eval_ep_method = self.__eval_ep_outer
        self.support_type = 'radial'
        self.params_kw = ['A', 'B']
        self.anisotropic = False

    def __eval_outer(self, x):
        d = np.linalg.norm(x - self.coordinates, axis=1)
        return generalized_exponential(self.__A, self.__B, d)

    def __integral_outer(self):
        return generalized_exponential_integral(self.__A, self.__B) * 4 * np.pi

    def __eval_ep_outer(self, x):
        d = np.linalg.norm(x - self.coordinates, axis=1)
        v = electrostatic_potential_exp(self.__A, self.__B, d)
        return v


class SupportAngular(Support):
    def __init__(self, **kwargs):
        """

        :param kwargs:
        """
        Support.__init__(self, **kwargs)
        try :
            self.__alpha = kwargs['alpha']
            self.__B = kwargs['B']
            self.__A = kwargs['A']
        except KeyError:
            raise IOError("missing parameters Alpha|G|U")
        try:  # the reason to handle the polynomial degree this way is to avoid
            self.__P = kwargs['P']
        except KeyError:
            self.__P = 1
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_trigo
        self.eval_ep_method = self.__eval_ep_ag
        self.params_kw = ['alpha','B', 'A']
        if 'P' in kwargs.keys():
            self.params_kw.append('P')
        self.support_type = 'angular'
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.coordinate_system is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_polar_rep(x, center=self.coordinates, ref_frame=self.coordinate_system)
            radial_component=generalized_exponential(self.__A, self.__B, d, self.__P)
            angular_component=gaussian(self.__alpha, 1.0, z)
            return radial_component * angular_component

    def __integral_trigo(self):
        return generalized_exponential_integral(
            self.__A, self.__B, self.__P
        ) * angular_gaussian_integral(self.__alpha) * 2 * np.pi

    def __eval_ep_ag(self, x):
        if self.__P != 1:
            raise NotImplementedError("still working in a generalized potential")
        if self.coordinate_system is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_polar_rep(x, center=self.coordinates, ref_frame=self.coordinate_system)
            u = electrostatic_potential_xexp_gaussian(self.__B, self.__alpha, d, z)
            return u * self.__A

class SupportHarmonic(Support):
    def __init__(self, **kwargs):
        """

        :param kwargs:
        """
        Support.__init__(self, **kwargs)
        try :
            self.__l = kwargs['l']
            self.__B = kwargs['B']
            self.__A = kwargs['A']
            self.__P = kwargs['P']
        except KeyError:
            raise IOError("missing parameters Alpha|G|U")
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_harmonic
        self.eval_ep_method = self.__eval_ep_harmonic
        self.params_kw = ['l','B', 'A', 'P']
        self.support_type = 'harmonic'
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.coordinate_system is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_polar_rep(x, center=self.coordinates, ref_frame=self.coordinate_system)
            radial_component=generalized_exponential(self.__A, self.__B, d, self.__P)
            angular_component = spherical_harmonic(z, self.__l)
            return radial_component * angular_component

    @staticmethod
    def __integral_harmonic():
        return 0.0

    def __eval_ep_harmonic(self, x):
        z, d = get_polar_rep(x, center=self.coordinates, ref_frame=self.coordinate_system)
        v = pe_harmonic(d, z, self.__l, self.__P, self.__B) * self.__A
        return v


class SupportEnsemble(A2MDBaseClass):

    def __init__(self, functions, name, map2atoms=None, map2bonds=None):
        A2MDBaseClass.__init__(self, name='ensamble of support functions', verbose=False)
        if isinstance(functions, list):
            self.fun = functions
        else:
            raise IOError("must provide at least an empty list, and fill it using append")
        self.name = name
        self.map2atoms = map2atoms
        self.map2bonds = map2bonds

    def append(self, fun):
        if issubclass(type(fun), Support):
            self.fun.append(fun)

    def eval(self, x):
        y = np.zeros(x.shape[0], dtype='float64')
        for fun in self.fun:
            y += fun.eval(x)
        return y

    def eval_ep(self, x):
        y = np.zeros(x.shape[0], dtype='float64')
        for fun in self.fun:
            y += fun.eval_ep(x)
        return y


    def integral(self):
        z = 0.0
        for fun in self.fun:
            z += fun.integral()
        return z