from a2md.baseclass import A2MD_basis
import numpy as np
from a2md.mathfunctions import expfun
from a2md.mathfunctions import xexpfun
from a2md.mathfunctions import x2expfun
from a2md.mathfunctions import expfun_integral
from a2md.mathfunctions import xexpfun_integral_trigo
from a2md.mathfunctions import expfun_integral_trigo
from a2md.mathfunctions import x2expfun_integral_trigo
from a2md.mathfunctions import get_angle
from a2md.mathfunctions import xexpfun_integral_vk
from a2md.mathfunctions import cosssin
from a2md.mathfunctions import gaussian, angular_gaussian_integral, dipole_gaussian
from a2md.mathfunctions import electrostatic_potential_exp, electrostatic_potential_xexp_gaussian

from a2md.mathfunctions import ANGULAR_PARAMS,ANGULAR_INTEGRAL_FUNCTIONS
from a2md.mathfunctions import RADIAL_PARAMS, RADIAL_INTEGRAL_FUNCTIONS

class support(A2MD_basis) :
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
        A2MD_basis.__init__(self, name ='support function ', verbose=False)
        self.support_coordinates = kwargs['coordinates']
        self.reference_frame = None
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
        B = np.zeros((3,3), dtype='float64')
        new_x = np.cross(w,np.array([0, 0, 1], dtype='float64'))
        new_x /= np.linalg.norm(new_x)
        new_y = np.cross(w, new_x)
        B[:, 0] = new_x
        B[:, 1] = new_y
        B[:, 2] = w
        try:
            self.reference_frame = np.linalg.inv(B)
        except np.linalg.LinAlgError:
            new_x = np.cross(w, np.array([1, 0, 0], dtype='float64'))
            new_y = np.cross(w, new_x)
            B[:, 0] = new_x
            B[:, 1] = new_y
            B[:, 2] = w
            self.reference_frame = np.linalg.inv(B)
        return True

class support_inner_radial(support):
    """
    support - inner radial function
        \ro(r) = A1 * e^-(B1*r) + A2 * e^-(B2*r)
    """
    def __init__(self, **kwargs):
        r"""
        Support functions. It requires:
        - coordinates : where function is centered
        - radial function type : 0 - sum of exp, 1 - exp, 2 - xexp
        - angular function type : 0 - isotropic, 1 - expansion of size 1, 2- expansion of size 2...
        - A1 : coefficient 1
        - A2 : coefficient 2
        - B1 : exponent 1
        - B2 : exponent 2
        :param \**kwargs['coordinates']:
            See below
        :keyword Arguments:
        * *coordinates* *(``np.ndarray``)
        * *radial_function_type* *(``int``)
        * *angular_function_type* *(``int``)
        * *A1* *(``float``)
        * *B1* *(``float``)
        * *A2* *(``float``)
        * *B2* *(``float``)

        """
        support.__init__(self, **kwargs)
        self.set_name(
            'inner function at %04.2f:%04.2f:%04.2f' % (
                kwargs['coordinates'][0],kwargs['coordinates'][1],kwargs['coordinates'][2]
            )
        )
        try :
            self.__A1 = kwargs['A1']
            self.__A2 = kwargs['A2']
            self.__B1 = kwargs['B1']
            self.__B2 = kwargs['B2']
        except KeyError:
            raise IOError("missing parameters A1|A2|B1|B2")
        self.eval_method = self.__eval_inner
        self.integral_method = self.__integral_inner
        self.eval_ep_method = self.__eval_ep_inner
        self.support_type = 'IR' # IR stands for inner radial
        self.params_kw = ['A1', 'A2', 'B1', 'B2']
        self.anisotropic = False

    def __eval_inner(self, x):
        d = np.linalg.norm(x - self.support_coordinates, axis=1)
        return expfun(self.__A1, self.__B1, d) + expfun(self.__A2, self.__B2, d)

    def __integral_inner(self):
        return 2*(expfun_integral(self.__A1, self.__B1) + expfun_integral(self.__A2, self.__B2))

    def __eval_ep_inner(self, x):
        d = np.linalg.norm(x - self.support_coordinates, axis=1)
        return electrostatic_potential_exp(
            self.__A1, self.__B1, d
        ) + electrostatic_potential_exp(self.__A2, self.__B2, d)

class support_outer_radial(support):
    """
    support - outer radial function
        \ro(r) = A3 * e^-(B3*r)
    """
    def __init__(self, **kwargs):
        r"""
        Support functions. It requires:
        - coordinates : where function is centered
        - radial function type : 0 - sum of exp, 1 - exp, 2 - xexp
        - angular function type : 0 - isotropic, 1 - expansion of size 1, 2- expansion of size 2...
        - A3 : coefficient 3
        - B3 : exponent 3
        :param \**kwargs['coordinates']:
            See below
        :keyword Arguments:
        * *coordinates* *(``np.ndarray``)
        * *radial_function_type* *(``int``)
        * *angular_function_type* *(``int``)
        * *A3* *(``float``)
        * *B3* *(``float``)

        """
        support.__init__(self, **kwargs)
        self.set_name(
            'outer function at %04.2f:%04.2f:%04.2f' % (
                kwargs['coordinates'][0], kwargs['coordinates'][1], kwargs['coordinates'][2]
            )
        )
        try:
            self.__A3 = kwargs['A3']
            self.__B3 = kwargs['B3']
        except KeyError:
            raise IOError("missing parameters A3|B3")
        self.eval_method = self.__eval_outer
        self.integral_method = self.__integral_outer
        self.eval_ep_method = self.__eval_ep_outer
        self.support_type = 'OR' # OR stands for outter radial
        self.params_kw = ['A3', 'B3']
        self.anisotropic = False

    def __eval_outer(self, x):
        d = np.linalg.norm(x - self.support_coordinates, axis=1)
        return expfun(self.__A3, self.__B3, d)

    def __integral_outer(self):
        return 2*(expfun_integral(self.__A3, self.__B3))

    def __eval_ep_outer(self, x):
        d = np.linalg.norm(x - self.support_coordinates, axis=1)
        return electrostatic_potential_exp(self.__A3, self.__B3, d)

class support_angular(support):
    """
    support - angular function
        \ro(r,z) = F*r*e^-(G*r) * \sum u_i*(1+cos(i*z) + v_i*(1+sin(i*z))
    """
    def __init__(self, **kwargs):
        r"""
        Support functions. It requires:
        - coordinates : where function is centered
        - radial function type : 0 - sum of exp, 1 - exp, 2 - xexp
        - angular function type : 0 - isotropic, 1 - expansion of size 1, 2- expansion of size 2...
        - Psi: expansion size
        - F: coefficient
        - G: exponent
        - U: cosine coefficients
        - V: sine coefficients
        :param \**kwargs['coordinates']:
            See below
        :keyword Arguments:
        * *coordinates* *(``np.ndarray``)
        * *radial_function_type* *(``int``)
        * *angular_function_type* *(``int``)
        * *F* *(``float``)
        * *G* *(``float``)
        * *Psi* *(``int``)
        * *U* *(``np.ndarray``)
        * *V* *(``np.ndarray``)

        """
        support.__init__(self, **kwargs)
        self.set_name(
            'angular function at %04.2f:%04.2f:%04.2f' % (
                kwargs['coordinates'][0], kwargs['coordinates'][1], kwargs['coordinates'][2]
            )
        )
        try :
            self.__Psi = kwargs['Psi']
            self.__F = kwargs['F']
            self.__G = kwargs['G']
            self.__U = np.array(kwargs['U'], dtype='float64')
            self.__V = np.array(kwargs['V'], dtype='float64')
        except KeyError:
            raise IOError("missing parameters Psi|F|G|U|V")
        self.eval_method = self.__eval_angular
        self.integral_method = self.__integral_angular
        self.support_type = 'ACS' # Angular cosine + sine
        self.params_kw = ['Psi', 'F', 'G', 'U', 'V',]
        self.anisotropic = True

    def __eval_angular(self, x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            angular_component = np.zeros(z.size, dtype='float64')
            radial_component = xexpfun(self.__F, self.__G, d)
            for p in range(self.__Psi):
                angular_component += cosssin(self.__U[p], self.__V[p], p+1, z)
            return angular_component * radial_component

    def __integral_angular(self):
        return xexpfun_integral_vk(self.__F, self.__G, self.__U, self.__V, self.__Psi)

class support_fund_trigonometric(support):
    def __init__(self, trigo = None, **kwargs):
        support.__init__(self, **kwargs)
        if trigo == np.sin :
            self.set_name(" sine function")
        elif trigo == np.cos:
            self.set_name(" cosine function")
        else:
            raise IOError("not implemented")
        self.__trigo = trigo
        try :
            self.__Psi = kwargs['Psi']
            self.__G = kwargs['G']
            self.__U = kwargs['U'] # U can be either U or V depending on whether it is on a sin or cos function
        except KeyError:
            raise IOError("missing parameters Psi|F|G|U|V")
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_trigo
        if trigo == np.sin :
            self.support_type = 'AS' # Angular sine
        elif trigo == np.cos :
            self.support_type = 'AC' # Angular cosine
        self.params_kw = ['Psi','G', 'U',]
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            radial_component=xexpfun(self.__U, self.__G, d)
            angular_component=self.__trigo(z* (self.__Psi +1))
            return radial_component * angular_component
    def __integral_trigo(self):
        return xexpfun_integral_trigo(self.__trigo, self.__G, self.__U, self.__Psi)
    def first_moment(self):
        moment = np.array([0.0, 0.0, 0.0], dtype='float64')
        if self.__trigo == np.sin :
            if self.__Psi == 0:
                moment[0] = 12 * (np.pi**2) / (self.__G**5)
        elif self.__trigo == np.cos:
            if self.__Psi == 0:
                moment[0] = 32 * np.pi / (self.__G**5)
        return moment

class support_fund_trigonometric_exp(support):
    def __init__(self, trigo = None, **kwargs):
        support.__init__(self, **kwargs)
        if trigo == np.sin :
            self.set_name(" sine function")
        elif trigo == np.cos:
            self.set_name(" cosine function")
        else:
            raise IOError("not implemented")
        self.__trigo = trigo
        try :
            self.__Psi = kwargs['Psi']
            self.__G = kwargs['G']
            self.__U = kwargs['U'] # U can be either U or V depending on whether it is on a sin or cos function
        except KeyError:
            raise IOError("missing parameters Psi|F|G|U|V")
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_trigo
        if trigo == np.sin :
            self.support_type = 'AES' # Angular exponential sine
        elif trigo == np.cos :
            self.support_type = 'AEC' # Angular exponential cosine
        self.params_kw = ['Psi','G', 'U',]
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            radial_component=expfun(self.__U, self.__G, d)
            angular_component=self.__trigo(z* (self.__Psi +1))
            return radial_component * angular_component
    def __integral_trigo(self):
        return expfun_integral_trigo(self.__trigo, self.__G, self.__U, self.__Psi)


class support_fund_trigonometric_x2exp(support):
    def __init__(self, trigo = None, **kwargs):
        support.__init__(self, **kwargs)
        if trigo == np.sin :
            self.set_name(" sine function")
        elif trigo == np.cos:
            self.set_name(" cosine function")
        else:
            raise IOError("not implemented")
        self.__trigo = trigo
        try :
            self.__Psi = kwargs['Psi']
            self.__G = kwargs['G']
            self.__U = kwargs['U'] # U can be either U or V depending on whether it is on a sin or cos function
        except KeyError:
            raise IOError("missing parameters Psi|F|G|U|V")
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_trigo
        if trigo == np.sin :
            self.support_type = 'A2S' # Angular exponential sine
        elif trigo == np.cos :
            self.support_type = 'A2C' # Angular exponential cosine
        self.params_kw = ['Psi','G', 'U',]
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            radial_component=x2expfun(self.__U, self.__G, d)
            angular_component=self.__trigo(z* (self.__Psi +1))
            return radial_component * angular_component
    def __integral_trigo(self):
        return x2expfun_integral_trigo(self.__trigo, self.__G, self.__U, self.__Psi)


class support_agnostic(support):
    def __init__(self, name, angular = None, radial = None, **kwargs):
        support.__init__(self, **kwargs)
        self.support_type = name
        try :
            general_radial_kw = RADIAL_PARAMS[radial]
        except KeyError:
            raise NotImplementedError("there is not a defined list of parameters for that function")
        try :
            general_angular_kw = ANGULAR_PARAMS[angular]
        except KeyError:
            raise NotImplementedError("there is not a defined list of parameters for that function")

        private_radial_params = dict()
        private_angular_params = dict()

        if len(general_radial_kw) > 0:
            for key in general_radial_kw :
                private_radial_params[key] = kwargs[key]
        if len(general_angular_kw) > 0:
            for key in general_angular_kw :
                private_angular_params[key] = kwargs[key]

        self.radial_kw = private_radial_params
        self.angular_kw = private_angular_params
        self.eval_method = self.agnostic_eval
        self.integral_method = self.agnostic_integral
        self.angular = lambda z,pars : angular(z = z, **pars)
        self.radial = lambda x,pars : radial(d = x, **pars)
        self.params_kw = RADIAL_PARAMS[radial] + ANGULAR_PARAMS[angular]

        try:
            self.integral_radial = RADIAL_INTEGRAL_FUNCTIONS[radial]
        except KeyError:
            raise NotImplementedError("there is not a defined integral function for that function")
        try:
            self.angular_radial = ANGULAR_INTEGRAL_FUNCTIONS[angular]
        except KeyError:
            raise NotImplementedError("there is not a defined integral function for that function")
        if angular is not None:
            self.anisotropic = True
        else:
            self.anisotropic = False

    def get_params(self):
        if self.params_kw is None:
            raise NotImplementedError
        else:
            params = dict()
            for key, item in self.radial_kw.items():
                params[key] = item
            for key, item in self.angular_kw.items():
                params[key] = item
            return params

    def agnostic_eval(self, x):
        z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
        radial_component  = self.radial(d, self.radial_kw)
        angular_component = self.angular(z,  self.angular_kw)
        return radial_component * angular_component
    def agnostic_integral(self):
        return self.integral_radial(**self.radial_kw) * self.angular_radial(**self.angular_kw)


class SupportEnsemble(A2MD_basis):

    def __init__(self, functions, name):
        A2MD_basis.__init__(self, name='ensamble of support functions', verbose=False)
        if isinstance(functions, list):
            self.fun = functions
        else:
            raise IOError("must provide at least an empty list, and fill it using append")
        self.name = name

    def append(self, fun):
        if issubclass(type(fun), support):
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

class support_angular_gaussian(support):
    def __init__(self, **kwargs):
        support.__init__(self, **kwargs)
        try :
            self.__Alpha = kwargs['Alpha']
            self.__G = kwargs['G']
            self.__U = kwargs['U']
            self.__Psi = kwargs['Psi']
        except KeyError:
            raise IOError("missing parameters Alpha|G|U")
        self.eval_method = self.__eval_trigo
        self.integral_method = self.__integral_trigo
        self.eval_ep_method = self.__eval_ep_ag
        self.params_kw = ['Alpha','G', 'U', 'Psi']
        self.support_type = 'AG'
        self.anisotropic = True

    def __eval_trigo(self,x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            radial_component=xexpfun(self.__U, self.__G, d)
            angular_component=gaussian(self.__Alpha, 1.0, z)
            return radial_component * angular_component

    def __integral_trigo(self):
        return angular_gaussian_integral(self.__G, self.__U, self.__Alpha)

    def first_moment(self):
        return dipole_gaussian(self.__G, self.__U, self.__Alpha)

    def __eval_ep_ag(self, x):
        if self.reference_frame is None:
            raise RuntimeError("reference frame was not created. can not calculate angular function")
        else:
            z,d = get_angle(x, center=self.support_coordinates, ref_frame=self.reference_frame)
            # radial_component=xexpfun(self.__U, self.__G, d)
            # angular_component=gaussian(self.__Alpha, 1.0, z)
            u = electrostatic_potential_xexp_gaussian(self.__G, self.__Alpha, d, z)
            return u * self.__U