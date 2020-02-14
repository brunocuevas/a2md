import numpy as np
import time

class A2MD_basis :
    def __init__(self, name , verbose = False):
        """
        aAMD basis class
        :param name: name of the class
        :param verbose: activates the log calls
        :type verbose: bool
        :type name: str
        """
        self.v = verbose
        self.__name = name
        self.__fun  = []
        self.__ang  = []

    def __str__(self):
        return 'A2MD | ' + self.__name

    def add_ang(self, ang):
        """
        adds angular function
        :param ang: function dependent of angle
        :type ang: function
        :return: None
        """
        self.__ang.append(ang)

    def add_fun(self, fun):
        """
        adds radial function
        :param fun: function dependent of radius
        :return: None
        """
        self.__fun.append(fun)

    def eval(self, x):
        """
        evaluates radial and angular function product by using pairwise
        multiplication.
        :param x: cartesian coordinates for evaluation
        :type x: np.ndarray
        :return: density value
        :rtype: np.ndarray
        """
        y = np.zeros(x.shape, dtype='float64')
        z = np.zeros(x.shape, dtype='float64')
        for f in self.__fun :
            y += f(x)
        for g in self.__ang :
            z += g(x)
        return y*z

    def get_fun(self):
        """
        get list of radial functions
        :return: list of radial functions
        :rtype: list
        """
        return self.__fun

    def get_ang(self):
        """
        get list of angular functions
        :return: list of angular functions
        :rtype: list
        """
        return self.__ang

    def get_name(self):
        """
        provides name of the instance
        :return: name of the instance
        :rtype: str
        """
        return self.__name

    def log(self, message):
        """
        logs a message with format: [name] - message - hh:min:ss
        only when self.__v is True
        :param message: message to log
        :return: None
        """
        if self.v :
            t = time.localtime()
            formatted_t = '%02d:%02d:%02d' % (t.tm_hour, t.tm_min, t.tm_sec)
            print("[{0}] - {1} - {2} ".format(self.__name, message, formatted_t))
        return True

    def remove_fun(self):
        """
        removes all the list of radial functions
        :return:
        """
        del self.__fun

    def set_name(self, name):
        """
        sets a new object name
        :param name: new name
        :type name: str
        :return:
        """
        self.__name = name