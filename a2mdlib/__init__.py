import os
import time

current_pwd = os.getcwd()
path2params = current_pwd + '/inp/amd_params_18'


class A2MDlib:
    def __init__(
            self,
            name='aAMDlib',
            verbose=False,
    ):
        """

        :param name:
        :param verbose:
        """
        self.__name = name
        self.__v = verbose

    def log(self, message):
        """

        :param message:
        :return:
        """
        if self.__v:
            ltm = time.localtime()
            mssg_head = "[{:02d}:{:02d}:{:02d}] ".format(ltm.tm_hour, ltm.tm_min, ltm.tm_sec)
            print(mssg_head + message)