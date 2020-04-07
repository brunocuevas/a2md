import time

class A2MDBaseClass :
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

    def set_name(self, name):
        """
        sets a new object name
        :param name: new name
        :type name: str
        :return:
        """
        self.__name = name