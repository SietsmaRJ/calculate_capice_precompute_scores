from PreComputeCapice import Utilities
import os
from datetime import datetime


class Logger:
    class __Logger:
        """
        Class to make a logfile on the progress being made.
        """
        def __init__(self):
            self.output_dir = None
            self.logfile = None
            self.utilities = Utilities()
            self.check_if_dir_exist()
            self.check_if_log_file_exist()

        def set_output_dir(self, output_loc):
            self.output_dir = output_loc

        def check_if_dir_exist(self):
            output_dir = os.path.join(self.output_dir, 'log_output')
            self.utilities.check_if_dir_exists(output_dir)
            self.output_dir = output_dir

        def check_if_log_file_exist(self):
            log_file_name = '{}_logfile.txt'.format(datetime.now().strftime(
                "%Y_%m_%d_%H%M%S_%f"))
            joined_path = os.path.join(self.output_dir, log_file_name)
            self.utilities.check_if_file_exists(joined_path)
            self.logfile = joined_path

        def get_output_dir(self):
            return self.output_dir

        def log(self, message):
            timestamp = datetime.now().strftime("%H:%M:%S_%f")
            timed_message = '[{}]: {}\n'.format(timestamp, message)
            with open(self.logfile, 'a') as logfile:
                logfile.write(timed_message)

    instance = None

    def __new__(cls):
        """
        Class method to set Logger instance
        :return: instance
        """
        if not Logger.instance:
            Logger.instance = Logger.__Logger()
        return Logger.instance

    def __init__(self):
        """
        __init__ method to set instance to Logger.__Logger()
        """
        if not Logger.instance:
            Logger.instance = Logger.__Logger()

    def __getattr__(self, name):
        """
        Method to return the value of the named attribute of name
        :param name: str
        :return: str
        """
        return getattr(self.instance, name)
