import psutil
import os
import gzip


class Utilities:
    @staticmethod
    def get_ram_usage():
        process = psutil.Process(os.getpid())
        memory_usage = process.memory_info().rss / 1000000  # Megabytes
        return memory_usage

    @staticmethod
    def check_if_dir_exists(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def check_if_file_exists(filepath, capice_ouput=False, return_value=False):
        if return_value:
            make_output = False
        else:
            make_output = True
        if not os.path.isfile(filepath):
            if return_value:
                return False
            if capice_ouput:
                if make_output:
                    with gzip.open(filepath, 'a+') as file:
                        file.close()
            else:
                if make_output:
                    with open(filepath, 'a+') as file:
                        file.close()
        elif return_value:
            return True
