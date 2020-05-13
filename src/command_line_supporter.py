import argparse


class ArgumentSupporter:
    """
    Class to handle the given command line input.
    Type python3 PreComputeCapice.py --help for more details.
    """

    def __init__(self):
        parser = self._create_argument_parser()
        self.arguments = parser.parse_args()

    @staticmethod
    def _create_argument_parser():
        parser = argparse.ArgumentParser(
            prog="PreComputeCapice.py",
            description="Python script to calculate Pre-computed"
                        " scores for CAPICE for every given CADD variant.")
        required = parser.add_argument_group("Required arguments")
        optional = parser.add_argument_group("Optional arguments")

        required.add_argument('-f',
                              '--file',
                              nargs=1,
                              type=str,
                              required=True,
                              help='The location of the CADD'
                                   ' annotated SNV file.')

        required.add_argument('-m',
                              '--model',
                              nargs=1,
                              type=str,
                              required=True,
                              help='The location of the CAPICE'
                                   ' model pickled file.')

        required.add_argument('-o',
                              '--output',
                              nargs=1,
                              type=str,
                              required=True,
                              help='The output directory to put the processed'
                                   'CADD variants in.')

        optional.add_argument('-s',
                              '--batchsize',
                              nargs=1,
                              type=int,
                              default=10000,
                              required=False,
                              help='The chunksize for the script to'
                                   ' read the gzipped archive.'
                                   ' (Default: 10000)')
        return parser

    def get_argument(self, argument_key):
        """
        Method to get a command line argument.
        :param argument_key: Command line argument.
        :return: List or string.
        """
        if self.arguments is not None and argument_key in self.arguments:
            value = getattr(self.arguments, argument_key)
        else:
            value = None

        return value
