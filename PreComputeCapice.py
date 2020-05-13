#!/usr/bin/env python3

from src.pre_compute_capice import CalculateCapiceScores
from src.logger import Logger
from src.command_line_supporter import ArgumentSupporter


def main():
    """
    Main method of the script. Will call the various classes.
    """
    arguments = ArgumentSupporter()
    cadd_loc = arguments.get_argument('file')
    model_loc = arguments.get_argument('model')
    output_loc = arguments.get_argument('output')
    batch_size = arguments.get_argument('batchsize')
    if isinstance(cadd_loc, list):
        cadd_loc = str(cadd_loc[0])
    if isinstance(model_loc, list):
        model_loc = str(model_loc[0])
    if isinstance(output_loc, list):
        output_loc = str(output_loc[0])
    if isinstance(batch_size, list):
        batch_size = int(batch_size[0])
    logger = Logger()
    logger.set_output_dir(output_loc)
    logger.log('CADD file location: {}'.format(cadd_loc))
    logger.log('Model file location: {}'.format(model_loc))
    logger.log('Output directory: {}'.format(output_loc))
    logger.log('Batch size set to: {}'.format(batch_size))
    precompute_capice = CalculateCapiceScores(filepath=cadd_loc,
                                              model_loc=model_loc,
                                              output_loc=output_loc,
                                              batch_size=batch_size)
    precompute_capice.calc_capice()


if __name__ == '__main__':
    main()
