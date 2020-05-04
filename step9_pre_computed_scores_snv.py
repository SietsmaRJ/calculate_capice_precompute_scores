#!/usr/bin/env python3

import pickle
import pandas as pd
from impute_preprocess import impute, preprocess
import gzip
import time
import os
import argparse
from datetime import datetime
import psutil


class CalculateCapiceScores:
    """
    Main class of the script to call all the various logger class functions and
    will process the iterative chunking processing of the CADD file.
    """
    def __init__(self, logger_instance, filepath, model_loc, output_loc,
                 batch_size):
        self.log = logger_instance
        self.filepath = filepath
        self.titles = None
        self.get_header()
        self.model = None
        self.model_feats = None
        self.load_model(model_loc)
        self.not_done = True
        self.batch_size = batch_size
        self.output_loc = output_loc
        self.utilities = Utilities()
        self.utilities.check_if_dir_exists(output_loc)

    def get_header(self):
        if not self.titles:
            with gzip.open(self.filepath, "r") as f:
                while True:
                    line = f.readline()
                    if line.decode('utf-8').startswith("#Chr"):
                        self.titles = line.decode('utf-8').strip().split("\t")
                        self.log.log('Title found: {}'.format(self.titles))
                        break
                    else:
                        continue

    def calculate_save_capice_score(self, skip_rows):
        variants_df = pd.read_csv(self.filepath, sep='\t', skiprows=skip_rows,
                                  nrows=self.batch_size, names=self.titles,
                                  comment='#', compression='gzip',
                                  low_memory=False)
        if variants_df.shape[0] < self.batch_size:
            self.not_done = False
            self.log.log('Processing the last entries! '
                         'Total variants processed:'
                         ' {}.'.format(skip_rows + variants_df.shape[0]))

        variants_df_preprocessed = preprocess(impute(variants_df),
                                              model_features=self.model_feats)
        variants_df['prediction'] = self.model.predict_proba(
            variants_df_preprocessed[self.model_feats])[:, 1]
        if variants_df['prediction'].isnull().any():
            self.log.log('NaN encounter in chunk: {}+'
                         '{}!'.format(skip_rows,
                                      self.batch_size))
        if variants_df[variants_df.duplicated()].shape[0] > 0:
            duplicate = variants_df[variants_df.duplicated()]
            self.log.log('Duplicate encountered in CADD dataset!: \nIndex:{},'
                         '\nEntry:{}'.format(duplicate.index, duplicate))
        for unique_chr in variants_df['#Chr'].unique():
            subset_variants_df = variants_df[variants_df['#Chr'] == unique_chr]
            output_dir = os.path.join(self.output_loc, 'chr{}'.format(
                unique_chr))
            self.utilities.check_if_dir_exists(output_dir)
            chunk = 'chr_{}'.format(unique_chr)
            output_filename = 'whole_genome_SNVs_{}.txt'.format(chunk)
            final_destination = os.path.join(output_dir, output_filename)
            self.utilities.check_if_file_exists(final_destination)
            features_of_interest = ['#Chr', 'Pos', 'Ref', 'Alt',
                                    'GeneID', 'CCDS', 'FeatureID', 'prediction']
            with open(final_destination, 'a') as f:
                subset_variants_df[features_of_interest].to_csv(f, sep="\t",
                                                                index=False,
                                                                header=None)

    def load_model(self, model_loc):
        self.model = pickle.load(open(model_loc, "rb")).best_estimator_
        self.model_feats = self.model.get_booster().feature_names

    def calc_capice(self):
        start = None
        first_iter = True
        start_time = time.time()
        reset_timer = time.time()
        while self.not_done:
            time_iwl = time.time()
            if time_iwl - reset_timer > (60 * 60) or first_iter:
                # Seconds times the amount of minutes.
                curr_time = time.time()
                time_difference = curr_time - start_time
                minutes, seconds = divmod(time_difference, 60)
                hours, minutes = divmod(minutes, 60)
                self.log.log(
                    'Still going for {} hours,'
                    ' {} minutes and {} seconds.'.format(
                        int(round(hours)),
                        int(round(minutes)),
                        int(round(seconds))
                    )
                )
                self.log.log('Memory usage: {} MB.'.format(
                    self.log.get_ram_usage()))
                if not start:
                    self.log.log('Currently working on rows {} -'
                                 ' {}.'.format(start, start + self.batch_size))
                reset_timer = time.time()

            self.calculate_save_capice_score(start)
            if first_iter:
                start = 2
                first_iter = False
            start += self.batch_size


class ArgumentSupporter:
    """
    Class to handle the given command line input.
    Type python3 step9_pre_computed_scores_snv.py --help for more details.
    """

    def __init__(self):
        parser = self._create_argument_parser()
        self.arguments = parser.parse_args()

    @staticmethod
    def _create_argument_parser():
        parser = argparse.ArgumentParser(
            prog="step9_pre_computed_scores_snv.py",
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


class Logger:
    """
    Class to make a logfile on the progress being made.
    """
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.logfile = None
        self.utilities = Utilities()
        self.check_if_dir_exist()
        self.check_if_log_file_exist()

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

    @staticmethod
    def get_ram_usage():
        process = psutil.Process(os.getpid())
        memory_usage = process.memory_info().rss / 1000000  # Megabytes
        return memory_usage

    def log(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S_%f")
        timed_message = '[{}]: {}\n'.format(timestamp, message)
        with open(self.logfile, 'a') as logfile:
            logfile.write(timed_message)


class Utilities:
    @staticmethod
    def check_if_dir_exists(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def check_if_file_exists(filepath):
        if not os.path.isfile(filepath):
            new_file = open(filepath, 'a+')
            new_file.close()


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
    logger = Logger(output_loc)
    logger.log('CADD file location: {}'.format(cadd_loc))
    logger.log('Model file location: {}'.format(model_loc))
    logger.log('Output directory: {}'.format(output_loc))
    logger.log('Batch size set to: {}'.format(batch_size))
    precompute_capice = CalculateCapiceScores(logger_instance=logger,
                                              filepath=cadd_loc,
                                              model_loc=model_loc,
                                              output_loc=output_loc,
                                              batch_size=batch_size)
    precompute_capice.calc_capice()


if __name__ == '__main__':
    main()
