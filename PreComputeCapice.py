#!/usr/bin/env python3

import pickle
import pandas as pd
from utilities.impute_preprocess import impute, preprocess
import gzip
import time
import os
import json
from pathlib import Path
from logger import Logger
from command_line_supporter import ArgumentSupporter
from utilities.utilities import Utilities


class CalculateCapiceScores:
    """
    Main class of the script to call all the various logger class functions and
    will process the iterative chunking processing of the CADD file.
    """
    def __init__(self, filepath, model_loc, output_loc,
                 batch_size):
        self.log = Logger()
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
        self.reinstance = OutputReInitializer(self.output_loc)
        self.previous_iteration_df = None

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

    def calculate_save_capice_score(self, skip_rows, batch_size):
        variants_df = pd.read_csv(self.filepath, sep='\t', skiprows=skip_rows,
                                  nrows=batch_size, names=self.titles,
                                  comment='#', compression='gzip',
                                  low_memory=False)
        if variants_df.shape[0] < batch_size:
            self.not_done = False
            self.log.log('Processing the last entries! '
                         'Total variants processed:'
                         ' {}.'.format(skip_rows + variants_df.shape[0] - 1))

        variants_df_preprocessed = preprocess(impute(variants_df),
                                              model_features=self.model_feats)
        variants_df['prediction'] = self.model.predict_proba(
            variants_df_preprocessed[self.model_feats])[:, 1]
        if variants_df['prediction'].isnull().any():
            self.log.log('NaN encounter in chunk: {}+'
                         '{}!'.format(skip_rows,
                                      batch_size))
        if variants_df[variants_df.duplicated()].shape[0] > 0:
            duplicate = variants_df[variants_df.duplicated()]
            self.log.log('Duplicate encountered in CADD dataset!: \nIndex:{},'
                         '\nEntry:{}'.format(duplicate.index, duplicate))
        for unique_chr in variants_df['#Chr'].unique():
            # //TODO: add a merge with previous dataframe and merge left only (
            # so subset_variants_df is the left one), by using how='left'
            # and indicator = True. Then select only those
            # with _merge == 'left_only' and boom no duplicates.

            # //TODO: add function to get a dataframe of the gzipped archive
            # of that chromosome, last start - x entries.
            subset_variants_df = variants_df[variants_df['#Chr'] == unique_chr]
            output_dir = os.path.join(self.output_loc, 'chr{}'.format(
                unique_chr))
            self.utilities.check_if_dir_exists(output_dir)
            chunk = 'chr_{}'.format(unique_chr)
            output_filename = 'whole_genome_SNVs_{}.tsv.gz'.format(chunk)
            final_destination = os.path.join(output_dir, output_filename)
            self.utilities.check_if_file_exists(final_destination,
                                                capice_ouput=True)
            features_of_interest = ['#Chr', 'Pos', 'Ref', 'Alt',
                                    'GeneID', 'CCDS', 'FeatureID', 'prediction']
            with gzip.open(final_destination, 'at') as f:
                subset_variants_df[features_of_interest].to_csv(f, sep="\t",
                                                                index=False,
                                                                header=None)

    def load_model(self, model_loc):
        self.model = pickle.load(open(model_loc, "rb")).best_estimator_
        self.model_feats = self.model.get_booster().feature_names

    def calc_capice(self):
        start, batch_size = self.reinstance.get_start_and_batchsize()
        if not batch_size:
            batch_size = self.batch_size
        self._calc_capice(start, batch_size)

    def _calc_capice(self, start, batch_size):
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
                    self.utilities.get_ram_usage()))
                if start:
                    self.log.log('Currently working on rows {} -'
                                 ' {}.'.format(start, start + batch_size))
                reset_timer = time.time()

            self.calculate_save_capice_score(start, batch_size)
            if first_iter:
                start = 2
                first_iter = False
            start += batch_size


class OutputReInitializer:
    """
    Class to check for existing files in terms of progress.
    """
    def __init__(self, output_loc):
        self.output = output_loc
        self.log = Logger()
        self.utilities = Utilities()
        self.progress_json = None
        self.start = None
        self.batch_size = None
        self.progress_json_loc = None
        self._check_for_progress_json()
        self._check_for_processed_files()

    def _check_for_progress_json(self):
        progress_json = os.path.join(self.output, 'log_output',
                                     'progression.json')
        self.progress_json_loc = progress_json
        is_json = self.utilities.check_if_file_exists(progress_json,
                                                      return_value=True)
        if is_json:
            with open(progress_json) as json_file:
                json_data = json.load(json_file)
            self.log.log('Progression json found! Continuing from {}'.format(
                json_data['start']))
            self.start = json_data['start']
            self.batch_size = json_data['batch_size']
        else:
            self.log.log('No progression json found,'
                         ' checking for processed files.')
            with open(progress_json, 'w+') as json_file:
                json.dump({'start': None, 'batch_size': None}, json_file)
        with open(progress_json) as p_json:
            self.progress_json = json.load(p_json)
            p_json.close()

    def _check_for_processed_files(self):
        need_to_process = []
        processed_file_nlines = self.progress_json
        for path in Path(self.output).rglob('whole_genome_SNVs_*.tsv.gz'):
            path = str(path)
            if path not in processed_file_nlines.keys():
                self.log.log('Found progress file!: {}'.format(path))
                processed_file_nlines[path] = 0
                need_to_process.append(path)
        if len(need_to_process) > 0:
            self.log.log('Attempting to track back'
                         ' the amount of processed lines.')
            for file in need_to_process:
                if os.path.isfile(file):
                    for line in gzip.open(file):
                        processed_file_nlines[file] += 1
            total_lines = 0
            for nline in processed_file_nlines.values():
                total_lines += nline
            self.log.log('Amount of lines found: {}'.format(total_lines))
            self.start = total_lines
            start_batchsize_json = {'start': self.start, 'batch_size': None}
            start_batchsize_json.update(processed_file_nlines)
            with open(self.progress_json, 'w') as progress_json:
                json.dump(start_batchsize_json,
                          progress_json)
            self.log.log('Start has been saved in: {}'.format(
                self.progress_json))
        with open(self.progress_json) as p_json:
            self.progress_json = json.load(p_json)
            p_json.close()

    def get_start_and_batchsize(self):
        return self.start, self.batch_size

    def get_progression_loc(self):
        return self.progress_json_loc

    def get_progression_json_value(self, key):
        return self.progress_json[key]


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
