import pickle
import pandas as pd
from src.utilities.impute_preprocess import impute, preprocess
import gzip
import time
import os
from src.logger import Logger
from src.utilities.utilities import Utilities
from src.progress_tracker import ProgressTracker


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
        self.progress_track = ProgressTracker(self.output_loc)
        self.features_of_interest = ['#Chr', 'Pos', 'Ref', 'Alt',
                                     'GeneID', 'CCDS', 'FeatureID',
                                     'prediction']
        self.previous_iteration_df = pd.DataFrame(
            columns=self.features_of_interest)

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
        variants_df = variants_df[self.features_of_interest]
        if variants_df['prediction'].isnull().any():
            self.log.log('NaN encounter in chunk: {}+'
                         '{}!'.format(skip_rows,
                                      batch_size))
        if variants_df[variants_df.duplicated()].shape[0] > 0:
            duplicate = variants_df[variants_df.duplicated()]
            self.log.log('Duplicate encountered in CADD dataset!: \nIndex:{},'
                         '\nEntry:{}'.format(duplicate.index, duplicate))
        for iteration, unique_chr in enumerate(variants_df['#Chr'].unique()):
            subset_variants_df = variants_df[variants_df['#Chr'] == unique_chr]
            output_dir = os.path.join(self.output_loc, 'chr{}'.format(
                unique_chr))
            self.utilities.check_if_dir_exists(output_dir)
            output_filename = 'whole_genome_SNVs_chr_{}.tsv.gz'.format(
                unique_chr)
            final_destination = os.path.join(output_dir, output_filename)
            self.utilities.check_if_file_exists(final_destination,
                                                capice_ouput=True)
            subset_variants_df = self._get_and_check_last_entries(
                final_destination, subset_variants_df)
            with gzip.open(final_destination, 'at') as f:
                subset_variants_df.to_csv(f, sep="\t",
                                          index=False,
                                          header=None)
            total_processed_rows = 0
            if self.progress_track.is_in_progression_json(final_destination):
                total_processed_rows = \
                    self.progress_track.get_progression_json_value(
                        final_destination)
            total_processed_rows += subset_variants_df.shape[0]
            self.progress_track.update_progression(final_destination,
                                                   total_processed_rows)
        self.previous_iteration_df = variants_df.tail(100)

    def load_model(self, model_loc):
        self.model = pickle.load(open(model_loc, "rb")).best_estimator_
        self.model_feats = self.model.get_booster().feature_names

    def _merge_and_remove_dupes(self, subset_df):
        nrows_before = subset_df.shape[0]
        self.previous_iteration_df['iter'] = 'p'
        subset_df['iter'] = 'c'

        merge = self.previous_iteration_df.append(subset_df,
                                                  ignore_index=True)

        merge.drop_duplicates(subset=self.features_of_interest[:-1],
                              inplace=True,
                              ignore_index=True)

        merge = merge[merge['iter'] == 'c']
        merge.drop('iter', axis=1, inplace=True)
        nrows_after = merge.shape[0]
        self.log.log('Removed {} duplicated rows.'.format(
            nrows_before - nrows_after))
        return merge

    def _get_and_check_last_entries(self, output_filename, subset_df):
        if self.progress_track.is_in_progression_json(output_filename):
            self.log.log('No progression dataframe found, loading in'
                         ' from file: {}.'.format(output_filename))
            lines_processed = \
                self.progress_track.get_progression_json_value(
                    output_filename
                )
            if lines_processed < 100:
                start = None
                get_nrows = lines_processed
            else:
                start = lines_processed - 99
                get_nrows = 100
            self.previous_iteration_df = pd.read_csv(
                output_filename,
                compression='gzip',
                sep='\t',
                names=self.features_of_interest,
                nrows=get_nrows,
                skiprows=start)
        else:
            self.log.log(
                'No progression dataframe found nor progression file.'
                ' Skipping duplicate check for file: {}.'.format(
                    output_filename
                )
            )
        return_df = self._merge_and_remove_dupes(subset_df)
        return return_df

    def calc_capice(self):
        start, batch_size = self.progress_track.get_start_and_batchsize()
        if batch_size is None:
            batch_size = self.batch_size
            self.progress_track.update_progression('batch_size', batch_size)
        elif batch_size != self.batch_size:
            batch_size = self.batch_size
            self.progress_track.update_progression('batch_size', batch_size)
        self._calc_capice(start, batch_size)

    def _calc_capice(self, start, batch_size):
        start_time = time.time()
        reset_timer = time.time()
        while self.not_done:
            time_iwl = time.time()
            if time_iwl - reset_timer > (60 * 60):
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
            if not start:
                start = 0
            start += batch_size - 1
            self.progress_track.update_progression('start', start)
