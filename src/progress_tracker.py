from src.logger import Logger
from src.utilities.utilities import Utilities
import json
from pathlib import Path
import gzip
import os


class ProgressTracker:
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
            self.progress_json = {'start': None, 'batch_size': None}
            with open(progress_json, 'w+') as json_file:
                json.dump(self.progress_json, json_file)

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
            total_lines = -10
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

    def is_in_progression_json(self, key):
        return_value = False
        if key in self.progress_json.keys():
            return_value = True
        return return_value

    def get_start_and_batchsize(self):
        return self.start, self.batch_size

    def get_progression_loc(self):
        return self.progress_json_loc

    def get_progression_json_value(self, key):
        return self.progress_json[key]

    def update_progression(self, key, value):
        if key in self.progress_json.keys():
            to_be_updated = self.progress_json[key]
            if to_be_updated != value:
                self.progress_json[key] = value
                self.log.log('Updating progression on: {}, value: {}'.format(
                    key, value
                ))
        else:
            self.progress_json[key] = value
            self.log.log('No progression found on {}, adding'
                         ' to progression json.'.format(key))
        with open(self.progress_json_loc, 'w') as json_file:
            json.dump(self.progress_json, json_file)
