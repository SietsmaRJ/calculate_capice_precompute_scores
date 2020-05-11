# PreComputeCapice
###### When you want the full capice experience, but don't want to use the model

### This is a program that works on the wonderfull [CAPICE](https://github.com/molgenis/capice) predictive model.

## Introduction

'PreComputeCapice' is a program developed by R.J. Sietsma and maintained by the Genomic Coordination Centre Groningen.
The program calculates capice scores for all entries in a given CADD file.

## Prerequisites

This program is developed in Pycharm 2020.1 (Professional Edition), performance on other systems is not guaranteed.

The program requires the following packages:
 
 * numpy ([v1.18.3](https://github.com/numpy/numpy/releases); [BSD 3-Clause License](https://www.numpy.org/license.html))
 * pandas ([v1.0.3](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))
 * psutil ([v5.7.0](https://github.com/giampaolo/psutil); [BSD 3-Clause License](https://github.com/giampaolo/psutil/blob/master/LICENSE))
 * scipy ([v1.4.1](https://github.com/scipy/scipy); [BSD 3-Clause License](https://github.com/giampaolo/psutil/blob/master/LICENSE))
 * scikit-learn ([v0.19.1](https://scikit-learn.org/stable/whats_new.html); [BSD 3-Clause License](https://github.com/scikit-learn/scikit-learn/blob/master/COPYING))
 * xgboost ([v0.72.1](https://github.com/dmlc/xgboost); [Apache 2 License](https://github.com/dmlc/xgboost/blob/master/LICENSE))
 
__Warning: this program works for python version 3.6, it does not work for python 3.7 or higher__

## Installing

**Step 1: acquire the source files**
Either [clone]() or [download]() the source files.

**Step 2: activate the virtual environment**
- Open a terminal in the cloned or downloaded folder.
- Make sure you have python3.6 installed by typing python3 --version.
-  Execute the following:
```console
mkdir ./venv
cd ./venv
python3 -m venv ./
cd ..
```
-  Activate the virtual environment:
```console
source ./venv/bin/activate
```
- Install the required packages by executing:
```console
pip install -r requirements.txt
```

_Note: if any package fails to install, please try to install the package using:_
```console
pip install package==version
```
__Make sure you have the virtual environment enabled before you install packages!__

## Usage

The program requires the following arguments:

- -f / --file: the cadd file in gzip (.gz) format.
- -m / --model: the pickled capice model in .dat format.
- -o / --output: the location where the program should place it's files.

Optional argument:

- -s / --batchsize: the amount of rows the program should read each iteration from the CADD file.

Example usage:

_with a batch size of 1 million_
```console
python3 PreComputeCapice.py -f path/to/cadd/file.gz -m path/to/model.dat -o path/to/output/folder -s 1000000
```

## Output

The program will output the following files:
- For each chromosome in the CADD file, it makes a folder named chrx (where x = chromosome) and places a gzipped tsv of all CADD entries for that chromosome.
__Note: The program continually adds entries to this file, do NOT remove or replace this file till the program is done!__
- Log_output: a file with timed messages on updates within the program. (Does not contain error messages or warnings).

## TODO:
- Add function to check if existing files are present and continue from where was left of (in case of crash).
    - Plan: export json with start and stop (start + batch size) each iteration to the output folder.