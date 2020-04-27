# Sirius B Data Reduction

PI: Michael Bottom, Institute for Astronomy

## Setup

To set up, I used [pipenv](https://github.com/pypa/pipenv) for python virtual environmenting. To loosely recreate my environment, simply

    $ pipenv install

To exactly replicate my environment from Pipfile.lock 

    $ pipenv install --deploy

## Usage

The processing and analysis code is all located in `notebooks`. The order is `preprocessing` > `adi` > `photometry`. Some helper library code is located in `src`. Of particular note, to make sure your paths are set up correctly, take a look at `src/paths.py`. 

