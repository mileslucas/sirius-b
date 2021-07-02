# Sirius B Data Reduction

PI: Michael Bottom, Institute for Astronomy


## Setup

To set up, I used [virtualenv](https://github.com/pypa/pipenv) for the python requirements.

    $ pip install -r requirements.txt

Then, to set up the Julia requirements

    $ julia --project=@. -e 'using Pkg; Pkg.instantiate()'

## Usage

The processing and analysis code is all located in `notebooks`. The order is `preprocessing` > `psf_subtraction` > `analysis`. Some helper library code is located in `src`.

