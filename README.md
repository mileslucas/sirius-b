# Sirius B Data Reduction

PI: Michael Bottom, Institute for Astronomy

## Setup

After cloning this repository, we can set up the python dependencies. I used [poetry](https://python-poetry.org) for the python requirements.

    $ poetry install

I had some issues getting scipy to play nicely, you'll have to make it work yourself. For me I had to set 

    $ poetry config experimental.new-installer false 

Pro tip: if you want a local `.venv` folder set 

    $ poetry config virtualenvs.in-project true

Then, to set up the Julia requirements

    # activate python environment first
    $ poetry shell
    # add environment variable so PyCall.jl links to virtual env
    $ PYTHON=$(which python) \
      julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'

(the final `Pkg.build()` is not necessarily required, but will re-link *PyCall.jl* in case something bad happened)

## Usage

The processing and analysis code is all located in `notebooks`. The order is `preprocessing` > `psf_subtraction` > `analysis`. Some helper library code is located in `src`.
