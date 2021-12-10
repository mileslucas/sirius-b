# Sirius B Data Reduction

**PI: Michael Bottom, Institute for Astronomy**

This work has been accepted for publication in the Astronomical Journal as of December 3, 2021

## Attribution


## Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5115225.svg)](https://doi.org/10.5281/zenodo.5115225)
[![License](https://img.shields.io/badge/license-CC--BY--4-orange.svg)](https://creativecommons.org/licenses/by/4.0/)

The data from these observations is available publicly under the CC-BY-4 open-source license. Provided are the pre-processed data cubes and parallactic angles from each epoch of observation. To download the data, visit the [zenodo page](https://doi.org/10.5281/zenodo.5115225) and please cite the appropriate DOI if you use it in your work.


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

The processing and analysis code is all located in `notebooks`. The order is `preprocessing` > `psf_subtraction`. Some helper library code is located in `src`. All of the paper figures were produced from scripts in `paper/figures`.

## Appendix figure sets

Many figures were excluded from the PDF version of the paper for brevity. The ADI processing reports are in `paper/figures/reports`, along with the PCA/NMF/GreeDS mosaics and reports.
