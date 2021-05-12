# Python4Python

**Python** is a radiative transfer and ionisation code written in *c*.
*Python* is a scripting language.
`py4py` is a module written in *Python* for reading, processing and visualising the input and output
files of the *c* code **Python**.

## Getting Started

Create a new virtual environnment and install the requirements using `pip install -r requirements.txt`
Then, use `pip install -e .` (or `Make install`) to install `py4py`.

If you would like to run the example Jupyter Notebook, use `pip install jupyter`.

## Examples

* A reverberation mapping example using a Jupyter Notebook can be found in [examples/reverb](../../examples/reverb/reverb.ipynb).

## Submodules

### reverb

This module handles the production of transfer functions
from the output data from **Python** simulation runs.

#### timeseries

This module handles producing artificial time series in *CARAMEL* and *MEMEcho*
format, to test deconvolution processes.