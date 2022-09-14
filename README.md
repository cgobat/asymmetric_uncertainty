# Asymmetric Uncertainty

[![GPLv3 License](https://img.shields.io/github/license/cgobat/asymmetric_uncertainty)](https://opensource.org/licenses/GPL-3.0) [![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) ![Python 3+](https://img.shields.io/badge/made%20with-Python%203-blue) [![ascl:2208.005](https://img.shields.io/badge/ascl-2208.005-blue.svg?colorB=262255)](https://ascl.net/2208.005)

In science and engineering, many values and quantities have uncertainties associated with them—for example, the expression $3.0\pm0.3$ denotes that a number whose measured or expected value is 3.0, but might reasonably be expected[^1] to have a true value anywhere between 2.7 and 3.3. This is a *symmetric* uncertainty, because the error is symmetric on either side of the central value. Not all numbers behave this way: for example, we could write $3.0_{-0.4}^{+0.2}$ to signify that although the expected value is still 3.0, we have inequal probabilities of the true value being above or below this mark.

## Usage

This repository contains the `asymmetric_uncertainty` Python package, which implements the `a_u` class for dealing with these kinds of numbers in practice. Usage is very simple:

```python
from asymmetric_uncertainty import a_u

x = a_u(3.0, 0.2, 0.4)

print(x + 1)
```

Here we see how to create an instance of the class representing the example number given before: $3.0_{-0.4}^{+0.2}$. The initialization arguments are the nominal value, positive error, and negative error—in that order. `x` will play nicely with other numeric objects under most mathematical operations, and its errors will propagate appropriately.

Further background information can be found in the [supporting materials](./supporting_matl.md).

## Installation

If you have a local installation of Git, you can install this package with a single command: run `pip install git+https://github.com/cgobat/asymmetric_uncertainty.git` to install the latest version.

Otherwise, clone/download this repository, then (from a command line running within the associated directory) run `pip install .` or `python setup.py install`.

`asymmetric_uncertainty` should then be available as a module that you can import like any other (see [§Usage](#usage) or [`example.ipynb`](./example.ipynb) for examples).

## [Documentation](../../wiki)

## [Citing this project](./CITATION.bib)

If you make use of this package in your own work, please cite it by referencing its entry in ASCL, the Astrophysics Source Code Library: [ascl:2208.005](http://ascl.net/2208.005). Citation information can also be exported from ADS ([2022ascl.soft08005G](https://ui.adsabs.harvard.edu/abs/2022ascl.soft08005G/exportcitation)) or the sidebar of this repository.

[^1]: Definitions and conventions for what this means vary. Common uses/interpretations include standard deviation, standard deviation divided by the square root of the number of samples (often called standard error), 90% confidence interval, etc. As long as usage is consistent and intentions are clear, any of these are valid.