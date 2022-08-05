# Asymmetric Uncertainty

[![GPLv3 License](https://img.shields.io/github/license/cgobat/asymmetric_uncertainty)](https://opensource.org/licenses/GPL-3.0) [![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) ![Python 3+](https://img.shields.io/badge/made%20with-Python%203-blue)

In science and engineering, many values and quantities have uncertainties associated with them—for example, the expression $3.0\pm0.3$ denotes that a number whose measured or expected value is 3.0, but might reasonably be expected to have a true value anywhere between 2.7 and 3.3. This is a *symmetric* uncertainty, because the error is symmetric to either side of the central value. Not all numbers behave this way: for example, we could write $3.0_{-0.4}^{+0.2}$ to signify that although the expected value is still 3.0, we have non-equal probabilities of the true value being above or below this mark.

## This package

This repository contains the `asymmetric_uncertainty` Python package, which implements the `AsymmetricUncertainty` class for dealing with these kinds of numbers in practice. Usage is very simple:

```python
from asymmetric_uncertainty import AsymmetricUncertainty

A = AsymmetricUncertainty(3.0,0.2,0.4)

print(A + 1)
```

Here we see how to create an instance of the class representing the example number given before: $3.0_{-0.4}^{+0.2}$. The arguments are the nominal value, positive error, and negative error—in that order. `A` will now play nicely with other numeric objects under most mathematical operations, and its errors will propagate appropriately.

Further background information can be found in the [supporting materials document](./supporting_matl.md).

## Installation

Clone this repository, then (from a command line running within the associated directory), run `pip install .` or `python setup.py install`.[^1] `asymmetric_uncertainty` will then be available as a module that you can import like any other.

[^1]: Alternatively, if you have Git installed, you can also remotely install using `pip install git+https://github.com/cgobat/asymmetric_uncertainty.git`