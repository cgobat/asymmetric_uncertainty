# Version history

Notable changes to `asymmetric_uncertainty` are documented here.

## 0.5.0 - 2026-10-13

[Compare 0.4.1...0.5.0](https://github.com/cgobat/asymmetric_uncertainty/compare/0.4.1...0.5.0)

### Added

- Added `propagate()` for Monte Carlo propagation using the split-normal
  interpretation of `a_u` uncertainties.
- Added recognition of simple self-correlation cases such as `x + x`,
  `x - x`, `x * x`, `x / x`, and `x ** x`.
- Added `median` and Fisher/excess `kurtosis` properties for the associated
  split-normal distribution.
- Allowed `pdfplot()` and `cdfplot()` to plot onto an existing Matplotlib
  `Axes`.
- Added NumPy 1.x and 2.x coverage to CI.
- Added a Google Colab link for the example notebook.
- Added `scinum` to the related-software documentation.

### Changed

- Restored addition and subtraction to the package's standard first-order
  asymmetric error-propagation model, reversing the moment-matched
  split-normal propagation introduced in 0.3.0. Nominal values therefore
  once again obey ordinary arithmetic exactly.
- Clarified the distinction between ordinary arithmetic propagation and the
  associated split-normal probability-distribution interpretation.
- Improved string representations and formatting.
- Expanded and standardized API docstrings.
- Reformatted `core.py` with Ruff.

### Fixed

- Replaced `numpy.concat` with `numpy.concatenate` to preserve NumPy 1.x
  compatibility.

### Documentation

- Updated the README and API documentation for the revised propagation model.
- Added a documentation-site favicon and various presentation improvements.


## 0.4.1 - 2026-10-12

[Compare 0.4.0...0.4.1](https://github.com/cgobat/asymmetric_uncertainty/compare/0.4.0...0.4.1)

### Added

- Added a documentation page comparing `asymmetric_uncertainty` with related
  software.

### Changed

- Moved the test suite to a top-level `tests/` directory.
- Updated the example notebook for the new addition and subtraction behavior
  introduced in 0.3.0.
- Improved documentation references and links.

### Fixed

- Fixed the `\pm` escape sequence in the LaTeX representation.
- Ensured comparison-style properties and operators return built-in `bool`
  values rather than `numpy.bool_` values.


## 0.4.0 - 2026-10-12

[Compare 0.3.0...0.4.0](https://github.com/cgobat/asymmetric_uncertainty/compare/0.3.0...0.4.0)

### Added

- Added inverse-CDF/quantile evaluation with `ppf()`.
- Added random sampling from the associated split-normal distribution with
  `rvs()`.
- Added mode-centered interval utilities.
- Added equal-tail confidence intervals.
- Added tests for the new distribution utilities.


## 0.3.0 - 2026-10-12

[Compare 0.2.10...0.3.0](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.10...0.3.0)

### Added

- Added analytical split-normal distribution properties including mean,
  variance, standard deviation, and skewness.
- Added moment-matched split-normal propagation for addition and subtraction.
  Note that this behavior was subsequently reverted in 0.5.0.
- Added `codemeta.json`.

### Changed

- Migrated package configuration from `setup.py` to `pyproject.toml`.
- Migrated project documentation from the GitHub wiki to Read the Docs.
- Reorganized and substantially expanded the documentation.
- Adopted the PyData Sphinx theme and improved site navigation.
- Updated legacy ADS references to SciX.


## 0.2.10 - 2026-08-12

[Compare 0.2.9...0.2.10](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.9...0.2.10)

### Fixed

- Ensured derived `a_u` state such as `maximum`, `minimum`, `sign`, and
  symmetry information remains correct after mutation.
- Ensured `UncertaintyArray` cached values and uncertainties remain
  synchronized after list mutations.


## 0.2.9 - 2026-08-12

[Compare 0.2.8...0.2.9](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.8...0.2.9)

### Changed

- Reworked `UncertaintyArray` to use actual `list` storage rather than
  maintaining a separate backing representation.
- Restored compatibility with NumPy 1.x after the NumPy 2-only requirement
  introduced in 0.2.5.
- Updated CI accordingly.


## 0.2.8 - 2026-08-12

[Compare 0.2.7...0.2.8](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.7...0.2.8)

### Added

- Added propagation support for uncertain exponents.
- Added reverse exponentiation support.

### Fixed

- Added explicit domain handling for exponentiation, including negative and
  zero bases.
- Made power propagation use the appropriate derivatives for both base and
  exponent uncertainties.


## 0.2.7 - 2026-08-12

[Compare 0.2.6...0.2.7](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.6...0.2.7)

### Fixed

- Changed asymmetric error propagation to account for the sign of each
  derivative.
- Corrected positive/negative uncertainty assignment for operations whose
  local derivative is negative, including multiplication and division by
  negative values.


## 0.2.6 - 2026-08-12

[Compare 0.2.5...0.2.6](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.5...0.2.6)

### Fixed

- Reworked multiplication and division propagation to avoid relative-error
  expressions that divide by zero.
- Corrected propagation for zero-valued operands and numerators.


## 0.2.5 - 2026-08-12

[Compare 0.2.4...0.2.5](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.4...0.2.5)

### Fixed

- Corrected normalization of the asymmetric Gaussian/split-normal CDF.

### Changed

- Temporarily raised the NumPy requirement to NumPy 2.x. NumPy 1.x
  compatibility was restored in 0.2.9.


## 0.2.4 - 2026-01-06

[Compare 0.2.3...0.2.4](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.3...0.2.4)

### Added

- Added a unit-test suite and GitHub Actions CI.
- Added richer string formatting and output options.

### Fixed

- Made propagation through negative powers consistent with the equivalent
  reciprocal/division expression.

### Changed

- Updated CI actions and workflow configuration.


## 0.2.3 - 2023-08-10

[Compare 0.2.2...0.2.3](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.2...0.2.3)

### Fixed

- Corrected problems in the CDF implementation.


## 0.2.2 - 2023-08-03

[Compare 0.2.1...0.2.2](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.1...0.2.2)

### Fixed

- Fixed `a_u.notna()`.
- Fixed a package dependency-installation ordering problem.

### Changed

- Added a deprecation warning for the legacy `AsymmetricUncertainty` class
  name.
- Expanded and clarified README and supporting documentation.


## 0.2.1 - 2022-09-13

[Compare 0.2.0...0.2.1](https://github.com/cgobat/asymmetric_uncertainty/compare/0.2.0...0.2.1)

### Added

- Added `exp()`.
- Added package metadata and citation information.
- Added documentation of package limitations.
- Added examples for the overloaded `<<` and `>>` comparisons.
- Exposed package metadata from the top-level module.

### Fixed

- Fixed a relative-import issue.

### Documentation

- Substantially expanded the README, documentation, supporting material, and
  citation information.


## 0.2.0 - 2022-08-26

[Compare 0.1.0...0.2.0](https://github.com/cgobat/asymmetric_uncertainty/compare/0.1.0...0.2.0)

### Changed

- Reorganized the package backend, including moving the implementation into
  `core.py`.
- Updated and expanded project documentation and supporting material.


## 0.1.0 - 2022-08-05

[Compare 0.0.2...0.1.0](https://github.com/cgobat/asymmetric_uncertainty/compare/0.0.2...0.1.0)

### Changed

- Renamed the primary `AsymmetricUncertainty` class to `a_u`.
- Updated examples and package exports for the new class name.
- Added README status badges.


## 0.0.2 - 2022-07-25

[Compare 0.0.1...0.0.2](https://github.com/cgobat/asymmetric_uncertainty/compare/0.0.1...0.0.2)

### Added

- Added an example Jupyter notebook demonstrating package use.

### Changed

- Removed the pandas dependency.
- Renamed the implementation module from `source.py` to `src.py`.
- Updated the README and repository configuration.


## 0.0.1 - 2022-07-09

Initial tagged release.

### Added

- Added the `AsymmetricUncertainty` class for representing numeric values with
  separate positive and negative uncertainties.
- Added arithmetic propagation for addition, subtraction, multiplication,
  division, and powers.
- Added PDF and CDF evaluation and plotting utilities for the asymmetric
  pseudo-Gaussian distribution.
- Added string and LaTeX representations.
- Added utilities for adding additional uncertainty and extracting positive
  and negative errors from collections of uncertain values.
