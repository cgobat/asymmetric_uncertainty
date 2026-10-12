Related software
================

Several packages address asymmetric uncertainties, but the notation
:math:`x_{-\sigma_-}^{+\sigma_+}` does not by itself specify a unique
statistical model. The quoted bounds may represent widths of an assumed
probability distribution, confidence or credible intervals, or intervals
inferred from a likelihood function. Packages that accept superficially
similar inputs therefore need not be statistically interchangeable.

``asymmetric_uncertainty`` treats ``value`` as the mode of a jointly normalized
split-normal distribution and ``minus`` and ``plus`` as its Gaussian widths on
the two sides of the mode. Addition and subtraction of independent values are
approximated by matching the first three moments of the convolution to another
split normal; nonlinear arithmetic uses first-order derivative propagation.
The packages below cover related use cases with different statistical models or
propagation strategies.

.. list-table::
   :header-rows: 1
   :widths: 20 24 25 31

   * - Package
     - Statistical representation
     - Propagation approach
     - Main difference from ``asymmetric_uncertainty``
   * - `asymmetric_uncertainties <https://github.com/muryelgp/asymmetric_uncertainties>`_
     - Barlow's Linear Variable Width Gaussian, with uncertainties specified at
       a chosen confidence level.
     - Monte Carlo sampling for arithmetic and user-defined functions.
     - Uses a Barlow variable-width model and numerical sampling rather than a
       split-normal object with primarily analytic/local propagation.
   * - `add_asym <https://github.com/anisotropela/add_asym>`_
     - Linear or quadratic transformed-Gaussian models following Barlow (2003).
     - Deterministic moment-based combination of independent asymmetric errors.
     - Specialized for addition rather than providing a general numeric type or
       arbitrary arithmetic. Its treatment of moments is the closest conceptual
       relative of this package's addition/subtraction model.
   * - `SOAD <https://github.com/kiyami/soad>`_
     - Barlow Variable Gaussian representation.
     - Monte Carlo propagation through arithmetic operations, with fitted output
       distributions and confidence-level utilities.
     - Keeps sampled distributions as the basis for propagation and includes
       plotting/confidence-level tools rather than using local analytic error
       propagation for nonlinear operations.
   * - `augment <https://github.com/travisthieme/augment>`_
     - Sampled split-normal, normal, or log-normal inputs, with optional physical
       bounds.
     - Monte Carlo evaluation of user-supplied model functions.
     - Model/workflow oriented rather than scalar-object oriented; additionally
       provides Astropy units, sample summaries, plots, and rank-correlation
       error budgets.
   * - `AsymmetricErrorsPy <https://github.com/igvgit/AsymmetricErrorsPy>`_
     - A large collection of probability-distribution and log-likelihood models
       accompanying `Barlow, Brazzale, and Volobouev (2026)
       <https://doi.org/10.1016/j.nima.2025.170857>`_.
     - Numerical routines for combining errors and measurements, backed by the
       C++ `AsymmetricErrors <https://github.com/igvgit/AsymmetricErrors>`_
       library.
     - Explicitly distinguishes PDF-based and likelihood-based asymmetric errors
       and supports many model families; it is a broader statistical framework
       rather than a lightweight scalar type built around one distribution.
   * - `BigDog <https://pypi.org/project/bigdog/>`_
     - Uncertain values are represented directly by sample arrays, including
       asymmetric values with physical units.
     - Element-wise Monte Carlo propagation through overloaded arithmetic.
     - Stores the sampled distribution itself, naturally retaining dependence
       when a sampled object is reused; also provides Astropy units and physical
       domain constraints, at greater sampling and memory cost.
   * - `split-normal <https://pypi.org/project/split-normal/>`_
     - The split-normal/two-piece-normal distribution, with NumPy and JAX
       implementations.
     - Provides distribution functions such as PDF, CDF, and PPF rather than an
       uncertainty-propagation system.
     - Uses the same general distribution family, but does not provide an
       uncertain numeric object or arithmetic error propagation.

The widely used `uncertainties <https://pythonhosted.org/uncertainties/>`_
package is also relevant as a symmetric-error counterpart. It propagates
standard deviations using derivatives and tracks correlations between derived
quantities, but it does not represent separate upper and lower uncertainties.

The moment-matching approach used here for independent addition and subtraction
is closest in motivation to ``add_asym`` and Barlow's work on combining
asymmetric errors. The implementation remains different: the result is
projected back into the split-normal family already used by
``asymmetric_uncertainty`` rather than adopting Barlow's transformed-Gaussian or
variable-width-Gaussian parameterizations.
