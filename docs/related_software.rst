Related software
================

Several packages address asymmetric uncertainties, but the notation
:math:`x_{-\sigma_-}^{+\sigma_+}` does not by itself specify a unique
statistical model or even a particular meaning for the quoted bounds. The
comparison below deals with both **representation** (an assumed distribution,
likelihood specification, error interval, etc.) and **propagation** (local
derivatives, moment-matching, Monte Carlo sampling, etc.) differences. See
:doc:`theory` for the statistical distinctions and :doc:`supporting_material`
for (some of) the relevant literature.

Comparison
----------

.. list-table::
   :header-rows: 1
   :widths: 20 24 25 31

   * - Package
     - Statistical representation
     - Propagation approach
     - Main difference from ``asymmetric_uncertainty``
   * - `asymmetric_uncertainties <https://github.com/muryelgp/asymmetric_uncertainties>`_
     - Barlow's Linear Variable Width Gaussian [2]_, with uncertainties specified
       at a chosen confidence level.
     - Monte Carlo sampling for arithmetic and user-defined functions.
     - Uses a Barlow variable-width model and numerical sampling rather than
       local asymmetric-error propagation with optional split-normal utilities.
   * - `add_asym <https://github.com/anisotropela/add_asym>`_
     - Linear or quadratic transformed-Gaussian models following Barlow (2003). [1]_
     - Deterministic moment-based combination of independent asymmetric errors.
     - Specialized for statistically motivated combination by addition rather
       than providing a general numeric type whose nominal values obey ordinary
       arithmetic semantics.
   * - `SOAD <https://github.com/kiyami/soad>`_
     - Barlow Variable Gaussian representation. [2]_
     - Monte Carlo propagation through arithmetic operations, with fitted output
       distributions and confidence-level utilities. The package and methodology
       are described by Erdim and Hüdaverdi (2019). [3]_
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
       accompanying Barlow, Brazzale, and Volobouev (2026). [4]_
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
   * - `scinum <https://github.com/riga/scinum>`_
     - Numeric values can carry multiple named symmetric or asymmetric
       uncertainties, specified as either absolute or relative effects.
     - Gaussian derivative-based propagation with configurable correlations
       between uncertainty sources, plus a library of uncertainty-aware math
       operations.
     - Tracks multiple uncertainty components and their correlations rather than
       interpreting a single pair of errors as one probability distribution; it
       also supports NumPy arrays and scientific formatting/rounding.
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

.. rubric:: References

.. [1] R. Barlow (2003), *Asymmetric Systematic Errors*, arXiv:physics/0306138.
       `bibcode:2003physics...6138B <https://scixplorer.org/abs/2003physics...6138B/abstract>`_.

.. [2] R. Barlow (2004), *Asymmetric Statistical Errors*, arXiv:physics/0406120.
       `bibcode:2004physics...6120B <https://scixplorer.org/abs/2004physics...6120B/abstract>`_.

.. [3] M. K. Erdim and M. Hüdaverdi (2019), *AIP Conference Proceedings*
       **2178**, 030023.
       `bibcode:2019AIPC.2178c0023E <https://scixplorer.org/abs/2019AIPC.2178c0023E/abstract>`_.

.. [4] R. Barlow, A. R. Brazzale, and I. Volobouev (2026), *Nuclear Instruments
       and Methods in Physics Research A* **1082**, 170857.
       `bibcode:2026NIMPA108270857B <https://scixplorer.org/abs/2026NIMPA108270857B/abstract>`_.
