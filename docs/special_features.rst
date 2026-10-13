Special features
================

In addition to the methods and attributes listed in the :doc:`api`, the
:class:`asymmetric_uncertainty.a_u` class provides several convenience
features.

Propagation model
-----------------

An ``a_u`` object is primarily a nominal numeric value with separate positive
and negative uncertainties. Arithmetic applies the ordinary operation to the
nominal values, so, for example, values of 2 and 3 always add to a value of 5.
Distinct operands are treated as independent when their uncertainty
contributions are combined in quadrature to first order. The sign of each local
derivative determines which side of an operand contributes to the positive and
negative uncertainty of the result.

For independent operands this gives, for example,

.. math::

   \sigma_{+,A+B}^2 = \sigma_{+,A}^2 + \sigma_{+,B}^2,
   \qquad
   \sigma_{-,A+B}^2 = \sigma_{-,A}^2 + \sigma_{-,B}^2,

and

.. math::

   \sigma_{+,A-B}^2 = \sigma_{+,A}^2 + \sigma_{-,B}^2,
   \qquad
   \sigma_{-,A-B}^2 = \sigma_{-,A}^2 + \sigma_{+,B}^2.

The same first-order convention is used for multiplication, division, powers,
and the logarithmic/exponential helpers. Arithmetic therefore propagates the
quoted local uncertainties; it does **not** claim to reproduce the convolution
or nonlinear transformation of a probability density.

The package recognizes the obvious fully correlated cases ``x + x``, ``x - x``,
``x * x``, ``x / x``, and ``x ** x`` when both operands are the same ``a_u``
object. These reduce to ``2*x``, an exact zero, ``x**2``, an exact one (except
exact zero divided by itself), and propagation through the single-variable
function ``x**x``, respectively.

Distribution utilities
----------------------

For distributional calculations, ``value`` is interpreted as the mode and
``minus`` and ``plus`` as the Gaussian widths of an associated jointly
normalized split-normal distribution. This is a choice for the interpretation
of the stored value and uncertainties rather than a propagation model used by
ordinary arithmetic.

The analytic :meth:`~asymmetric_uncertainty.a_u.cdf` is complemented by
:meth:`~asymmetric_uncertainty.a_u.ppf`, which evaluates the inverse CDF, and
:meth:`~asymmetric_uncertainty.a_u.rvs`, which draws random values from that
associated split-normal distribution. Sampling uses the unequal probability
mass on either side of the mode, so the fraction of samples below ``value`` is
``minus / (minus + plus)`` rather than one half for an asymmetric value.

``rvs`` accepts either a seed or a :class:`numpy.random.Generator` through its
``rng`` argument, allowing reproducible simulations without introducing a
SciPy dependency.

:meth:`~asymmetric_uncertainty.a_u.interval` returns intervals centered on the
mode and scaled independently by the lower and upper Gaussian widths. The
interval can be specified either by a number of sigma or by the corresponding
contained probability. For example, ``interval(n_sigma=1)`` contains the same
:math:`\operatorname{erf}(1/\sqrt{2})` probability as a one-sigma Gaussian
interval, even when the two widths differ.

This is distinct from an equal-tail credible interval. Use
:meth:`~asymmetric_uncertainty.a_u.equal_tail_interval` when the probability
excluded below the lower bound should equal the probability excluded above the
upper bound. With no confidence specified, it uses the same
:math:`\operatorname{erf}(1/\sqrt{2})` probability as a one-sigma Gaussian
interval. For an asymmetric distribution, an equal-tail interval is not in
general centered on ``value``.

Monte Carlo propagation
-----------------------

The :func:`~asymmetric_uncertainty.propagate` function provides an explicit
Monte Carlo route when the split-normal interpretation should be propagated as
a probability distribution rather than as local asymmetric uncertainties. It
replaces each ``a_u`` argument with samples drawn from its split-normal
distribution and evaluates a supplied callable on those arrays. For example:

.. code-block:: python

   samples = propagate(lambda x, y: x / y, A, B, size=10000, rng=42)

The raw return value from the callable is returned directly; ``propagate`` does
not fit the Monte Carlo result back to an ``a_u`` object. This makes the
function useful for diagnostics and direct comparison with the package's
first-order arithmetic without imposing another approximation on the sampled
output.

Each distinct ``a_u`` object is sampled once per call. Passing the same object
multiple times therefore reuses the same sample realization and preserves
perfect correlation, while separate objects are sampled independently.

Intelligent display formatting
------------------------------

In environments that support rich display, such as Jupyter notebooks, an
``a_u`` instance renders itself using LaTeX-style mathematical notation.
For example:

.. code-block:: python

   A = a_u(3.0, 0.2, 0.4)
   A

is displayed as

.. math::

   3.0_{-0.4}^{+0.2}

while a symmetric uncertainty:

.. code-block:: python

   B = a_u(3.0, 0.3, 0.3)
   B

is displayed as

.. math::

   3.0 \pm 0.3

"Strict" inequalities
----------------------

The ``<<`` and ``>>`` operators are overloaded for ``a_u`` instances to
represent comparisons in which the quoted error intervals do not overlap.
This mirrors the common mathematical notation :math:`\ll` and :math:`\gg`.

For two quantities
:math:`A_{-\sigma_{A_-}}^{+\sigma_{A_+}}` and
:math:`B_{-\sigma_{B_-}}^{+\sigma_{B_+}}`, the package evaluates

.. math::

   A \ll B
   \iff
   A + \sigma_{A_+} < B - \sigma_{B_-},

and similarly

.. math::

   A \gg B
   \iff
   A - \sigma_{A_-} > B + \sigma_{B_+}.

For example:

.. code-block:: python

   A = a_u(3.0, 0.2, 0.4)
   B = a_u(2.5, 0.3, 0.3)

   A > B
   # True

   A >> B
   # False

The ordinary comparison operators compare nominal values only; ``<<`` and
``>>`` additionally require separation between the quoted uncertainty
intervals.
