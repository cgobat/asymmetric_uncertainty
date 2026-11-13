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
Separately constructed quantities start as independent uncertainty sources.
Arithmetic retains signed derivatives with respect to these sources, combines
contributions from each shared source, and then combines independent sources in
quadrature to first order. The sign of each combined derivative determines
which side of a source operand contributes to each side of the result.

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
or nonlinear transformation of a probability density. The statistical model
associated with an ``a_u`` value, and the distinction between local uncertainty
propagation and PDF propagation, are discussed in :doc:`theory`.

The package recognizes the obvious fully correlated cases ``x + x``, ``x - x``,
``x * x``, ``x / x``, and ``x ** x`` when both operands are the same ``a_u``
object. These reduce to ``2*x``, an exact zero, ``x**2``, an exact one (except
exact zero divided by itself), and propagation through the single-variable
function ``x**x``, respectively.

Distribution utilities
----------------------

For distributional calculations, ``value`` is interpreted as the mode and
``minus`` and ``plus`` as the Gaussian widths of an associated continuous
split-normal distribution. The PDF, parameter interpretation, and differences
from other asymmetric-Gaussian conventions are described in :doc:`theory`.

The analytic :meth:`~asymmetric_uncertainty.a_u.cdf` is complemented by
:meth:`~asymmetric_uncertainty.a_u.ppf`, which evaluates the inverse CDF, and
:meth:`~asymmetric_uncertainty.a_u.rvs`, which draws random values from the
associated distribution. ``rvs`` accepts either a seed or a
:class:`numpy.random.Generator` through its ``rng`` argument.

The :meth:`~asymmetric_uncertainty.a_u.fit` class method estimates the same
split-normal distribution from a one-dimensional sample by maximum likelihood,
returning an ``a_u`` with the fitted mode and Gaussian widths. This is a
parametric approximation to the sample distribution; the returned uncertainties
are widths of that distribution, not errors on the fitted parameters.

:meth:`~asymmetric_uncertainty.a_u.interval` returns mode-centered intervals
scaled by the lower and upper Gaussian widths, while
:meth:`~asymmetric_uncertainty.a_u.equal_tail_interval` returns quantile-based
equal-tail intervals. These coincide for symmetric uncertainties but generally
differ for an asymmetric distribution.

Monte Carlo propagation
-----------------------

The :func:`~asymmetric_uncertainty.propagate` function provides an explicit
Monte Carlo route when the split-normal interpretation should be propagated as
a probability distribution rather than as local asymmetric uncertainties. It
replaces each ``a_u`` argument with samples drawn from its split-normal
distribution and evaluates a supplied callable on those arrays. For example:

.. code-block:: python

   samples = propagate(lambda x, y: x / y, A, B, size=10000, rng=42)

By default, ``propagate`` returns the callable's output directly, preserving
the Monte Carlo samples for analysis without imposing a parametric model. Set
``fit=True`` to instead return an ``a_u`` obtained by fitting a split-normal
distribution to a one-dimensional sample:

.. code-block:: python

   result = propagate(lambda x, y: x + y, A, B, size=10000, rng=42, fit=True)

This fit approximates the output PDF and can change the nominal value as well
as the uncertainties. It may poorly represent distributions that are not
approximately split-normal. For multi-step calculations, propagate the complete
expression and fit only the final result to retain intermediate correlations.

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

Upper and lower limits
----------------------

The ``a_u`` class has no dedicated statistical model for handling upper or
lower limits. The table below shows several possible ways to encode a limit as
a nominal value with asymmetric errors. There are tradeoffs to each and none
are uniquely correct representations.

.. list-table::
   :header-rows: 1
   :widths: 24 24 52

   * - Upper limit
     - Lower limit
     - Notes
   * - :math:`x_{-\infty}^{+0}`
     - :math:`x_{-0}^{+\infty}`
     - Infinities are not always well-behaved.
   * - :math:`0_{-0}^{+x}` if :math:`x > 0`
     - :math:`0_{-|x|}^{+0}` if :math:`x < 0`
     - Be careful of division-by-zero errors. This makes sense only when the
       quantity has a true physical limit at zero, such as a non-negative flux.
   * - :math:`x_{-x}^{+0}` if :math:`x > 0`
     - :math:`x_{-0}^{+|x|}` if :math:`x < 0`
     - Poor representation of a true limit when the uncertainty is of the same
       order as the nominal value.
