Special features
================

In addition to the methods and attributes listed in the :doc:`api`, the
:class:`asymmetric_uncertainty.a_u` class provides several convenience
features.

Propagation model
-----------------

An ``a_u`` value represents a jointly normalized split-normal distribution.
The stored ``value`` is the mode of that distribution, and ``minus`` and
``plus`` are the Gaussian widths below and above the mode.

Addition and subtraction treat distinct operands as independent random
variables. Their exact convolution is generally not another split normal, so
the result is approximated by a split normal with the same mean, variance, and
third central moment as the convolution. Consequently, the returned ``value``
need not equal the arithmetic sum or difference of the two input modes when the
operands are asymmetric.

The package recognizes the obvious fully correlated cases ``x + x``, ``x - x``,
``x * x``, ``x / x``, and ``x ** x`` when both operands are the same ``a_u``
object. These reduce to ``2*x``, an exact zero, ``x**2``, an exact one (except
exact zero divided by itself), and propagation through the single-variable
function ``x**x``, respectively.

Multiplication, division, powers, and the logarithmic/exponential helpers use
first-order derivative propagation about the nominal values instead. These
operations therefore propagate the quoted local uncertainty rather than
attempting to reproduce the full transformed probability distribution.

Distribution utilities
----------------------

The analytic :meth:`~asymmetric_uncertainty.a_u.cdf` is complemented by
:meth:`~asymmetric_uncertainty.a_u.ppf`, which evaluates the inverse CDF, and
:meth:`~asymmetric_uncertainty.a_u.rvs`, which draws random values from the
same jointly normalized split-normal distribution. Sampling uses the unequal
probability mass on either side of the mode, so the fraction of samples below
``value`` is ``minus / (minus + plus)`` rather than one half for an asymmetric
value.

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
