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

Shared sources in derived quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No correlation argument is needed when expressions share an input:

.. code-block:: python

   x = a_u(3, 0.2, 0.4)
   y = 2 * x
   y - x             # a_u(3.0, 0.2, 0.4), to within machine precision
   x - y             # a_u(-3.0, 0.4, 0.2), i.e. -x
   (x + 5) - x       # exactly 5
   (2 * x) / x       # exact 2 at the nominal point

For a source :math:`s_i`, the chain rule combines the sensitivities of the
operands :math:`A` and :math:`B` to a result :math:`R=f(A,B)`:

.. math::

   D_i = \frac{\partial f}{\partial A}\frac{\partial A}{\partial s_i}
       + \frac{\partial f}{\partial B}\frac{\partial B}{\partial s_i}.

The positive variance contribution is :math:`(D_i\sigma_{i,+})^2` if
:math:`D_i\ge0`, or :math:`(D_i\sigma_{i,-})^2` otherwise; the negative
contribution uses the opposite width. Contributions from independent sources
are summed. Combining the signed derivatives **before** choosing the widths
allows shared contributions to cancel even for asymmetric errors. With disjoint
sources, this reproduces the existing directional quadrature convention.
Negation, arithmetic, powers (including uncertain exponents), and the
logarithmic/exponential helpers all retain this information.

The existing same-object special cases remain: ``x - x`` is exact zero,
``x / x`` is exact one unless ``x`` is exact zero, and ``x ** x`` uses its
single-variable derivative. Distinct derived expressions still obey ordinary
nominal-value domain checks: ``x / (x + 0)`` raises at nominal zero. First-order
limitations also remain: ``x**2`` has zero propagated uncertainty at nominal
zero, and no higher-order curvature or transformed PDF is retained.

Source identity and boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Separate constructor calls are independent even if all three numbers agree.
  ``a_u(*x.items())`` explicitly makes an independent quantity.
* Shallow copies, deep copies, and pickle round trips preserve source identity.
  Equality still compares only ``value``, ``plus``, and ``minus``; it does not
  indicate independence or correlation.
* Assigning to ``value``, ``plus``, or ``minus`` starts a new independent source
  from the object's current state, even when assigning the same number. Earlier
  results retain snapshots of their original sources and are not recomputed.
  Use arithmetic such as ``x + offset`` to retain correlations under a shift.
* ``add_error()`` starts a new independent source in both modes. Its straight
  and split adjustments have no unique source decomposition, so this limited
  implementation does not infer a correlation with the original quantity.
* Distribution utilities and Monte Carlo ``propagate()`` retain their existing
  semantics. To sample a correlated nonlinear expression, pass its original
  independent inputs and put the full expression inside the callable, e.g.
  ``propagate(lambda x: 2*x - x, x)``. Passing ``y`` and ``x`` as distinct
  arguments does not reuse the arithmetic source information.

This is local derivative tracking for shared inputs, not a general covariance
API. Correlations between separately measured inputs cannot be specified.
Non-finite values and nondifferentiable expressions do not acquire a new
statistical interpretation. Source metadata is private, and serialization of
that metadata is not a stable cross-version interchange format.

Implementation tradeoffs
~~~~~~~~~~~~~~~~~~~~~~~~

A source label alone can detect overlap but cannot determine its magnitude or
sign; a single correlation coefficient per result cannot represent several
partially shared inputs. Retaining a full expression graph could support later
reevaluation and joint nonlinear sampling, but would require a larger change to
mutation, serialization, and distribution semantics.

The limited implementation instead stores a sparse mapping of source identities
to signed sensitivities and immutable pairs of input widths. It retains no
operand objects or expression graph. UUID identities survive copying and pickle
serialization without relying on reusable memory addresses. Storage per result
and propagation work grow linearly with the number of contributing sources;
exactly zero sensitivities are dropped, without tolerance-based clipping. A
large reduction over many independent sources therefore costs more than the
previous constant-size representation. Explicit covariance input and joint
sampling of derived results remain separate future work.

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
