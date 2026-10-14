Statistical model and interpretation
====================================

The notation :math:`x_{-\sigma_-}^{+\sigma_+}` does not uniquely specify a
probability distribution. The quoted central value may be a mode, mean, median,
or best-fit estimator, while the two uncertainties may represent local scale
parameters, confidence limits, posterior quantiles, likelihood crossings, or
something else entirely. Different statistical models can therefore assign
meaningfully different probability distributions to the same three reported
numbers.

``asymmetric_uncertainty`` makes this distinction explicit. An
:class:`~asymmetric_uncertainty.a_u` object is first of all a nominal numeric
value with separate positive and negative uncertainties. When a probability
distribution is required, the package associates those three values with a
continuous two-piece normal (split-normal, or Fechner) distribution. Ordinary
arithmetic, by contrast, propagates the quoted local uncertainties and does not
claim to perform an exact convolution or nonlinear transformation of that PDF.

This page defines that statistical interpretation, contrasts it with other
common asymmetric-Gaussian conventions, and states the assumptions under which
the model is intended to be used.

The associated split-normal distribution
----------------------------------------

For an ``a_u`` value

.. math::

   X = c_{-\sigma_-}^{+\sigma_+},

``asymmetric_uncertainty`` interprets ``value`` as the mode :math:`c`, ``minus``
as the lower Gaussian scale :math:`\sigma_-`, and ``plus`` as the upper
Gaussian scale :math:`\sigma_+`. The associated PDF is

.. math::

   f(x) =
   \frac{\sqrt{2/\pi}}{\sigma_-+\sigma_+}
   \begin{cases}
      \exp\!\left[-\dfrac{(x-c)^2}{2\sigma_-^2}\right], & x<c, \\
      \exp\!\left[-\dfrac{(x-c)^2}{2\sigma_+^2}\right], & x\ge c.
   \end{cases}

This continuous two-piece normal distribution originates in Fechner's
posthumously published *Kollektivmasslehre* (1897). [1]_ A modern historical
account is given by Wallis (2014). [2]_ Possolo, Merkatas, and Bodnar use the
same Fechner distribution in a modern metrology treatment of asymmetric
measurement uncertainties, describing it as two half-normal densities rescaled
so that the resulting PDF is continuous. [3]_ The common normalization on the
two sides makes the density continuous at :math:`c`, which is its mode
("central" peak). The probability masses on the two sides are consequently

.. math::

   P(X<c) = \frac{\sigma_-}{\sigma_-+\sigma_+},
   \qquad
   P(X\ge c) = \frac{\sigma_+}{\sigma_-+\sigma_+}.

Thus the wider side of the distribution contains more probability mass. When
:math:`\sigma_-=\sigma_+`, the model reduces exactly to an ordinary normal
distribution.

The corresponding CDF is

.. math::

   F(x) =
   \begin{cases}
      \dfrac{2\sigma_-}{\sigma_-+\sigma_+}
      \Phi\!\left(\dfrac{x-c}{\sigma_-}\right), & x<c, \\
      \dfrac{\sigma_- - \sigma_+
      + 2\sigma_+\Phi\!\left(\dfrac{x-c}{\sigma_+}\right)}
      {\sigma_-+\sigma_+}, & x\ge c,
   \end{cases}

where :math:`\Phi` is the standard-normal CDF. In particular,

.. math::

   F(c) = \frac{\sigma_-}{\sigma_-+\sigma_+},

so the mode is not generally the median.

The mean and variance are [2]_

.. math::

   E[X] = c + \sqrt{\frac{2}{\pi}}(\sigma_+-\sigma_-),

and

.. math::

   \operatorname{Var}(X)
   = \sigma_-\sigma_+
   + \left(1-\frac{2}{\pi}\right)(\sigma_+-\sigma_-)^2.

For positive skew (:math:`\sigma_+>\sigma_-`), the mean and median lie above the
mode; for negative skew the ordering is reversed. The PDF is continuous at the
mode, although its derivatives are not all smooth there when the two scales
differ. [2]_ This is a limitation of the model, not an unacknowledged property.

What ``one sigma`` means in this model
--------------------------------------

The two scale parameters retain an exact Gaussian-width interpretation on their
respective sides. For any :math:`n\ge 0`,

.. math::

   P\left(c-n\sigma_- \le X \le c+n\sigma_+\right)
   = \operatorname{erf}\!\left(\frac{n}{\sqrt{2}}\right).

In particular,

.. math::

   [c-\sigma_-,\ c+\sigma_+]

contains the usual Gaussian one-sigma probability,
:math:`\operatorname{erf}(1/\sqrt{2})\simeq0.6827`, even for an asymmetric
value. This is the interval returned by
:meth:`~asymmetric_uncertainty.a_u.interval` with its default arguments.

This interval is **not** generally an equal-tail interval. Because
:math:`F(c)\ne 1/2` when the scales differ, :math:`c-\sigma_-` and
:math:`c+\sigma_+` are not generally the 15.87th and 84.13th percentiles.
:meth:`~asymmetric_uncertainty.a_u.equal_tail_interval` and
:meth:`~asymmetric_uncertainty.a_u.ppf` are provided when quantile-based
intervals are required.

This distinction is important when interpreting published error bars. If a
source explicitly defines its quoted lower and upper errors as distances from a
median to the 16th and 84th percentiles, those numbers are not, in general,
the split-normal scale parameters used here. Treating them as such is an
additional modeling assumption.

Possolo et al. make this distinction between a reported interval and a model's
parameters explicit. In their treatment, reported values such as
:math:`x-u_L`, :math:`x`, and :math:`x+u_R` are assigned specified percentile
meanings, and the parameters of a Fechner, skew-normal, or generalized
extreme-value distribution are fitted to reproduce those percentiles. [3]_
Thus even choosing the same Fechner family does not imply that reported
asymmetric error bars should automatically be identified with
:math:`\sigma_-` and :math:`\sigma_+`; that mapping depends on what the
reported numbers are intended to mean.

Terminology
-----------

Names for asymmetric Gaussian constructions are not used consistently across
the literature or software. ``split normal``, ``two-piece normal``,
``bifurcated Gaussian``, and ``double Gaussian`` have all been used for related
but not always identical constructions. [2]_ This documentation uses
*split normal* and *two-piece normal* specifically for the continuous Fechner
form defined above, and *dimidiated Gaussian* for the equal-area construction
defined below. The formula and parameter interpretation are therefore more
important than the name alone.

An alternative: the dimidiated Gaussian
---------------------------------------

A different convention is common in the asymmetric-error literature associated
with Barlow and in software descended from that approach. [4]_ [5]_ Starting
from a standard-normal variate and applying one scale to negative deviations
and another to positive deviations gives the *dimidiated* (or bifurcated)
Gaussian

.. math::

   f_{\mathrm{dim}}(x) =
   \begin{cases}
      \dfrac{1}{\sqrt{2\pi}\sigma_-}
      \exp\!\left[-\dfrac{(x-c)^2}{2\sigma_-^2}\right], & x<c, \\
      \dfrac{1}{\sqrt{2\pi}\sigma_+}
      \exp\!\left[-\dfrac{(x-c)^2}{2\sigma_+^2}\right], & x\ge c.
   \end{cases}

Each half contains exactly one half of the probability, so :math:`c` is the
median and

.. math::

   c-\sigma_- = F^{-1}(0.1587\ldots),
   \qquad
   c+\sigma_+ = F^{-1}(0.8413\ldots).

That convention therefore preserves an equal-tail interpretation of the quoted
``minus`` and ``plus`` values.

It preserves a different property at the center, however. The one-sided limits
of the PDF are

.. math::

   f_{\mathrm{dim}}(c^-)
   = \frac{1}{\sqrt{2\pi}\sigma_-},
   \qquad
   f_{\mathrm{dim}}(c^+)
   = \frac{1}{\sqrt{2\pi}\sigma_+}.

Unless :math:`\sigma_-=\sigma_+`, the PDF is therefore discontinuous at the
central value. Barlow, Brazzale, and Volobouev explicitly note this
characteristic of the dimidiated model and also emphasize that such
three-parameter forms are models rather than uniquely ``correct``
reconstructions of an unknown distribution. [5]_

The two conventions can be summarized as follows:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Property
     - Continuous split normal used here
     - Dimidiated Gaussian
   * - Meaning of :math:`c`
     - Mode
     - Median / splice point
   * - Density at :math:`c`
     - Continuous
     - Discontinuous when :math:`\sigma_-\ne\sigma_+`
   * - Probability below :math:`c`
     - :math:`\sigma_-/(\sigma_-+\sigma_+)`
     - :math:`1/2`
   * - Probability above :math:`c`
     - :math:`\sigma_+/(\sigma_-+\sigma_+)`
     - :math:`1/2`
   * - Meaning of :math:`c\pm\sigma_\pm`
     - Mode-centered interval containing 68.27% in total
     - 15.87th and 84.13th percentiles
   * - Property deliberately preserved
     - Common density and common mode
     - Equal probability on either side

Neither set of properties follows uniquely from the notation
:math:`c_{-\sigma_-}^{+\sigma_+}`. They represent different choices about what
information those three numbers are intended to encode.

Why use the continuous two-piece normal?
----------------------------------------

The split-normal interpretation used by ``asymmetric_uncertainty`` is intended
for situations in which the reported central value is naturally treated as the
most probable or best-fit location and the two uncertainties are naturally
interpreted as local Gaussian scales on either side.

Under those assumptions, requiring the two Gaussian pieces to meet at a common
height is a natural modeling choice. It produces a continuous, unimodal density
with a single mode at the reported central value while retaining Gaussian
falloff with the quoted scale on each side. This is an established probability
distribution with a history extending back to Fechner's nineteenth-century
work, rather than a package-specific construction. [1]_ [2]_ The Fechner
distribution also continues to be used in contemporary metrology as one of
several defensible models for asymmetric measurement uncertainty. [3]_

Continuity is not a theorem that makes this model uniquely correct. It is a
modeling constraint. The dimidiated model instead imposes equal probability on
the two sides and accepts a discontinuity at the center. Other three-parameter
models preserve still other quantities. Barlow, Brazzale, and Volobouev survey
several such alternatives and explicitly note that no chosen form can be
guaranteed to be the true underlying distribution. [5]_

No three-parameter asymmetric model is closed under propagation
---------------------------------------------------------------

Choosing an input PDF does not, by itself, solve the propagation problem. If
:math:`X` and :math:`Y` are independent random variables, the exact PDF of
:math:`Z=X+Y` is the convolution

.. math::

   f_Z(z) = \int_{-\infty}^{\infty} f_X(x)f_Y(z-x)\,dx.

Apart from special families such as the Gaussian, the resulting function does
not generally belong to the same parametric family as the inputs. Barlow,
Brazzale, and Volobouev state this explicitly for the asymmetric models they
consider and use moment matching to project a convolution back onto a chosen
three-parameter family when a compact representation is required. [5]_ That
projection is a coherent approximation, but it is still a projection rather
than an exact closure property.

The point is especially clear for the dimidiated Gaussian. An asymmetric
input PDF has a jump at its central value, whereas the convolution of two such
integrable densities is continuous. The exact sum of two asymmetric dimidiated
Gaussians therefore cannot, in general, itself be a dimidiated Gaussian.
Barlow, Brazzale, and Volobouev derive this convolution explicitly and show
that fitting another dimidiated Gaussian to it requires a new approximation;
the fitted central value and asymmetry change. [5]_

The continuous split-normal family used by ``asymmetric_uncertainty`` is also
not closed under convolution. Nor is either family closed under arbitrary
nonlinear transformations. For a monotonic transformation :math:`Y=g(X)`, the
exact transformed density is

.. math::

   f_Y(y)
   = f_X\!\left(g^{-1}(y)\right)
     \left|\frac{d g^{-1}}{dy}\right|,

which will not generally retain a Gaussian-on-each-side form.

Possolo et al. adopt a complementary strategy for this reason: fit an
appropriate asymmetric distribution to each input and propagate draws from
those distributions through the measurement model with Monte Carlo sampling.
[3]_ This approach propagates the input distributions without requiring the
output to remain in the same parametric family.

The relevant distinction is therefore not between one three-parameter model
that is ``correct under propagation`` and another that is not. In general, no
such asymmetric three-parameter family is preserved by the operations of
interest. The substantive choices are instead

* what the three input numbers mean;
* which properties the input PDF is intended to preserve;
* how uncertainties are propagated; and
* what approximation, if any, is made when the propagated result is compressed
  back to a small number of parameters.

Arithmetic and PDF propagation are deliberately separate
----------------------------------------------------------

Ordinary ``a_u`` arithmetic is designed to preserve ordinary numeric semantics.
For example, uncertain values with nominal values 2 and 3 add to a result with
nominal value 5. Positive and negative uncertainties are propagated locally
using first-order derivatives and quadrature, with the sign of each derivative
determining which side of an input contributes to each side of the result. See
:doc:`special_features` for the explicit propagation rules.

This convention is **not** a claim that the split-normal PDF associated with an
arithmetic result is the exact convolution or nonlinear transformation of the
input PDFs. It is a local asymmetric-error propagation rule. The associated
split-normal distribution of the result is a compact interpretation of the
resulting nominal value and directional scales.

This distinction also clarifies the usual Central Limit Theorem objection to
combining positive uncertainties in quadrature and negative uncertainties in
quadrature. Barlow and later Barlow, Brazzale, and Volobouev correctly point out
that repeatedly applying that rule while simultaneously claiming that the
result is the exact PDF convolution would retain a fixed asymmetry and conflict
with the tendency of sums of independent variables toward a Gaussian. [4]_
[5]_ ``asymmetric_uncertainty`` does not make that claim about its ordinary
arithmetic.

When propagation of the associated PDFs themselves is required, use
:func:`~asymmetric_uncertainty.propagate`. It draws from each input's
split-normal distribution and evaluates the requested transformation on the
samples directly. The output is left as the sampled result rather than being
forced back into another split-normal ``a_u`` object. Repeated independent sums
therefore exhibit the expected reduction in asymmetry without requiring the
split-normal family to be closed under convolution.

When is the split-normal interpretation appropriate?
----------------------------------------------------

The associated split-normal PDF is a reasonable model when all of the following
are approximately true:

* the nominal value can reasonably be interpreted as a mode or local best-fit
  location;
* the positive and negative uncertainties can reasonably be interpreted as
  local Gaussian scale parameters on the two sides;
* the uncertainty is unimodal and approximately Gaussian in each direction
  near the central value;
* detailed tail behavior is not important to the conclusion; and
* a three-parameter approximation is appropriate for the information actually
  available.

Greater care is required when the reported numbers have a different explicit
meaning. Examples include a median with 16th/84th-percentile errors, a profile-
likelihood interval, a frequentist confidence interval with specified coverage,
a hard bound, a systematic envelope, or a multimodal or strongly heavy-tailed
posterior. Such values can still be represented numerically by ``a_u`` when its
arithmetic semantics are useful, but invoking the distribution utilities then
adds a split-normal modeling assumption that was not implied by the original
notation alone.

For strongly asymmetric problems, tail-sensitive calculations, or cases where
the full likelihood or posterior is available, the full distribution should be
preferred over any reconstruction from only three numbers. The existence of
multiple plausible asymmetric models is itself evidence that three reported
numbers cannot uniquely recover all of the lost distributional information.

Summary
-------

``asymmetric_uncertainty`` uses a specific and established statistical model:
a continuous two-piece normal distribution whose central value is the mode and
whose lower and upper uncertainties are the two Gaussian scales. This choice
preserves continuity at the mode and an exact mode-centered Gaussian
``n-sigma`` probability interpretation.

Other conventions, notably the dimidiated Gaussian, make different choices:
they preserve equal probability on either side and a quantile interpretation
of the two errors, at the cost of a discontinuous PDF when the uncertainties
are asymmetric. Neither convention is implied uniquely by asymmetric-error
notation, and neither family is generally preserved by convolution or nonlinear
transformation.

The package therefore separates the two tasks. ``a_u`` arithmetic propagates
nominal values and local directional uncertainties; distribution utilities give
those values a clearly defined split-normal interpretation; and
:func:`~asymmetric_uncertainty.propagate` is available when the probability
distributions themselves should be propagated numerically.

References
----------

.. [1] G. T. Fechner (1897), *Kollektivmasslehre*, edited by G. F. Lipps,
       Leipzig: Wilhelm Engelmann.

.. [2] K. F. Wallis (2014), *Statistical Science* **29**\ (1), 106.
       `bibcode:2014arXiv1405.4995W <https://scixplorer.org/abs/2014arXiv1405.4995W/abstract>`_.

.. [3] A. Possolo, C. Merkatas, and O. Bodnar (2019), *Metrologia* **56**\ (4),
       045009. `bibcode:2019Metro..56d5009P
       <https://scixplorer.org/abs/2019Metro..56d5009P/abstract>`_.

.. [4] R. Barlow (2003), arXiv:physics/0306138. `bibcode:2003physics...6138B
       <https://scixplorer.org/abs/2003physics...6138B/abstract>`_.

.. [5] R. Barlow, A. R. Brazzale, and I. Volobouev (2026), *Nuclear Instruments
       and Methods in Physics Research A* **1082**, 170857.
       `bibcode:2026NIMPA108270857B
       <https://scixplorer.org/abs/2026NIMPA108270857B/abstract>`_.
