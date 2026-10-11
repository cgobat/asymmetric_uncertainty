Special features
================

In addition to the methods and attributes listed in the :doc:`api`, the
:class:`asymmetric_uncertainty.a_u` class provides several convenience
features.

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

   3.0_{-0.4}^{+0.2},

while a symmetric uncertainty:

.. code-block:: python

   B = a_u(3.0, 0.3, 0.3)
   B

is displayed as

.. math::

   3.0 \pm 0.3.

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
