API reference
=============

.. currentmodule:: asymmetric_uncertainty

``a_u``
-------

.. autoclass:: a_u
   :members:
   :undoc-members:

The principal scalar attributes are:

``value``
   Nominal value of the represented quantity.

``plus``
   Positive-side uncertainty.

``minus``
   Negative-side uncertainty.

``maximum``
   ``value + plus``.

``minimum``
   ``value - minus``.

``sign``
   ``1`` when ``value`` is positive or zero, otherwise ``-1``.

``is_symmetric``
   ``True`` when ``plus`` and ``minus`` are numerically equal within
   :func:`numpy.isclose` tolerance.

``mean``
   Mean of the represented split-normal distribution. For asymmetric errors,
   this generally differs from ``value``, which is the mode.

``variance``
   Variance of the represented split-normal distribution.

``std``
   Standard deviation of the represented split-normal distribution.

``skewness``
   Standardized third central moment of the represented split-normal
   distribution. This is ``nan`` for an exact, zero-variance value.

``UncertaintyArray``
--------------------

.. autoclass:: UncertaintyArray
   :members:
   :undoc-members:

Convenience functions
---------------------

.. autofunction:: propagate

.. autofunction:: pos_errors

.. autofunction:: neg_errors
