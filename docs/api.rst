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

``UncertaintyArray``
--------------------

.. autoclass:: UncertaintyArray
   :members:
   :undoc-members:

Convenience functions
---------------------

.. autofunction:: pos_errors

.. autofunction:: neg_errors
