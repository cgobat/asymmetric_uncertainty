Getting started
===============

Installation
------------

Install the latest version directly from GitHub with pip:

.. code-block:: console

   pip install git+https://github.com/cgobat/asymmetric_uncertainty.git

Alternatively, clone the repository and install it locally:

.. code-block:: console

   git clone https://github.com/cgobat/asymmetric_uncertainty
   cd asymmetric_uncertainty
   pip install .

Basic usage
-----------

The package's main utility is the :class:`asymmetric_uncertainty.a_u` class,
which represents and manipulates numbers with asymmetric errors:

.. code-block:: python

   from asymmetric_uncertainty import a_u

Create an instance as ``a_u(nominal, pos_err, neg_err)``, where ``nominal`` is
the central value and ``pos_err`` and ``neg_err`` are the positive- and
negative-side uncertainties:

.. code-block:: python

   x = a_u(3.0, 0.2, 0.4)
   print(x + 1)

.. code-block:: text

   4.0 (+0.2, -0.4)

Instances can also be created from strings in either asymmetric form:

.. code-block:: python

   x = a_u("3.0 (+0.2, -0.4)")

or symmetric form:

.. code-block:: python

   x = a_u("3.0 +/- 0.3")

The symmetric string form accepted by the package uses the Unicode plus-minus
symbol (``±``); for example, ``a_u("3.0 ± 0.3")``.

Arithmetic applies the ordinary operation to the nominal values and propagates
the positive and negative uncertainties locally with first-order derivatives
and quadrature. Distribution methods such as ``pdf()``, ``cdf()``, and
``rvs()`` additionally interpret the same three stored values as the mode and
two Gaussian widths of a split-normal distribution.

See the :doc:`api` for the complete public interface and
:doc:`special_features` for additional behavior.
