Asymmetric Uncertainty
======================

``asymmetric_uncertainty`` provides objects for representing and propagating
numeric quantities with asymmetric uncertainties. Its primary interface is the
:class:`asymmetric_uncertainty.a_u` class.

Getting started
---------------

Installation
~~~~~~~~~~~~

Install the latest version directly from GitHub with pip:

.. code-block:: console

   pip install git+https://github.com/cgobat/asymmetric_uncertainty.git

Alternatively, clone the repository and install it locally:

.. code-block:: console

   git clone https://github.com/cgobat/asymmetric_uncertainty
   cd asymmetric_uncertainty
   pip install .

Basic usage
~~~~~~~~~~~

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

In the case of symmetric errors, the Unicode plus-minus symbol (``±``) is also
accepted; for example, ``a_u("3.0 ± 0.3")``.

Arithmetic applies the ordinary operation to the nominal values and propagates
the positive and negative uncertainties locally with first-order derivatives
and quadrature. Distribution methods such as ``pdf()``, ``cdf()``, and
``rvs()`` additionally interpret the same three stored values as the mode and
two Gaussian widths of a split-normal distribution.

More examples are available in the `example Jupyter notebook
<https://github.com/cgobat/asymmetric_uncertainty/blob/main/example.ipynb>`_.

See the :doc:`api` for the complete public interface and
:doc:`special_features` for additional behavior.

.. toctree::
   :hidden:

   Home <self>
   special_features
   theory
   related_software
   api

Documentation
-------------

* :doc:`Special features <special_features>` — uncertainty propagation,
  distribution utilities, formatting, and edge cases.
* :doc:`Statistical model and interpretation <theory>` — the split-normal
  model, its assumptions, and alternative statistical conventions.
* :doc:`Related software <related_software>` — comparison with other
  uncertainty-handling packages.
* :doc:`API reference <api>` — classes, functions, and attributes.

Project links and citations
---------------------------

This package is free software, and its source code is available on `GitHub
<https://github.com/cgobat/asymmetric_uncertainty/>`_. If you use this package
in published work, please cite its `ASCL entry <https://ascl.net/2208.005>`_
(`bibcode:2022ascl.soft08005G
<https://scixplorer.org/abs/2022ascl.soft08005G/abstract>`_).

Applications
~~~~~~~~~~~~

``asymmetric_uncertainty`` was originally developed for and subsequently used
in the following projects:

* C. Gobat et al. (2023), *MNRAS* **523**\ (1), 775.
  `bibcode:2023MNRAS.523..775G <https://scixplorer.org/abs/2023MNRAS.523..775G/abstract>`_.

  * GitHub repo: `cgobat/dark-GRBs <https://github.com/cgobat/dark-GRBs>`_

* C. Gobat et al. (2022), *RNAAS* **6**\ (8), 163. `bibcode:2022RNAAS...6..163G
  <https://scixplorer.org/abs/2022RNAAS...6..163G/abstract>`_.

  * GitHub repo: `cgobat/XDBS <https://github.com/cgobat/XDBS>`_

A (potentially non-exhaustive) list of other works that have used this software
can be retrieved from SciX/ADS:

* https://scixplorer.org/abs/2022ascl.soft08005G/citations
