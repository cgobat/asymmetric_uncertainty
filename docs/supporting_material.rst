Further reading and applications
================================

The :doc:`statistical model and interpretation <theory>` page defines the
split-normal distribution used by the package and explains its assumptions.
This page collects the literature behind that choice, broader approaches to
asymmetric uncertainties, and examples of the package's use. For comparisons
with other implementations, see :doc:`related_software`.

* The probability distribution used by this package to model asymmetrically
  uncertain quantities is given by `Eq. (2) <https://iopscience.iop.org/article/10.1086/521975/fulltext/71979.text.html#sc5.1>`_
  of Starling et al. (2008). [1]_ The methodology presented there partially
  inspired the creation of this package.

* A dedicated statistical treatment of the same three-parameter two-piece
  normal family is given by John (1982), [2]_ including maximum-likelihood and
  method-of-moments fitting and a test of symmetry and normality.

* Broader treatments of asymmetric errors and their statistical combination
  are given by Barlow (2003) [3]_ and Barlow (2004). [4]_ These approaches are
  useful context for the problem of asymmetric errors, but this package's
  ordinary arithmetic instead preserves nominal-value operations and uses
  first-order propagation of the quoted upper and lower uncertainties.

* See :doc:`related_software` for a comparison with other packages that model
  or propagate asymmetric uncertainties, and :doc:`theory` for a more detailed
  discussion of the mathematical theory surrounding various interpretations of
  such quantities.

* This package was originally developed for and subsequently used in
  `cgobat/dark-GRBs <https://github.com/cgobat/dark-GRBs>`_ [5]_ and
  `cgobat/XDBS <https://github.com/cgobat/XDBS>`_. [6]_

.. rubric:: References

.. [1] R. L. C. Starling et al. (2008), *ApJ* **672**\ (1), 433.
       `bibcode:2008ApJ...672..433S <https://scixplorer.org/abs/2008ApJ...672..433S/abstract>`_.

.. [2] S. John (1982), *Comm. in Stat. - Theory & Methods* **11**\ (8), 879.
       `doi:10.1080/03610928208828279 <https://doi.org/10.1080/03610928208828279>`_.

.. [3] R. Barlow (2003), arXiv:physics/0306138. `bibcode:2003physics...6138B
       <https://scixplorer.org/abs/2003physics...6138B/abstract>`_.

.. [4] R. Barlow (2004), arXiv:physics/0406120. `bibcode:2004physics...6120B
       <https://scixplorer.org/abs/2004physics...6120B/abstract>`_.

.. [5] C. Gobat et al. (2023), *MNRAS* **523**\ (1), 775.
       `bibcode:2023MNRAS.523..775G <https://scixplorer.org/abs/2023MNRAS.523..775G/abstract>`_.

.. [6] C. Gobat et al. (2022), *RNAAS* **6**\ (8), 163.
       `bibcode:2022RNAAS...6..163G <https://scixplorer.org/abs/2022RNAAS...6..163G/abstract>`_.
