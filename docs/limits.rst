Representing upper and lower limits
===================================

The table below records several possible conventions for representing a number
:math:`x` that is an upper or lower limit. These are possibilities rather than
a single prescribed representation.

.. list-table::
   :header-rows: 1
   :widths: 24 24 52

   * - Upper limit
     - Lower limit
     - Commentary
   * - :math:`x_{-\infty}^{+0}`
     - :math:`x_{-0}^{+\infty}`
     - Infinities are not always well behaved.
   * - :math:`0_{-0}^{+x}` if :math:`x > 0`
     - :math:`0_{-|x|}^{+0}` if :math:`x < 0`
     - Be careful of divide-by-zero errors. This makes sense only when the
       quantity has a true physical limit at zero, such as a non-negative flux.
   * - :math:`x_{-x}^{+0}` if :math:`x > 0`
     - :math:`x_{-0}^{+|x|}` if :math:`x < 0`
     - This is a poor representation of a true limit when the uncertainty is
       of the same order as the nominal value.
