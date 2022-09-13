# Representing upper/lower limits

Below are some possibilities for different ways one might represent a number $x$ that is an upper limit or lower limit.

| Upper limit        | Lower limit        | Commentary |
|--------------------|--------------------|---------|
| $x_{-\infty}^{+0}$ | $x_{-0}^{+\infty}$ | Infinities are not always well-behaved |
| $0_{-0}^{+x}$ (if $x$>0) | $0_{-\|x\|}^{+0}$ (if $x$<0) | Be careful of divide-by-zero errors; makes sense only if quantity has a true limit at zero (for instance, we expect fluxes to be strictly ≥0 physically) |
| $x_{-x}^{+0}$ (if $x$>0) | $x_{-0}^{+\|x\|}$ (if $x$<0) | Poor representation of a true limit if errors are on the order of the nominal value anyway. |