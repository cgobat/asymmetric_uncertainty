# References and further reading

- The probability distribution used by this package to model asymmetrically uncertain quantities is given in Eq. (2) of Starling, et al. (2008).[^1] The methodology presented in this paper partially inspired the creation of this package in the first place.
- [@muryelgp](https://github.com/muryelgp)'s [`asymmetric_uncertainties`](https://github.com/muryelgp/asymmetric_uncertainties) package, which (as the name suggests) has some similar functionality to this, but uses [Monte Carlo methods](https://en.wikipedia.org/wiki/Monte_Carlo_method) rather than an empirical/analytical function to model the error distributions. The mathematical basis behind this alternative method is described in Barlow (2004).[^2]
- This package has been used in two of my own works: [cgobat/dark-GRBs](https://github.com/cgobat/dark-GRBs)[^3] and [cgobat/XDBS](https://github.com/cgobat/XDBS).[^4]


[^1]: R.L.C. Starling, *et al.* (2008). *ApJ* **672**(1), 433. [bibcode:2008ApJ...672..433S](https://ui.adsabs.harvard.edu/abs/2008ApJ...672..433S/abstract).
[^2]: R. Barlow (2004). arXiv:physics/0406120 [bibcode:2004physics...6120B](https://ui.adsabs.harvard.edu/abs/2004physics...6120B/abstract)
[^3]: C. Gobat, *et al.* (2023). *MNRAS* **523**(1), 775. [bibcode:2023MNRAS.523..775G](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523..775G/abstract)
[^4]: C. Gobat, *et al.* (2022). *RNAAS* **6**(8), 163. [bibcode:2022RNAAS...6..163G](https://ui.adsabs.harvard.edu/abs/2022RNAAS...6..163G/abstract)