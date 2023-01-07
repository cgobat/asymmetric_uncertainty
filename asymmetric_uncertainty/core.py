"""Asymmetric Uncertainty: A package for handling non-standard numerical uncertainties."""

__author__ = "Caden Gobat"
__author_affiliation__ = ["George Washington University",
                          "Southwest Research Institute"]
__contact__ = "<cgobat@gwu.edu>"
__deprecated__ = False
__version__ = "0.2.2"

import warnings
from numbers import Number
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


class a_u(u.Quantity):
    """
    Class for representing and handling propagation of asymmetric uncertainties assuming a pseudo-
    Gaussian probability distribution where the errors on either side of the nominal value are like
    modified 1-sigma standard deviations.

    Parameters
    ----------
    nominal : numeric
        the nominal value of the represented quantity
    pos_err : numeric
        the plus error on the value
    neg_err : numeric
        the minus error on the value
    unit    : astropy.units.Unit
        physical unit with which to initialize the quantity

    Attributes
    ----------
    value : numeric
        the nominal value of the represented quantity
    plus  : numeric
        the positive error on the value
    minus : numeric
        the negative error on the value
    unit  : astropy.units.Unit
        the physical units of the quantity
    """
    value: Number
    plus : Number
    minus: Number
    unit : u.Unit
    
    def __new__(cls, nominal: "Number|u.Quantity",
                     pos_err: "Number|u.Quantity"=0.,
                     neg_err: "Number|u.Quantity"=0.,
                     unit   : "str|u.UnitBase"=u.dimensionless_unscaled,
                     **quantity_kwargs) -> "a_u":
        
        for check_val in (nominal, pos_err, neg_err):
            if isinstance(check_val, u.Quantity):
                assert check_val.isscalar, "Input values must be scalars, not arrays."
        
        if unit == u.dimensionless_unscaled: # unit unspecified
            for check_val in (nominal, pos_err, neg_err):
                if hasattr(check_val, "unit"): # does input already have units?
                    unit = check_val.unit
                    break
                else: # if it doesn't, move on
                    pass
        else: # unit is specified outright,
            pass # so just use that
        
        if hasattr(nominal, "unit"):
            nominal = nominal.to(unit).value
        if hasattr(pos_err, "unit"):
            pos_err = pos_err.to(unit).value
        if hasattr(neg_err, "unit"):
            neg_err = neg_err.to(unit).value
        
        obj: cls = super().__new__(cls, value=nominal, unit=unit, **quantity_kwargs)
        # astropy Quantity initializes self.value and self.unit
        try:
            obj._initialize_err(pos_err, neg_err)
        except:
            warnings.warn("Failed to initialize errors.", RuntimeWarning)
        return obj

    # **************** INSTANCE PROPERTIES ****************
    @property
    def maximum(self) -> Number:
        return (self.value + self.plus)*self.unit
    
    @property
    def minimum(self) -> u.Quantity:
        return (self.value - self.minus)*self.unit
    
    @property
    def issymmetric(self) -> bool:
        if self.isscalar:
            return np.isclose(self.plus, self.minus)
        else:
            return [np.isclose(pos, neg) for pos, neg in zip(self.plus, self.minus)]
    
    # **************** ERROR INITIALIZATION/MANIPULATION ****************
    def _initialize_err(self, pos_err: Number, neg_err: Number) -> None:
        self.__dict__["value"] = self.value
        if self.isscalar:
            self.plus = abs(float(pos_err))
            self.minus = abs(float(neg_err))
            self._sign = 1. if self.value > 0 else -1. if self.value < 0 else 0.
        else:
            self.plus = np.array([abs(float(x)) for x in pos_err])
            self.minus = np.array([abs(float(x)) for x in neg_err])
            self._sign = np.array([1. if x > 0 else -1. if x < 0 else 0. for x in self.value])

    def add_error(self, delta, how="quadrature", inplace=False):
        """
        Adds `delta` to an instance's existing error. Possible `how`s are `quadrature`, `straight`, or `split`.
        If `inplace` is `True`, the existing object's errors are modified in place. If it is `False`, a new instance is returned.
        """
        if how=="quadrature":
            new_pos = np.sqrt(self.plus**2 + delta**2)
            new_neg = np.sqrt(self.minus**2 + delta**2)
        elif how=="straight":
            new_pos = self.plus + delta
            new_neg = self.minus + delta
        elif how=="split":
            new_pos = self.plus + delta/2
            new_neg = self.minus + delta/2
        else:
            raise ValueError(f"'how' should be one of {{'quadrature', 'straight', 'split'}}")
        if inplace:
            self.plus = new_pos
            self.minus = new_neg
        else:
            return a_u(self.value, new_pos, new_neg, unit=self.unit)
        
    # **************** OUTPUT FORMATTING/DISPLAY MAGIC METHODS ****************
    def __str__(self) -> str:
        return self.__repr__()
            
    def __repr__(self) -> str:
        if self.isscalar:
            if np.isclose(self.plus, self.minus):
                str_repr = f"<{self.value} ± {self.plus}{self._unitstr}>"
            else:
                str_repr = f"<{self.value} (+{self.plus}, -{self.minus}){self._unitstr}>"
            return str_repr
        else:
            return "[" + ", ".join([element.__repr__() for element in self]) + "]"
    
    def _repr_latex_(self) -> str:
        if self.isscalar:
            val_str = format(self.value, "f")
            while val_str.endswith("0"):
                if val_str.endswith(".0"):
                    break
                val_str = val_str[:-1]
            pos_str = format(self.plus, "f")
            while pos_str.endswith("0"):
                if pos_str.endswith(".0"):
                    break
                pos_str = pos_str[:-1]
            neg_str = format(self.minus, "f")
            while neg_str.endswith("0"):
                if neg_str.endswith(".0"):
                    break
                neg_str = neg_str[:-1]
            if np.isclose(self.plus, self.minus):
                str_repr = f"${val_str} \pm {pos_str}$"
            else:
                str_repr = f"${val_str}_{{-{neg_str}}}^{{+{pos_str}}}$"
            if self._unitstr:
                str_repr += " " + self.unit._repr_latex_()
            return str_repr
        else:
            return "[" + ", ".join([element._repr_latex_() for element in self]) + "]"

    def __format__(self, format_spec) -> str:
        try:
            val = format(self.value, format_spec)
            pos = format(self.plus, format_spec)
            neg = format(self.minus, format_spec)
            sym_flag = (pos == neg)
            overall_fmt = "s"
        except ValueError:
            val = self.value
            pos = self.plus
            neg = self.minus
            sym_flag = self.issymmetric
            overall_fmt = format_spec
        if sym_flag:
            return format(f"{val} ± {pos}{self._unitstr}", overall_fmt)
        else:
            return format(f"{val} (+{pos}, -{neg}){self._unitstr}", overall_fmt)
    
    # **************** MISCELLANEOUS ****************
    def __reduce__(self):
        return NotImplemented
    
    def items(self, units=True) -> "tuple[float]|tuple[u.Quantity]":
        """
        Returns a tuple of `(value, plus, minus)`, either with units or without.
        """
        if units:
            return (self.value*self.unit, self.plus*self.unit, self.minus*self.unit)
        else:
            return (self.value, self.plus, self.minus)
    
    # **************** CONVERSIONS AND CASTING METHODS ****************
    def as_Quantity(self) -> u.Quantity:
        return self.value * self.unit

    def to(self, target: u.UnitBase) -> "a_u":
        '''Unit conversion method'''
        converted = [q.to(target) for q in self.items(units=True)]
        return a_u(*converted)

    def __int__(self) -> int:
        return int(self.as_Quantity())
    
    def __float__(self) -> float:
        return float(self.as_Quantity())
    
    def __neg__(self):
        return a_u(-self.value, self.minus, self.plus, self.unit)

    def __abs__(self):
        '''Absolute value'''
        if self < 0:
            return -self
        else:
            return self
    
    # **************** ARITHMETIC OPERATORS AND MATH OPERATIONS ****************
    def __array_ufunc__(self, function, method, *inputs, **kwargs):
        fname = function.__name__
        if fname in ("sin", "cos", "tan"):
            if self.unit == u.dimensionless_unscaled:
                x = self.value
            elif self.unit in u.radian.find_equivalent_units():
                x = self.as_Quantity()
            else:
                raise u.UnitTypeError(f"Inputs to trigonometric functions must either be angles or dimensionless, not '{self.unit}'.")
            result = function(x)
            deriv = (np.cos if fname=="sin" else
                     lambda x: -np.sin(x) if fname=="cos" else
                     lambda x: 1/(np.cos(x)**2) if fname=="tan" else
                     lambda  : None)
            pos = self.plus*deriv(x)
            neg = self.minus*deriv(x)
            return self.__class__(result, pos, neg)
        elif fname in ("log", "log2", "log10", "exp", "sqrt"):
            return getattr(self, fname)()
        elif fname == "sign":
            return self._sign
        elif fname == "absolute":
            return abs(self)
        elif fname == "multiply":
            return [self.__mul__(x) for x in inputs[-1]]
        else:
            warnings.warn(f"ufunc {fname!r} is not (yet) implemented for quantities with asymmetric uncertainties.", Warning)
            return NotImplemented

    def __add__(self, other): # self + other
        if isinstance(other, type(self)):
            other = other.to(self.unit)
            result = self.as_Quantity() + other.as_Quantity()
            pos = np.sqrt(self.plus**2 + other.plus**2)
            neg = np.sqrt(self.minus**2 + other.minus**2)
            #print("added", self, "+", other, "=", a_u(result, pos, neg))
            return a_u(result, pos, neg, unit=result.unit)
        elif isinstance(other, u.Quantity):
            result = self.as_Quantity() + other
            pos = (self.plus*self.unit).to(result.unit)
            neg = (self.minus*self.unit).to(result.unit)
            return a_u(result.value, pos.value, neg.value, result.unit)
        elif isinstance(other, Number):
            result: u.Quantity = self.as_Quantity() + other
            return a_u(result.value, self.plus, self.minus, result.unit)
        else:
            return NotImplemented
    
    def __radd__(self, other): # other + self
        return self + other # addition is commutative
    
    def __sub__(self, other): # self - other
        if isinstance(other, type(self)):
            result = self.as_Quantity() - other.as_Quantity()
            pos = np.sqrt((self.plus*self.unit)**2 + (other.minus*other.unit)**2)
            neg = np.sqrt((self.minus*self.unit)**2 + (other.plus*other.unit)**2)
            #print("subtracted", other, "from", self, "=", a_u(result, pos, neg))
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        elif isinstance(other, (u.Quantity, Number)):
            result = self.as_Quantity() - other
            pos = (self.plus*self.unit).to(result.unit)
            neg = (self.minus*self.unit).to(result.unit)
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        else:
            return NotImplemented
    
    def __rsub__(self, other): # other - self
        return -(self - other) # leverage existing __sub__() and __neg__() implementations
    
    def __mul__(self, other): # self * other
        if isinstance(other, type(self)):
            result = self.as_Quantity() * other.as_Quantity()
            pos = np.sqrt((self.plus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
            neg = np.sqrt((self.minus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
            return a_u(result.value, pos, neg, unit=result.unit)
        elif isinstance(other, (u.Quantity, Number)):
            result = self.as_Quantity() * other
            pos = (self.plus/self.value) * np.abs(result)
            neg = (self.minus/self.value) * np.abs(result)
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        elif isinstance(other, u.UnitBase):
            return a_u(self.value, self.plus, self.minus, unit=self.unit*other)
        else:
            return NotImplemented
    
    def __rmul__(self, other): # other * self
        return self*other # multiplication is commutative
    
    def __truediv__(self, other): # self divided by other
        if isinstance(other, type(self)):
            result = self.as_Quantity() / other.as_Quantity()
            pos = np.sqrt((self.plus  / self.value)**2 + \
                          (other.minus/other.value)**2) * np.abs(result)
            neg = np.sqrt((self.minus / self.value)**2 + \
                          (other.plus /other.value)**2) * np.abs(result)
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        elif isinstance(other, (u.Quantity, Number)):
            result = self.as_Quantity() / other
            pos = (self.plus/self.value) * np.abs(result)
            neg = (self.minus/self.value) * np.abs(result)
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        elif isinstance(other, u.UnitBase):
            return self * (1/other)
        else:
            return NotImplemented
    
    def __rtruediv__(self, other): # other divided by self
        # no need to check if type(other) is a_u because that case would have invoked __truediv__()
        if isinstance(other, (u.Quantity, Number)):
            result = other / self.as_Quantity()
            pos = (self.minus/self.value) * np.abs(result)
            neg = (self.plus/self.value) * np.abs(result)
            return a_u(result.value, pos.value, neg.value, unit=result.unit)
        elif isinstance(other, u.UnitBase):
            return other * (1/self)
        else:
            return NotImplemented
    
    def __pow__(self, other): # self to the something power
        if isinstance(other, type(self)):
            pass
        elif isinstance(other, u.Quantity):
            other = a_u(other, 0, 0, unit=other.unit)
        elif isinstance(other, Number):
            other = a_u(other, 0, 0, unit=self.unit)
        else:
            return NotImplemented
        result = self.as_Quantity()**other.as_Quantity()
        pos = np.abs(result)*np.sqrt((self.plus*other.value/self.value)**2 + (other.plus*np.log(self.value))**2)
        neg = np.abs(result)*np.sqrt((self.minus*other.value/self.value)**2 + (other.minus*np.log(self.value))**2)
        #print("raised", self, "to", other, "=", a_u(result, pos, neg))        
        return a_u(result, pos, neg)
    
    def __rpow__(self, other): # something to the self power
        if isinstance(other, type(self)):
            pass
        elif isinstance(other, u.Quantity):
            other = a_u(other, 0, 0, unit=other.unit)
        elif isinstance(other, Number):
            other = a_u(other, 0, 0, unit=self.unit)
        else:
            return NotImplemented
        result = other.as_Quantity()**self.as_Quantity()
        pos = np.abs(result)*np.sqrt((other.plus*self.value/other.value)**2 + (self.plus*np.log(other.value))**2)
        neg = np.abs(result)*np.sqrt((other.minus*self.value/other.value)**2 + (self.minus*np.log(other.value))**2)
        #print("raised", other, "to", self, "=", a_u(result, pos, neg))                
        return a_u(result.value, pos, neg, unit=result.unit)
    
    def log10(self):
        result = np.log10(self.as_Quantity())
        pos = (self.plus*self.unit)/(self.as_Quantity()*np.log(10))
        neg = (self.minus*self.unit)/(self.as_Quantity()*np.log(10))
        #print("logged", self, "=", a_u(result, pos, neg))
        return a_u(result, pos, neg)

    def log2(self):
        result = np.log2(self.as_Quantity())
        pos = (self.plus*self.unit)/(self.as_Quantity()*np.log(10))
        neg = (self.minus*self.unit)/(self.as_Quantity()*np.log(10))
        return a_u(result, pos, neg)
    
    def log(self):
        result = np.log(self.as_Quantity())
        pos = (self.plus*self.unit)/self.as_Quantity()
        neg = (self.minus*self.unit)/self.as_Quantity()
        return a_u(result, pos, neg)

    def exp(self):
        e = float(np.exp(1))
        return e**self

    def sqrt(self):
        return self**0.5

    # **************** COMPARISON OPERATIONS AND CHECKS ****************
    def __eq__(self, other) -> bool:
        if isinstance(other, type(self)):
            val_match = (self.as_Quantity() == other.as_Quantity())
            pos_match = (self.plus*self.unit == other.plus*other.unit)
            neg_match = (self.minus*self.unit == other.minus*other.unit)
            return val_match and pos_match and neg_match
        elif isinstance(other, (u.Quantity, Number)):
            return (self.as_Quantity() == other) and np.isclose(self.plus, 0.) and np.isclose(self.minus, 0.)
        else:
            return False
    
    def __gt__(self, other) -> bool:
        if isinstance(other, type(self)):
            return self.as_Quantity() > other.as_Quantity()
        elif isinstance(other, (u.Quantity, Number)):
            return self.as_Quantity() > other
        else:
            return False
    
    def __lt__(self, other) -> bool:
        if isinstance(other, type(self)):
            return self.as_Quantity() < other.as_Quantity()
        elif isinstance(other, (u.Quantity, Number)):
            return self.as_Quantity() < other
        else:
            return False
    
    def __lshift__(self, other) -> bool: # overloaded <<; definitively less than
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        return self.maximum < other.minimum
    
    def __rshift__(self, other) -> bool: # overloaded >>; definitively greater than
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        return self.minimum > other.maximum
    
    def __le__(self, other) -> bool:
        return not (self > other)
    
    def __ge__(self, other) -> bool:
        return not (self < other)
    
    def __isfinite__(self) -> bool:
        return all(np.isfinite(self.items()))
    
    def isna(self) -> bool:
        """
        `pandas`-style NaN checker. Returns True if value is NaN or None, and False if neither.        
        """
        return (np.isnan(self.as_Quantity()) or (self.value is None))
    
    def notna(self) -> bool:
        """
        Inverse of `isna()`. Returns True if value is neither NaN nor None, and False if it is.
        """
        return not self.isna()

    # **************** DISTRIBUTION FUNCTIONS AND PLOTTING ****************
    def pdf(self, x) -> np.ndarray:
        """
        Computes and returns the values of the probability distribution function for the specified input.
        """
        return np.piecewise(x, [x<self.value, x>=self.value],
                            [lambda x : np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * \
                                        np.exp(-1*(x-self.value)**2 / (2*self.minus**2)),
                             lambda x : np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * \
                                        np.exp(-1*(x-self.value)**2 / (2*self.plus**2))])
    
    def cdf(self, x) -> np.ndarray:
        """
        Computes and returns the values of the cumulative distribution function for the specified input.
        """
        return np.cumsum(self.pdf(x))/np.sum(self.pdf(x))
    
    def pdfplot(self, num_sigma=5, discretization=100, **kwargs) -> None:
        """
        Plots the associated PDF over the specified number of sigma, using 2*`discretization` points.
        `**kwargs` are passed on to `matplotlib` for configuration of the resulting plot.
        """
        neg_x = np.linspace(self.value-(num_sigma*self.minus), self.value, discretization)
        pos_x = np.linspace(self.value, self.value+(num_sigma*self.minus), discretization)
        x = np.hstack([neg_x, pos_x])
        pdf = self.pdf(x)
        plt.plot(x, pdf, **kwargs)
        return None
    
    def cdfplot(self, num_sigma=5, discretization=100, **kwargs) -> None:
        """
        Plots the associated CDF over the specified number of sigma, using 2*`discretization` points.
        `**kwargs` are passed on to `matplotlib` for configuration of the resulting plot.
        """
        neg_x = np.linspace(self.value-(num_sigma*self.minus), self.value, discretization)
        pos_x = np.linspace(self.value, self.value+(num_sigma*self.minus), discretization)
        x = np.hstack([neg_x, pos_x])
        pdf = self.pdf(x)
        cdf = np.cumsum(pdf)/np.sum(pdf)
        plt.plot(x, cdf, **kwargs)
        return None


def AsymmetricUncertainty(*args, **kwargs) -> a_u:
    '''Factory function for legacy namespace support.'''
    warnings.warn("AsymmetricUncertainty has been renamed. Use a_u instead.",
                  DeprecationWarning, stacklevel=2)
    return a_u(*args, **kwargs)


class UncertaintyArray(list):
    """
    Class for representing an array of `a_u` objects.
    Mostly provides utilities for slicing and accessing various attributes.
    """
    def refresh(self):
        for i in range(len(self)):
            try:
                self[i].value
                self[i].plus
                self[i].minus
            except AttributeError:
                self[i] = a_u(self[i], 0, 0)
                
        self.as_numpy = np.array(self.as_list)
        self.flattened = self.as_numpy.flatten()
        self.shape = self.as_numpy.shape
        self.ndim = self.as_numpy.ndim

        self.minus = [v.minus for v in self.as_list]
        self.plus = [v.plus for v in self.as_list]
        self.values = [v.value for v in self.as_list]
    
    def __init__(self, array=[]):

        self.as_list = list(array)
        self.refresh()

    def __len__(self):
        return len(self.as_list)

    def __iter__(self):
        return iter(self.as_list)

    def __getitem__(self, key):
        return self.as_list[key]

    def __setitem__(self, key, val):
        self.as_list[key] = val

    def __str__(self):
        return str([str(entry) for entry in self.as_list])

    def __repr__(self):
        return str(self)

    def __contains__(self, item):
        return item in self.as_list

    def append(self, entry):
        self.as_list.append(entry)
        self.refresh()
    
    def pdf(self, x):
        return np.sum([entry.pdf(x) for entry in self], axis=0)

def pos_errors(array):
    """
    Stand-alone function to return an array of the positive errors of an array of `a_u` objects.
    Functional equivalent to `UncertaintyArray(array).plus`.
    """
    return [v.plus for v in array]

def neg_errors(array):
    """
    Stand-alone function to return an array of the negative errors of an array of `a_u` objects.
    Functional equivalent to `UncertaintyArray(array).minus`.
    """
    return [v.minus for v in array]