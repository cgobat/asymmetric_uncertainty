"""Asymmetric Uncertainty: A package for handling non-standard numerical uncertainties.

    Copyright (C) 2022 Caden Gobat

    This program is free software: you can redistribute and/or modify it under
    the terms of the GNU General Public License (either version 3 or, at your
    option, any later version).

    This program is distributed in the hope that it will be useful, but comes
    with ABSOLUTELY NO WARRANTY. See <https://www.gnu.org/licenses/> or the
    LICENSE file distributed with this program for more details.
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings
from math import erf
from numbers import Number
try:
    from astropy import units as u
except ImportError:  # Astropy is installed by the package, but keep source-tree use graceful.
    u = None


np_erf = np.vectorize(erf, doc="Vectorized error function at x.")  # prevents need for SciPy dependency

def _has_astropy():
    return u is not None

def _is_quantity(value):
    return _has_astropy() and isinstance(value, u.Quantity)

def _is_unit(value):
    return _has_astropy() and isinstance(value, u.UnitBase)

def _unit_is_dimensionless(unit):
    return not _has_astropy() or unit is None or unit == u.dimensionless_unscaled

def _sum_in_quadrature(terms):
    if not terms:
        return 0.0
    total = terms[0]
    for term in terms[1:]:
        if _is_quantity(total) and _is_quantity(term):
            term = term.to(total.unit)
        total = total + term
    return np.sqrt(total)

def _derivative_sign(derivative):
    value = derivative.value if _is_quantity(derivative) else derivative
    return 1 if value >= 0 else -1

def _propagated_errors(operands, derivatives, result_unit=None):
    """Return first-order asymmetric errors for the supplied partial derivatives."""
    pos_terms = []
    neg_terms = []
    for operand, derivative in zip(operands, derivatives):
        if _derivative_sign(derivative) >= 0:
            pos_error = operand.plus_quantity
            neg_error = operand.minus_quantity
        else:
            pos_error = operand.minus_quantity
            neg_error = operand.plus_quantity

        pos_term = derivative*pos_error
        neg_term = derivative*neg_error
        if _is_quantity(pos_term) and result_unit is not None:
            pos_term = pos_term.to(result_unit)
            neg_term = neg_term.to(result_unit)
        pos_terms.append(pos_term**2)
        neg_terms.append(neg_term**2)

    pos = _sum_in_quadrature(pos_terms)
    neg = _sum_in_quadrature(neg_terms)
    if _is_quantity(pos):
        return pos.to_value(result_unit), neg.to_value(result_unit)
    return pos, neg


class a_u:
    """
    Class for representing and handling propagation of asymmetric uncertainties assuming a pseudo-Gaussian
    probability distribution where the plus and minus errors in each direction of the nominal value are like
    modified 1-sigma standard deviations.

    Parameters
    ----------
    numeric or astropy.units.Quantity : nominal
        the nominal value of the represented quantity
    numeric or astropy.units.Quantity : pos_err
        the plus error on the value
    numeric or astropy.units.Quantity : neg_err
        the minus error on the value
    str or astropy.units.UnitBase, optional : unit
        physical unit for the value and its uncertainties. If omitted, the first
        constructor argument with units supplies the unit for the object.

    Attributes
    ----------
    numeric : value
        the nominal value of the represented quantity, expressed in `unit`
    numeric : plus
        the positive error on the value, expressed in `unit`
    numeric : minus
        the negative error on the value, expressed in `unit`
    astropy.units.UnitBase : unit
        the physical unit, or `dimensionless_unscaled` for ordinary numbers
    """

    def __init__(self, nominal, pos_err=0, neg_err=0, unit=None):
        if isinstance(nominal, str):
            stripped = nominal.replace(" ", "")
            if "±" in stripped:
                nominal = float(stripped.split("±")[0])
                pos_err = neg_err = float(stripped.split("±")[1])
            elif "+" in stripped and "-" in stripped:
                nominal = float(stripped.split("(")[0])
                err_str = stripped.split("(")[1].replace(")", "")
                pos_err = float(err_str.split(",")[0][1:])
                neg_err = float(err_str.split(",")[1][1:])
            elif _has_astropy():
                try:
                    nominal = u.Quantity(nominal)
                except (TypeError, ValueError) as exc:
                    raise ValueError("Failed to parse string, likely due to improper formatting.") from exc
            else:
                raise ValueError("Failed to parse string, likely due to improper formatting.")

        if not _has_astropy():
            if unit is not None:
                raise ImportError("Astropy is required to construct an a_u with physical units.")
            self.unit = None
            self.value = float(nominal)
            self.plus = np.abs(float(pos_err))
            self.minus = np.abs(float(neg_err))
            return

        if unit is None:
            target_unit = u.dimensionless_unscaled
            for candidate in (nominal, pos_err, neg_err):
                if _is_quantity(candidate):
                    target_unit = candidate.unit
                    break
        else:
            target_unit = u.Unit(unit)

        self.unit = target_unit
        self.value = self._value_in_unit(nominal, target_unit)
        self.plus = np.abs(self._value_in_unit(pos_err, target_unit))
        self.minus = np.abs(self._value_in_unit(neg_err, target_unit))

    @staticmethod
    def _value_in_unit(value, unit):
        if _is_quantity(value):
            return float(value.to_value(unit))
        return float(value)

    @property
    def quantity(self):
        """Nominal value as an Astropy `Quantity` (if Astropy is available)."""
        if not _has_astropy():
            return self.value
        return self.value*self.unit

    @property
    def plus_quantity(self):
        if not _has_astropy():
            return self.plus
        return self.plus*self.unit

    @property
    def minus_quantity(self):
        if not _has_astropy():
            return self.minus
        return self.minus*self.unit

    def to(self, target, equivalencies=None):
        """Return a copy converted to `target` units."""
        if not _has_astropy():
            raise ImportError("Astropy is required for unit conversion.")
        if equivalencies is None:
            equivalencies = []
        target = u.Unit(target)
        value = self.quantity.to_value(target, equivalencies=equivalencies)
        plus = self.plus_quantity.to_value(target, equivalencies=equivalencies)
        minus = self.minus_quantity.to_value(target, equivalencies=equivalencies)
        return a_u(value, plus, minus, unit=target)

    def to_value(self, target=None, equivalencies=None):
        """Return the nominal value, optionally converted to `target` units."""
        if target is None or not _has_astropy():
            return self.value
        if equivalencies is None:
            equivalencies = []
        return self.quantity.to_value(target, equivalencies=equivalencies)

    @property
    def maximum(self):
        result = self.value + self.plus
        if _unit_is_dimensionless(self.unit):
            return result
        return result*self.unit

    @property
    def minimum(self):
        result = self.value - self.minus
        if _unit_is_dimensionless(self.unit):
            return result
        return result*self.unit

    @property
    def sign(self):
        return 1 if self.value >= 0 else -1

    @property
    def is_symmetric(self):
        return np.isclose(self.plus, self.minus)

    def __str__(self):
        return f"{self}"

    def __format__(self, fmt_spec):
        unit_suffix = "" if _unit_is_dimensionless(self.unit) else f" {self.unit}"
        if np.isclose(self.plus, self.minus):
            return f"{format(self.value, fmt_spec)} ± {format(self.plus, fmt_spec)}{unit_suffix}"
        return "%s (+%s, -%s)%s" % (*(format(x, fmt_spec) for x in self.items()), unit_suffix)

    def _repr_latex_(self, float_fmt="f"):
        if np.isclose(self.plus, self.minus):
            representation = f"${format(self.value, float_fmt)} \\pm {format(self.plus, float_fmt)}"
        else:
            representation = "$%s^{+%s}_{-%s}" % tuple(format(x, float_fmt) for x in self.items())
        if not _unit_is_dimensionless(self.unit):
            unit_latex = self.unit.to_string("latex_inline").strip("$")
            representation += f"\\,{unit_latex}"
        return representation + "$"

    def _x_values(self, x):
        if _is_quantity(x):
            return np.asanyarray(x.to_value(self.unit), dtype=float)
        return np.asanyarray(x, dtype=float)

    def pdf(self, x):
        """
        Computes and returns the values of the probability distribution function for the specified input.

        For values with units, the density returned by this function has inverse units.
        """
        x_arr = self._x_values(x)
        pdf_arr = np.piecewise(x_arr, [x_arr<self.value, x_arr>=self.value],
                               [lambda values: np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(values - self.value)**2 / (2*self.minus**2)),
                                lambda values: np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(values - self.value)**2 / (2*self.plus**2))])
        result = pdf_arr.item() if x_arr.ndim == 0 else pdf_arr
        if _unit_is_dimensionless(self.unit):
            return result
        return result/self.unit

    def cdf(self, x):
        """
        Computes and returns the values of the cumulative distribution function for the specified input.
        """
        x_arr = self._x_values(x)
        width_sum = self.plus + self.minus
        cdf_arr = np.piecewise(x_arr, [x_arr<self.value, x_arr>=self.value],
                               [lambda values: self.minus/width_sum * (1 + np_erf((values-self.value)/(self.minus*np.sqrt(2)))),
                                lambda values: (self.minus + self.plus*np_erf((values-self.value)/(self.plus*np.sqrt(2))))/width_sum])
        return cdf_arr.item() if x_arr.ndim == 0 else cdf_arr

    def pdfplot(self, num_sigma=5, discretization=100, **kwargs):
        """
        Plots the associated PDF over the specified number of sigma, using 2*`discretization` points.
        `**kwargs` are passed on to `matplotlib` for configuration of the resulting plot.
        """
        neg_x = np.linspace(self.value-(num_sigma*self.minus), self.value, discretization)
        pos_x = np.linspace(self.value, self.value+(num_sigma*self.plus), discretization)
        x = np.array(list(neg_x)+list(pos_x))
        pdf = self.pdf(x)
        if _is_quantity(pdf):
            pdf = pdf.value
        plt.plot(x, pdf, **kwargs)
        plt.show()

    def cdfplot(self, num_sigma=5, discretization=100, **kwargs):
        """
        Plots the associated CDF over the specified number of sigma, using 2*`discretization` points.
        `**kwargs` are passed on to `matplotlib` for configuration of the resulting plot.
        """
        neg_x = np.linspace(self.value-(num_sigma*self.minus), self.value, discretization)
        pos_x = np.linspace(self.value, self.value+(num_sigma*self.plus), discretization)
        x = np.array(list(neg_x)+list(pos_x))
        cdf = self.cdf(x)
        plt.plot(x, cdf, **kwargs)
        plt.show()

    def add_error(self, delta, method="quadrature", inplace=False):
        """
        Adds `delta` to an instance's existing error. Possible `method`s are `quadrature`, `straight`, or `split`.
        If `inplace` is `True`, the existing object's errors are modified in place. If it is `False`, a new instance is returned.
        """
        if _is_quantity(delta):
            delta = delta.to_value(self.unit)
        if method=="quadrature":
            new_pos = np.sqrt(self.plus**2 + delta**2)
            new_neg = np.sqrt(self.minus**2 + delta**2)
        elif method=="straight":
            new_pos = self.plus + delta
            new_neg = self.minus + delta
        elif method=="split":
            new_pos = self.plus + delta/2
            new_neg = self.minus + delta/2
        else:
            raise ValueError
        if inplace:
            self.plus = new_pos
            self.minus = new_neg
        else:
            return a_u(self.value, new_pos, new_neg, unit=self.unit)

    def items(self, units=False):
        """
        Returns a tuple of `(value, plus, minus)`. Set `units=True` to return
        Astropy `Quantity` objects.
        """
        if units and _has_astropy():
            return (self.quantity, self.plus_quantity, self.minus_quantity)
        return (self.value, self.plus, self.minus)

    def __int__(self):
        if _unit_is_dimensionless(self.unit):
            return int(self.value)
        return int(self.quantity)

    def __float__(self):
        if _unit_is_dimensionless(self.unit):
            return float(self.value)
        return float(self.quantity)

    def __neg__(self):
        return a_u(-self.value, self.minus, self.plus, unit=self.unit)

    def __abs__(self):
        if self.value < 0:
            return -self
        return a_u(self.value, self.plus, self.minus, unit=self.unit)

    def _coerce_additive_operand(self, other):
        if isinstance(other, type(self)):
            if _has_astropy():
                return other.to(self.unit)
            return other
        if _is_quantity(other):
            return a_u(other).to(self.unit)
        if isinstance(other, Number):
            other = a_u(other, 0, 0)
            if _has_astropy():
                return other.to(self.unit)
            return other
        return NotImplemented

    @staticmethod
    def _coerce_multiplicative_operand(other):
        if isinstance(other, a_u):
            return other
        if _is_quantity(other):
            return a_u(other)
        if isinstance(other, Number):
            return a_u(other, 0, 0)
        return NotImplemented

    def __add__(self, other):
        if _has_astropy() and isinstance(other, Number):
            result_quantity = self.quantity + other
            pos = self.plus_quantity.to_value(result_quantity.unit)
            neg = self.minus_quantity.to_value(result_quantity.unit)
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        other = self._coerce_additive_operand(other)
        if other is NotImplemented:
            return NotImplemented
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        return a_u(result, pos, neg, unit=self.unit)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if _has_astropy() and isinstance(other, Number):
            result_quantity = self.quantity - other
            pos = self.plus_quantity.to_value(result_quantity.unit)
            neg = self.minus_quantity.to_value(result_quantity.unit)
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        other = self._coerce_additive_operand(other)
        if other is NotImplemented:
            return NotImplemented
        result = self.value - other.value
        pos = np.sqrt(self.plus**2 + other.minus**2)
        neg = np.sqrt(self.minus**2 + other.plus**2)
        return a_u(result, pos, neg, unit=self.unit)

    def __rsub__(self, other):
        if _has_astropy() and isinstance(other, Number):
            result_quantity = other - self.quantity
            pos = self.minus_quantity.to_value(result_quantity.unit)
            neg = self.plus_quantity.to_value(result_quantity.unit)
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        other = self._coerce_additive_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return other - self

    def __mul__(self, other):
        if _is_unit(other):
            return a_u(self.value, self.plus, self.minus, unit=self.unit*other)

        other = self._coerce_multiplicative_operand(other)
        if other is NotImplemented:
            return NotImplemented

        if _has_astropy():
            result_quantity = self.quantity*other.quantity
            derivatives = (other.quantity, self.quantity)
            pos, neg = _propagated_errors(
                (self, other),
                derivatives,
                result_unit=result_quantity.unit,
            )
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        result = self.value*other.value
        pos, neg = _propagated_errors((self, other), (other.value, self.value))
        return a_u(result, pos, neg)

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        if _is_unit(other):
            return a_u(self.value, self.plus, self.minus, unit=self.unit/other)

        other = self._coerce_multiplicative_operand(other)
        if other is NotImplemented:
            return NotImplemented

        if _has_astropy():
            result_quantity = self.quantity/other.quantity
            derivatives = (
                1/other.quantity,
                -self.quantity/other.quantity**2,
            )
            pos, neg = _propagated_errors((self, other), derivatives, result_unit=result_quantity.unit)
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        result = self.value / other.value
        pos, neg = _propagated_errors((self, other), (1/other.value, -self.value/other.value**2))
        return a_u(result, pos, neg)

    def __rtruediv__(self, other):
        if _is_unit(other):
            other = a_u(1, 0, 0, unit=other)
        else:
            other = self._coerce_multiplicative_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return other/self

    @staticmethod
    def _dimensionless_exponent(exponent):
        if isinstance(exponent, a_u):
            if _has_astropy():
                return exponent.to(u.dimensionless_unscaled)
            return exponent
        if _is_quantity(exponent):
            value = exponent.to_value(u.dimensionless_unscaled)
            return a_u(value, 0, 0)
        if isinstance(exponent, Number):
            return a_u(exponent, 0, 0)
        return NotImplemented

    def __pow__(self, other): # self to the something power
        exponent = self._dimensionless_exponent(other)
        if exponent is NotImplemented:
            return NotImplemented

        exponent_is_uncertain = exponent.plus != 0 or exponent.minus != 0
        exponent_is_integer = exponent.value.is_integer()

        if self.value < 0 and (exponent_is_uncertain or not exponent_is_integer):
            raise ValueError("negative bases require an exact integer exponent")
        if self.value == 0 and exponent_is_uncertain:
            raise ValueError("uncertain exponents require a positive base")
        if self.value == 0 and not exponent_is_integer:
            raise ValueError("zero bases require an integer exponent")
        if self.value == 0 and exponent.value < 0:
            raise ValueError("zero cannot be raised to a negative power")
        if exponent_is_uncertain and not _unit_is_dimensionless(self.unit):
            raise ValueError("uncertain exponents require a dimensionless base")

        if _has_astropy():
            result_quantity = self.quantity**exponent.value
            if self.value == 0:
                if exponent.value == 1:
                    base_derivative = 1*self.unit**0
                else:
                    base_derivative = 0*(result_quantity.unit/self.unit)
            else:
                base_derivative = exponent.value*self.quantity**(exponent.value - 1)

            if exponent_is_uncertain:
                exponent_derivative = result_quantity*np.log(self.value)
            else:
                exponent_derivative = 0*result_quantity.unit

            pos, neg = _propagated_errors(
                (self, exponent),
                (base_derivative, exponent_derivative),
                result_unit=result_quantity.unit,
            )
            return a_u(result_quantity.value, pos, neg, unit=result_quantity.unit)

        result = self.value**exponent.value

        if self.value == 0:
            base_derivative = 1.0 if exponent.value == 1 else 0.0
        else:
            base_derivative = exponent.value*self.value**(exponent.value - 1)
        exponent_derivative = result*np.log(self.value) if exponent_is_uncertain else 0.0
        pos, neg = _propagated_errors(
            (self, exponent),
            (base_derivative, exponent_derivative),
        )
        return a_u(result, pos, neg)

    def __rpow__(self, other): # something to the self power
        base = self._coerce_multiplicative_operand(other)
        if base is NotImplemented:
            return NotImplemented
        return base**self

    def _dimensionless_components(self):
        if not _has_astropy():
            return self.value, self.plus, self.minus
        value = self.quantity.to_value(u.dimensionless_unscaled)
        plus = self.plus_quantity.to_value(u.dimensionless_unscaled)
        minus = self.minus_quantity.to_value(u.dimensionless_unscaled)
        return value, plus, minus

    def log10(self):
        value, plus, minus = self._dimensionless_components()
        result = np.log10(value)
        pos = plus/(value*np.log(10))
        neg = minus/(value*np.log(10))
        return a_u(result, pos, neg)

    def log2(self):
        value, plus, minus = self._dimensionless_components()
        result = np.log2(value)
        pos = plus/(value*np.log(2))
        neg = minus/(value*np.log(2))
        return a_u(result, pos, neg)

    def log(self):
        value, plus, minus = self._dimensionless_components()
        result = np.log(value)
        pos = plus/value
        neg = minus/value
        return a_u(result, pos, neg)

    def exp(self):
        value, plus, minus = self._dimensionless_components()
        result = np.exp(value)
        return a_u(result, result*plus, result*minus)

    def sqrt(self):
        return self**0.5

    def _comparison_operand(self, other):
        if isinstance(other, type(self)):
            if _has_astropy():
                return other.to(self.unit)
            return other
        if _is_quantity(other):
            return a_u(other).to(self.unit)
        if isinstance(other, Number):
            other = a_u(other)
            if _has_astropy():
                return other.to(self.unit)
            return other
        return NotImplemented

    def __eq__(self, other):
        try:
            other = self._comparison_operand(other)
        except Exception:
            return False
        if other is NotImplemented:
            return False
        return self.value == other.value and self.plus == other.plus and self.minus == other.minus

    def __gt__(self, other):
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value > other.value

    def __lt__(self, other):
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value < other.value

    def __lshift__(self, other): # overloaded <<; definitively less than
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.maximum < other.minimum

    def __rshift__(self, other): # overloaded >>; definitively greater than
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.minimum > other.maximum

    def __le__(self, other):
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value <= other.value

    def __ge__(self, other):
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value >= other.value

    def conjugate(self):
        return self.value

    def __isfinite__(self):
        return all(np.isfinite(self.items()))

    def isna(self):
        """
        `pandas`-style NaN checker. Returns True if value is NaN or None, and False if neither.
        """
        return np.isnan(self.value) or self.value is None

    def notna(self):
        """
        Inverse of `isna()`. Returns True if value is neither NaN nor None, and False if it is.
        """
        return not self.isna()


class AsymmetricUncertainty(a_u):
    """Alias for legacy namespace support."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        warnings.warn("class AsymmetricUncertainty has been renamed. Use a_u instead.",
                      DeprecationWarning, stacklevel=2)


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
                list.__setitem__(self, i, a_u(self[i], 0, 0))

        self.minus = [v.minus for v in self]
        self.plus = [v.plus for v in self]
        self.values = [v.value for v in self]
        self.units = [v.unit for v in self]

    def __init__(self, array=[]):
        super().__init__(array)
        self.refresh()

    @property
    def as_list(self):
        return self

    @property
    def as_numpy(self):
        return np.asarray(self)

    @property
    def quantities(self):
        """Nominal values as Astropy `Quantity` objects when available."""
        return [entry.quantity for entry in self]

    @property
    def plus_quantities(self):
        return [entry.plus_quantity for entry in self]

    @property
    def minus_quantities(self):
        return [entry.minus_quantity for entry in self]

    @property
    def shape(self) -> tuple:
        return self.as_numpy.shape

    @property
    def ndim(self) -> int:
        return self.as_numpy.ndim

    def flatten(self):
        return self.as_numpy.flatten()

    def to(self, target, equivalencies=None):
        """Return a copy with every entry converted to `target` units."""
        return UncertaintyArray([entry.to(target, equivalencies=equivalencies) for entry in self])

    def __str__(self):
        return str([str(entry) for entry in self])

    def __repr__(self):
        return str(self)

    def append(self, entry):
        super().append(entry)
        self.refresh()

    def extend(self, iterable):
        super().extend(iterable)
        self.refresh()

    def insert(self, index, entry):
        super().insert(index, entry)
        self.refresh()

    def pop(self, index=-1):
        entry = super().pop(index)
        self.refresh()
        return entry

    def remove(self, entry):
        super().remove(entry)
        self.refresh()

    def clear(self):
        super().clear()
        self.refresh()

    def reverse(self):
        super().reverse()
        self.refresh()

    def sort(self, *args, **kwargs):
        super().sort(*args, **kwargs)
        self.refresh()

    def __setitem__(self, key, val):
        super().__setitem__(key, val)
        self.refresh()

    def __delitem__(self, key):
        super().__delitem__(key)
        self.refresh()

    def __iadd__(self, other):
        super().__iadd__(other)
        self.refresh()
        return self

    def __imul__(self, other):
        super().__imul__(other)
        self.refresh()
        return self

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
