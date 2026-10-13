"""Asymmetric Uncertainty: A package for handling non-standard numerical uncertainties.

Copyright (C) 2022 Caden Gobat

This program is free software: you can redistribute and/or modify it under the
terms of the GNU General Public License (either version 3 or, at your option,
any later version).

This program is distributed in the hope that it will be useful, but comes with
ABSOLUTELY NO WARRANTY. See <https://www.gnu.org/licenses/> or the LICENSE file
distributed with this program for more details.
"""

import warnings
from math import erf
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np

np_erf = np.vectorize(erf, doc="Vectorized error function at x.")  # prevents need for SciPy dependency
np_ppf = np.vectorize(NormalDist().inv_cdf, otypes=[float])

_SPLIT_NORMAL_MEAN_FACTOR = np.sqrt(2 / np.pi)
_SPLIT_NORMAL_VARIANCE_FACTOR = 3 / 4 - 2 / np.pi
_SPLIT_NORMAL_THIRD_MOMENT_FACTOR = 6 / np.pi - 2


def _propagated_errors(operands, derivatives):
    """Return first-order asymmetric errors for the supplied partial derivatives."""
    pos_terms = []
    neg_terms = []
    for operand, derivative in zip(operands, derivatives):
        if derivative >= 0:
            pos_error = operand.plus
            neg_error = operand.minus
        else:
            pos_error = operand.minus
            neg_error = operand.plus
        pos_terms.append((derivative * pos_error)**2)
        neg_terms.append((derivative * neg_error)**2)
    return np.sqrt(np.sum(pos_terms)), np.sqrt(np.sum(neg_terms))


class a_u:
    """
    Class for representing and propagating numeric values with asymmetric
    uncertainties. Arithmetic operates on the nominal values and propagates the
    positive and negative uncertainties independently using first-order
    derivatives and quadrature. Reuse of the same object is recognized for
    simple self-addition, subtraction, multiplication, division, and
    exponentiation.

    Distribution utilities interpret the nominal value as the mode and the
    positive and negative uncertainties as the two Gaussian widths of an
    associated jointly normalized split-normal distribution. Arithmetic does
    not attempt to convolve or transform that probability distribution.

    Parameters
    ----------
    nominal : numeric
        the nominal value of the represented quantity
    pos_err : numeric
        the plus error on the value
    neg_err : numeric
        the minus error on the value

    Attributes
    ----------
    value : numeric
        the nominal value of the represented quantity
    plus : numeric
        the positive error on the value
    minus : numeric
        the negative error on the value
    """

    def __init__(self, nominal, pos_err=0, neg_err=0):
        if isinstance(nominal, str):
            stripped = nominal.replace(" ", "")
            if "±" in stripped:
                self.value = float(stripped.split("±")[0])
                self.plus = self.minus = float(nominal.split("±")[1])
            elif "+" in stripped and "-" in stripped:
                self.value = float(stripped.split("(")[0])
                err_str = stripped.split("(")[1].replace(")", "")
                self.plus = float(err_str.split(",")[0][1:])
                self.minus = float(err_str.split(",")[1][1:])
            else:
                raise ValueError(
                    "Failed to parse string, likely due to improper formatting."
                )
        else:
            self.value = float(nominal)
            self.plus = abs(float(pos_err))
            self.minus = abs(float(neg_err))

    @property
    def maximum(self):
        """Return the nominal value plus the positive uncertainty."""
        return self.value + self.plus

    @property
    def minimum(self):
        """Return the nominal value minus the negative uncertainty."""
        return self.value - self.minus

    @property
    def sign(self):
        """Return 1 for a non-negative nominal value and -1 otherwise."""
        return 1 if self.value >= 0 else -1

    @property
    def is_symmetric(self) -> bool:
        """Return whether the positive and negative uncertainties are equal within tolerance."""
        return bool(np.isclose(self.plus, self.minus))

    @property
    def mean(self):
        """Mean of the associated split-normal distribution."""
        return self.value + _SPLIT_NORMAL_MEAN_FACTOR * (self.plus - self.minus)

    @property
    def median(self):
        """Median of the associated split-normal distribution."""
        return self.ppf(0.5)

    @property
    def variance(self):
        """Variance of the associated split-normal distribution."""
        width_sum = self.plus + self.minus
        width_diff = self.plus - self.minus
        return (width_sum**2) / 4 + _SPLIT_NORMAL_VARIANCE_FACTOR * width_diff**2

    @property
    def std(self):
        """Standard deviation of the associated split-normal distribution."""
        return np.sqrt(self.variance)

    @property
    def _third_central_moment(self):
        """Third central moment of the associated split-normal distribution."""
        width_diff = self.plus - self.minus
        return (
            _SPLIT_NORMAL_MEAN_FACTOR
            * width_diff
            * (self.variance + _SPLIT_NORMAL_THIRD_MOMENT_FACTOR * width_diff**2)
        )

    @property
    def skewness(self):
        """Standardized third central moment of the associated split-normal distribution."""
        variance = self.variance
        if variance == 0:
            return np.nan
        return self._third_central_moment / variance**1.5

    @property
    def kurtosis(self):
        """Fisher (excess) kurtosis of the associated split-normal distribution."""
        variance = self.variance
        if variance == 0:
            return np.nan

        width_sum = self.plus + self.minus
        mean_shift = self.mean - self.value
        second_about_mode = (self.minus**3 + self.plus**3) / width_sum
        third_about_mode = (
            2 * _SPLIT_NORMAL_MEAN_FACTOR * (self.plus**4 - self.minus**4) / width_sum
        )
        fourth_about_mode = 3 * (self.minus**5 + self.plus**5) / width_sum
        fourth_central_moment = (
            fourth_about_mode
            - 4 * mean_shift * third_about_mode
            + 6 * mean_shift**2 * second_about_mode
            - 3 * mean_shift**4
        )
        return fourth_central_moment / variance**2 - 3

    def __str__(self) -> str:
        return format(self)

    def __format__(self, fmt_spec) -> str:
        if self.is_symmetric:
            return f"{format(self.value, fmt_spec)} ± {format(self.plus, fmt_spec)}"
        else:
            return "{} (+{}, -{})".format(*(format(x, fmt_spec) for x in self.items()))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.value!r}, {self.plus!r}, {self.minus!r})"

    def _repr_latex_(self, float_fmt="f"):
        if self.is_symmetric:
            return (
                f"${format(self.value, float_fmt)}"
                r" \pm "
                f"{format(self.plus, float_fmt)}$"
            )
        else:
            return "${}^{{+{}}}_{{-{}}}$".format(*(format(x, float_fmt) for x in self.items()))

    def pdf(self, x):
        """Evaluate the associated split-normal probability density function."""
        x_arr = np.asanyarray(x).astype(float)
        pdf_arr = np.piecewise(
            x_arr,
            [x_arr < self.value, x_arr >= self.value],
            [
                lambda x: np.sqrt(2) / np.sqrt(np.pi) / (self.plus + self.minus) * np.exp(-1 * (x - self.value)**2 / (2 * self.minus**2)),
                lambda x: np.sqrt(2) / np.sqrt(np.pi) / (self.plus + self.minus) * np.exp(-1 * (x - self.value)**2 / (2 * self.plus**2)),
            ],
        )
        return pdf_arr.item() if np.isscalar(x) else pdf_arr

    def cdf(self, x):
        """Evaluate the associated split-normal cumulative distribution function."""
        x_arr = np.asanyarray(x).astype(float)
        width_sum = self.plus + self.minus
        lower_cdf = (
            0.0 if self.minus == 0
            else lambda x: self.minus / width_sum * (1 + np_erf((x - self.value) / (self.minus * np.sqrt(2))))
        )
        upper_cdf = (
            1.0 if self.plus == 0
            else lambda x: (self.minus + self.plus * np_erf((x - self.value) / (self.plus * np.sqrt(2)))) / width_sum
        )
        cdf_arr = np.piecewise(
            x_arr,
            [x_arr < self.value, x_arr >= self.value],
            [lower_cdf, upper_cdf]
        )
        return cdf_arr.item() if np.isscalar(x) else cdf_arr

    def ppf(self, q):
        """Return quantiles of the associated split-normal distribution."""
        q_arr = np.asanyarray(q).astype(float)
        finite = np.isfinite(q_arr)
        if np.any((q_arr[finite] < 0) | (q_arr[finite] > 1)):
            raise ValueError("quantiles must lie between 0 and 1")

        result = np.full(q_arr.shape, np.nan, dtype=float)
        width_sum = self.plus + self.minus
        if width_sum == 0:
            result[finite] = self.value
            return result.item() if np.isscalar(q) else result

        lower_endpoint = finite & (q_arr == 0)
        upper_endpoint = finite & (q_arr == 1)
        result[lower_endpoint] = self.value if self.minus == 0 else -np.inf
        result[upper_endpoint] = self.value if self.plus == 0 else np.inf

        interior = finite & (q_arr > 0) & (q_arr < 1)
        mode_probability = self.minus / width_sum
        lower = interior & (q_arr < mode_probability)
        upper = interior & ~lower

        if np.any(lower):
            probability = q_arr[lower] * width_sum / (2 * self.minus)
            result[lower] = self.value + self.minus * np_ppf(probability)
        if np.any(upper):
            tail_probability = (1 - q_arr[upper]) * width_sum / (2 * self.plus)
            result[upper] = self.value - self.plus * np_ppf(tail_probability)

        return result.item() if np.isscalar(q) else result

    def rvs(self, size=None, rng=None):
        """Draw random values from the associated split-normal distribution."""
        if isinstance(rng, np.random.Generator):
            generator = rng
        else:
            generator = np.random.default_rng(rng)

        width_sum = self.plus + self.minus
        if width_sum == 0:
            if size is None:
                return self.value
            return np.full(size, self.value, dtype=float)

        negative_side = generator.random(size=size) < self.minus / width_sum
        magnitude = np.abs(generator.standard_normal(size=size))
        samples = np.where(
            negative_side,
            self.value - self.minus * magnitude,
            self.value + self.plus * magnitude,
        )
        return samples.item() if size is None else samples

    def interval(self, n_sigma=None, confidence=None):
        """Return a mode-centered interval in units of the two Gaussian widths."""
        if n_sigma is not None and confidence is not None:
            raise ValueError("specify either n_sigma or confidence, not both")

        if confidence is not None:
            confidence = float(confidence)
            if not 0 < confidence < 1:
                raise ValueError("confidence must lie strictly between 0 and 1")
            n_sigma = -NormalDist().inv_cdf((1 - confidence) / 2)
        elif n_sigma is None:
            n_sigma = 1.0

        n_sigma = float(n_sigma)
        if not np.isfinite(n_sigma) or n_sigma < 0:
            raise ValueError("n_sigma must be a finite non-negative number")

        return (
            self.value - n_sigma * self.minus,
            self.value + n_sigma * self.plus,
        )

    def equal_tail_interval(self, confidence=None):
        """Return an equal-tail interval containing the requested probability mass.

        The default confidence is the conventional one-sigma Gaussian probability.
        """
        if confidence is None:
            confidence = erf(1 / np.sqrt(2))
        confidence = float(confidence)
        if not 0 < confidence < 1:
            raise ValueError("confidence must lie strictly between 0 and 1")

        tail_probability = (1 - confidence) / 2
        lower = self.ppf(tail_probability)
        upper = -(-self).ppf(tail_probability)
        return lower, upper

    def pdfplot(self, num_sigma=5, discretization=100, ax=None, **kwargs):
        """Plot the associated PDF over the specified number of sigma, using ``2*discretization`` points.

        ``**kwargs`` are passed on to ``matplotlib`` for configuration of the resulting plot.
        """
        neg_x = np.linspace(
            self.value - (num_sigma * self.minus),
            self.value,
            discretization
        )
        pos_x = np.linspace(
            self.value,
            self.value + (num_sigma * self.plus),
            discretization
        )
        x = np.concatenate([neg_x, pos_x])
        pdf = self.pdf(x)
        if ax is None:
            ax = plt.axes()
        ax.plot(x, pdf, **kwargs)

    def cdfplot(self, num_sigma=5, discretization=100, ax=None, **kwargs):
        """Plot the associated CDF over the specified number of sigma, using ``2*discretization`` points.

        ``**kwargs`` are passed on to ``matplotlib`` for configuration of the resulting plot.
        """
        neg_x = np.linspace(self.value - (num_sigma * self.minus), self.value, discretization)
        pos_x = np.linspace(self.value, self.value + (num_sigma * self.plus), discretization)
        x = np.concatenate([neg_x, pos_x])
        cdf = self.cdf(x)
        if ax is None:
            ax = plt.axes()
        ax.plot(x, cdf, **kwargs)

    def add_error(self, delta, method="quadrature", inplace=False):
        """Add ``delta`` to an instance's existing error.

        Possible ``method``s are ``quadrature``, ``straight``, or ``split``. If
        ``inplace`` is ``True``, the existing object's errors are modified in
        place. If it is ``False`` (default), a new instance is returned.
        """
        if method == "quadrature":
            new_pos = np.hypot(self.plus, delta)
            new_neg = np.hypot(self.minus, delta)
        elif method == "straight":
            new_pos = self.plus + delta
            new_neg = self.minus + delta
        elif method == "split":
            new_pos = self.plus + delta / 2
            new_neg = self.minus + delta / 2
        else:
            raise ValueError
        if inplace:
            self.plus = float(new_pos)
            self.minus = float(new_neg)
        else:
            return a_u(self.value, new_pos, new_neg)

    def items(self):
        """Return a tuple of ``(value, plus, minus)``."""
        return (self.value, self.plus, self.minus)

    def __int__(self):
        """Convert the nominal value to an integer."""
        return int(self.value)

    def __float__(self):
        """Convert the nominal value to a float."""
        return float(self.value)

    def __neg__(self):
        """Return the additive inverse, swapping the positive and negative uncertainties."""
        return a_u(-self.value, self.minus, self.plus)

    def __add__(self, other):
        """Add another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return self * 2
        if not isinstance(other, type(self)):
            other = a_u(other, 0, 0)
        result = self.value + other.value
        pos, neg = _propagated_errors((self, other), (1, 1))
        return a_u(result, pos, neg)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """Subtract another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return a_u(0.0)
        if not isinstance(other, type(self)):
            other = a_u(other, 0, 0)
        result = self.value - other.value
        pos, neg = _propagated_errors((self, other), (1, -1))
        return a_u(result, pos, neg)

    def __rsub__(self, other):
        if not isinstance(other, type(self)):
            other = a_u(other, 0, 0)
        return other - self

    def __mul__(self, other):
        """Multiply by another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return self**2
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)

        result = self.value * other.value
        pos, neg = _propagated_errors((self, other), (other.value, self.value))
        return a_u(result, pos, neg)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):  # self divided by something
        """Divide by another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            if self.variance == 0 and self.value == 0:
                raise ZeroDivisionError("division by zero")
            return a_u(1.0)
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        result = self.value / other.value
        pos, neg = _propagated_errors(
            (self, other),
            (1 / other.value, -self.value / other.value**2)
        )
        return a_u(result, pos, neg)

    def __rtruediv__(self, other):  # something divided by self
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        result = other.value / self.value
        pos, neg = _propagated_errors(
            (other, self),
            (1 / self.value, -other.value / self.value**2)
        )
        return a_u(result, pos, neg)

    def __pow__(self, other):  # self to the something power
        """Raise this value to a power with first-order asymmetric uncertainty propagation."""
        if isinstance(other, type(self)):
            exponent = other
        else:
            exponent = a_u(other, 0, 0)

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

        if exponent is self and exponent_is_uncertain:
            result = self.value**self.value
            derivative = result * (np.log(self.value) + 1)
            pos, neg = _propagated_errors((self,), (derivative,))
            return a_u(result, pos, neg)

        result = self.value**exponent.value

        if self.value == 0:
            if exponent.value == 1:
                base_derivative = 1.0
            else:
                base_derivative = 0.0
        else:
            base_derivative = exponent.value * self.value**(exponent.value - 1)

        if exponent_is_uncertain:
            exponent_derivative = result * np.log(self.value)
        else:
            exponent_derivative = 0.0

        pos, neg = _propagated_errors(
            (self, exponent),
            (base_derivative, exponent_derivative),
        )
        return a_u(result, pos, neg)

    def __rpow__(self, other):  # something to the self power
        if isinstance(other, type(self)):
            base = other
        else:
            base = a_u(other, 0, 0)
        return base**self

    def log10(self):
        """Return the base-10 logarithm with first-order propagated uncertainties."""
        result = np.log10(self.value)
        pos = self.plus / (self.value * np.log(10))
        neg = self.minus / (self.value * np.log(10))
        return a_u(result, pos, neg)

    def log(self):
        """Return the natural logarithm with first-order propagated uncertainties."""
        result = np.log(self.value)
        pos = self.plus / self.value
        neg = self.minus / self.value
        return a_u(result, pos, neg)

    def exp(self):
        """Return the exponential with first-order propagated uncertainties."""
        return np.exp(1)**self

    def sqrt(self):
        """Return the square root with first-order propagated uncertainties."""
        return self**0.5

    def __eq__(self, other):
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        return (
            self.value == other.value
            and self.plus == other.plus
            and self.minus == other.minus
        )

    def __gt__(self, other):  # self > other ?
        if isinstance(other, type(self)):
            return self.value > other.value
        return self.value > other

    def __lt__(self, other):  # self < other ?
        if isinstance(other, type(self)):
            return self.value < other.value
        return self.value < other

    def __lshift__(self, other):  # overloaded <<; definitively less than
        """Return whether this uncertainty interval lies entirely below the other."""
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        return self.maximum < other.minimum

    def __rshift__(self, other):  # overloaded >>; definitively greater than
        """Return whether this uncertainty interval lies entirely above the other."""
        if isinstance(other, type(self)):
            pass
        else:
            other = a_u(other, 0, 0)
        return self.minimum > other.maximum

    def __le__(self, other):  # self <= other ?
        if isinstance(other, type(self)):
            return self.value <= other.value
        return self.value <= other

    def __ge__(self, other):  # self >= other ?
        if isinstance(other, type(self)):
            return self.value >= other.value
        return self.value >= other

    def conjugate(self):
        """Return the real-valued nominal value unchanged by conjugation."""
        return self.value

    def __isfinite__(self):
        """Return whether the nominal value and both uncertainties are finite."""
        return all(np.isfinite(self.items()))

    def isna(self):
        """Return whether the nominal value is NaN or ``None``."""
        return np.isnan(self.value) or (self.value is None)

    def notna(self):
        """Return whether the nominal value is neither NaN nor ``None``."""
        return not self.isna()


def propagate(function, *args, size=10000, rng=None):
    """Evaluate a callable on split-normal samples drawn from ``a_u`` arguments.

    Each distinct ``a_u`` object is sampled once, so passing the same object
    more than once preserves perfect correlation. Separate objects are sampled
    independently, even if they have identical values and uncertainties.
    Non-``a_u`` arguments are passed to ``function`` unchanged.

    Parameters
    ----------
    function : callable
        Function to evaluate on the sampled arguments.
    *args
        Positional arguments for ``function``. ``a_u`` instances are replaced
        by sample arrays; other arguments are passed through unchanged.
    size : int, optional
        Number of Monte Carlo samples to draw. Default is 10000.
    rng : numpy.random.Generator or int or None, optional
        Random generator or seed passed to the split-normal sampler.

    Returns
    -------
    object : array-like
        The return value from ``function``. No split-normal approximation is
        fitted to the Monte Carlo result.
    """
    if isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size <= 0:
        raise ValueError("size must be a positive integer")

    if isinstance(rng, np.random.Generator):
        generator = rng
    else:
        generator = np.random.default_rng(rng)

    samples = {}
    sampled_args = []
    for arg in args:
        if isinstance(arg, a_u):
            key = id(arg)
            if key not in samples:
                samples[key] = arg.rvs(size=size, rng=generator)
            sampled_args.append(samples[key])
        else:
            sampled_args.append(arg)

    return function(*sampled_args)


class AsymmetricUncertainty(a_u):  # alias for legacy namespace support
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        warnings.warn(
            "class AsymmetricUncertainty has been renamed. Use a_u instead.",
            DeprecationWarning,
            stacklevel=2,
        )


class UncertaintyArray(list):
    """Class for representing an array of ``a_u`` objects.

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

    def __init__(self, array=None):
        if array is None:
            array = []
        super().__init__(array)
        self.refresh()

    @property
    def as_list(self):
        return self

    @property
    def as_numpy(self):
        return np.asarray(self)

    @property
    def shape(self) -> tuple:
        return self.as_numpy.shape

    @property
    def ndim(self) -> int:
        return self.as_numpy.ndim

    def flatten(self):
        return self.as_numpy.flatten()

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
    """Standalone function to return an array of the positive errors of an array of ``a_u`` objects.

    Functionally equivalent to ``UncertaintyArray(array).plus``.
    """
    return [v.plus for v in array]


def neg_errors(array):
    """Standalone function to return an array of the negative errors of an array of ``a_u`` objects.

    Functionally equivalent to ``UncertaintyArray(array).minus``.
    """
    return [v.minus for v in array]
