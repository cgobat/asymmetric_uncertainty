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
from numbers import Number
try:
    from astropy import units as u
except ImportError:  # Astropy is installed by the package, but keep source-tree use graceful.
    u = None


np_erf = np.vectorize(erf, doc="Vectorized error function at x.")  # prevents need for SciPy dependency
np_ppf = np.vectorize(NormalDist().inv_cdf, otypes=[float])

_SPLIT_NORMAL_MEAN_FACTOR = np.sqrt(2 / np.pi)
_SPLIT_NORMAL_VARIANCE_FACTOR = 3 / 4 - 2 / np.pi
_SPLIT_NORMAL_THIRD_MOMENT_FACTOR = 6 / np.pi - 2


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

def _propagated_errors(operands, derivatives, result_unit=None) -> tuple[float, float]:
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
    nominal : numeric or astropy.units.Quantity
        the nominal value of the represented quantity
    pos_err : numeric or astropy.units.Quantity
        the plus error on the value
    neg_err : numeric or astropy.units.Quantity
        the minus error on the value
    str or astropy.units.UnitBase, optional : unit
        physical unit for the value and its uncertainties. If omitted, the first
        constructor argument with units supplies the unit for the object.

    Attributes
    ----------
    value : numeric
        the nominal value of the represented quantity, expressed in `unit`
    plus : numeric
        the positive error on the value, expressed in `unit`
    minus : numeric
        the negative error on the value, expressed in `unit`
    unit : astropy.units.UnitBase
        the physical unit, or `dimensionless_unscaled` for ordinary numbers
    """

    def __init__(self, nominal, pos_err=0, neg_err=0, unit=None):
        if isinstance(nominal, str):
            stripped = nominal.replace(" ", "")
            if "+/-" in stripped:
                stripped = stripped.replace("+/-", "±")
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
                raise ValueError(
                    "Failed to parse string, likely due to improper formatting."
                )

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
    def _value_in_unit(value, unit) -> float:
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
        """Return a copy converted to `target` units.

        For nonlinear equivalencies, transform the interval endpoints and
        approximate the result with new asymmetric widths. The resulting PDF
        need not remain split-normal.
        """
        if not _has_astropy():
            raise ImportError("Astropy is required for unit conversion.")
        if equivalencies is None:
            equivalencies = []
        target = u.Unit(target)
        value = self.quantity.to_value(target, equivalencies=equivalencies)
        lower = (self.value - self.minus)*self.unit
        upper = (self.value + self.plus)*self.unit
        lower = lower.to_value(target, equivalencies=equivalencies)
        upper = upper.to_value(target, equivalencies=equivalencies)
        return a_u(value, max(lower, upper) - value,
                   value - min(lower, upper), unit=target)

    def to_value(self, target=None, equivalencies=None):
        """Return the nominal value, optionally converted to `target` units."""
        if target is None or not _has_astropy():
            return self.value
        if equivalencies is None:
            equivalencies = []
        return self.quantity.to_value(target, equivalencies=equivalencies)

    @property
    def maximum(self) -> float:
        """Return the nominal value plus the positive uncertainty."""
        result = self.value + self.plus
        if _unit_is_dimensionless(self.unit):
            return result
        return result*self.unit

    @property
    def minimum(self) -> float:
        """Return the nominal value minus the negative uncertainty."""
        result = self.value - self.minus
        if _unit_is_dimensionless(self.unit):
            return result
        return result*self.unit

    @property
    def sign(self) -> int:
        """Return 1 for a non-negative nominal value and -1 otherwise."""
        return 1 if self.value >= 0 else -1

    @property
    def is_symmetric(self) -> bool:
        """Return whether the positive and negative uncertainties are equal within tolerance."""
        return bool(np.isclose(self.plus, self.minus))

    @property
    def mean(self) -> float:
        """Mean of the associated split-normal distribution."""
        mean = self.value + _SPLIT_NORMAL_MEAN_FACTOR * (self.plus - self.minus)
        return mean if _unit_is_dimensionless(self.unit) else mean*self.unit

    @property
    def median(self) -> float:
        """Median of the associated split-normal distribution."""
        return self.ppf(0.5)

    @property
    def variance(self) -> float:
        """Variance of the associated split-normal distribution."""
        width_sum = self.plus + self.minus
        width_diff = self.plus - self.minus
        variance = (width_sum**2) / 4 + _SPLIT_NORMAL_VARIANCE_FACTOR * width_diff**2
        return variance if _unit_is_dimensionless(self.unit) else variance*self.unit**2

    @property
    def std(self) -> float:
        """Standard deviation of the associated split-normal distribution."""
        return np.sqrt(self.variance)

    @property
    def _third_central_moment(self) -> float:
        """Third central moment of the associated split-normal distribution."""
        width_diff = self.plus - self.minus
        variance = self.variance
        if _is_quantity(variance):
            variance = variance.to_value(self.unit**2)
        return (
            _SPLIT_NORMAL_MEAN_FACTOR
            * width_diff
            * (variance + _SPLIT_NORMAL_THIRD_MOMENT_FACTOR * width_diff**2)
        )

    @property
    def skewness(self) -> float:
        """Standardized third central moment of the associated split-normal distribution."""
        variance = self.variance
        if _is_quantity(variance):
            variance = variance.to_value(self.unit**2)
        if variance == 0:
            return np.nan
        return self._third_central_moment / variance**1.5

    @property
    def kurtosis(self) -> float:
        """Fisher (excess) kurtosis of the associated split-normal distribution."""
        variance = self.variance
        if _is_quantity(variance):
            variance = variance.to_value(self.unit**2)
        if variance == 0:
            return np.nan

        width_sum = self.plus + self.minus
        mean = self.mean
        if _is_quantity(mean):
            mean = mean.to_value(self.unit)
        mean_shift = mean - self.value
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
        unit_suffix = "" if _unit_is_dimensionless(self.unit) else f" {self.unit}"
        if self.is_symmetric:
            return f"{format(self.value, fmt_spec)} ± {format(self.plus, fmt_spec)}{unit_suffix}"
        return "%s (+%s, -%s)%s" % (*(format(x, fmt_spec) for x in self.items()), unit_suffix)

    def __repr__(self) -> str:
        params = f"{self.value!r}, {self.plus!r}, {self.minus!r}"
        if not _unit_is_dimensionless(self.unit):
            params += f", unit={self.unit.to_string()!r}"
        return f"{type(self).__name__}({params})"

    def _repr_latex_(self, float_fmt="f"):
        if self.is_symmetric:
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
        """Evaluate the associated split-normal probability density function.

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
        """Evaluate the associated split-normal cumulative distribution function."""
        x_arr = self._x_values(x)
        width_sum = self.plus + self.minus
        lower_cdf = (
            0.0 if self.minus == 0
            else lambda values: self.minus/width_sum * (1 + np_erf((values-self.value)/(self.minus*np.sqrt(2))))
        )
        upper_cdf = (
            1.0 if self.plus == 0
            else lambda values: (self.minus + self.plus*np_erf((values-self.value)/(self.plus*np.sqrt(2))))/width_sum
        )
        cdf_arr = np.piecewise(x_arr, [x_arr<self.value, x_arr>=self.value],
                               [lower_cdf, upper_cdf])
        return cdf_arr.item() if x_arr.ndim == 0 else cdf_arr

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
            values = result.item() if np.isscalar(q) else result
            return values if _unit_is_dimensionless(self.unit) else values*self.unit

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

        values = result.item() if np.isscalar(q) else result
        return values if _unit_is_dimensionless(self.unit) else values*self.unit

    def rvs(self, size=None, rng=None):
        """Draw random values from the associated split-normal distribution."""
        if isinstance(rng, np.random.Generator):
            generator = rng
        else:
            generator = np.random.default_rng(rng)

        width_sum = self.plus + self.minus
        if width_sum == 0:
            samples = self.value if size is None else np.full(size, self.value, dtype=float)
            return samples if _unit_is_dimensionless(self.unit) else samples*self.unit

        negative_side = generator.random(size=size) < self.minus / width_sum
        magnitude = np.abs(generator.standard_normal(size=size))
        samples = np.where(
            negative_side,
            self.value - self.minus * magnitude,
            self.value + self.plus * magnitude,
        )
        samples = samples.item() if size is None else samples
        return samples if _unit_is_dimensionless(self.unit) else samples*self.unit

    @classmethod
    def fit(cls, samples):
        """Fit the split-normal distribution to a sample by maximum likelihood.

        Parameters
        ----------
        samples : array-like or astropy.units.Quantity, shape (n,)
            At least two finite observations from the distribution. Quantity
            samples retain their physical unit in the fitted result.

        Returns
        -------
        a_u
            A fitted distribution whose ``value`` is the mode and whose
            ``plus`` and ``minus`` are its upper and lower Gaussian widths.

        Notes
        -----
        The fitted widths describe the sample distribution, not the statistical
        uncertainty on the estimated mode. The fit assumes independent samples
        of the continuous split-normal distribution used by :meth:`pdf` and
        :meth:`rvs`. Identical samples give an exact value with zero widths.
        For small or strongly one-sided samples, one fitted width may also be
        zero.
        """
        unit = samples.unit if _is_quantity(samples) else None
        data = np.asarray(samples.to_value(unit) if unit is not None else samples,
                          dtype=float)
        if data.ndim != 1 or data.size < 2:
            raise ValueError("samples must be a 1-D array with at least two points")
        if not np.all(np.isfinite(data)):
            raise ValueError("samples must contain only finite values")

        center = float(np.median(data))
        residuals = data - center
        scale = float(np.max(np.abs(residuals)))
        if scale == 0:
            return cls(center, 0, 0, unit=unit)
        scaled = residuals / scale
        count = scaled.size

        # Profile out the two widths: for a fixed mode, the MLE widths are
        # proportional to the cube roots of the sums of squared deviations
        # on their respective sides. The remaining objective is 1D.
        def profile(mode):
            offsets = scaled - mode
            lower = np.sum(np.square(offsets[offsets < 0]))
            upper = np.sum(np.square(offsets[offsets >= 0]))
            return np.cbrt(lower) + np.cbrt(upper)

        grid = np.unique(np.quantile(scaled, np.linspace(0, 1, 33)))
        scores = np.array([profile(value) for value in grid])
        index = int(np.argmin(scores))
        left = grid[max(0, index - 1)]
        right = grid[min(grid.size - 1, index + 1)]

        # Refine the best grid point with a bounded golden-section search.
        ratio = (np.sqrt(5) - 1) / 2
        a = right - ratio * (right - left)
        b = left + ratio * (right - left)
        fa, fb = profile(a), profile(b)
        for _ in range(64):
            if fa <= fb:
                right, b, fb = b, a, fa
                a = right - ratio * (right - left)
                fa = profile(a)
            else:
                left, a, fa = a, b, fb
                b = left + ratio * (right - left)
                fb = profile(b)

        mode = (left + right) / 2
        if scores[index] <= profile(mode):
            mode = grid[index]

        offsets = scaled - mode
        lower = np.cbrt(np.sum(np.square(offsets[offsets < 0])))
        upper = np.cbrt(np.sum(np.square(offsets[offsets >= 0])))
        factor = scale * np.sqrt((lower + upper) / count)
        return cls(center + scale * mode, upper * factor, lower * factor, unit=unit)

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

        lower = self.value - n_sigma * self.minus
        upper = self.value + n_sigma * self.plus
        if not _unit_is_dimensionless(self.unit):
            return lower*self.unit, upper*self.unit
        return lower, upper

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
        if _is_quantity(pdf):
            pdf = pdf.value
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
        if _is_quantity(delta):
            delta = delta.to_value(self.unit)
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
            return a_u(self.value, new_pos, new_neg, unit=self.unit)

    def items(self, units=False):
        """
        Return a tuple of `(value, plus, minus)`. Set `units=True` to return
        Astropy `Quantity` objects.
        """
        if units and _has_astropy():
            return (self.quantity, self.plus_quantity, self.minus_quantity)
        return (self.value, self.plus, self.minus)

    def __int__(self):
        """Convert the nominal value to an integer."""
        if _unit_is_dimensionless(self.unit):
            return int(self.value)
        return int(self.quantity)

    def __float__(self):
        """Convert the nominal value to a float."""
        if _unit_is_dimensionless(self.unit):
            return float(self.value)
        return float(self.quantity)

    def __neg__(self):
        """Return the additive inverse, swapping the positive and negative uncertainties."""
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
        """Add another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return self * 2
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
        """Subtract another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return a_u(0, unit=self.unit)
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
        """Multiply by another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            return self**2
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

        result = self.value * other.value
        pos, neg = _propagated_errors((self, other), (other.value, self.value))
        return a_u(result, pos, neg)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):  # self divided by something
        """Divide by another value with first-order asymmetric uncertainty propagation."""
        if other is self:
            if self.value == 0 and self.plus == 0 and self.minus == 0:
                raise ZeroDivisionError("division by zero")
            return a_u(1)
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

    def __rtruediv__(self, other):  # something divided by self
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

    def __pow__(self, other):  # self to the something power
        """Raise this value to a power with first-order asymmetric uncertainty propagation."""
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

        if other is self and exponent_is_uncertain:
            result = self.value**self.value
            derivative = result * (np.log(self.value) + 1)
            pos, neg = _propagated_errors(
                (self,), (derivative,),
                result_unit=u.dimensionless_unscaled if _has_astropy() else None,
            )
            return a_u(result, pos, neg)

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

    def __rpow__(self, other):  # something to the self power
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
        """Return the base-10 logarithm with first-order propagated uncertainties."""
        value, plus, minus = self._dimensionless_components()
        result = np.log10(value)
        pos = plus/(value*np.log(10))
        neg = minus/(value*np.log(10))
        return a_u(result, pos, neg)

    def log2(self):
        """Return the base-2 logarithm with first-order propagated uncertainties."""
        value, plus, minus = self._dimensionless_components()
        result = np.log2(value)
        pos = plus/(value*np.log(2))
        neg = minus/(value*np.log(2))
        return a_u(result, pos, neg)

    def log(self):
        """Return the natural logarithm with first-order propagated uncertainties."""
        value, plus, minus = self._dimensionless_components()
        result = np.log(value)
        pos = plus/value
        neg = minus/value
        return a_u(result, pos, neg)

    def exp(self):
        """Return the exponential with first-order propagated uncertainties."""
        return np.exp(1)**self

    def sqrt(self):
        """Return the square root with first-order propagated uncertainties."""
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

    def __gt__(self, other) -> bool:  # self > other ?
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value > other.value

    def __lt__(self, other) -> bool:  # self < other ?
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value < other.value

    def __lshift__(self, other) -> bool:  # overloaded <<; definitively less than
        """Return whether this uncertainty interval lies entirely below the other."""
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.maximum < other.minimum

    def __rshift__(self, other) -> bool:  # overloaded >>; definitively greater than
        """Return whether this uncertainty interval lies entirely above the other."""
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.minimum > other.maximum

    def __le__(self, other) -> bool:  # self <= other ?
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value <= other.value

    def __ge__(self, other) -> bool:  # self >= other ?
        other = self._comparison_operand(other)
        if other is NotImplemented:
            return NotImplemented
        return self.value >= other.value

    def conjugate(self):
        """Return the real-valued nominal value unchanged by conjugation."""
        return self.value

    def __isfinite__(self) -> bool:
        """Return whether the nominal value and both uncertainties are finite."""
        return all(np.isfinite(self.items()))

    def isna(self) -> bool:
        """Return whether the nominal value is NaN or ``None``."""
        return np.isnan(self.value) or self.value is None

    def notna(self) -> bool:
        """Return whether the nominal value is neither NaN nor ``None``."""
        return not self.isna()

    def _trig(self, function_name):
        if _has_astropy() and self.unit.is_equivalent(u.dimensionless_unscaled):
            angle = self.quantity.to_value(u.dimensionless_unscaled)
            plus = self.plus_quantity.to_value(u.dimensionless_unscaled)
            minus = self.minus_quantity.to_value(u.dimensionless_unscaled)
        elif _has_astropy():
            angle = self.quantity.to_value(u.rad)
            plus = self.plus_quantity.to_value(u.rad)
            minus = self.minus_quantity.to_value(u.rad)
        else:
            angle = self.value
            plus = self.plus
            minus = self.minus

        if function_name == "sin":
            result = np.sin(angle)
            derivative = np.cos(angle)
        elif function_name == "cos":
            result = np.cos(angle)
            derivative = -np.sin(angle)
        elif function_name == "tan":
            result = np.tan(angle)
            derivative = 1/np.cos(angle)**2
        else:
            raise ValueError(function_name)

        if derivative >= 0:
            pos = derivative*plus
            neg = derivative*minus
        else:
            pos = derivative*minus
            neg = derivative*plus
        return a_u(result, np.abs(pos), np.abs(neg))

    __array_priority__ = 50000 # astropy Quantities are 10000 and FunctionQuantities are 40000

    def __array_ufunc__(self, function, method, *inputs, **kwargs):
        if method != "__call__" or kwargs.get("out") is not None:
            return NotImplemented

        name = function.__name__
        if name in ("add", "subtract", "multiply", "true_divide", "divide", "power"):
            left, right = inputs
            if name == "add":
                return self.__add__(right) if left is self else self.__radd__(left)
            if name == "subtract":
                return self.__sub__(right) if left is self else self.__rsub__(left)
            if name == "multiply":
                return self.__mul__(right) if left is self else self.__rmul__(left)
            if name in ("true_divide", "divide"):
                return self.__truediv__(right) if left is self else self.__rtruediv__(left)
            return self.__pow__(right) if left is self else self.__rpow__(left)

        if name in ("sin", "cos", "tan"):
            return self._trig(name)
        if name == "log":
            return self.log()
        if name == "log2":
            return self.log2()
        if name == "log10":
            return self.log10()
        if name == "exp":
            return self.exp()
        if name == "sqrt":
            return self.sqrt()
        if name in ("absolute", "fabs"):
            return abs(self)
        if name == "negative":
            return -self
        if name == "positive":
            return a_u(self.value, self.plus, self.minus, unit=self.unit)
        if name == "isfinite":
            return all(np.isfinite(self.items()))
        if name == "isnan":
            return self.isna()
        if name == "sign":
            return self.sign
        return NotImplemented


def propagate(function, *args, size=10000, rng=None, fit=False):
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
    fit : bool, optional
        If True, fit an ``a_u`` split-normal distribution to the output samples.
        Default is False, returning the callable's output unchanged.

    Returns
    -------
    object : array-like or a_u
        The callable's output when ``fit=False``, or an ``a_u`` fitted to its
        one-dimensional, finite output sample when ``fit=True``. Fitting is a
        parametric approximation and may change the nominal value.
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

    result = function(*sampled_args)
    if fit:
        return a_u.fit(result)
    return result


class AsymmetricUncertainty(a_u):  # alias for legacy namespace support
    """Alias for legacy namespace support."""

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
        self.units = [v.unit for v in self]

    def __init__(self, array=None):
        if array is None:
            array = []
        super().__init__(array)
        self.refresh()

    @property
    def as_list(self) -> list:
        return self

    @property
    def as_numpy(self) -> np.ndarray:
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

    def flatten(self) -> np.ndarray:
        return self.as_numpy.flatten()

    def to(self, target, equivalencies=None):
        """Return a copy with every entry converted to `target` units."""
        return UncertaintyArray([entry.to(target, equivalencies=equivalencies) for entry in self])

    def __str__(self) -> str:
        return str([str(entry) for entry in self])

    def __repr__(self) -> str:
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
