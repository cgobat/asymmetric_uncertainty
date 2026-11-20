"""Regression tests for unit-aware uncertainty propagation and fitting."""

import unittest

import numpy as np
from astropy import units as u

from asymmetric_uncertainty import a_u, propagate


class TestUnitRegressions(unittest.TestCase):

    def test_self_correlation(self):
        x = a_u(3 * u.m, 0.2 * u.m, 0.4 * u.m)

        self.assertEqual(x + x, x * 2)
        self.assertEqual((x - x).items(), (0.0, 0.0, 0.0))
        self.assertEqual((x - x).unit, u.m)
        self.assertEqual(x * x, x**2)
        self.assertEqual((x / x).items(), (1.0, 0.0, 0.0))
        self.assertEqual((x / x).unit, u.dimensionless_unscaled)

        y = a_u(3.0, 0.2, 0.4)
        expected = y.value**y.value * (np.log(y.value) + 1)
        self.assertAlmostEqual((y**y).plus, expected * y.plus)
        self.assertAlmostEqual((y**y).minus, expected * y.minus)

        zero = a_u(0 * u.m)
        with self.assertRaises(ZeroDivisionError):
            _ = zero / zero

    def test_cdf_one_sided_and_exact(self):
        exact = a_u(2 * u.m)
        np.testing.assert_array_equal(exact.cdf([1, 2, 3] * u.m), [0, 1, 1])
        self.assertEqual(a_u(2 * u.m, 0 * u.m, 1 * u.m).cdf(2 * u.m), 1)
        self.assertEqual(a_u(2 * u.m, 1 * u.m, 0 * u.m).cdf(1 * u.m), 0)

    def test_sampling_fitting_and_propagation_preserve_units(self):
        x = a_u(3 * u.m, 0.2 * u.m, 0.4 * u.m)
        samples = x.rvs(size=20000, rng=123)
        self.assertEqual(samples.unit, u.m)
        self.assertEqual(x.rvs(rng=123).unit, u.m)
        self.assertEqual(a_u(2 * u.cm).rvs(size=3).unit, u.cm)

        fitted = a_u.fit(samples)
        self.assertEqual(fitted.unit, u.m)
        self.assertAlmostEqual(fitted.value, x.value, delta=0.04)
        self.assertAlmostEqual(fitted.plus, x.plus, delta=0.04)
        self.assertAlmostEqual(fitted.minus, x.minus, delta=0.04)

        exact = a_u.fit([3, 3, 3] * u.cm)
        self.assertEqual(exact.unit, u.cm)
        self.assertEqual(exact.items(), (3.0, 0.0, 0.0))

        y = a_u(50 * u.cm, 10 * u.cm, 20 * u.cm)
        raw = propagate(lambda a, b: a + b, x, y, size=10000, rng=123)
        self.assertEqual(raw.unit, u.m)
        result = propagate(lambda a, b: a + b, x, y, size=10000, rng=123, fit=True)
        self.assertEqual(result.unit, u.m)
        np.testing.assert_allclose(result.items(), a_u.fit(raw).items())

        shared = propagate(lambda a, b: a - b, x, x, size=100, rng=123, fit=True)
        self.assertEqual(shared.unit, u.m)
        self.assertEqual(shared.items(), (0.0, 0.0, 0.0))

    def test_unit_aware_distribution_properties(self):
        x = a_u(3 * u.m, 0.2 * u.m, 0.4 * u.m)
        bare = a_u(3, 0.2, 0.4)

        self.assertEqual(x.mean.unit, u.m)
        self.assertEqual(x.median.unit, u.m)
        self.assertEqual(x.variance.unit, u.m**2)
        self.assertEqual(x.std.unit, u.m)
        self.assertAlmostEqual(x.mean.to_value(u.m), bare.mean)
        self.assertAlmostEqual(x.variance.to_value(u.m**2), bare.variance)
        self.assertAlmostEqual(x.skewness, bare.skewness)
        self.assertAlmostEqual(x.kurtosis, bare.kurtosis)
        self.assertEqual(x.ppf(0.5).unit, u.m)
        self.assertEqual(x.ppf([0.1, 0.9]).unit, u.m)
        self.assertEqual(x.interval()[0].unit, u.m)
        self.assertEqual(x.equal_tail_interval()[0].unit, u.m)
        self.assertEqual(eval(repr(x), {"a_u": a_u}), x)

    def test_nonlinear_spectral_equivalency(self):
        wavelength = a_u(500 * u.nm, 10 * u.nm, 20 * u.nm)
        frequency = wavelength.to(u.THz, equivalencies=u.spectral())
        self.assertEqual(frequency.unit, u.THz)
        expected_mode = (500 * u.nm).to_value(u.THz, equivalencies=u.spectral())
        expected_upper = (480 * u.nm).to_value(u.THz, equivalencies=u.spectral())
        expected_lower = (510 * u.nm).to_value(u.THz, equivalencies=u.spectral())
        self.assertAlmostEqual(frequency.value, expected_mode)
        self.assertAlmostEqual(frequency.plus, expected_upper - expected_mode)
        self.assertAlmostEqual(frequency.minus, expected_mode - expected_lower)


if __name__ == "__main__":
    unittest.main()
