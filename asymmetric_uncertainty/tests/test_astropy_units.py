"""Tests for interoperability with ``astropy.units``."""

import unittest

import numpy as np

try:
    from astropy import units as u
except ImportError:
    u = None

from asymmetric_uncertainty import UncertaintyArray, a_u


@unittest.skipIf(u is None, "Astropy is not installed")
class TestAstropyUnits(unittest.TestCase):

    def test_initialization_and_conversion(self):
        value = a_u(3*u.m, 20*u.cm, 40*u.cm)

        self.assertEqual(value.unit, u.m)
        self.assertAlmostEqual(value.value, 3.0)
        self.assertAlmostEqual(value.plus, 0.2)
        self.assertAlmostEqual(value.minus, 0.4)
        self.assertEqual(value.maximum, 3.2*u.m)
        self.assertEqual(value.minimum, 2.6*u.m)

        converted = value.to(u.cm)
        self.assertEqual(converted.unit, u.cm)
        self.assertAlmostEqual(converted.value, 300.0)
        self.assertAlmostEqual(converted.plus, 20.0)
        self.assertAlmostEqual(converted.minus, 40.0)

        explicit_unit = a_u(3.0, 0.2, 0.4, unit=u.m)
        self.assertEqual(explicit_unit, value)

    def test_addition_and_subtraction_convert_units(self):
        left = a_u(3*u.m, 20*u.cm, 40*u.cm)
        right = a_u(50*u.cm, 10*u.cm, 20*u.cm)

        added = left + right
        self.assertEqual(added.unit, u.m)
        self.assertAlmostEqual(added.value, 3.5)
        self.assertAlmostEqual(added.plus, np.sqrt(0.2**2 + 0.1**2))
        self.assertAlmostEqual(added.minus, np.sqrt(0.4**2 + 0.2**2))

        subtracted = left - right
        self.assertEqual(subtracted.unit, u.m)
        self.assertAlmostEqual(subtracted.value, 2.5)
        self.assertAlmostEqual(subtracted.plus, np.sqrt(0.2**2 + 0.2**2))
        self.assertAlmostEqual(subtracted.minus, np.sqrt(0.4**2 + 0.1**2))

        exact_added = left + 50*u.cm
        self.assertAlmostEqual(exact_added.value, 3.5)
        self.assertAlmostEqual(exact_added.plus, 0.2)
        self.assertAlmostEqual(exact_added.minus, 0.4)

        reverse_added = 50*u.cm + left
        self.assertEqual(reverse_added, exact_added)

        with self.assertRaises(u.UnitConversionError):
            _ = left + 1*u.s

        with self.assertRaises(u.UnitConversionError):
            _ = left + 1

    def test_multiplication_and_division_propagate_units(self):
        length = a_u(2*u.m, 0.1*u.m, 0.2*u.m)
        time = a_u(3*u.s, 0.3*u.s, 0.4*u.s)

        product = length*time
        self.assertEqual(product.unit, u.m*u.s)
        self.assertAlmostEqual(product.value, 6.0)
        self.assertAlmostEqual(product.plus, np.sqrt(0.3**2 + 0.6**2))
        self.assertAlmostEqual(product.minus, np.sqrt(0.6**2 + 0.8**2))

        scaled = length*(2*u.s)
        self.assertEqual(scaled.unit, u.m*u.s)
        self.assertAlmostEqual(scaled.plus, 0.2)
        self.assertAlmostEqual(scaled.minus, 0.4)

        reverse_scaled = (2*u.s)*length
        self.assertEqual(reverse_scaled, scaled)

        unit_scaled = length*u.s
        self.assertEqual(unit_scaled.unit, u.m*u.s)
        self.assertEqual(u.s*length, unit_scaled)

        numpy_scaled = np.multiply(2*u.s, length)
        self.assertEqual(numpy_scaled, scaled)

        quotient = length/time
        self.assertEqual(quotient.unit, u.m/u.s)
        self.assertAlmostEqual(quotient.value, 2/3)

        inverse = (6*u.m)/time
        self.assertEqual(inverse.unit, u.m/u.s)
        self.assertAlmostEqual(inverse.value, 2.0)
        self.assertAlmostEqual(inverse.plus, 6/9*0.4)
        self.assertAlmostEqual(inverse.minus, 6/9*0.3)

    def test_power_and_math_operations(self):
        area = a_u(4*u.m**2, 0.4*u.m**2, 0.8*u.m**2)
        length = area.sqrt()

        self.assertEqual(length.unit, u.m)
        self.assertAlmostEqual(length.value, 2.0)
        self.assertAlmostEqual(length.plus, 0.1)
        self.assertAlmostEqual(length.minus, 0.2)

        self.assertEqual(np.sqrt(area), length)

        with self.assertRaises(ValueError):
            _ = a_u(2*u.m, 0.1*u.m, 0.2*u.m)**a_u(2.0, 0.1, 0.1)

        with self.assertRaises(u.UnitConversionError):
            _ = a_u(2*u.m, 0.1*u.m, 0.2*u.m).log()

        dimensionless = a_u(2*u.dimensionless_unscaled, 0.1, 0.2)
        self.assertAlmostEqual(dimensionless.log().value, np.log(2.0))
        self.assertAlmostEqual(dimensionless.exp().value, np.exp(2.0))

    def test_trigonometric_ufuncs(self):
        angle = a_u(30*u.deg, 1*u.deg, 2*u.deg)

        sine = np.sin(angle)
        self.assertAlmostEqual(sine.value, 0.5)
        self.assertAlmostEqual(sine.plus, np.cos(np.pi/6)*np.deg2rad(1))
        self.assertAlmostEqual(sine.minus, np.cos(np.pi/6)*np.deg2rad(2))

        cosine = np.cos(angle)
        self.assertAlmostEqual(cosine.value, np.cos(np.pi/6))
        self.assertAlmostEqual(cosine.plus, np.sin(np.pi/6)*np.deg2rad(2))
        self.assertAlmostEqual(cosine.minus, np.sin(np.pi/6)*np.deg2rad(1))

    def test_distribution_accepts_quantities(self):
        value = a_u(3*u.m, 0.2*u.m, 0.4*u.m)

        self.assertAlmostEqual(value.cdf(300*u.cm), 2/3)
        density = value.pdf(3*u.m)
        self.assertTrue(density.unit.is_equivalent(1/u.m))
        self.assertAlmostEqual(
            density.to_value(1/u.m),
            np.sqrt(2)/np.sqrt(np.pi)/0.6,
        )

    def test_add_error_accepts_quantity(self):
        value = a_u(3*u.m, 0.2*u.m, 0.4*u.m)
        updated = value.add_error(10*u.cm, method="straight")

        self.assertEqual(updated.unit, u.m)
        self.assertAlmostEqual(updated.plus, 0.3)
        self.assertAlmostEqual(updated.minus, 0.5)

    def test_uncertainty_array_preserves_unit_aware_entries(self):
        values = UncertaintyArray([
            1*u.m,
            a_u(2*u.m, 0.1*u.m, 0.2*u.m),
        ])

        self.assertEqual(values.units, [u.m, u.m])
        self.assertEqual(values.quantities, [1*u.m, 2*u.m])

        converted = values.to(u.cm)
        self.assertEqual(converted.values, [100.0, 200.0])
        self.assertEqual(converted.plus, [0.0, 10.0])
        self.assertEqual(converted.minus, [0.0, 20.0])
        self.assertEqual(converted.units, [u.cm, u.cm])


if __name__ == "__main__":
    unittest.main()
