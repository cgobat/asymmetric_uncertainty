"""First-order propagation through shared independent uncertainty sources."""

import copy
import pickle
import unittest

import numpy as np

from asymmetric_uncertainty import a_u, propagate


class TestCorrelations(unittest.TestCase):
    def assertQuantity(self, actual, expected):
        np.testing.assert_allclose(actual.items(), expected, rtol=1e-12, atol=1e-14)

    def test_derived_linear_expressions(self):
        x = a_u(3, 0.2, 0.4)
        for result in (2 * x - x, (x + 7) - 7, (x / 2) + (x / 2)):
            with self.subTest(result=result):
                self.assertQuantity(result, x.items())
        for result in ((x + 7) - x, x + (-x), (2 * x) - (x + x)):
            self.assertQuantity(result, (result.value, 0, 0))

    def test_asymmetric_sign_after_combining_derivatives(self):
        x = a_u(3, 0.2, 0.4)
        self.assertQuantity(x - 2 * x, (-3, 0.4, 0.2))
        self.assertQuantity((-2 * x) + x, (-3, 0.4, 0.2))
        self.assertQuantity(-(2 * x - x), (-3, 0.4, 0.2))
        self.assertQuantity((2 * x - x) * -3, (-9, 1.2, 0.6))

    def test_partial_overlap(self):
        x = a_u(3, 0.2, 0.4)
        z = a_u(5, 0.3, 0.7)
        self.assertQuantity((x + z) - x, z.items())
        self.assertQuantity((x - z) - x, (-5, 0.7, 0.3))
        self.assertQuantity((2 * x + z) - x,
                            (8, np.hypot(0.2, 0.3), np.hypot(0.4, 0.7)))
        self.assertQuantity((x + z) - (x - z), (10, 0.6, 1.4))

    def test_independent_operands_keep_directional_quadrature(self):
        x = a_u(-3, 0.2, 0.4)
        z = a_u(2, 0.3, 0.7)
        cases = (
            (x + z, -1, (1, 1)),
            (x - z, -5, (1, -1)),
            (x * z, -6, (2, -3)),
            (x / z, -1.5, (0.5, 0.75)),
            (2 / x, -2 / 3, (-2 / 9, 0)),
        )
        for actual, value, derivatives in cases:
            with self.subTest(value=value):
                plus = []
                minus = []
                for operand, derivative in zip((x, z), derivatives):
                    upper, lower = (operand.plus, operand.minus) if derivative >= 0 else (operand.minus, operand.plus)
                    plus.append((derivative * upper)**2)
                    minus.append((derivative * lower)**2)
                self.assertQuantity(actual, (value, np.sqrt(sum(plus)), np.sqrt(sum(minus))))
        a = a_u(3, 0.2, 0.4)
        b = a_u(3, 0.2, 0.4)
        self.assertQuantity(a - b, (0, np.hypot(0.2, 0.4), np.hypot(0.2, 0.4)))

    def test_nonlinear_shared_sources(self):
        x = a_u(3, 0.2, 0.4)
        self.assertQuantity((x * x) / x, x.items())
        self.assertQuantity((2 * x) / x, (2, 0, 0))
        self.assertQuantity((x + 0) * x, (9, 1.2, 2.4))
        self.assertQuantity(x / (2 * x), (0.5, 0, 0))
        self.assertQuantity(x / (x + 1), (0.75, 0.2 / 16, 0.4 / 16))
        self.assertQuantity(x ** (x + 0), (27, 27 * (np.log(3) + 1) * 0.2,
                                            27 * (np.log(3) + 1) * 0.4))
        self.assertQuantity((x ** x) - (x ** (x + 0)), (0, 0, 0))
        self.assertQuantity((x ** 2) ** 0.5, x.items())

    def test_independent_uncertain_powers(self):
        for base in (0.2, 3):
            x = a_u(base, 0.02, 0.04)
            z = a_u(2, 0.3, 0.7)
            value = base**2
            base_derivative = 2 * base
            exponent_derivative = value * np.log(base)
            upper, lower = (z.plus, z.minus) if exponent_derivative >= 0 else (z.minus, z.plus)
            self.assertQuantity(x**z, (value,
                                       np.hypot(base_derivative * x.plus, exponent_derivative * upper),
                                       np.hypot(base_derivative * x.minus, exponent_derivative * lower)))
            self.assertQuantity((x**z) - (x**(z + 0)), (0, 0, 0))
        x = a_u(0.1, 0.02, 0.04)
        derivative = x.value**x.value * (1 + np.log(x.value))
        self.assertQuantity(x ** (x + 0), (x.value**x.value,
                                          -derivative * x.minus, -derivative * x.plus))

    def test_unary_helpers_preserve_sources(self):
        x = a_u(3, 0.2, 0.4)
        for result in (x.log().exp(), 10 ** x.log10(), x.sqrt() ** 2):
            self.assertQuantity(result, x.items())
            self.assertQuantity(result - x, (0, 0, 0))

    def test_symmetric_and_one_sided_sources(self):
        for plus, minus in ((0.2, 0.2), (0, 0.4), (0.2, 0)):
            with self.subTest(plus=plus, minus=minus):
                x = a_u(3, plus, minus)
                self.assertQuantity(2 * x - x, x.items())
                self.assertQuantity(x - 2 * x, (-3, minus, plus))

    def test_public_assignment_starts_new_source(self):
        for attribute, replacement in (("value", 4), ("plus", 0.5), ("minus", 0.8)):
            with self.subTest(attribute=attribute):
                x = a_u(3, 0.2, 0.4)
                y = 2 * x
                setattr(x, attribute, replacement)
                self.assertQuantity(y - x, (6 - x.value,
                                            np.hypot(0.4, x.minus),
                                            np.hypot(0.8, x.plus)))
                self.assertQuantity(2 * x - x, x.items())
                # Restore the field before doing arithmetic: assignment still rebases.
                setattr(x, attribute, {"value": 3, "plus": 0.2, "minus": 0.4}[attribute])
                self.assertQuantity(y - x, (3, np.hypot(0.4, 0.4), np.hypot(0.8, 0.2)))
        x = a_u(3, 0.2, 0.4)
        y = 2 * x
        y.plus = 1
        self.assertQuantity(y - x, (3, np.hypot(1, 0.4), np.hypot(0.8, 0.2)))

    def test_add_error_is_explicit_source_boundary(self):
        for method in ("quadrature", "straight", "split"):
            for inplace in (False, True):
                with self.subTest(method=method, inplace=inplace):
                    x = a_u(3, 0.2, 0.4)
                    old = x + 0
                    adjusted = x.add_error(0.1, method=method, inplace=inplace)
                    if inplace:
                        adjusted = x
                    self.assertQuantity(adjusted - old,
                                        (0, np.hypot(adjusted.plus, old.minus),
                                         np.hypot(adjusted.minus, old.plus)))
                    if not inplace:
                        self.assertQuantity(x - old, (0, 0, 0))

    def test_copy_and_pickle_preserve_sources(self):
        x = a_u(3, 0.2, 0.4)
        for clone in (copy.copy(x), copy.deepcopy(x), pickle.loads(pickle.dumps(x))):
            self.assertQuantity(clone - x, (0, 0, 0))
            clone.plus = 1
            self.assertQuantity(x, (3, 0.2, 0.4))
        restored_x, restored_y = pickle.loads(pickle.dumps((x, 2 * x)))
        self.assertQuantity(restored_y - restored_x, x.items())
        self.assertQuantity(a_u(*x.items()) - x, (0, np.hypot(0.2, 0.4), np.hypot(0.2, 0.4)))

    def test_zero_derivative_remains_first_order(self):
        x = a_u(0, 0.2, 0.4)
        self.assertQuantity(x * x, (0, 0, 0))
        self.assertQuantity((x + 0) * x, (0, 0, 0))
        self.assertQuantity((x + 2) / (x + 2), (1, 0, 0))
        self.assertQuantity(x / x, (1, 0, 0))  # Existing identity special case.
        with self.assertRaises(ZeroDivisionError):
            x / (x + 0)

    def test_monte_carlo_keeps_existing_sampling_semantics(self):
        x = a_u(3, 0.2, 0.4)
        y = 2 * x
        actual = propagate(lambda a, b: a - b, y, x, size=20, rng=42)
        generator = np.random.default_rng(42)
        expected = y.rvs(size=20, rng=generator) - x.rvs(size=20, rng=generator)
        np.testing.assert_array_equal(actual, expected)
