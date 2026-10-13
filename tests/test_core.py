"""Asymmetric Uncertainty: A package for handling non-standard numerical uncertainties.

    Copyright (C) 2022 Caden Gobat

    This program is free software: you can redistribute and/or modify it under
    the terms of the GNU General Public License (either version 3 or, at your
    option, any later version).

    This program is distributed in the hope that it will be useful, but comes
    with ABSOLUTELY NO WARRANTY. See <https://www.gnu.org/licenses/> or the
    LICENSE file distributed with this program for more details.
"""

from math import erf
import unittest
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from asymmetric_uncertainty import UncertaintyArray, a_u, propagate

if hasattr(np, "trapz") and not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz


class TestCore(unittest.TestCase):

    @staticmethod
    def _distribution_moments(value):
        width_sum = value.plus + value.minus
        width_diff = value.plus - value.minus
        mean = value.value + np.sqrt(2/np.pi)*width_diff
        variance = width_sum**2/4 + (3/4 - 2/np.pi)*width_diff**2
        third = (
            np.sqrt(2/np.pi)
            * width_diff
            * (variance + (6/np.pi - 2)*width_diff**2)
        )
        return mean, variance, third

    def setUp(self): #We are only going to load the workbook once for each test case
        self.testValue = a_u(3.0, 0.2, 0.4)

        self.testValue_a = a_u(10.0, 1.0, 2.0)
        self.testValue_b = a_u(20.0, 1.0, 2.0)

        self.testValue_one = a_u(1.0, 0.0, 0.0)

        self.testValue_d = a_u(1.0, 1.0, 2.0)

    def test_Initialization(self):

        self.assertIsInstance(self.testValue, a_u)

        self.assertAlmostEqual(self.testValue.value, 3.0)
        self.assertAlmostEqual(self.testValue.plus, 0.2)
        self.assertAlmostEqual(self.testValue.minus, 0.4)

        self.assertAlmostEqual(self.testValue.maximum, 3.2)
        self.assertAlmostEqual(self.testValue.minimum, 2.6)


        self.assertFalse(self.testValue.is_symmetric)
        self.assertIs(self.testValue.sign, 1)

        self.assertAlmostEqual(self.testValue.exp().value, np.exp(3.0)) 

        self.assertAlmostEqual(self.testValue.log10().value, np.log10(3.0)) 

        self.assertAlmostEqual(self.testValue.log().value, np.log(3.0)) 

        self.assertAlmostEqual(self.testValue.sqrt().value, np.sqrt(3.0)) 

    def test_Native_Scalar_Types(self):

        value = a_u(3.0, 0.2, 0.4)
        other = a_u(2.5, 0.3, 0.3)

        self.assertIs(type(value.value), float)
        self.assertIs(type(value.plus), float)
        self.assertIs(type(value.minus), float)
        self.assertIs(type(value.maximum), float)
        self.assertIs(type(value.minimum), float)
        self.assertIs(type(value.is_symmetric), bool)
        self.assertIs(type(value > other), bool)
        self.assertIs(type(value >> other), bool)

        value.add_error(0.1, inplace=True)
        self.assertIs(type(value.plus), float)
        self.assertIs(type(value.minus), float)
        self.assertIs(type(value >> other), bool)

    def test_Distribution_Moments(self):

        mean, variance, third = self._distribution_moments(self.testValue)

        self.assertAlmostEqual(self.testValue.mean, mean)
        self.assertAlmostEqual(self.testValue.median, self.testValue.ppf(0.5))
        self.assertAlmostEqual(self.testValue.variance, variance)
        self.assertAlmostEqual(self.testValue.std, np.sqrt(variance))
        self.assertAlmostEqual(self.testValue.skewness, third/variance**1.5)

        samples = self.testValue.rvs(size=500000, rng=12345)
        sample_mean = np.mean(samples)
        sample_variance = np.var(samples)
        sample_kurtosis = np.mean((samples - sample_mean)**4)/sample_variance**2 - 3
        self.assertAlmostEqual(self.testValue.kurtosis, sample_kurtosis, delta=0.03)

        symmetric = a_u(3.0, 0.4, 0.4)
        self.assertAlmostEqual(symmetric.mean, symmetric.value)
        self.assertAlmostEqual(symmetric.variance, 0.4**2)
        self.assertAlmostEqual(symmetric.std, 0.4)
        self.assertAlmostEqual(symmetric.skewness, 0.0)
        self.assertAlmostEqual(symmetric.kurtosis, 0.0)

        exact = a_u(3.0)
        self.assertEqual(exact.variance, 0.0)
        self.assertEqual(exact.median, 3.0)
        self.assertTrue(np.isnan(exact.skewness))
        self.assertTrue(np.isnan(exact.kurtosis))

    def test_Distribution(self):

        width_sum = self.testValue.plus + self.testValue.minus

        self.assertAlmostEqual(
            self.testValue.cdf(self.testValue.value),
            self.testValue.minus/width_sum,
        )
        self.assertAlmostEqual(self.testValue.cdf(-np.inf), 0.0)
        self.assertAlmostEqual(self.testValue.cdf(np.inf), 1.0)

        neg_x = np.linspace(
            self.testValue.value - 10*self.testValue.minus,
            self.testValue.value,
            50001,
        )
        pos_x = np.linspace(
            self.testValue.value,
            self.testValue.value + 10*self.testValue.plus,
            50001,
        )
        x = np.concatenate((neg_x[:-1], pos_x))
        pdf = self.testValue.pdf(x)
        norm = np.trapezoid(pdf, x)
        self.assertAlmostEqual(norm, 1.0, places=6)

        numerical_mean = np.trapezoid(x*pdf, x)/norm
        numerical_variance = np.trapezoid((x - numerical_mean)**2*pdf, x)/norm
        numerical_third = np.trapezoid((x - numerical_mean)**3*pdf, x)/norm
        self.assertAlmostEqual(self.testValue.mean, numerical_mean, places=8)
        self.assertAlmostEqual(self.testValue.variance, numerical_variance, places=8)
        self.assertAlmostEqual(
            self.testValue.skewness,
            numerical_third/numerical_variance**1.5,
            places=8,
        )

        self.assertAlmostEqual(
            np.trapezoid(self.testValue.pdf(neg_x), neg_x),
            self.testValue.cdf(self.testValue.value),
            places=6,
        )

        exact = a_u(2.0)
        np.testing.assert_array_equal(exact.cdf([1.0, 2.0, 3.0]), [0.0, 1.0, 1.0])

        right_sided = a_u(2.0, 1.0, 0.0)
        left_sided = a_u(2.0, 0.0, 1.0)
        quantiles = np.array([0.1, 0.5, 0.9])
        for one_sided in (right_sided, left_sided):
            np.testing.assert_allclose(
                one_sided.cdf(one_sided.ppf(quantiles)),
                quantiles,
                rtol=0,
                atol=1e-12,
            )
        self.assertEqual(right_sided.cdf(right_sided.ppf(0.0)), 0.0)
        self.assertEqual(left_sided.cdf(left_sided.ppf(1.0)), 1.0)

    def test_Plotting(self):
        for plot in (self.testValue.pdfplot, self.testValue.cdfplot):
            plot(discretization=10)
            figure = plt.gcf()
            self.assertEqual(len(figure.axes), 1)
            self.assertEqual(len(figure.axes[0].lines), 1)
            plt.close(figure)

            figure, ax = plt.subplots()
            plot(discretization=10, ax=ax)
            self.assertEqual(len(ax.lines), 1)
            plt.close(figure)

    def test_PPF(self):

        mode_probability = self.testValue.minus/(self.testValue.plus + self.testValue.minus)
        quantiles = np.array([0.01, 0.1, mode_probability, 0.9, 0.99])
        values = self.testValue.ppf(quantiles)

        np.testing.assert_allclose(self.testValue.cdf(values), quantiles, rtol=0, atol=1e-12)
        self.assertAlmostEqual(self.testValue.ppf(mode_probability), self.testValue.value)
        self.assertEqual(self.testValue.ppf(0), -np.inf)
        self.assertEqual(self.testValue.ppf(1), np.inf)
        self.assertTrue(np.isnan(self.testValue.ppf(np.nan)))

        with self.assertRaises(ValueError):
            self.testValue.ppf(-0.1)
        with self.assertRaises(ValueError):
            self.testValue.ppf(1.1)

        exact = a_u(2.0)
        self.assertEqual(exact.ppf(0.5), 2.0)

        near_one = np.nextafter(1.0, 0.0)
        upper_tail = a_u(0.0, 3.0, 1.0)
        self.assertTrue(np.isfinite(upper_tail.ppf(near_one)))

    def test_Random_Sampling(self):

        samples = self.testValue.rvs(size=50000, rng=12345)
        repeated = self.testValue.rvs(size=50000, rng=12345)

        np.testing.assert_array_equal(samples, repeated)
        self.assertEqual(samples.shape, (50000,))
        self.assertAlmostEqual(np.mean(samples), self.testValue.mean, delta=0.01)
        self.assertAlmostEqual(
            np.mean(samples < self.testValue.value),
            self.testValue.minus/(self.testValue.plus + self.testValue.minus),
            delta=0.01,
        )

        exact = a_u(2.0)
        np.testing.assert_array_equal(exact.rvs(size=3, rng=123), np.full(3, 2.0))
        self.assertEqual(exact.rvs(rng=123), 2.0)

    def test_Intervals(self):

        self.assertEqual(self.testValue.interval(), (2.6, 3.2))
        self.assertEqual(self.testValue.interval(n_sigma=2), (2.2, 3.4))

        one_sigma_confidence = erf(1/np.sqrt(2))
        confidence_interval = self.testValue.interval(confidence=one_sigma_confidence)
        np.testing.assert_allclose(confidence_interval, self.testValue.interval(), rtol=0, atol=1e-12)

        equal_tail = self.testValue.equal_tail_interval(one_sigma_confidence)
        lower_probability = (1 - one_sigma_confidence)/2
        self.assertAlmostEqual(self.testValue.cdf(equal_tail[0]), lower_probability)
        self.assertAlmostEqual(self.testValue.cdf(equal_tail[1]), 1 - lower_probability)
        self.assertNotEqual(equal_tail, self.testValue.interval())

        symmetric = a_u(3.0, 0.4, 0.4)
        np.testing.assert_allclose(
            symmetric.equal_tail_interval(one_sigma_confidence),
            symmetric.interval(),
            rtol=0,
            atol=1e-12,
        )

        near_one = np.nextafter(1.0, 0.0)
        self.assertTrue(np.all(np.isfinite(self.testValue.interval(confidence=near_one))))
        self.assertTrue(
            np.all(np.isfinite(self.testValue.equal_tail_interval(near_one)))
        )

        with self.assertRaises(ValueError):
            self.testValue.interval(n_sigma=1, confidence=0.68)
        with self.assertRaises(ValueError):
            self.testValue.interval(n_sigma=-1)
        with self.assertRaises(ValueError):
            self.testValue.interval(confidence=1)
        with self.assertRaises(ValueError):
            self.testValue.equal_tail_interval(0)

    def test_Addition(self):

        testValue_c : a_u = self.testValue_a + self.testValue_b

        self.assertAlmostEqual(testValue_c.value, 30.0)
        self.assertAlmostEqual(testValue_c.plus, np.sqrt(2.0))
        self.assertAlmostEqual(testValue_c.minus, np.sqrt(8.0))

        shifted = self.testValue_a + 2
        self.assertAlmostEqual(shifted.value, 12.0)
        self.assertAlmostEqual(shifted.plus, self.testValue_a.plus)
        self.assertAlmostEqual(shifted.minus, self.testValue_a.minus)

        extreme = a_u(10.0, 1e-12, 2.0)
        shifted = extreme + 2
        self.assertEqual(shifted.value, 12.0)
        self.assertEqual(shifted.plus, extreme.plus)
        self.assertEqual(shifted.minus, extreme.minus)

        shifted = a_u(2.0) + extreme
        self.assertEqual(shifted.value, 12.0)
        self.assertEqual(shifted.plus, extreme.plus)
        self.assertEqual(shifted.minus, extreme.minus)

    def test_Addition_Associativity(self):

        first = a_u(1.0, 0.1, 0.4)
        second = a_u(2.0, 0.3, 0.2)
        third = a_u(-0.5, 0.2, 0.5)

        left = (first + second) + third
        right = first + (second + third)

        self.assertAlmostEqual(left.value, right.value)
        self.assertAlmostEqual(left.plus, right.plus)
        self.assertAlmostEqual(left.minus, right.minus)

    def test_Self_Correlation(self):

        value = a_u(3.0, 0.2, 0.4)

        summed = value + value
        scaled = value * 2
        self.assertEqual(summed.items(), scaled.items())
        self.assertAlmostEqual(summed.variance, 4*value.variance)

        independent = value + a_u(*value.items())
        self.assertAlmostEqual(independent.variance, 2*value.variance)

        difference = value - value
        self.assertEqual(difference.items(), (0.0, 0.0, 0.0))

        product = value * value
        squared = value**2
        self.assertEqual(product.items(), squared.items())

        ratio = value / value
        self.assertEqual(ratio.items(), (1.0, 0.0, 0.0))

        zero = a_u(0.0)
        with self.assertRaises(ZeroDivisionError):
            zero / zero

        power = value**value
        derivative = value.value**value.value * (np.log(value.value) + 1)

        self.assertAlmostEqual(power.value, value.value**value.value)
        self.assertAlmostEqual(power.plus, derivative*value.plus)
        self.assertAlmostEqual(power.minus, derivative*value.minus)

    def test_Monte_Carlo_Propagation(self):

        value = a_u(3.0, 0.2, 0.4)

        samples = propagate(lambda x: x**2, value, size=5000, rng=12345)
        repeated = propagate(lambda x: x**2, value, size=5000, rng=12345)
        np.testing.assert_array_equal(samples, repeated)
        self.assertEqual(samples.shape, (5000,))

        shared = propagate(lambda x, y: x - y, value, value, size=1000, rng=12345)
        np.testing.assert_array_equal(shared, np.zeros(1000))

        independent = propagate(
            lambda x, y: x - y,
            value,
            a_u(*value.items()),
            size=1000,
            rng=12345,
        )
        self.assertGreater(np.std(independent), 0)

        with self.assertRaises(ValueError):
            propagate(lambda x: x, value, size=0)

    def test_Subtraction(self):

        testValue_c : a_u = self.testValue_a - self.testValue_b

        self.assertAlmostEqual(testValue_c.value, -10.0)
        self.assertAlmostEqual(testValue_c.plus, np.sqrt(5.0))
        self.assertAlmostEqual(testValue_c.minus, np.sqrt(5.0))

        shifted = 2 - self.testValue_a
        self.assertAlmostEqual(shifted.value, -8.0)
        self.assertAlmostEqual(shifted.plus, self.testValue_a.minus)
        self.assertAlmostEqual(shifted.minus, self.testValue_a.plus)


    def test_Multiply(self):

        testValue_c : a_u = self.testValue_a * self.testValue_one

        self.assertAlmostEqual(testValue_c.plus, self.testValue_a.plus) 
        self.assertAlmostEqual(testValue_c.minus, self.testValue_a.minus)  

    def test_Multiply_Zero(self):

        zero_value = a_u(0.0, 0.1, 0.2)

        result = zero_value * 2
        self.assertAlmostEqual(result.value, 0.0)
        self.assertAlmostEqual(result.plus, 0.2)
        self.assertAlmostEqual(result.minus, 0.4)

        result = self.testValue_a * 0
        self.assertAlmostEqual(result.value, 0.0)
        self.assertAlmostEqual(result.plus, 0.0)
        self.assertAlmostEqual(result.minus, 0.0)

    def test_Multiply_Negative(self):

        result = self.testValue_a * -3
        self.assertAlmostEqual(result.value, -30.0)
        self.assertAlmostEqual(result.plus, 6.0)
        self.assertAlmostEqual(result.minus, 3.0)

        result = self.testValue_a * a_u(-3.0, 0.5, 1.0)
        self.assertAlmostEqual(result.plus, np.sqrt(61.0))
        self.assertAlmostEqual(result.minus, np.sqrt(109.0))

    def test_Divison_OverOne(self):

        testValue_c : a_u = self.testValue_a / self.testValue_one

        self.assertAlmostEqual(testValue_c.plus, self.testValue_a.plus) 
        self.assertAlmostEqual(testValue_c.minus, self.testValue_a.minus)  

    def test_Division_Zero_Numerator(self):

        zero_value = a_u(0.0, 0.1, 0.2)
        result = zero_value / 2

        self.assertAlmostEqual(result.value, 0.0)
        self.assertAlmostEqual(result.plus, 0.05)
        self.assertAlmostEqual(result.minus, 0.1)

        result = 0 / self.testValue_a
        self.assertAlmostEqual(result.value, 0.0)
        self.assertAlmostEqual(result.plus, 0.0)
        self.assertAlmostEqual(result.minus, 0.0)

    def test_Division_Negative(self):

        result = self.testValue_a / -2
        self.assertAlmostEqual(result.value, -5.0)
        self.assertAlmostEqual(result.plus, 1.0)
        self.assertAlmostEqual(result.minus, 0.5)

        result = -20 / self.testValue_a
        self.assertAlmostEqual(result.value, -2.0)
        self.assertAlmostEqual(result.plus, 0.2)
        self.assertAlmostEqual(result.minus, 0.4)

    def test_Divison_UnderOne(self):
        testValue_e : a_u = self.testValue_one/self.testValue_d

        self.assertAlmostEqual(testValue_e.value, 1.0) 

        self.assertAlmostEqual(testValue_e.plus, self.testValue_a.minus) 
        self.assertAlmostEqual(testValue_e.minus, self.testValue_a.plus)  

    def test_Pow(self):
        d_squared : a_u = self.testValue_d**2

        d_times_d : a_u = self.testValue_d * self.testValue_d

        self.assertAlmostEqual(d_squared.value, 1.0) 
        self.assertAlmostEqual(d_squared.plus, 2.0) 
        self.assertAlmostEqual(d_squared.minus, 4.0)

        self.assertAlmostEqual(d_squared.value, d_times_d.value)
        #self.assertAlmostEqual(d_squared.plus, d_times_d.plus) 
        #self.assertAlmostEqual(d_squared.minus, d_times_d.minus)  

        minus_two_power: a_u = self.testValue**(-2)
        reciprocal_squared: a_u = 1/(self.testValue**2)
        self.assertAlmostEqual(minus_two_power.value, reciprocal_squared.value)
        self.assertAlmostEqual(minus_two_power.plus, reciprocal_squared.plus)
        self.assertAlmostEqual(minus_two_power.minus, reciprocal_squared.minus)

    def test_Pow_Uncertain_Exponent(self):

        base = a_u(2.0, 0.1, 0.2)
        exponent = a_u(2.0, 0.1, 0.1)
        result = base**exponent

        exponent_term = 4*np.log(2)*0.1
        self.assertAlmostEqual(result.value, 4.0)
        self.assertAlmostEqual(result.plus, np.sqrt(0.4**2 + exponent_term**2))
        self.assertAlmostEqual(result.minus, np.sqrt(0.8**2 + exponent_term**2))

        reverse_result = 2**exponent
        self.assertAlmostEqual(reverse_result.value, 4.0)
        self.assertAlmostEqual(reverse_result.plus, exponent_term)
        self.assertAlmostEqual(reverse_result.minus, exponent_term)

    def test_Pow_Negative_Base(self):

        result = a_u(-2.0, 0.1, 0.2)**2
        self.assertAlmostEqual(result.value, 4.0)
        self.assertAlmostEqual(result.plus, 0.8)
        self.assertAlmostEqual(result.minus, 0.4)

        with self.assertRaises(ValueError):
            a_u(-2.0, 0.1, 0.2)**0.5

        with self.assertRaises(ValueError):
            a_u(-2.0, 0.1, 0.2)**a_u(2.0, 0.1, 0.1)

    def test_UncertaintyArray_List_Storage(self):

        first = a_u(1.0, 0.1, 0.2)
        second = a_u(2.0, 0.2, 0.3)
        array = UncertaintyArray([first, second])

        self.assertEqual(list.__len__(array), 2)
        self.assertEqual(array, [first, second])
        self.assertEqual(array.index(second), 1)
        self.assertIs(array.as_list, array)

        popped = array.pop()
        self.assertIs(popped, second)
        self.assertEqual(len(array), 1)
        self.assertEqual(array, [first])

    def test_Derived_State_After_Mutation(self):

        value = a_u(3.0, 0.2, 0.4)
        value.add_error(0.3, method="straight", inplace=True)

        self.assertAlmostEqual(value.maximum, 3.5)
        self.assertAlmostEqual(value.minimum, 2.3)

        value.minus = 0.5
        self.assertTrue(value.is_symmetric)

        value.value = -3.0
        self.assertEqual(value.sign, -1)
        self.assertAlmostEqual(value.maximum, -2.5)
        self.assertAlmostEqual(value.minimum, -3.5)

    def test_UncertaintyArray_State_After_Mutation(self):

        array = UncertaintyArray([a_u(1.0, 0.1, 0.2)])

        array[0] = 2.0
        self.assertIsInstance(array[0], a_u)
        self.assertEqual(array.values, [2.0])
        self.assertEqual(array.plus, [0.0])
        self.assertEqual(array.minus, [0.0])

        array.extend([3.0, a_u(4.0, 0.4, 0.5)])
        self.assertEqual(array.values, [2.0, 3.0, 4.0])
        self.assertEqual(array.plus, [0.0, 0.0, 0.4])
        self.assertEqual(array.minus, [0.0, 0.0, 0.5])

        array.reverse()
        self.assertEqual(array.values, [4.0, 3.0, 2.0])

        array.pop()
        self.assertEqual(array.values, [4.0, 3.0])
        self.assertEqual(array.shape, (2,))
