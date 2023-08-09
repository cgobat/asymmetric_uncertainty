"""Asymmetric Uncertainty: A package for handling non-standard numerical uncertainties.

    Copyright (C) 2022 Caden Gobat

    This program is free software: you can redistribute and/or modify it under
    the terms of the GNU General Public License (either version 3 or, at your
    option, any later version).

    This program is distributed in the hope that it will be useful, but comes
    with ABSOLUTELY NO WARRANTY. See <https://www.gnu.org/licenses/> or the
    LICENSE file distributed with this program for more details.
"""

import unittest
import numpy as np
from asymmetric_uncertainty import a_u

class TestCore(unittest.TestCase):

    def setUp(self): #We are only going to load the workbook once for each test case
        self.testValue = a_u(3.0,0.2,0.4)

        self.testValue_a = a_u(10.0,1.0,2.0)
        self.testValue_b = a_u(20.0,1.0,2.0)

        self.testValue_one = a_u(1.0,0.0,0.0)

        self.testValue_d = a_u(1.0,1.0,2.0)

    def test_Initialization(self):

        self.assertIsInstance(self.testValue,a_u)

        self.assertAlmostEqual(self.testValue.value, 3.0)
        self.assertAlmostEqual(self.testValue.plus, 0.2)
        self.assertAlmostEqual(self.testValue.minus, 0.4)

        self.assertAlmostEqual(self.testValue.maximum, 3.2)
        self.assertAlmostEqual(self.testValue.minimum, 2.6)


        self.assertFalse(self.testValue.is_symmetric)
        self.assertIs(self.testValue.sign, +1)

        # exp

        self.assertAlmostEqual(self.testValue.exp().value, np.exp(3.0)) 

        self.assertAlmostEqual(self.testValue.log10().value, np.log10(3.0)) 

        self.assertAlmostEqual(self.testValue.log().value, np.log(3.0)) 

        self.assertAlmostEqual(self.testValue.sqrt().value, np.sqrt(3.0)) 

    def test_Summation(self):

        testValue_c : a_u = self.testValue_a + self.testValue_b

        self.assertAlmostEqual(testValue_c.plus, np.sqrt(2)) 
        self.assertAlmostEqual(testValue_c.minus, np.sqrt(8))       


    def test_Subtraction(self):

        testValue_c : a_u = self.testValue_a - self.testValue_b

        self.assertAlmostEqual(testValue_c.plus, np.sqrt(2)) 
        self.assertAlmostEqual(testValue_c.minus, np.sqrt(8))  


    def test_Multiply(self):

        testValue_c : a_u = self.testValue_a * self.testValue_one

        self.assertAlmostEqual(testValue_c.plus, self.testValue_a.plus) 
        self.assertAlmostEqual(testValue_c.minus, self.testValue_a.minus)  

    def test_Divison_OverOne(self):

        testValue_c : a_u = self.testValue_a / self.testValue_one

        self.assertAlmostEqual(testValue_c.plus, self.testValue_a.plus) 
        self.assertAlmostEqual(testValue_c.minus, self.testValue_a.minus)  

    def test_Divison_UnderOne(self):
        testValue_e : a_u = self.testValue_one/self.testValue_d

        self.assertAlmostEqual(testValue_e.value, 1.0) 

        self.assertAlmostEqual(testValue_e.plus, self.testValue_a.plus) 
        self.assertAlmostEqual(testValue_e.minus, self.testValue_a.minus)  

    def test_Exp(self):
        testValue_c : a_u = self.testValue_d**2

        #testValue_e : a_u = self.testValue_d * self.testValue_d

        self.assertAlmostEqual(testValue_c.value, 1.0) 

        #self.assertAlmostEqual(testValue_c.plus, testValue_e.plus) 
        #self.assertAlmostEqual(testValue_c.minus, testValue_e.minus)  


