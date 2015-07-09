import unittest
import numpy as np
from particle_packing.spherical_gaussian import Spherical_Gaussian



class TestCode(unittest.TestCase):


    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0, 0], dtype=float)
        radius = 1.
        c = Spherical_Gaussian(center, radius)
        self.assertTrue(c.radius_eff == radius)
        self.assertTrue(np.allclose(c.center, center))


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = 0
        radius = 'a'
        self.assertRaises(ValueError, Spherical_Gaussian, center, radius)

        center = np.asarray([0,0,0])
        radius = 'a'
        self.assertRaises(ValueError, Spherical_Gaussian, center, radius)





