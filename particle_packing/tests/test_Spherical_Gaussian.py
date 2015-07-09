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
        self.assertRaises(ValueError, Spherical_Gaussian, center, radius_eff)

        center = np.asarray([0,0,0])
        radius = 'a'
        self.assertRaises(ValueError, Spherical_Gaussian, center, radius_eff)






    def test1_generate_volume_spherical_gaussian(self):
        """

        Test generate_volume_spherical_gaussian() method for Sphere class

        """

        radius = 0.2
        center = 0.5 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)
        y_ax = np.linspace(0, 1, 3)
        z_ax = np.linspace(0, 1, 3)

        c = Spherical_Gaussian(center, radius)
        subvol = c.generate_volume_(x_ax, y_ax, z_ax)

        arr = np.zeros((3,3,3))
        arr[0,1,1] = 1
        arr[1,1,1] = 1
        arr[1,1,0] = 1
        arr[1,0,1] = 1
        arr[1,1,2] = 1
        arr[1,2,1] = 1
        arr[2,1,1] = 1
        self.assertTrue(np.allclose(arr, subvol))