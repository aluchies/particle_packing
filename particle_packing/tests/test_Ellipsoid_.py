import unittest
import numpy as np
from particle_packing.ellipsoid import Ellipsoid

class TestCode(unittest.TestCase):




    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0, 0], dtype=float)
        radii = np.ones(3)
        alpha, beta, gamma = 0, 0, 0
        c = Ellipsoid(center, radii, alpha, beta, gamma)
        self.assertTrue(np.allclose(c.radii, radii))
        self.assertTrue(np.allclose(c.center, center))
        self.assertTrue(c.alpha == alpha)
        self.assertTrue(c.beta == beta)
        self.assertTrue(c.gamma == gamma)




    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = 0
        radii = [1, 1, 1]
        alpha, beta, gamma = 0, 0, 0
        self.assertRaises(ValueError, Ellipsoid, center, radii,
            alpha, beta, gamma)

        center = [0, 0, 0]
        radii = 'a'
        alpha, beta, gamma = 0, 0, 0
        self.assertRaises(ValueError, Ellipsoid, center, radii,
            alpha, beta, gamma)










    def test1_gen_ellipsoid(self):
        """Sphere at origin
        """
        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        alpha = 0.
        beta = 0.
        gamma = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)
        vol = gen_ellipsoid(x, y, z, center, radii, alpha, beta, gamma)
        arr = np.array([[[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]]])
        npt.assert_array_equal(arr, vol)

    def test2_gen_ellipsoid(self):
        """Sphere that is narrow along z-axis
        """
        a = 1.
        b = 1.
        c = 0.5
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        alpha = 0.
        beta = 0.
        gamma = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)
        vol = gen_ellipsoid(x, y, z, center, radii, alpha, beta, gamma)
        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]]])
        npt.assert_array_equal(arr, vol)

    def test3_gen_ellipsoid(self):
        """Sphere offset along z-axis
        """
        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 0., 1.])
        radii = np.array([a, b, c])
        alpha = 0.
        beta = 0.
        gamma = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)
        vol = gen_ellipsoid(x, y, z, center, radii, alpha, beta, gamma)
        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]]])
        npt.assert_array_equal(arr, vol)

    def test4_gen_ellipsoid(self):
        """Sphere offset along y-axis
        """
        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 1., 0.])
        radii = np.array([a, b, c])
        alpha = 0.
        beta = 0.
        gamma = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)
        vol = gen_ellipsoid(x, y, z, center, radii, alpha, beta, gamma)
        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [1, 1, 1]],

                        [[0, 0, 0],
                         [0, 0, 0],
                         [0, 1, 0]]])
        npt.assert_array_equal(arr, vol)

    def test5_gen_ellipsoid(self):
        """Rotated ellipsoid
        """
        a = 0.5
        b = np.sqrt(2)
        c = 1
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        alpha = np.pi / 4.
        beta = 0.
        gamma = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)
        vol = gen_ellipsoid(x, y, z, center, radii, alpha, beta, gamma)
        arr = np.array([[[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]]])
        npt.assert_array_equal(arr, vol)


if __name__ == '__main__':
    print 'Running unit tests for Ellipsoid.py'
    unittest.main()