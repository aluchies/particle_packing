import unittest
import numpy as np
from particle_packing.ellipsoid import Ellipsoid

from particle_packing.ellipsoid.Ellipsoid \
    import _quaternion_to_rotation_matrix, rotation_matrix

from random import random

class TestCode(unittest.TestCase):




    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0, 0], dtype=float)
        radii = np.ones(3)
        rt_ax = [1., 0., 0.]
        phi = 0.
        c = Ellipsoid(center, radii, rt_ax, phi)
        self.assertTrue(np.allclose(c.radii, radii))
        self.assertTrue(np.allclose(c.center, center))
        self.assertTrue(c.rt_ax[0] == rt_ax[0])
        self.assertTrue(c.rt_ax[1] == rt_ax[1])
        self.assertTrue(c.rt_ax[2] == rt_ax[2])
        self.assertTrue(c.phi == phi )




    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = 0
        radii = [1, 1, 1]
        rt_ax = [1., 0., 0.]
        phi = 0.
        self.assertRaises(ValueError, Ellipsoid, center, radii,
            rt_ax, phi)


        center = [0, 0, 0]
        rt_ax = [1., 0., 0.]
        phi = 'a'
        self.assertRaises(ValueError, Ellipsoid, center, radii,
            rt_ax, phi)


        center = [0, 0, 0]
        rt_ax = [1., 1., 0.]
        phi = np.pi / 2.
        self.assertRaises(ValueError, Ellipsoid, center, radii,
            rt_ax, phi)











    def test1_generate_volume(self):
        """

        Sphere at origin

        """

        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        rt_ax = [1., 0., 0.]
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)

        c = Ellipsoid(center, radii, rt_ax, phi)
        vol = c.generate_volume(x, y, z)

        arr = np.array([[[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]]])

        self.assertTrue(np.allclose(arr, vol))


    def test2_generate_volume(self):
        """

        Sphere that is narrow along z-axis

        """
        
        a = 1.
        b = 1.
        c = 0.5
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        rt_ax = [1., 0., 0.]
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)

        c = Ellipsoid(center, radii, rt_ax, phi)
        vol = c.generate_volume(x, y, z)

        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]]])

        self.assertTrue(np.allclose(arr, vol))


    def test3_generate_volume(self):
        """

        Sphere offset along z-axis

        """

        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 0., 1.])
        radii = np.array([a, b, c])
        rt_ax = [1., 0., 0.]
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)

        c = Ellipsoid(center, radii, rt_ax, phi)
        vol = c.generate_volume(x, y, z)

        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]]])

        self.assertTrue(np.allclose(arr, vol))

    def test4_generate_volume(self):
        """

        Sphere offset along y-axis

        """

        a = 1.
        b = 1.
        c = 1.
        center = np.array([0., 1., 0.])
        radii = np.array([a, b, c])
        rt_ax = [1., 0., 0.]
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)

        c = Ellipsoid(center, radii, rt_ax, phi)
        vol = c.generate_volume(x, y, z)

        arr = np.array([[[0, 0, 0],
                         [0, 0, 0],
                         [0, 1, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [1, 1, 1]],

                        [[0, 0, 0],
                         [0, 0, 0],
                         [0, 1, 0]]])

        self.assertTrue(np.allclose(arr, vol))

    def test5_generate_volume(self):
        """

        Rotate ellipsoid

        """

        a = 0.5
        b = np.sqrt(2)
        c = 1.
        center = np.array([0., 0., 0.])
        radii = np.array([a, b, c])
        rt_ax = [0., 0., 1.]
        phi = np.pi / 4.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)
        z = np.linspace(-1., 1., 3)

        c = Ellipsoid(center, radii, rt_ax, phi)
        vol = c.generate_volume(x, y, z)

        arr = np.array([[[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]],

                        [[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]],

                        [[0, 0, 0],
                         [0, 1, 0],
                         [0, 0, 0]]])


        self.assertTrue(np.allclose(arr, vol))






    def test1__quaternion_to_rotation_matrix(self):
        """

        Test simple axial rotations

        """


        # x-axis
        for i in xrange(100):
            phi = np.pi * random()
            rot_ax = np.asarray([1., 0., 0.])
            R1 = _quaternion_to_rotation_matrix(rot_ax, phi)

            alpha, beta, gamma, = -phi, 0., 0.
            R2 = rotation_matrix(alpha, beta, gamma)

            self.assertTrue(np.allclose(R1, R2))


        # y-axis
        for i in xrange(100):
            phi = np.pi * random()
            rot_ax = np.asarray([0., 1., 0.])
            R1 = _quaternion_to_rotation_matrix(rot_ax, phi)

            alpha, beta, gamma, = 0., -phi, 0.
            R2 = rotation_matrix(alpha, beta, gamma)

            self.assertTrue(np.allclose(R1, R2))


        # z-axis
        for i in xrange(100):
            phi = np.pi * random()
            rot_ax = np.asarray([0., 0., 1.])
            R1 = _quaternion_to_rotation_matrix(rot_ax, phi)

            alpha, beta, gamma, = 0., 0., -phi
            R2 = rotation_matrix(alpha, beta, gamma)

            self.assertTrue(np.allclose(R1, R2))








if __name__ == '__main__':
    print 'Running unit tests for Ellipsoid.py'
    unittest.main()