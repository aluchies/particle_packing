import unittest
import numpy as np
from particle_packing.ellipse import Ellipse

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0 ], dtype=float)
        radii = [1, 1]
        phi = 0.

        c = Ellipse(center, radii, phi)
        self.assertTrue(np.allclose(c.radii, np.asarray(radii)))
        self.assertTrue(np.allclose(c.center, center))
        self.assertTrue(c.phi == phi)


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = [0, 0, 0]
        radii = [1,1]
        phi = 0.
        self.assertRaises(ValueError, Ellipse, center, radii, phi)

        center = np.asarray([0,0])
        radii = [1,1,1]
        phi = 0.
        self.assertRaises(ValueError, Ellipse, center, radii, phi)

        center = np.asarray([0,0])
        radii = [1,1,1]
        phi = 'a'
        self.assertRaises(ValueError, Ellipse, center, radii, phi)









    def test1_generate_volume(self):
        """

        Sphere at origin with radius 1

        """

        a = 1.
        b = 1.
        center = np.array([0., 0.])
        radii = np.array([a, b])
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Ellipse(center, radii, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 1, 0],
                        [1, 1, 1],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))


    def test2_generate_volume(self):
        """

        Ellipse, major y-axis, minor x-axis

        """

        a = 0.5
        b = 1.
        center = np.array([0., 0.])
        radii = np.array([a, b])
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Ellipse(center, radii, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 1, 0],
                        [0, 1, 0],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))



    def test3_generate_volume(self):
        """

        Sphere at origin, radius 1, rotated 90

        """

        a = 1.
        b = 1.
        center = np.array([0., 0.])
        radii = np.array([a, b])
        phi = np.pi / 2.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Ellipse(center, radii, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 1, 0],
                        [1, 1, 1],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))


    def test4_generate_volume(self):
        """

        Narrow ellipse rotated so 45

        """

        a = 0.5
        b = np.sqrt(2)
        center = np.array([0., 0.])
        radii = np.array([a, b])
        phi = np.pi / 4.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Ellipse(center, radii, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]])

        self.assertTrue(np.allclose(arr, vol))


    def test5_generate_volume(self):
        """

        Sphere not at origin, radius 1

        """

        a = 1.
        b = 1.
        center = np.array([0., 0.5])
        radii = np.array([a, b])
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Ellipse(center, radii, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 0, 0],
                        [0, 1, 0],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))






    def test1_find_ellipse_subvolume(self):
        """

        Test find_subvolume() method for Ellipse class.

        """


        radii = [0.5, 0.5]
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol = \
        c.find_subvolume(x_ax, y_ax )

        self.assertTrue(np.allclose(x_ax, x_ax_subvol))
        self.assertTrue(np.allclose(y_ax, y_ax_subvol))


    def test2_find_ellipse_subvolume(self):
        """

        Test find_subvolume() method for Ellipse class

        """


        radii = [0.5, 0.5]
        center = 2. * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol = \
        c.find_subvolume(x_ax, y_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol))
        self.assertTrue(np.allclose(arr, y_ax_subvol))



    def test3_find_ellipse_subvolume(self):
        """

        Test find_subvolume() method for Ellipse class

        """


        radii = 0.4
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol= \
        c.find_subvolume(x_ax, y_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax_subvol))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax_subvol))







    def test1_find_ellipse_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Ellipse class.

        """


        radii = [0.5, 0.5]
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax )

        self.assertTrue(np.allclose(x_ax, x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax, y_ax[y_ax_subvol_ix]))


    def test2_find_ellipse_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Ellipse class

        """


        radii = [0.5, 0.5]
        center = 2. * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol_ix))
        self.assertTrue(np.allclose(arr, y_ax_subvol_ix))



    def test3_find_ellipse_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Ellipse class

        """


        radii = 0.4
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Ellipse(center, radii, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax[y_ax_subvol_ix]))















    def test1_overlap_potential(self):
        """

        Test overlap_potential method for Ellipse class

        """

        center = np.zeros(2)
        radii = 0.5 * np.ones(2)
        phi = 0.
        c1 = Ellipse(center, radii, phi)
        c2 = Ellipse(center, radii, phi)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)




if __name__ == '__main__':
    print 'Running unit tests for Ellipse.py'
    unittest.main()