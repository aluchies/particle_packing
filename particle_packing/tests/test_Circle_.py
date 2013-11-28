import unittest
import numpy as np
from particle_packing.circle import Circle

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0 ], dtype=float)
        radius = 1
        phi = 0.

        c = Circle(center, radius, phi)
        self.assertTrue(c.radius, radius)
        self.assertTrue(np.allclose(c.center, center))
        self.assertTrue(c.phi == phi)


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = [0, 0, 0]
        radius = 1
        phi = 0.
        self.assertRaises(ValueError, Circle, center, radius, phi)

        center = np.asarray([0,0])
        radius = 'a'
        phi = 0.
        self.assertRaises(ValueError, Circle, center, radius, phi)

        center = np.asarray([0,0])
        radius = 1
        phi = 'a'
        self.assertRaises(ValueError, Circle, center, radius, phi)









    def test1_generate_volume(self):
        """

        Sphere at origin with radius 1

        """


        center = np.array([0., 0.])
        radius = 1
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Circle(center, radius, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 1, 0],
                        [1, 1, 1],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))




    def test3_generate_volume(self):
        """

        Sphere at origin, radius 1, rotated 90

        """


        center = np.array([0., 0.])
        radius = 1
        phi = np.pi / 2.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Circle(center, radius, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 1, 0],
                        [1, 1, 1],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))





    def test5_generate_volume(self):
        """

        Sphere not at origin, radius 1

        """


        center = np.array([0., 0.5])
        radius = 1.
        phi = 0.
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Circle(center, radius, phi)

        vol = c.generate_volume(x, y)

        arr = np.array([[0, 0, 0],
                        [0, 1, 0],
                        [0, 1, 0]])

        self.assertTrue(np.allclose(arr, vol))






    def test1_find_circle_subvolume(self):
        """

        Test find_subvolume() method for Circle class.

        """


        radius = 0.5
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol = \
        c.find_subvolume(x_ax, y_ax )

        self.assertTrue(np.allclose(x_ax, x_ax_subvol))
        self.assertTrue(np.allclose(y_ax, y_ax_subvol))


    def test2_find_circle_subvolume(self):
        """

        Test find_subvolume() method for Circle class

        """


        radius = 0.5
        center = 2. * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol = \
        c.find_subvolume(x_ax, y_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol))
        self.assertTrue(np.allclose(arr, y_ax_subvol))



    def test3_find_circle_subvolume(self):
        """

        Test find_subvolume() method for Circle class

        """


        radius = 0.4
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol= \
        c.find_subvolume(x_ax, y_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax_subvol))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax_subvol))







    def test1_find_circle_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Circle class.

        """


        radius = 0.5
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax )

        self.assertTrue(np.allclose(x_ax, x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax, y_ax[y_ax_subvol_ix]))


    def test2_find_circle_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Circle class

        """


        radius = 0.5
        center = 2. * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol_ix))
        self.assertTrue(np.allclose(arr, y_ax_subvol_ix))



    def test3_find_circle_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Circle class

        """


        radius = 0.4
        center = 0.5 * np.ones(2)
        phi = 0.
        c = Circle(center, radius, phi)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax[y_ax_subvol_ix]))




if __name__ == '__main__':
    print 'Running unit tests for Circle.py'
    unittest.main()