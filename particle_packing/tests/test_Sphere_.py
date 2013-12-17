import unittest
import numpy as np
from particle_packing.sphere import Sphere

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0, 0], dtype=float)
        radius = 1.
        c = Sphere(center, radius)
        self.assertTrue(c.radius == radius)
        self.assertTrue(np.allclose(c.center, center))


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = 0
        radius = 'a'
        self.assertRaises(ValueError, Sphere, center, radius)

        center = np.asarray([0,0,0])
        radius = 'a'
        self.assertRaises(ValueError, Sphere, center, radius)







    def test1_generate_volume(self):
        """

        Test generate_volume() method for Sphere class

        """

        radius = 0.5
        center = 0.5 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)
        y_ax = np.linspace(0, 1, 3)
        z_ax = np.linspace(0, 1, 3)

        c = Sphere(center, radius)
        subvol = c.generate_volume(x_ax, y_ax, z_ax)

        arr = np.zeros((3,3,3))
        arr[0,1,1] = 1
        arr[1,1,1] = 1
        arr[1,1,0] = 1
        arr[1,0,1] = 1
        arr[1,1,2] = 1
        arr[1,2,1] = 1
        arr[2,1,1] = 1
        self.assertTrue(np.allclose(arr, subvol))

    def test2_generate_sphere(self):
        """

        Test generate_volume() method for Sphere class

        """

        radius = 0.5
        center = 2 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)
        y_ax = np.linspace(0, 1, 3)
        z_ax = np.linspace(0, 1, 3)

        c = Sphere(center, radius)
        subvol = c.generate_volume(x_ax, y_ax, z_ax)

        self.assertTrue(np.allclose(np.zeros((3,3,3)), subvol))

    def test3_generate_sphere(self):
        """

        Test generate_volume() method for Sphere class

        """

        radius = 0.4
        center = 0.5 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)
        y_ax = np.linspace(0, 1, 3)
        z_ax = np.linspace(0, 1, 3)

        c = Sphere(center, radius)
        subvol = c.generate_volume(x_ax, y_ax, z_ax)

        arr = np.zeros((3,3,3))
        arr[1,1,1] = 1

        self.assertTrue(np.allclose(arr, subvol))











    def test1_find_sphere_subvolume(self):
        """

        Test find_subvolume() method for Sphere class.

        """


        radius = 0.5
        center = 0.5 * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol, z_ax_subvol = \
        c.find_subvolume(x_ax, y_ax, z_ax)

        self.assertTrue(np.allclose(x_ax, x_ax_subvol))
        self.assertTrue(np.allclose(y_ax, y_ax_subvol))
        self.assertTrue(np.allclose(z_ax, z_ax_subvol))


    def test2_find_sphere_subvolume(self):
        """

        Test find_subvolume() method for Sphere class

        """


        radius = 0.5
        center = 2. * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol, z_ax_subvol = \
        c.find_subvolume(x_ax, y_ax, z_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol))
        self.assertTrue(np.allclose(arr, y_ax_subvol))
        self.assertTrue(np.allclose(arr, z_ax_subvol))



    def test3_find_sphere_subvolume(self):
        """

        Test find_subvolume() method for Sphere class

        """


        radius = 0.4
        center = 0.5 * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol, y_ax_subvol, z_ax_subvol = \
        c.find_subvolume(x_ax, y_ax, z_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax_subvol))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax_subvol))
        self.assertTrue(np.allclose(z_ax[1:-1], z_ax_subvol))









    def test1_find_sphere_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Sphere class.

        """


        radius = 0.5
        center = 0.5 * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax, z_ax)

        self.assertTrue(np.allclose(x_ax, x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax, y_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(z_ax, z_ax[x_ax_subvol_ix]))


    def test2_find_sphere_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Sphere class

        """


        radius = 0.5
        center = 2. * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax, z_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol_ix))
        self.assertTrue(np.allclose(arr, y_ax_subvol_ix))
        self.assertTrue(np.allclose(arr, z_ax_subvol_ix))



    def test3_find_sphere_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Sphere class

        """


        radius = 0.4
        center = 0.5 * np.ones(3)
        c = Sphere(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)
        z_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax, z_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax[y_ax_subvol_ix]))
        self.assertTrue(np.allclose(z_ax[1:-1], z_ax[z_ax_subvol_ix]))








    def test1_overlap_potential(self):
        """

        Test overlap_potential method for Sphere class

        """

        center = np.zeros(3)
        radius = 0.5
        c1 = Sphere(center, radius)
        c2 = Sphere(center, radius)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)





if __name__ == '__main__':
    print 'Running unit tests for Sphere.py'
    unittest.main()