import unittest
import numpy as np
from particle_packing.boxcar import Boxcar

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = 0.
        radius = 1.
        c = Boxcar(center, radius)
        self.assertTrue(c.radius == radius)
        self.assertTrue(c.center == center


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = 0
        radius = 'a'
        self.assertRaises(ValueError, Boxcar, center, radius)

        center = 0.
        radius = 'a'
        self.assertRaises(ValueError, Boxcar, center, radius)







    def test1_generate_volume(self):
        """

        Test generate_volume() method for Boxcar class

        """

        radius = 0.5
        center = 0.5


        x_ax = np.linspace(0, 1, 3)


        c = Boxcar(center, radius)
        subvol = c.generate_volume(x_ax)

        arr = np.ones(3)
        self.assertTrue(np.allclose(arr, subvol))

    def test2_generate_boxcar(self):
        """

        Test generate_volume() method for Boxcar class

        """

        radius = 0.5
        center = 2.


        x_ax = np.linspace(0, 1, 3)

        c = Boxcar(center, radius)
        subvol = c.generate_volume(x_ax)

        self.assertTrue(np.allclose(np.zeros((3)), subvol))

    def test3_generate_boxcar(self):
        """

        Test generate_volume() method for Boxcar class

        """

        radius = 0.4
        center = 0.5 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)

        c = Boxcar(center, radius)
        subvol = c.generate_volume(x_ax, y_ax, z_ax)

        arr = np.zeros(3)
        arr[1] = 1

        self.assertTrue(np.allclose(arr, subvol))











    def test1_find_boxcar_subvolume(self):
        """

        Test find_subvolume() method for Boxcar class.

        """


        radius = 0.5
        center = 0.5
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol = \
        c.find_subvolume(x_ax)

        self.assertTrue(np.allclose(x_ax, x_ax_subvol))


    def test2_find_boxcar_subvolume(self):
        """

        Test find_subvolume() method for Boxcar class

        """


        radius = 0.5
        center = 2.
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol = \
        c.find_subvolume(x_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol))



    def test3_find_boxcar_subvolume(self):
        """

        Test find_subvolume() method for Boxcar class

        """


        radius = 0.4
        center = 0.5
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol = \
        c.find_subvolume(x_ax, y_ax, z_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax_subvol))









    def test1_find_boxcar_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Boxcar class.

        """


        radius = 0.5
        center = 0.5
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax)

        self.assertTrue(np.allclose(x_ax, x_ax[x_ax_subvol_ix]))


    def test2_find_boxcar_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Boxcar class

        """


        radius = 0.5
        center = 2.
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax)

        arr = np.array([])

        self.assertTrue(np.allclose(arr, x_ax_subvol_ix))



    def test3_find_boxcar_subvolume_ix(self):
        """

        Test find_subvolume_ix() method for Boxcar class

        """


        radius = 0.4
        center = 0.5
        c = Boxcar(center, radius)

        x_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax[x_ax_subvol_ix]))





if __name__ == '__main__':
    print 'Running unit tests for Boxcar.py'
    unittest.main()