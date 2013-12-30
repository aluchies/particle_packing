import unittest
import numpy as np
from particle_packing.boxcar import Boxcar, \
    overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = 0.
        radius = 1.
        c = Boxcar(center, radius)
        self.assertTrue(c.radius == radius)
        self.assertTrue(c.center == center)


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
        center = 0.5


        x_ax = np.linspace(0, 1, 3)

        c = Boxcar(center, radius)
        subvol = c.generate_volume(x_ax)

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
        c.find_subvolume(x_ax)


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





    def test1_boxcar_overlap(self):
        """

        Full overlap, identical ellipses.

        """

        # A
        rA = np.array([0.])
        radiiA = np.array([1.])

        # B
        rB = np.array([0.])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)

        self.assertTrue(np.allclose(F_py, 0.))

        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))


    def test2_boxcar_overlap(self):
        """

        Tangent

        """

        # Ellipse A
        rA = np.array([[0.]])
        radiiA = np.array([1.])

        # Ellipse B
        rB = np.array([[2.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))


    def test3_boxcar_overlap(self):
        """

        Overlap

        """

        # Ellipse A
        rA = np.array([[0.]])
        radiiA = np.array([1.])

        # Ellipse B
        rB = np.array([[2.]])
        radiiB = np.array([1.1])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))


    def test4_boxcar_overlap(self):
        """

        No overlap

        """

        # Ellipse A
        rA = np.array([[0.]])
        radiiA = np.array([1.])

        # Ellipse B
        rB = np.array([[2.]])
        radiiB = np.array([0.9])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))











    def test5_overlap_potential(self):
        """

        Test overlap_potential method for Boxcar class

        """

        center = 0.
        radii = 0.5
        phi = 0.
        c1 = Boxcar(center, radii)
        c2 = Boxcar(center, radii)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)





if __name__ == '__main__':
    print 'Running unit tests for Boxcar.py'
    unittest.main()