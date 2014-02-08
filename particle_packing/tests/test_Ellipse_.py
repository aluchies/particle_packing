import unittest
import numpy as np
from particle_packing.ellipse import Ellipse, \
    overlap_potential, overlap_potential_py, \
    square_container_potential_py, square_container_potential

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

        arr = np.array([[0, 0, 1],
                        [0, 1, 0],
                        [1, 0, 0]])

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


















    def test1_ellipse_overlap(self):
        """

        Full overlap, identical ellipses.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F_py, 0.))


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))










    def test2_ellipse_overlap(self):
        """

        Tangent, x-axis.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))



    def test3_ellipse_overlap(self):
        """

        Tangent, y-axis.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[0., 2.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))


    def test4_ellipse_overlap(self):
        """

        Tangent, all.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))











    def test5_ellipse_overlap(self):
        """

        Overlap, x-axis.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1.1, 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))



    def test6_ellipse_overlap(self):
        """

        Overlap, y-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[0., 2.]])
        radiiB = np.array([1., 1.1])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))


    def test7_ellipse_overlap(self):
        """

        Overlap, all.

        """
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([1.1, 1.1])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))















    def test8_ellipse_overlap(self):
        """

        No overlap, x-axis.

        """

        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([0.9, 1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))



    def test9_ellipse_overlap(self):
        """

        No overlap, y-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[0., 2.]])
        radiiB = np.array([1., 0.9])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))


    def test10_ellipse_overlap(self):
        """

        No overlap, all.

        """
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([0.9, 0.9])
        phiB = 0.



        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))




    def test11_ellipse_overlap(self):
        """

        Overlap, rotated ellipses

        """

        # Ellipse A
        rA = np.array([[0.46728551, 0.47801381]])
        radiiA = np.array([0.1, 0.2])
        phiA = -np.pi / 4.

        # Ellipse B
        rB = np.array([[0.53806873, 0.71002836]])
        radiiB = np.array([0.1, 0.2])
        phiB = -np.pi / 4.



        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        self.assertTrue(np.allclose(F, F_py))





    def test12_overlap_potential(self):
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

















    def test1_container_potential(self):
        """

        completely inside

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.5] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(F_py > 1.)

        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))

 









    def test2_container_potential(self):
        """

        tangent, left side

        """

        # Ellipse A
        rA = np.array( [0.4, 0.5] )
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, 1.))


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))





    def test3_container_potential(self):
        """

        tangent, right side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.6, 0.5] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, 1.))


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))




    def test4_container_potential(self):
        """

        tangent, top side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.6] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, 1.))


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))





    def test5_container_potential(self):
        """

        tangent, bottom side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.4] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, 1.))


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))











    def test6_container_potential(self):
        """

        outside, left side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.3, 0.5] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(F_py < 1.)


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))






    def test7_container_potential(self):
        """

        outside, right side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.7, 0.5] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(F_py < 1.)


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))






    def test8_container_potential(self):
        """

        outside, top side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.7] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(F_py < 1.)


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))







    def test9_container_potential(self):
        """

        outside, bottom side

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.3] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        F_py = square_container_potential_py(rA, radiiA, phiA)
        self.assertTrue(F_py < 1.)


        F = square_container_potential(rA, radiiA, phiA)
        self.assertTrue(np.allclose(F_py, F))



    def test10_container_potential(self):
        """

        Test square_container_potential method for Ellipse class

        """

        # Ellipse A
        rA = np.matrix( np.array( [0.5, 0.5] ))
        radiiA = 0.4 * np.ones(2)
        phiA = 0.

        c1 = Ellipse(rA, radiiA, phiA)

        F = c1.square_container_potential()

        self.assertTrue(F > 1.)










if __name__ == '__main__':
    print 'Running unit tests for Ellipse.py'
    unittest.main()