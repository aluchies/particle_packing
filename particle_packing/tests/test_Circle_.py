import unittest
import numpy as np
from particle_packing.circle import Circle, \
    overlap_potential, overlap_potential_py, \
    contain_potential_py, container_potential_square_py

class TestCode(unittest.TestCase):





    def test1_constuctor(self):
        """

        Test simple constructor example

        """

        center = np.asarray([0, 0 ], dtype=float)
        radius = 1

        c = Circle(center, radius)
        self.assertTrue(c.radius, radius)
        self.assertTrue(np.allclose(c.center, center))


    def test2_constuctor(self):
        """

        Test constructor failure

        """

        center = [0, 0, 0]
        radius = 1
        self.assertRaises(ValueError, Circle, center, radius)

        center = np.asarray([0,0])
        radius = 'a'
        self.assertRaises(ValueError, Circle, center, radius)









    def test1_generate_volume(self):
        """

        Sphere at origin with radius 1

        """


        center = np.array([0., 0.])
        radius = 1
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Circle(center, radius)

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
        x = np.linspace(-1., 1., 3)
        y = np.linspace(-1., 1., 3)

        c = Circle(center, radius)

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
        c = Circle(center, radius)

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
        c = Circle(center, radius)

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
        c = Circle(center, radius)

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
        c = Circle(center, radius)

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
        c = Circle(center, radius)

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
        c = Circle(center, radius)

        x_ax = np.linspace(0, 1, 10)
        y_ax = np.linspace(0, 1, 10)

        x_ax_subvol_ix, y_ax_subvol_ix = \
        c.find_subvolume_ix(x_ax, y_ax)


        self.assertTrue(np.allclose(x_ax[1:-1], x_ax[x_ax_subvol_ix]))
        self.assertTrue(np.allclose(y_ax[1:-1], y_ax[y_ax_subvol_ix]))











    def test1_overlap_potential(self):
        """

        Full overlap, identical circles.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 0.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))












    def test2_overlap_potential(self):
        """

        Tangent, x-axis.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test3_overlap_potential(self):
        """

        Tangent, y-axis.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[0., 2.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))




    def test4_overlap_potential(self):
        """

        Tangent, all.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))













    def test5_overlap_potential(self):
        """

        Overlap, x-axis.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1.1])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test6_overlap_potential(self):
        """

        Overlap, y-axis.

        """
        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[0., 2.]])
        radiiB = np.array([1.1])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))




    def test7_overlap_potential(self):
        """

        Overlap, all.

        """
        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([1.1])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))
















    def test8_overlap_potential(self):
        """

        No overlap, x-axis.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[2., 0.]])
        radiiB = np.array([0.9])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)



        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test9_overlap_potential(self):
        """

        No overlap, y-axis.

        """
        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.array([[0., 2.]])
        radiiB = np.array([0.9])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))



    def test10_overlap_potential(self):
        """

        No overlap, all.

        """
        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])

        # Circle B
        rB = np.sqrt(2) * np.array([[1., 1.]])
        radiiB = np.array([0.9])



        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))










    def test11_overlap_potential(self):
        """

        Test overlap_potential method for Circle class

        """

        center = np.zeros(2)
        radius = 0.5
        c1 = Circle(center, radius)
        c2 = Circle(center, radius)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)











    def test1_contain_potential(self):
        """

        Full containment

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([0.9])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py > 1.)



    def test2_contain_potential(self):
        """

        Full containment, tangent

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py == 1.)



    def test3_contain_potential(self):
        """

        Partially oustide due to size

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.1])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)


    def test4_contain_potential(self):
        """

        Partially oustide due to location

        """

        # Circle A
        rA = np.array([[0.1, 0.]])
        radiiA = np.array([1.])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)




    def test5_contain_potential(self):
        """

        Completely oustide due to location, x-axis

        """

        # Circle A
        rA = np.array([[2., 0.]])
        radiiA = np.array([1.])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)





    def test5_contain_potential(self):
        """

        Completely oustide due to location, y-axis

        """

        # Circle A
        rA = np.array([[0., 2.]])
        radiiA = np.array([1.])


        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])


        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)



    def test6_contain_potential(self):
        """

        Test contain_potential_py called as class method

        """

        center = np.array([[0., 0.]])
        radius = np.array([0.9])
        c1 = Circle(center, radius)

        center = np.array([[0., 0.]])
        radius = np.array([1.])
        c2 = Circle(center, radius)


        G_py = c1.contain_potential(c2)
        self.assertTrue(G_py > 1.)














    def test1_container_potential(self):
        """

        Completely inside square container

        """

        # Circle A
        rA = np.array([[0.5, 0.5]])
        radiiA = np.array([0.1])

        H = container_potential_square_py(rA, radiiA)

        self.assertTrue(H > 1.)



    def test2_container_potential(self):
        """

        Completely outside square container

        """

        # Circle A
        rA = np.array([[2., 2.]])
        radiiA = np.array([0.1])

        H = container_potential_square_py(rA, radiiA)

        self.assertTrue(H < 1.)


    def test3_container_potential(self):
        """

        Completely inside and tange to square container

        """

        # Circle A
        rA = np.array([[0.5, 0.5]])
        radiiA = np.array([0.5])

        H = container_potential_square_py(rA, radiiA)

        self.assertTrue(np.allclose(H, 1.))






    def test4_container_potential(self):
        """

        Test calling container_potential class method

        """

        # Circle A
        rA = np.array([[0.5, 0.5]])
        radiiA = np.array([0.1])

        c = Circle(rA, radiiA)

        H = c.container_potential('square')


        self.assertTrue(H > 1.)










if __name__ == '__main__':
    print 'Running unit tests for Circle.py'
    unittest.main()