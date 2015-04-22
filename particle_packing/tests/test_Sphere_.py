import unittest
import numpy as np
from particle_packing.sphere import Sphere, \
    overlap_potential, overlap_potential_py, \
    contain_potential_py, container_potential_cube_py

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



    def test2_generate_volume(self):
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

    def test3_generate_volume(self):
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






    def test1_generate_volume_spherical_gaussian(self):
        """

        Test generate_volume_spherical_gaussian() method for Sphere class

        """

        radius = 0.2
        center = 0.5 * np.ones(3)


        x_ax = np.linspace(0, 1, 3)
        y_ax = np.linspace(0, 1, 3)
        z_ax = np.linspace(0, 1, 3)

        c = Sphere(center, radius)
        subvol = c.generate_volume_spherical_gaussian(x_ax, y_ax, z_ax)

        arr = np.zeros((3,3,3))
        arr[0,1,1] = 1
        arr[1,1,1] = 1
        arr[1,1,0] = 1
        arr[1,0,1] = 1
        arr[1,1,2] = 1
        arr[1,2,1] = 1
        arr[2,1,1] = 1
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

        Full overlap.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 0.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))








    def test2_overlap_potential(self):
        """

        Tangent, x-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])

        # Sphere B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test3_overlap_potential(self):
        """

        Tangent, y-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])


        # Sphere B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))







    def test4_overlap_potential(self):
        """

        Tangent, z-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 2.]])
        radiiA = np.array([1.])


        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))






    def test5_overlap_potential(self):
        """

        Tangent, all.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])


        # Sphere B
        rB = np.sqrt(4. / 3. ) * np.array([[1., 1., 1.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))















    def test6_overlap_potential(self):
        """

        Overlap, x-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.1])

        # Sphere B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))




    def test7_overlap_potential(self):
        """

        Overlap, y-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.1])


        # Sphere B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1.])

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))




    def test8_overlap_potential(self):
        """

        Overlap, z-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.1])


        # Sphere B
        rB = np.array([[0., 0., 2.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))




    def test9_overlap_potential(self):
        """

        Overlap, all.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])


        # Sphere B
        rB = np.sqrt(4. / 3.) * np.array([[1., 1., 1.]]) - 0.1
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))
















    def test10_overlap_potential(self):
        """

        No overlap, x-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9])


        # Sphere B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test11_overlap_potential(self):
        """

        No overlap, y-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9])


        # Sphere B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test12_overlap_potential(self):
        """

        No overlap, z-axis.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9])


        # Sphere B
        rB = np.array([[0., 0., 2.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))





    def test13_overlap_potential(self):
        """

        No overlap, all.

        """
        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9])


        # Sphere B
        rB = np.sqrt(4. / 3.) * np.array([[1., 1., 1.]])
        radiiB = np.array([1.])


        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(F_py > 1.)










    def test14_overlap_potential(self):
        """

        Test overlap_potential method for Sphere class

        """

        center = np.zeros(3)
        radius = 0.5
        c1 = Sphere(center, radius)
        c2 = Sphere(center, radius)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)














    def test1_contain_potential(self):
        """

        Full containment

        """

        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py > 1.)



    def test2_contain_potential(self):
        """

        Full containment, tangent

        """

        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py == 1.)



    def test3_contain_potential(self):
        """

        Partially oustide due to size

        """

        # Sphere A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.1])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)


    def test4_contain_potential(self):
        """

        Partially oustide due to location

        """

        # Sphere A
        rA = np.array([[0.1, 0., 0.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)




    def test5_contain_potential(self):
        """

        Completely oustide due to location, x-axis

        """

        # Sphere A
        rA = np.array([[2., 0., 0.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)





    def test5_contain_potential(self):
        """

        Completely oustide due to location, y-axis

        """

        # Sphere A
        rA = np.array([[0., 2., 0.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)




    def test6_contain_potential(self):
        """

        Completely oustide due to location, z-axis

        """

        # Sphere A
        rA = np.array([[0., 0., 2.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Sphere B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)





    def test7_contain_potential(self):
        """

        Completely inside

        """

        # Sphere A
        rA = np.array([[0.3, 0.5, 0.5]])
        radiiA = np.array([0.25])

        # Sphere B
        rB = np.array([[0.5, 0.5, 0.5]])
        radiiB = np.array([0.5])

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py > 1.)




    def test8_contain_potential(self):
        """

        Completely inside, tangent

        """

        # Sphere A
        rA = np.array([[0.25, 0.5, 0.5]])
        radiiA = np.array([0.25])

        # Sphere B
        rB = np.array([[0.5, 0.5, 0.5]])
        radiiB = np.array([0.5])

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py == 1.)




    def test9_contain_potential(self):
        """

        Completely inside, tangent

        """

        # Sphere A
        rA = np.array([[0.2, 0.5, 0.5]])
        radiiA = np.array([0.25])

        # Sphere B
        rB = np.array([[0.5, 0.5, 0.5]])
        radiiB = np.array([0.5])

        G_py = contain_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(G_py < 1.)












    def test10_contain_potential(self):
        """

        Test contain_potential_py called as class method

        """

        center = np.array([[0., 0., 0.]])
        radius = np.array([0.9])
        c1 = Sphere(center, radius)

        center = np.array([[0., 0., 0.]])
        radius = np.array([1.])
        c2 = Sphere(center, radius)


        G_py = c1.contain_potential(c2)
        self.assertTrue(G_py > 1.)













    def test1_container_potential(self):
        """

        Completely inside square container

        """

        # Sphere A
        rA = np.array([[0.5, 0.5, 0.5]])
        radiiA = np.array([0.1])

        H = container_potential_cube_py(rA, radiiA)

        self.assertTrue(H > 1.)



    def test2_container_potential(self):
        """

        Completely outside container

        """

        # Sphere A
        rA = np.array([[2., 2., 2.]])
        radiiA = np.array([0.1])

        H = container_potential_cube_py(rA, radiiA)

        self.assertTrue(H < 1.)


    def test3_container_potential(self):
        """

        Completely inside and tangent to container

        """

        # Sphere A
        rA = np.array([[0.5, 0.5, 0.5]])
        radiiA = np.array([0.5])

        H = container_potential_cube_py(rA, radiiA)

        self.assertTrue(np.allclose(H, 1.))





    def test4_container_potential(self):
        """

        Test calling container_potential class method

        """

        # Sphere A
        rA = np.array([[0.5, 0.5, 0.5]])
        radiiA = np.array([0.1])

        c = Sphere(rA, radiiA)

        H = c.container_potential('cube')


        self.assertTrue(H > 1.)







if __name__ == '__main__':
    print 'Running unit tests for Sphere.py'
    unittest.main()