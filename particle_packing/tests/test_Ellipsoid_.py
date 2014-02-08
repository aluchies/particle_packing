import unittest
import numpy as np
from particle_packing.ellipsoid import Ellipsoid, \
    overlap_potential, overlap_potential_py, \
    cube_container_potential_py, cube_container_potential

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









    def test1_overlap_potential(self):
        """

        Full overlap.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, 0.))


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))






    def test2_overlap_potential(self):
        """

        Tangent, x-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))



    def test3_overlap_potential(self):
        """

        Tangent, y-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))





    def test4_overlap_potential(self):
        """

        Tangent, z-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 2.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 0., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))



    def test5_overlap_potential(self):
        """

        Tangent, all.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.sqrt(4. / 3. ) * np.array([[1., 1., 1.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, 1.))


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))













    def test6_overlap_potential(self):
        """

        Overlap, x-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1.1, 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))


    def test7_overlap_potential(self):
        """

        Overlap, y-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1.1, 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))


    def test8_overlap_potential(self):
        """

        Overlap, z-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.1])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 0., 2.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))


    def test9_overlap_potential(self):
        """

        Overlap, all.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.sqrt(4. / 3.) * np.array([[1., 1., 1.]]) - 0.1
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py < 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))














    def test10_overlap_potential(self):
        """

        No overlap, x-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9, 1., 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[2., 0., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))



    def test11_overlap_potential(self):
        """

        No overlap, y-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 0.9, 1.])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 2., 0.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))


    def test12_overlap_potential(self):
        """

        No overlap, z-axis.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([1., 1., 0.9])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.array([[0., 0., 2.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))



    def test13_overlap_potential(self):
        """

        No overlap, all.

        """
        # Ellipse A
        rA = np.array([[0., 0., 0.]])
        radiiA = np.array([0.9, 0.9, 0.9])
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        # Ellipse B
        rB = np.sqrt(4. / 3.) * np.array([[1., 1., 1.]])
        radiiB = np.array([1., 1., 1.])
        phiB = 0.
        rotaxB = np.array([1., 0., 0.])

        F_py = overlap_potential_py(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(F_py > 1.)


        F = overlap_potential(rA, radiiA, phiA, rotaxA,
            rB, radiiB, phiB, rotaxB)
        self.assertTrue(np.allclose(F_py, F))











    def test14_overlap_potential(self):
        """

        Test overlap_potential method for Ellipsoid class

        """

        center = np.zeros(3)
        radii = 0.5 * np.ones(3)
        phi = 0.
        rot_ax = np.array([1., 0., 0.])
        c1 = Ellipsoid(center, radii, rot_ax, phi)
        c2 = Ellipsoid(center, radii, rot_ax, phi)

        F = c1.overlap_potential(c2)

        self.assertTrue(F == 0.)
















    def test1_container_potential(self):
        """

        completely contained

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py > 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))






    def test2_container_potential(self):
        """

        tangent top

        """

        # Ellipsoid A
        rA = [0.4, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))



    def test3_container_potential(self):
        """

        tangent bottom

        """

        # Ellipsoid A
        rA = [0.6, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))





    def test4_container_potential(self):
        """

        tangent right

        """

        # Ellipsoid A
        rA = [0.5, 0.4, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test5_container_potential(self):
        """

        tangent left

        """

        # Ellipsoid A
        rA = [0.5, 0.6, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))



    def test6_container_potential(self):
        """

        tangent front

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.4]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))



    def test7_container_potential(self):
        """

        tangent back

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.6]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, 1.))

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test8_container_potential(self):
        """

        outside, top

        """

        # Ellipsoid A
        rA = [0.3, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test9_container_potential(self):
        """

        outside, bottom

        """

        # Ellipsoid A
        rA = [0.7, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test10_container_potential(self):
        """

        outside, left

        """

        # Ellipsoid A
        rA = [0.5, 0.3, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))



    def test11_container_potential(self):
        """

        outside, right

        """

        # Ellipsoid A
        rA = [0.5, 0.7, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))



    def test12_container_potential(self):
        """

        outside, front

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.3]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test13_container_potential(self):
        """

        outside, back

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.7]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        F_py = cube_container_potential_py(rA, radiiA, phiA, rotaxA)
        self.assertTrue(F_py < 1.)

        F = cube_container_potential(rA, radiiA, phiA, rotaxA)
        self.assertTrue(np.allclose(F_py, F))




    def test14_container_potential(self):
        """

        Test cube_container_potential method for Ellipsoid class

        """

        # Ellipsoid A
        rA = [0.5, 0.5, 0.5]
        radiiA = 0.4 * np.ones(3)
        phiA = 0.
        rotaxA = np.array([1., 0., 0.])

        c1 = Ellipsoid(rA, radiiA, rotaxA, phiA)

        F = c1.cube_container_potential()

        self.assertTrue(F > 1.)







if __name__ == '__main__':
    print 'Running unit tests for Ellipsoid.py'
    unittest.main()