import unittest
import numpy as np
from particle_packing.sphere import overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):



    def test1_sphere_overlap_py(self):
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








    def test2_sphere_overlap_py(self):
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





    def test3_sphere_overlap_py(self):
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







    def test4_sphere_overlap_py(self):
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






    def test5_sphere_overlap_py(self):
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















    def test6_sphere_overlap_py(self):
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




    def test7_sphere_overlap_py(self):
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




    def test8_sphere_overlap_py(self):
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




    def test9_sphere_overlap_py(self):
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
















    def test10_sphere_overlap_py(self):
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





    def test11_sphere_overlap_py(self):
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





    def test12_sphere_overlap_py(self):
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





    def test13_sphere_overlap_py(self):
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












if __name__ == '__main__':
    print 'Running unit tests for ellipsoid_overlap.py'
    unittest.main()