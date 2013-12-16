import unittest
import numpy as np
from particle_packing.circle import overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):



    def test1_circle_overlap(self):
        """

        Full overlap, identical circles.

        """

        # Circle A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1.])
        phiA = 0.

        # Circle B
        rB = np.array([[0., 0.]])
        radiiB = np.array([1.])
        phiB = 0.

        F_py = overlap_potential_py(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F_py, 0.))


        F = overlap_potential(rA, radiiA, rB, radiiB)
        self.assertTrue(np.allclose(F, F_py))












    def test2_circle_overlap(self):
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





    def test3_circle_overlap(self):
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




    def test4_circle_overlap(self):
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













    def test5_circle_overlap(self):
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





    def test6_circle_overlap(self):
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




    def test7_circle_overlap(self):
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
















    def test8_circle_overlap(self):
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





    def test9_circle_overlap(self):
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



    def test10_circle_overlap(self):
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









if __name__ == '__main__':
    print 'Running unit tests for circle_overlap'
    unittest.main()