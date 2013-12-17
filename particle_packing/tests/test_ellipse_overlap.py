import unittest
import numpy as np
from particle_packing.ellipse import overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):



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










if __name__ == '__main__':
    print 'Running unit tests for ellipse_overlap'
    unittest.main()