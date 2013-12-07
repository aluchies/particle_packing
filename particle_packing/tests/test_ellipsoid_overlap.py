import unittest
import numpy as np
from particle_packing.ellipsoid import overlap_potential_py

class TestCode(unittest.TestCase):



    def test1_ellipse_overlap_py(self):
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






    def test2_ellipse_overlap_py(self):
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


    def test3_ellipse_overlap_py(self):
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





    def test4_ellipse_overlap_py(self):
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


    def test5_ellipse_overlap_py(self):
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













    def test6_ellipse_overlap_py(self):
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


    def test7_ellipse_overlap_py(self):
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


    def test8_ellipse_overlap_py(self):
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


    def test9_ellipse_overlap_py(self):
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














    def test10_ellipse_overlap_py(self):
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



    def test11_ellipse_overlap_py(self):
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


    def test12_ellipse_overlap_py(self):
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



    def test13_ellipse_overlap_py(self):
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










if __name__ == '__main__':
    print 'Running unit tests for ellipsoid_overlap.py'
    unittest.main()