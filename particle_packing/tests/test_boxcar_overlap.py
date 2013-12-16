import unittest
import numpy as np
from particle_packing.boxcar import overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):



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




if __name__ == '__main__':
    print 'Running unit tests for boxcar_overlap'
    unittest.main()