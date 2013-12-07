import unittest
import numpy as np
from particle_packing.ellipse import overlap_potential, overlap_potential_py

class TestCode(unittest.TestCase):



    def test1_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertEqual(F, F_py)





    def test2_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([1.1, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)






    def test3_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[2., 0.]])
        radiiB = np.array([0.9, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)




    def test4_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[3., 0.]])
        radiiB = np.array([2., 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)



    def test5_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[3., 0.]])
        radiiB = np.array([2.1, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)




    def test6_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[3., 0.]])
        radiiB = np.array([1.9, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)










    def test7_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[np.sqrt(2)], [np.sqrt(2)]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)



    def test8_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[np.sqrt(2)], [np.sqrt(2)]])
        radiiB = np.array([0.9, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)





    def test8_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[np.sqrt(2)], [np.sqrt(2)]])
        radiiB = np.array([1.1, 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)


    def test9_ellipse_overlap(self):
        # Ellipse A
        rA = np.array([[0., 0.]])
        radiiA = np.array([1., 1.])
        phiA = 0.

        # Ellipse B
        rB = np.array([[0.], [0.]])
        radiiB = np.array([1., 1.])
        phiB = 0.

        F = overlap_potential(rA, radiiA, phiA, rB, radiiB, phiB)
        F_py = overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB)

        self.assertAlmostEqual(F, F_py)







if __name__ == '__main__':
    print 'Running unit tests for ellipse_overlap'
    unittest.main()