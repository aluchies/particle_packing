import unittest
import numpy as np
from particle_packing.ellipsoid import overlap_potential_py

class TestCode(unittest.TestCase):



    def test1_ellipse_overlap(self):
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

        import code; code.interact(local=dict(globals(), **locals()))










if __name__ == '__main__':
    print 'Running unit tests for ellipse_overlap'
    unittest.main()