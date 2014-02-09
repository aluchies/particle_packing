from particle_packing import ellipsoid
from particle_packing.ellipsoid import Ellipsoid
from scipy.spatial.distance import pdist
import unittest
import numpy as np

class TestCode(unittest.TestCase):


    def test1_pack_rsa_mda(self):
        """

        sphere, npoints small

        """

        npoints = 5
        radius = 0.05 * np.ones(3)
        rotax = np.array([1., 0., 0.])
        phi = 0.
        step_limit = 10 ** 2

        x, y, z = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i], z[i]])
            c = Ellipsoid(center, radius, rotax, phi)
            F = c.cube_container_potential()
            self.assertTrue(F >= 1.)




        xy = np.vstack([x, y, z]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))




    def test2_pack_rsa_mda(self):
        """

        sphere, npoints large

        """

        npoints = 25
        radius = 0.05 * np.ones(3)
        rotax = np.array([1., 0., 0.])
        phi = 0.
        step_limit = 10 ** 3

        x, y, z = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i], z[i]])
            c = Ellipsoid(center, radius, rotax, phi)
            F = c.cube_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y, z]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))






    def test3_pack_rsa_mda(self):
        """

        sphere, test random seed

        """

        npoints = 5
        radius = 0.05 * np.ones(3)
        rotax = np.array([1., 0., 0.])
        phi = 0.
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0, z0 = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit, randSeed)
        x1, y1, z1 = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit, randSeed)



        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))
        self.assertTrue(np.allclose(z0, z1))





    def test4_pack_rsa_mda(self):
        """

        ellipsoid, npoints small

        """

        npoints = 5
        radius = 0.05 * np.array([2., 1., 1.])
        rotax = np.array([1., 0., 0.])
        phi = 0.
        step_limit = 10 ** 2

        x, y, z = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i], z[i]])
            c = Ellipsoid(center, radius, rotax, phi)
            F = c.cube_container_potential()
            self.assertTrue(F >= 1.)




        xyz = np.vstack([x, y, z]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i], z[i]])
                ci = Ellipsoid(center, radius, rotax, phi)

                center = np.asarray([x[k], y[k], z[k]])
                ck = Ellipsoid(center, radius, rotax, phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)

        self.assertTrue(npoints == len(x))





    def test5_pack_rsa_mda(self):
        """

        ellipsoid, npoints small

        """

        npoints = 25
        radius = 0.05 * np.array([2., 1., 1.])
        rotax = np.array([1., 0., 0.])
        phi = 0.
        step_limit = 10 ** 2

        x, y, z = ellipsoid.pack.rsa_mda(npoints, radius, rotax, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i], z[i]])
            c = Ellipsoid(center, radius, rotax, phi)
            F = c.cube_container_potential()
            self.assertTrue(F >= 1.)




        xyz = np.vstack([x, y, z]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i], z[i]])
                ci = Ellipsoid(center, radius, rotax, phi)

                center = np.asarray([x[k], y[k], z[k]])
                ck = Ellipsoid(center, radius, rotax, phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)

        self.assertTrue(npoints == len(x))







if __name__ == '__main__':
    print 'Running unit tests for ellipsoid.so'
    unittest.main()