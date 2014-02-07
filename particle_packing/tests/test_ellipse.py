from particle_packing import ellipse
from particle_packing.ellipse import Ellipse
from scipy.spatial.distance import pdist
import unittest
import numpy as np

class TestCode(unittest.TestCase):


    def test1_pack_rsa_mda(self):
        """

        circle, npoints small

        """

        npoints = 5
        radius = 0.05 * np.ones(2)
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi)
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)




        xy = np.vstack([x, y]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))


    def test2_pack_rsa_mda(self):
        """

        circle, npoints large

        """

        npoints = 50
        radius = 0.05 * np.ones(2)
        phi = 0.
        step_limit = 10 ** 4

        x, y = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit)


        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi)
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))




    def test3_pack_rsa_mda(self):
        """

        circle, random seed test

        """


        npoints = 5
        radius = 0.05 * np.ones(2)
        phi = 0.
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0 = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit,
            randSeed)
        x1, y1 = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))


    def test4_pack_rsa_mda(self):
        """

        ellipse, npoints small

        """

        npoints = 5
        radius = 0.05 * np.array([1., 2.])
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi)
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = Ellipse(center=center, radii=radius, phi=phi)

                center = np.asarray([x[k], y[k]])
                ck = Ellipse(center=center, radii=radius, phi=phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)


        self.assertTrue(npoints == len(x))



    def test5_pack_rsa_mda(self):
        """

        ellipse, npoints larger

        """

        npoints = 15
        radius = 0.05 * np.array([1., 2.])
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack.rsa_mda(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi)
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = Ellipse(center=center, radii=radius, phi=phi)

                center = np.asarray([x[k], y[k]])
                ck = Ellipse(center=center, radii=radius, phi=phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)

        self.assertTrue(npoints == len(x))







    def test1_pack_rsa_md(self):
        """

        circle, npoints small

        """

        npoints = 5
        radius = 0.05 * np.ones(2)
        phi = 0.
        step_limit = 10 ** 2

        x, y, phi = ellipse.pack.rsa_md(npoints, radius, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi[i])
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)




        xy = np.vstack([x, y]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))




    def test2_pack_rsa_md(self):
        """

        circle, npoints large

        """

        npoints = 50
        radius = 0.05 * np.ones(2)
        step_limit = 10 ** 4

        x, y, phi = ellipse.pack.rsa_md(npoints, radius, step_limit)


        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi[i])
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        d = pdist(xy)
        self.assertTrue(d.min() > 2. * radius[0])

        self.assertTrue(npoints == len(x))




    def test3_pack_rsa_md(self):
        """

        circle, random seed test

        """


        npoints = 5
        radius = 0.05 * np.ones(2)
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0, phi0 = ellipse.pack.rsa_md(npoints, radius, step_limit,
            randSeed)
        x1, y1, phi1 = ellipse.pack.rsa_md(npoints, radius, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))
        self.assertTrue(np.allclose(phi0, phi1))






    def test4_pack_rsa_md(self):
        """

        ellipse, npoints small

        """

        npoints = 5
        radius = 0.05 * np.array([1., 2.])
        step_limit = 10 ** 2

        x, y, phi = ellipse.pack.rsa_md(npoints, radius, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi[i])
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = Ellipse(center=center, radii=radius, phi=phi[i])

                center = np.asarray([x[k], y[k]])
                ck = Ellipse(center=center, radii=radius, phi=phi[k])

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)


        self.assertTrue(npoints == len(x))






    def test5_pack_rsa_md(self):
        """

        ellipse, npoints larger

        """

        npoints = 15
        radius = 0.05 * np.array([1., 2.])
        step_limit = 10 ** 2

        x, y, phi = ellipse.pack.rsa_md(npoints, radius, step_limit)



        for i in xrange(len(x)):
            center = np.array([x[i], y[i]])
            c = Ellipse(center, radius, phi[i])
            F = c.square_container_potential()
            self.assertTrue(F >= 1.)


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = Ellipse(center=center, radii=radius, phi=phi[i])

                center = np.asarray([x[k], y[k]])
                ck = Ellipse(center=center, radii=radius, phi=phi[k])

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)

        self.assertTrue(npoints == len(x))






if __name__ == '__main__':
    print 'Running unit tests for ellipse.so'
    unittest.main()