from particle_packing import ellipse
import unittest
import numpy as np
from scipy.spatial.distance import pdist

class TestCode(unittest.TestCase):



    """ pack_rsa_md() """

    def test1_pack_rsa_md(self):
        """

        circle, npoints small

        """

        npoints = 5
        radius = 0.05 * np.ones(2)
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[0])
            self.assertTrue(x[i] < 1. - radius[0])
            self.assertTrue(y[i] > radius[0])
            self.assertTrue(y[i] < 1. - radius[0])


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
        phi = 0.
        step_limit = 10 ** 4

        x, y = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit)


        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[0])
            self.assertTrue(x[i] < 1. - radius[0])
            self.assertTrue(y[i] > radius[0])
            self.assertTrue(y[i] < 1. - radius[0])


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
        phi = 0.
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0 = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit,
            randSeed)
        x1, y1 = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))


    def test4_pack_rsa_md(self):
        """

        ellipse, npoints small

        """

        npoints = 5
        radius = 0.05 * np.array([1., 2.])
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[1])
            self.assertTrue(x[i] < 1. - radius[1])
            self.assertTrue(y[i] > radius[0])
            self.assertTrue(y[i] < 1. - radius[0])


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = ellipse.Ellipse(center=center, radii=radius, phi=phi)

                center = np.asarray([x[k], y[k]])
                ck = ellipse.Ellipse(center=center, radii=radius, phi=phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)


        self.assertTrue(npoints == len(x))



    def test5_pack_rsa_md(self):
        """

        ellipse, npoints larger

        """

        npoints = 15
        radius = 0.05 * np.array([1., 2.])
        phi = 0.
        step_limit = 10 ** 2

        x, y = ellipse.pack_rsa_md_aligned(npoints, radius, phi, step_limit)



        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[1])
            self.assertTrue(x[i] < 1. - radius[1])
            self.assertTrue(y[i] > radius[0])
            self.assertTrue(y[i] < 1. - radius[0])


        xy = np.vstack([x, y]).transpose()
        for i in xrange(len(x)):
            for k in xrange(i):

                center = np.asarray([x[i], y[i]])
                ci = ellipse.Ellipse(center=center, radii=radius, phi=phi)

                center = np.asarray([x[k], y[k]])
                ck = ellipse.Ellipse(center=center, radii=radius, phi=phi)

                F = ci.overlap_potential(ck)
                self.assertTrue(F >= 1.)

        self.assertTrue(npoints == len(x))






if __name__ == '__main__':
    print 'Running unit tests for ellipse.so'
    unittest.main()