from particle_packing import spheres
import unittest
import numpy as np
from scipy.spatial.distance import pdist

class TestCode(unittest.TestCase):

    """ pack_grid_md() """

    def test1_pack_grid_md(self):
        """

        Test case with zero points.

        """

        x, y, z = spheres.pack_grid_md(npoints=0, radius=0.05)
        self.assertTrue(x.size == 0)
        self.assertTrue(y.size == 0)
        self.assertTrue(z.size == 0)

    def test2_pack_grid_md(self):
        """

        Test case with default arguments.

        """

        npoints = 5
        radius = 0.05
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)
        self.assertTrue(y.size == npoints)
        self.assertTrue(z.size == npoints)

    def test3_pack_grid_md(self):
        """

        Test case with npoints large

        """

        npoints = 500
        radius = 0.05
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)
        self.assertTrue(y.size == npoints)
        self.assertTrue(z.size == npoints)

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius)


    def test4_pack_grid_md(self):
        """

        Test case with npoints too large

        """

        npoints = 1000
        radius = 0.05
        self.assertRaises(ValueError, spheres.pack_grid_md, npoints, 0.05)









    """ pack_metro_md() """


    def test1_pack_metro_md(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        success_steps = spheres.pack_metro_md(x, y, z, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)
            self.assertTrue(y[i] > radius)
            self.assertTrue(y[i] < 1. - radius)
            self.assertTrue(z[i] > radius)
            self.assertTrue(z[i] < 1. - radius)

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius)

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_md(self):
        """

        Test case with npoints small

        """

        npoints = 500
        radius = 0.05
        step_limit = 10 ** 3
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        success_steps = spheres.pack_metro_md(x, y, z, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)
            self.assertTrue(y[i] > radius)
            self.assertTrue(y[i] < 1. - radius)
            self.assertTrue(z[i] > radius)
            self.assertTrue(z[i] < 1. - radius)

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius)

        self.assertTrue(success_steps > 0)


    def test3_pack_metro_md(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        z0 = np.ascontiguousarray([0.1, 0.3, 0.5])

        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        z1 = np.ascontiguousarray([0.1, 0.3, 0.5])

        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100

        success_steps0 = spheres.pack_metro_md(x0, y0, z0, radius, step_limit,
            randSeed)
        success_steps1 = spheres.pack_metro_md(x1, y1, z1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))
        self.assertTrue(np.allclose(z0, z1))






    """ pack_metro_pd() """


    def test1_pack_metro_pd(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = spheres.pack_metro_pd(x, y, z, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[i])
            self.assertTrue(x[i] < 1. - radius[i])
            self.assertTrue(y[i] > radius[i])
            self.assertTrue(y[i] < 1. - radius[i])
            self.assertTrue(z[i] > radius[i])
            self.assertTrue(z[i] < 1. - radius[i])

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius.min())

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_pd(self):
        """

        Test case with npoints small

        """

        npoints = 500
        radius = 0.05
        step_limit = 10 ** 3
        x, y, z = spheres.pack_grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = spheres.pack_metro_pd(x, y, z, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius[i])
            self.assertTrue(x[i] < 1. - radius[i])
            self.assertTrue(y[i] > radius[i])
            self.assertTrue(y[i] < 1. - radius[i])
            self.assertTrue(z[i] > radius[i])
            self.assertTrue(z[i] < 1. - radius[i])

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius.min())

        self.assertTrue(success_steps > 0)


    def test3_pack_metro_pd(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        z0 = np.ascontiguousarray([0.1, 0.3, 0.5])

        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        z1 = np.ascontiguousarray([0.1, 0.3, 0.5])

        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100
        npoints = 3
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))

        success_steps0 = spheres.pack_metro_pd(x0, y0, z0, radius, step_limit,
            randSeed)
        success_steps1 = spheres.pack_metro_pd(x1, y1, z1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))
        self.assertTrue(np.allclose(z0, z1))







    """ pack_rsa_md() """

    def test1_pack_rsa_md(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2

        x, y, z = spheres.pack_rsa_md(npoints, radius, step_limit)


        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)
            self.assertTrue(y[i] > radius)
            self.assertTrue(y[i] < 1. - radius)
            self.assertTrue(z[i] > radius)
            self.assertTrue(z[i] < 1. - radius)

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius)

        self.assertTrue(npoints == len(x))


    def test2_pack_rsa_md(self):
        """

        Test case with npoints large

        """

        npoints = 250
        radius = 0.05
        step_limit = 10 ** 4

        x, y, z = spheres.pack_rsa_md(npoints, radius, step_limit)


        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)
            self.assertTrue(y[i] > radius)
            self.assertTrue(y[i] < 1. - radius)
            self.assertTrue(z[i] > radius)
            self.assertTrue(z[i] < 1. - radius)

        xyz = np.vstack([x, y, z]).transpose()
        d = pdist(xyz)
        self.assertTrue(d.min() > 2. * radius)

        self.assertTrue(npoints == len(x))


    def test3_pack_rsa_md(self):
        """

        Test case random seed

        """


        npoints = 5
        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0, z0 = spheres.pack_rsa_md(npoints, radius, step_limit,
            randSeed)
        x1, y1, z1 = spheres.pack_rsa_md(npoints, radius, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))
        self.assertTrue(np.allclose(z0, z1))




if __name__ == '__main__':
    print 'Running unit tests for spheres.so'
    unittest.main()