from particle_packing import circle
from particle_packing.circle import Circle
import unittest
import numpy as np
from scipy.spatial.distance import pdist

class TestCode(unittest.TestCase):

    """ pack.grid_md() """

    def test1_pack_grid_md(self):
        """

        Test case with zero points.

        """

        x, y = circle.pack.grid_md(npoints=0, radius=0.05)
        self.assertTrue(x.size == 0)
        self.assertTrue(y.size == 0)

    def test2_pack_grid_md(self):
        """

        Test case with default arguments.

        """

        npoints = 5
        radius = 0.05
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)
        self.assertTrue(y.size == npoints)

    def test3_pack_grid_md(self):
        """

        Test case with npoints large

        """

        npoints = 50
        radius = 0.05
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)
        self.assertTrue(y.size == npoints)

        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius)

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))




    def test4_pack_grid_md(self):
        """

        Test case with npoints too large

        """

        npoints = 1000
        radius = 0.05
        self.assertRaises(ValueError, circle.pack.grid_md, npoints, 0.05)









    """ pack.metro_md() """


    def test1_pack_metro_md(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        success_steps = circle.pack.metro_md(x, y, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius)

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_md(self):
        """

        Test case with npoints small

        """

        npoints = 50
        radius = 0.05
        step_limit = 10 ** 3
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        success_steps = circle.pack.metro_md(x, y, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius)

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))


        self.assertTrue(success_steps > 0)




    def test3_pack_metro_md(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y0 = np.ascontiguousarray([0.1, 0.3, 0.5])

        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y1 = np.ascontiguousarray([0.1, 0.3, 0.5])

        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100

        success_steps0 = circle.pack.metro_md(x0, y0, radius, step_limit,
            randSeed)
        success_steps1 = circle.pack.metro_md(x1, y1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))



    def test4_pack_metro_md(self):
        """

        Test case when all steps are successful

        """

        npoints = 50
        radius = 0.0
        step_limit = 10 ** 3
        x, y = circle.pack.poisson_point(npoints=npoints)
        success_steps = circle.pack.metro_md(x, y, radius, step_limit)

        self.assertTrue(success_steps == step_limit)









    def test1_pack_metro_pd(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = circle.pack.metro_pd(x, y, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius[i])

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius[k])
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_pd(self):
        """

        Test case with npoints small

        """

        npoints = 50
        radius = 0.05
        step_limit = 10 ** 3
        x, y = circle.pack.grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = circle.pack.metro_pd(x, y, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius[i])

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius[k])
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test3_pack_metro_pd(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y0 = np.ascontiguousarray([0.1, 0.3, 0.5])

        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])
        y1 = np.ascontiguousarray([0.1, 0.3, 0.5])


        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100
        npoints = 3
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))

        success_steps0 = circle.pack.metro_pd(x0, y0, radius, step_limit,
            randSeed)
        success_steps1 = circle.pack.metro_pd(x1, y1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))



    def test4_pack_metro_pd(self):
        """

        Test case when all steps are successful

        """

        npoints = 500
        radius = 0.0
        radius = np.ascontiguousarray(radius * np.ones(npoints))
        step_limit = 10 ** 3
        x, y = circle.pack.poisson_point(npoints=npoints)
        success_steps = circle.pack.metro_pd(x, y, radius, step_limit)

        self.assertTrue(success_steps == step_limit)







    """ pack.rsa_md() """

    def test1_pack_rsa_md(self):
        """

        Test case with npoints small

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 2

        x, y = circle.pack.rsa_md(npoints, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius)

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(npoints == len(x))


    def test2_pack_rsa_md(self):
        """

        Test case with npoints large

        """

        npoints = 50
        radius = 0.05
        step_limit = 10 ** 4

        x, y = circle.pack.rsa_md(npoints, radius, step_limit)


        for i in xrange(npoints):
            ci = Circle([x[i], y[i]], radius)

            # inside the container
            H = ci.container_potential('square')
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Circle([x[k], y[k]], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(npoints == len(x))


    def test3_pack_rsa_md(self):
        """

        Test case random seed

        """


        npoints = 5
        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100

        x0, y0 = circle.pack.rsa_md(npoints, radius, step_limit,
            randSeed)
        x1, y1 = circle.pack.rsa_md(npoints, radius, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))
        self.assertTrue(np.allclose(y0, y1))




if __name__ == '__main__':
    print 'Running unit tests for circle.so'
    unittest.main()