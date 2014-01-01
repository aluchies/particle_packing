from particle_packing import boxcar
from particle_packing.boxcar import Boxcar
import unittest
import numpy as np
from scipy.spatial.distance import pdist

class TestCode(unittest.TestCase):

    """ pack_grid_md() """

    def test1_pack_grid_md(self):
        """

        Test case with zero points.

        """

        x = boxcar.pack.grid_md(npoints=0, radius=0.05)
        self.assertTrue(x.size == 0)


    def test2_pack_grid_md(self):
        """

        Test case with default arguments.

        """

        npoints = 1
        radius = 0.01
        x = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)



    def test3_pack_grid_md(self):
        """

        Test case with npoints large

        """

        npoints = 50
        radius = 0.01
        x = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        self.assertTrue(x.size == npoints)

        for i in xrange(npoints):
            ci = Boxcar(x[i], radius)

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))



    def test4_pack_grid_md(self):
        """

        Test case with npoints too large

        """

        npoints = 10 ** 4
        radius = 0.01
        self.assertRaises(ValueError, boxcar.pack.grid_md, npoints, 0.05)









    """ pack_metro_md() """


    def test1_pack_metro_md(self):
        """

        Test case with npoints small

        """

        npoints = 2
        radius = 0.05
        step_limit = 10 ** 2
        x = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        success_steps = boxcar.pack.metro_md(x, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)

        for i in xrange(npoints):
            ci = Boxcar(x[i], radius)

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_md(self):
        """

        Test case with npoints large

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 3
        x = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        success_steps = boxcar.pack.metro_md(x, radius, step_limit)
        for i in xrange(len(x)):
            self.assertTrue(x[i] > radius)
            self.assertTrue(x[i] < 1. - radius)


        for i in xrange(npoints):
            ci = Boxcar(x[i], radius)

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)




    def test3_pack_metro_md(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])


        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])


        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100

        success_steps0 = boxcar.pack.metro_md(x0, radius, step_limit,
            randSeed)
        success_steps1 = boxcar.pack.metro_md(x1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))




    def test4_pack_metro_md(self):
        """

        Test case when all steps are successful

        """

        npoints = 500
        radius = 0.0
        step_limit = 10 ** 3
        x = boxcar.pack.poisson_point(npoints=npoints)
        success_steps = boxcar.pack.metro_md(x, radius, step_limit)

        self.assertTrue(success_steps == step_limit)






    """ pack_metro_pd() """


    def test1_pack_metro_pd(self):
        """

        Test case with npoints small

        """

        npoints = 2
        radius = 0.05
        step_limit = 10 ** 2
        x  = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = boxcar.pack.metro_pd(x, radius, step_limit)

        for i in xrange(npoints):
            ci = Boxcar(x[i], radius[i])

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius[k])
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test2_pack_metro_pd(self):
        """

        Test case with npoints large

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 3
        x = boxcar.pack.grid_md(npoints=npoints, radius=radius)
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))
        success_steps = boxcar.pack.metro_pd(x, radius, step_limit)

        for i in xrange(npoints):
            ci = Boxcar(x[i], radius[i])

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius[k])
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(success_steps > 0)


    def test3_pack_metro_pd(self):
        """

        Test case random seed

        """

        x0 = np.ascontiguousarray([0.1, 0.3, 0.5])


        x1 = np.ascontiguousarray([0.1, 0.3, 0.5])


        radius = 0.05
        step_limit = 10 ** 3
        randSeed = 100
        npoints = 3
        radius = np.ascontiguousarray(0.05 * np.ones(npoints))

        success_steps0 = boxcar.pack.metro_pd(x0, radius, step_limit,
            randSeed)
        success_steps1 = boxcar.pack.metro_pd(x1, radius, step_limit,
            randSeed )

        self.assertTrue(np.allclose(x0, x1))




    def test4_pack_metro_pd(self):
        """

        Test case when all steps are successful

        """

        npoints = 500
        radius = 0.0
        radius = np.ascontiguousarray(radius * np.ones(npoints))
        step_limit = 10 ** 3
        x = boxcar.pack.poisson_point(npoints=npoints)
        success_steps = boxcar.pack.metro_pd(x, radius, step_limit)

        self.assertTrue(success_steps == step_limit)







    """ pack_rsa_md() """

    def test1_pack_rsa_md(self):
        """

        Test case with npoints small

        """

        npoints = 2
        radius = 0.05
        step_limit = 10 ** 2

        x = boxcar.pack.rsa_md(npoints, radius, step_limit)


        for i in xrange(npoints):
            ci = Boxcar(x[i], radius)

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius)
                F = ci.overlap_potential(ck)
                if i != k:
                    self.assertTrue(F >= 1. or np.allclose(F, 1.))

        self.assertTrue(npoints == len(x))


    def test2_pack_rsa_md(self):
        """

        Test case with npoints large

        """

        npoints = 5
        radius = 0.05
        step_limit = 10 ** 4

        x = boxcar.pack.rsa_md(npoints, radius, step_limit)


        for i in xrange(npoints):
            ci = Boxcar(x[i], radius)

            # inside the container
            H = ci.container_potential()
            self.assertTrue(H >= 1.)

            # overlap with others
            for k in xrange(npoints):
                ck = Boxcar(x[k], radius)
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

        x0 = boxcar.pack.rsa_md(npoints, radius, step_limit,
            randSeed)
        x1= boxcar.pack.rsa_md(npoints, radius, step_limit,
            randSeed)

        self.assertTrue(np.allclose(x0, x1))





if __name__ == '__main__':
    print 'Running unit tests for boxcar.so'
    unittest.main()