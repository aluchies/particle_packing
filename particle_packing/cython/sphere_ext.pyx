import numpy as np
cimport numpy as np
from itertools import product
import random
import sys



cdef extern from "c/sphere.c":
    double sphere_overlap(double *rA, double radiiA, double *rB, double radiiB)

    void sphere_collection_overlap(double *x, double *y, double *z, double *radii,
    size_t npoints, double *stats)

    size_t gen_pts_rsa_3d(double *x, double *y, double *z,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)

    size_t gen_pts_rsa_3d_2(double *x, double *y, double *z,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)

    unsigned int metro_md_3d(double *x, double *y, double *z,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_md_3d_2(double *x, double *y, double *z,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_pd_3d(double *x, double *y, double *z,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)








def overlap_potential(r1, radii1, r2, radii2):
    """

    Overlap potential function provides a distance measure
    for spheres A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of sphere A
    radiiA -- radii of sphere A
    rB -- center of sphere B
    radiiB -- radii of sphere B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """

    # Input argument checking.

    cdef double radiiA, radiiB, F
    cdef np.ndarray[double, ndim=1, mode="c"] rA
    cdef np.ndarray[double, ndim=1, mode="c"] rB


    r1 = np.asarray(r1).flatten()
    if len(r1) != 3:
        raise ValueError('input error for r1')
    rA = np.ascontiguousarray(r1, dtype=np.float64)


    r2 = np.asarray(r2).flatten()
    if len(r2) != 3:
        raise ValueError('input error for r2')
    rB = np.ascontiguousarray(r2, dtype=np.float64)


    radii1 = np.asarray(radii1).flatten()
    if len(radii1) != 1:
        raise ValueError('input error for radii1')
    radiiA = radii1[0]


    radii2 = np.asarray(radii2).flatten()
    if len(radii2) != 1:
        raise ValueError('input error for radii2')
    radiiB = radii2[0]



    F = sphere_overlap(&rA[0], radiiA, &rB[0], radiiB)


    return F





def collection_overlap_potential(x1, y1, z1, radii1):
    """

    Minimum overlap potential for a collection of spheres.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    x --
    y --
    z --
    radii --

    Return values:
    F -- minimum overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """

    # Input argument checking.

    cdef np.ndarray[double, ndim=1, mode="c"] xA
    cdef np.ndarray[double, ndim=1, mode="c"] yA
    cdef np.ndarray[double, ndim=1, mode="c"] zA
    cdef np.ndarray[double, ndim=1, mode="c"] radiiA
    cdef np.ndarray[double, ndim=1, mode="c"] statsA

    cdef size_t npoints
    cdef double F


    x1 = np.asarray(x1).flatten()
    xA = np.ascontiguousarray(x1, dtype=np.float64)
    npoints = len(xA)


    y1 = np.asarray(y1).flatten()
    if len(y1) != npoints:
        raise ValueError('input error for y')
    yA = np.ascontiguousarray(y1, dtype=np.float64)

    z1 = np.asarray(z1).flatten()
    if len(z1) != npoints:
        raise ValueError('input error for z')
    zA = np.ascontiguousarray(z1, dtype=np.float64)

    radii1 = np.asarray(radii1).flatten()
    if len(radii1) != npoints:
        raise ValueError('input error for radii')
    radiiA = np.ascontiguousarray(radii1, dtype=np.float64)

    statsA = np.ascontiguousarray(np.asarray([0, 0, 0]), dtype=np.float64)


    sphere_collection_overlap(&xA[0], &yA[0], &zA[0], &radiiA[0], npoints, &statsA[0])


    return statsA






def metro_md(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] z not None,
    double radius, int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse size spheres.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
    z -- array of z coordinates
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    success_steps -- number of metropolis steps that were successful

    """

    cdef unsigned int success_steps
    cdef size_t npoints 


    npoints = len(x)


    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)



    success_steps =  metro_md_3d(&x[0], &y[0], &z[0], radius, npoints, step_limit, rseed)


    return success_steps



def metro_md_2(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] z not None,
    double radius, int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse size spheres.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
    z -- array of z coordinates
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    success_steps -- number of metropolis steps that were successful

    """

    cdef unsigned int success_steps
    cdef size_t npoints 


    npoints = len(x)


    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)



    success_steps =  metro_md_3d_2(&x[0], &y[0], &z[0], radius, npoints, step_limit, rseed)


    return success_steps







def metro_pd(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] z not None,
    np.ndarray[double, ndim=1, mode="c"] radius not None,
    int step_limit, rand_seed=None):

    """Metropolis algorithm for poly-disperse size spheres.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
    z -- array of z coordinates
    radius -- array of sphere radii
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    success_steps -- number of metropolis steps that were successful

    """

    cdef unsigned int success_steps
    cdef size_t npoints 


    npoints = len(x)


    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    success_steps =  metro_pd_3d(&x[0], &y[0], &z[0], &radius[0], npoints, step_limit, rseed)


    return success_steps










def rsa_md(int npoints, double radius, int step_limit, rand_seed=None):
    """RSA algorithm for mono-disperse size spheres.

    Keyword arguments:
    npoints -- number of spheres positions to generate
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates
    z -- array of z-coordinates

    """


    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] z = np.zeros((npoints, ))



    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    cdef int valid_pts
    valid_pts = gen_pts_rsa_3d(&x[0], &y[0], &z[0], npoints, radius, step_limit, rseed)

    return x[:valid_pts], y[:valid_pts], z[:valid_pts]





def rsa_md_2(int npoints, double radius, int step_limit, rand_seed=None):
    """RSA algorithm for mono-disperse size spheres.

    Keyword arguments:
    npoints -- number of spheres positions to generate
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates
    z -- array of z-coordinates

    """


    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] z = np.zeros((npoints, ))



    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    cdef int valid_pts
    valid_pts = gen_pts_rsa_3d_2(&x[0], &y[0], &z[0], npoints, radius, step_limit, rseed)

    return x[:valid_pts], y[:valid_pts], z[:valid_pts]







def grid_md(int npoints=5, double radius=0.05):
    """Algorithm for placing mono-disperse size spheres on a square grid. May be used to generate initial positions for metropolis algorithm.

    Keyword arguments:
    npoints -- number of sphere positions to generate
    radius -- sphere radius

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates
    z -- array of z-coordinates


    """

    space = 1.01 * 2. * radius
    ppdim = np.floor(1. / space)
    grid = np.arange(ppdim) * space + space / 2.
    xlist = product(grid, grid, grid)
    xlist = list(xlist)

    np.random.shuffle(xlist)

    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints,))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints,))
    cdef np.ndarray[double, ndim=1, mode="c"] z = np.zeros((npoints,))

    if npoints == 0:

        return x, y, z


    elif len(xlist) >= npoints:

        X = np.asarray(xlist[0:npoints])
        x = X[:,0].copy()
        y = X[:,1].copy()
        z = X[:,2].copy()

        return x, y, z


    else:
        raise ValueError('Implement case when grid provides too few points')








def poisson_point(int npoints=5):
    """Generate i.i.d. uniformly distributed positions.

    Keyword arguements:
    npoints -- number of independent points to generate

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates
    z -- array of z-coordinates

    """

    x = np.ascontiguousarray(np.random.rand(npoints))
    y = np.ascontiguousarray(np.random.rand(npoints))
    z = np.ascontiguousarray(np.random.rand(npoints))

    return x, y, z