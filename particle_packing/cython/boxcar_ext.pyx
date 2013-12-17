import numpy as np
cimport numpy as np
from itertools import product
import random
import sys




cdef extern from "c/boxcar.c":
    double boxcar_overlap(double rA, double radiiA, double rB, double radiiB)

    size_t gen_pts_rsa_1d(double *x,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)

    unsigned int metro_md_1d(double *x,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_pd_1d(double *x,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)







def overlap_potential(r1, radii1, r2, radii2):
    """

    Overlap potential function provides a distance measure
    for boxcars A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of boxcar A
    radiiA -- radii of boxcar A
    rB -- center of boxcar B
    radiiB -- radii of boxcar B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """

    # Input argument checking.

    cdef double rA, radiiA, rB, radiiB, F

    r1 = np.asarray(r1).flatten()
    if len(r1) != 1:
        raise ValueError('input error for r1')
    rA = r1[0]

    r2 = np.asarray(r2).flatten()
    if len(r2) != 1:
        raise ValueError('input error for r2')
    rB = r2[0]


    radii1 = np.asarray(radii1).flatten()
    if len(radii1) != 1:
        raise ValueError('input error for radii1')
    radiiA = radii1[0]

    radii2 = np.asarray(radii2).flatten()
    if len(radii2) != 1:
        raise ValueError('input error for radii2')
    radiiB = radii2[0]

    F = boxcar_overlap(rA, radiiA, rB, radiiB)

    return F








def pack_metro_md(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    double radius, int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse size boxcar functions.

    Keyword arguments:
    x -- array of x coordinates
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



    success_steps = metro_md_1d(&x[0], radius, npoints, step_limit, rseed)


    return success_steps







def pack_metro_pd(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] radius not None,
    int step_limit, rand_seed=None):

    """Metropolis algorithm for poly-disperse size spheres.

    Keyword arguments:
    x -- array of x coordinates
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

    success_steps = metro_pd_1d(&x[0], &radius[0], npoints, step_limit, rseed)


    return success_steps










def pack_rsa_md(int npoints, double radius, int step_limit, rand_seed=None):
    """RSA algorithm for mono-disperse size spheres.

    Keyword arguments:
    npoints -- number of spheres positions to generate
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    x -- array of x-coordinates

    """


    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints, ))



    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    cdef int valid_pts
    valid_pts = gen_pts_rsa_1d(&x[0], npoints, radius, step_limit, rseed)

    return x[:valid_pts]







def pack_grid_md(int npoints=5, double radius=0.05):
    """Algorithm for placing mono-disperse size boxcars on grid. May be used to generate
    initial positions for metropolis algorithm.

    Keyword arguments:
    npoints -- number of sphere positions to generate
    radius -- sphere radius

    Return values:
    x -- array of x-coordinates


    """

    space = 1.1 * 2. * radius
    ppdim = np.floor(1. / space)
    xlist = np.arange(ppdim) * space + space / 2.

    np.random.shuffle(xlist)

    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints,))

    if npoints == 0:

        return x


    elif len(xlist) >= npoints:

        X = np.asarray(xlist[0:npoints])
        x = X.copy()

        return x


    else:
        raise ValueError('Implement case when grid provides too few points')








def pack_uniform(int npoints=5):
    """Generate i.i.d. uniformly distributed positions.

    Keyword arguements:
    npoints -- number of independent points to generate

    Return values:
    x -- array of x-coordinates


    """

    x = np.ascontiguousarray(np.random.rand(npoints))

    return x