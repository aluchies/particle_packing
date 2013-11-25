import numpy as np
cimport numpy as np
from itertools import product
import random
import sys



cdef extern from "c/circle_metro.c":
    unsigned int metro_md_2d(double *x, double *y,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_pd_2d(double *x, double *y,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)


cdef extern from "c/circle_rsa.c":
    size_t gen_pts_rsa_2d(double *x, double *y,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)







def pack_metro_md(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    double radius, int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse size hard disks.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
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



    success_steps =  metro_md_2d(&x[0], &y[0], radius, npoints, step_limit, rseed)


    return success_steps







def pack_metro_pd(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] radius not None,
    int step_limit, rand_seed=None):

    """Metropolis algorithm for poly-disperse size hard disks.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
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

    success_steps =  metro_pd_2d(&x[0], &y[0], &radius[0], npoints, step_limit, rseed)


    return success_steps










def pack_rsa_md(int npoints, double radius, int step_limit, rand_seed=None):
    """RSA algorithm for mono-disperse size hard disks.

    Keyword arguments:
    npoints -- number of spheres positions to generate
    radius -- sphere radius
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates

    """


    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints, ))



    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    cdef int valid_pts
    valid_pts = gen_pts_rsa_2d(&x[0], &y[0], npoints, radius, step_limit, rseed)

    return x[:valid_pts], y[:valid_pts]







def pack_grid_md(int npoints=5, double radius=0.05):
    """Algorithm for placing mono-disperse size hard disks on a square grid.
    May be used to generate initial positions for metropolis algorithm.

    Keyword arguments:
    npoints -- number of sphere positions to generate
    radius -- sphere radius

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates


    """

    space = 1.1 * 2. * radius
    ppdim = np.floor(1. / space)
    grid = np.arange(ppdim) * space + space / 2.
    xlist = product(grid, grid)
    xlist = list(xlist)

    np.random.shuffle(xlist)

    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints,))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints,))

    if npoints == 0:

        return x, y


    elif len(xlist) >= npoints:

        X = np.asarray(xlist[0:npoints])
        x = X[:,0].copy()
        y = X[:,1].copy()

        return x, y


    else:
        raise ValueError('Implement case when grid provides too few points')








def pack_uniform(int npoints=5):
    """Generate i.i.d. uniformly distributed positions.

    Keyword arguements:
    npoints -- number of independent points to generate

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates

    """

    x = np.ascontiguousarray(np.random.rand(npoints))
    y = np.ascontiguousarray(np.random.rand(npoints))

    return x, y