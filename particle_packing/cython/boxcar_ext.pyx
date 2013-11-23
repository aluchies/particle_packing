import numpy as np
cimport numpy as np
from itertools import product
import random
import sys



cdef extern from "c/boxcar_metro.c":
    unsigned int metro_md_1d(double *x,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_pd_1d(double *x,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)


cdef extern from "c/boxcar_rsa.c":
    size_t gen_pts_rsa_1d(double *x,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)







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



    success_steps =  metro_md_1d(&x[0], radius, npoints, step_limit, rseed)


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

    success_steps =  metro_pd_1d(&x[0], &radius[0], npoints, step_limit, rseed)


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
    """Generate mono-disperse size boxcar positions.

    Keyword arguements:

    Return values:
    x -- array of x-coordinates


    """

    x = np.ascontiguousarray(np.random.rand(npoints))

    return x