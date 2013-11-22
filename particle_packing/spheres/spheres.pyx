import numpy as np
cimport numpy as np
from itertools import product
import random
import sys



cdef extern from "src/spheres_metro.c":
    unsigned int metro_md_3d(double *x, double *y, double *z,
    double radius, size_t npoints, int step_limit,
    unsigned long randSeed)

    unsigned int metro_pd_3d(double *x, double *y, double *z,
    double *radius, size_t npoints, int step_limit,
    unsigned long randSeed)


cdef extern from "src/spheres_rsa.c":
    size_t gen_pts_rsa_3d(double *x, double *y, double *z,
    size_t npoints, double radius, int step_limit,
    unsigned long randSeed)







def pack_metro_md(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] z not None,
    double radius, int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse spheres.

    Keyword arguments:
    x --
    y --
    z --

    Return values:
    success_steps --

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







def pack_metro_pd(
    np.ndarray[double, ndim=1, mode="c"] x not None,
    np.ndarray[double, ndim=1, mode="c"] y not None,
    np.ndarray[double, ndim=1, mode="c"] z not None,
    np.ndarray[double, ndim=1, mode="c"] radius not None,
    int step_limit, rand_seed=None):

    """Metropolis algorithm for mono-disperse spheres.

    Keyword arguments:
    x -- array of x coordinates
    y -- array of y coordinates
    z -- array of z coordinates
    radius -- array of radii

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










def pack_rsa_md(int npoints, double radius, int step_limit, rand_seed=None):
    """
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







def pack_grid_md(int npoints=5, double radius=0.05):
    """
    """

    space = 1.1 * 2. * radius
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








def pack_uniform(int npoints=5):
    """Generate monodisperse sphere locations in n-dimensional cube volume.

    Keyword arguements:
    ndim -- number of dimensions
    L -- container length (accepts numeric value, list, tuple, ndarray)
    npoints -- numper of points
    nsim -- number of center point configurations to generate

    Return values:
    x -- center point locations

    """

    x = np.ascontiguousarray(np.random.rand(npoints))
    y = np.ascontiguousarray(np.random.rand(npoints))
    z = np.ascontiguousarray(np.random.rand(npoints))

    return x, y, z