import numpy as np
cimport numpy as np
import random
import sys

cdef extern from "c/ellipse.c":
    double ellipse_overlap(double *rA, double *radiiA, double phiA, double *rB, double *radiiB, double phiB)
    size_t gen_pts_rsa_aligned(double *x, double *y,
    size_t npoints, double *radius, double phi, int step_limit,
    unsigned long randSeed)


def overlap_potential(r1, radii1, phi1, r2, radii2, phi2):
    """

    Overlap potential function provides a distance measure
    for ellipses A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of ellipse A
    radiiA -- radii of ellipse A
    phiA -- rotation angle of ellipse A
    rB -- center of ellipse B
    radiiB -- radii of ellipse B
    phiB -- rotation angle of ellipse B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """





    # Input argument checking.
    r1 = np.asarray(r1.flatten())
    if len(r1) != 2:
        raise ValueError('input error for r1')

    r2 = np.asarray(r2.flatten())
    if len(r2) != 2:
        raise ValueError('input error for r2')


    radii1 = np.asarray(radii1.flatten())
    if len(radii1) != 2:
        raise ValueError('input error for radii1')

    radii2 = np.asarray(radii2.flatten())
    if len(radii2) != 2:
        raise ValueError('input error for radii2')

    phi1 = float(phi1)
    phi2 = float(phi2)



    cdef np.ndarray[double, ndim=1, mode="c"] rA = np.ascontiguousarray(r1, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiA = np.ascontiguousarray(radii1, dtype=np.float64)
    cdef double phiA = float(phi1)
    
    cdef np.ndarray[double, ndim=1, mode="c"] rB = np.ascontiguousarray(r2, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiB = np.ascontiguousarray(radii2, dtype=np.float64)
    cdef double phiB = float(phi2)

    cdef double F


    F = ellipse_overlap(&rA[0], &radiiA[0], phiA, &rB[0], &radiiB[0], phiB)

    return F










def pack_rsa_md_aligned(npoints, radius, phi, step_limit, rand_seed=None):
    """RSA algorithm for mono-disperse size hard ellipses.

    Keyword arguments:
    npoints -- number of spheres positions to generate
    radius -- ellipse radii
    phi --
    step_limit -- number of steps in metropolis algorithm
    rand_seed -- seed for the random number generator

    Return values:
    x -- array of x-coordinates
    y -- array of y-coordinates

    """

    # Input checking
    radius = np.asarray(radius).flatten()
    if len(radius) != 2:
        raise ValueError('input error for radius')
    cdef np.ndarray[double, ndim=1, mode="c"] radii = np.ascontiguousarray(radius)

    cdef double alpha = float(phi)




    cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros((npoints, ))
    cdef np.ndarray[double, ndim=1, mode="c"] y = np.zeros((npoints, ))



    # take care of random seed
    cdef unsigned long rseed
    if rand_seed is None:
        rseed = random.randint(0, sys.maxint)
    else:
        rseed = long(rand_seed)

    cdef int valid_pts
    valid_pts = gen_pts_rsa_aligned(&x[0], &y[0], npoints, &radii[0], alpha, step_limit, rseed)

    return x[:valid_pts], y[:valid_pts]