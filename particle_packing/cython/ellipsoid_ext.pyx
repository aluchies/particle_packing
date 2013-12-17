import numpy as np
cimport numpy as np

cdef extern from "c/ellipsoid.c":
    double ellipsoid_overlap(double *rA, double *radiiA, double phiA, double *rotaxA, double *rB, double *radiiB, double phiB, double *rotaxB)



def overlap_potential(r1, radii1, phi1, rotax1, r2, radii2, phi2, rotax2):
    """

    Overlap potential function provides a distance measure
    for ellipsoids A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of A
    radiiA -- radii of A
    phiA -- rotation angle of A
    rB -- center of B
    radiiB -- radii of B
    phiB -- rotation angle of B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """




    """Input argument checking."""
    r1 = np.asarray(r1.flatten())
    if len(r1) != 3:
        raise ValueError('input error for r1')

    r2 = np.asarray(r2.flatten())
    if len(r2) != 3:
        raise ValueError('input error for r2')


    radii1 = np.asarray(radii1.flatten())
    if len(radii1) != 3:
        raise ValueError('input error for radii1')

    radii2 = np.asarray(radii2.flatten())
    if len(radii2) != 3:
        raise ValueError('input error for radii2')

    phi1 = float(phi1)
    phi2 = float(phi2)

    rotax1 = np.asarray(rotax1.flatten())
    if len(rotax1) != 3:
        raise ValueError('input error for rotax1')
    if np.allclose(rotax1, 1):
        raise ValueError('input error for rotax1')

    rotax2 = np.asarray(rotax2.flatten())
    if len(rotax2) != 3:
        raise ValueError('input error for rotax2')
    if np.allclose(rotax2, 1):
        raise ValueError('input error for rotax2')





    cdef np.ndarray[double, ndim=1, mode="c"] rA = np.ascontiguousarray(r1, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiA = np.ascontiguousarray(radii1, dtype=np.float64)
    cdef double phiA = float(phi1)
    cdef np.ndarray[double, ndim=1, mode="c"] rotaxA = np.ascontiguousarray(rotax1, dtype=np.float64)

    cdef np.ndarray[double, ndim=1, mode="c"] rB = np.ascontiguousarray(r2, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiB = np.ascontiguousarray(radii2, dtype=np.float64)
    cdef double phiB = float(phi2)
    cdef np.ndarray[double, ndim=1, mode="c"] rotaxB = np.ascontiguousarray(rotax2, dtype=np.float64)


    cdef double F


    F = ellipsoid_overlap(&rA[0], &radiiA[0], phiA, &rotaxA[0], &rB[0], &radiiB[0], phiB, &rotaxB[0])

    return F


