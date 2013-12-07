import numpy as np
cimport numpy as np

cdef extern from "c/ellipse/ellipse_overlap.c":
    double ellipse_overlap(double *rA, double *radiiA, double phiA, double *rB, double *radiiB, double phiB)


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





    """Input argument checking."""
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