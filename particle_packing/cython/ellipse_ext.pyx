import numpy as np
cimport numpy as np

cdef extern from "c/ellipse_overlap.c":
    double ellipse_overlap(double *rA, double *radiiA, double phiA, double *rB, double *radiiB, double phiB)
    void characteristic_ellipse_matrix(double *X, double *R, double phi, double exponent)


def overlap_potential(r1, radii1, phi1, r2, radii2, phi2):
    """Metropolis algorithm for mono-disperse size hard disks.

    Keyword arguments:
    r1 --
    radii1 --
    phi1 --
    r2 --
    radii2 --
    phi2 --

    Return values:
    F --

    """
    cdef np.ndarray[double, ndim=1, mode="c"] rA = np.ascontiguousarray(r1.flatten(), dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiA = np.ascontiguousarray(radii1.flatten(), dtype=np.float64)
    cdef double phiA = float(phi1)
    
    cdef np.ndarray[double, ndim=1, mode="c"] rB = np.ascontiguousarray(r2.flatten(), dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] radiiB = np.ascontiguousarray(radii2.flatten(), dtype=np.float64)
    cdef double phiB = float(phi2)

    cdef double F


    F = ellipse_overlap(&rA[0], &radiiA[0], phiA, &rB[0], &radiiB[0], phiB)

    return F


def char_mat(radii, angle):

    cdef np.ndarray[double, ndim=1, mode="c"] rad
    cdef np.ndarray[double, ndim=2, mode="c"] X
    cdef double phi
    cdef exponent

    rad = np.ascontiguousarray(radii.flatten(), dtype=np.float64)
    X = np.ascontiguousarray(np.zeros((2,2)), dtype=np.float64)
    phi = float(angle)
    exponent = 1

    characteristic_ellipse_matrix(&X[0,0], &rad[0], phi, exponent)

    return X