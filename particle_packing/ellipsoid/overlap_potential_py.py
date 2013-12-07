import numpy as np

def overlap_potential_py(rA, radiiA, phiA, rotaxA, rB, radiiB, phiB, rotaxB):
    """

    Overlap potential function (Python version) provides a distance measure
    for ellipsoids A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of ellipsoid A
    radiiA -- radii of ellipsoid A
    phiA -- rotation angle of ellipsoid A
    rotaxA --
    rB -- center of ellipsoid B
    radiiB -- radii of ellipsoid B
    phiB -- rotation angle of ellipsoid B
    rotaxB --

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """






    """Input argument checking."""

    rA = np.asarray(rA.flatten())
    if len(rA) != 3:
        raise ValueError('input error for rA')
    rA = np.matrix(rA).T

    rB = np.asarray(rB.flatten())
    if len(rB) != 3:
        raise ValueError('input error for rB')
    rB = np.matrix(rB).T


    radiiA = np.asarray(radiiA.flatten())
    if len(radiiA) != 3:
        raise ValueError('input error for radiiA')

    radiiB = np.asarray(radiiB.flatten())
    if len(radiiB) != 3:
        raise ValueError('input error for radiiB')

    phiA = float(phiA)
    phiB = float(phiB)

    rotaxA = np.asarray(rotaxA.flatten())
    if len(rotaxA) != 3:
        raise ValueError('input error for rotaxA')
    if np.allclose(rotaxA, 1):
        raise ValueError('input error for rotaxA')

    rotaxB = np.asarray(rotaxB.flatten())
    if len(rotaxB) != 3:
        raise ValueError('input error for rotaxB')
    if np.allclose(rotaxB, 1):
        raise ValueError('input error for rotaxB')







    """

    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and 
    rotation axis rotax, and O is diagonal matrix with radii values along
    main diagonal.

    q = [s, p] = [cos(phi / 2), sin(phi / 2) phi^]
    ||q||^2 = s^2 + ||p||^2 = 1
    Q = 2 [ p p^T - s p p^T + (s^2 - 1 / 2) I ]
    
    What's needed is XA^-1 = Q^T * O^2 * Q

    """

    O = np.matrix(np.diag(radiiA ** 2))

    rotaxA = np.matrix(rotaxA).T
    I = np.matrix(np.eye(3))
    s = np.cos(phiA)
    p = np.sin(phiA) * rotaxA
    Q = 2. * (p * p.T * (1. - s) + (s ** 2 - 0.5) * I)


    XA = Q.T * O * Q





    """

    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and 
    rotation axis rotax, and O is diagonal matrix with radii values along
    main diagonal.

    q = [s, p] = [cos(phi / 2), sin(phi / 2) phi^]
    ||q||^2 = s^2 + ||p||^2 = 1
    Q = 2 [ p p^T - s p p^T + (s^2 - 1 / 2) I ]
    
    What's needed is XB^1/2 = Q^T * O^2 * Q

    """

    O = np.matrix(np.diag(radiiB ** -1))

    rotaxB = np.matrix(rotaxB).T
    I = np.matrix(np.eye(3))
    s = np.cos(phiB)
    p = np.sin(phiB) * rotaxB
    Q = 2. * (p * p.T * (1. - s) + (s ** 2 - 0.5) * I)


    XB = Q.T * O * Q











    # Find rAB
    rAB = rB - rA

    # Find A_AB
    A_AB = XB * XA * XB

    # Find a_AB
    a_AB = XB * rAB



    coeffs = h_coeffs(A_AB, a_AB)
    roots = np.roots(coeffs)

    F = None
    for r in roots:
        if np.isreal(r) and (0. <= r <= 1.):
            #print 'max in (0, 1)'
            F =  f_AB(A_AB, a_AB, np.real(r))

    if F is None:
        if f_AB(A_AB, a_AB, 0) > f_AB(A_AB, a_AB, 1):
            #print 'max at lambda = 0'
            F = f_AB(A_AB, a_AB, 0)
        else:
            #print 'max at lambda = 1'
            F = f_AB(A_AB, a_AB, 1)



    return F



def p_AB(A, b, L):
    """

    Numerator of function f.

    """


    a11, a12, a13, a21, a22, a23, a31, a32, a33 = \
    A[0, 0], A[0, 1], A[0, 2], \
    A[1, 0], A[1, 1], A[1, 2], \
    A[2, 0], A[2, 1], A[2, 2]

    b1, b2, b3 = b[0, 0], b[1, 0], b[2, 0]


    return L**4*(-a11*a22*b3**2 + a11*a23*b2*b3 + a11*a32*b2*b3 - a11*a33*b2**2 + a11*b2**2 + a11*b3**2 + a12*a21*b3**2 - a12*a23*b1*b3 - a12*a31*b2*b3 + a12*a33*b1*b2 - a12*b1*b2 - a13*a21*b2*b3 + a13*a22*b1*b3 + a13*a31*b2**2 - a13*a32*b1*b2 - a13*b1*b3 - a21*a32*b1*b3 + a21*a33*b1*b2 - a21*b1*b2 + a22*a31*b1*b3 - a22*a33*b1**2 + a22*b1**2 + a22*b3**2 - a23*a31*b1*b2 + a23*a32*b1**2 - a23*b2*b3 - a31*b1*b3 - a32*b2*b3 + a33*b1**2 + a33*b2**2 - b1**2 - b2**2 - b3**2) + L**3*(3*a11*a22*b3**2 - 3*a11*a23*b2*b3 - 3*a11*a32*b2*b3 + 3*a11*a33*b2**2 - 2*a11*b2**2 - 2*a11*b3**2 - 3*a12*a21*b3**2 + 3*a12*a23*b1*b3 + 3*a12*a31*b2*b3 - 3*a12*a33*b1*b2 + 2*a12*b1*b2 + 3*a13*a21*b2*b3 - 3*a13*a22*b1*b3 - 3*a13*a31*b2**2 + 3*a13*a32*b1*b2 + 2*a13*b1*b3 + 3*a21*a32*b1*b3 - 3*a21*a33*b1*b2 + 2*a21*b1*b2 - 3*a22*a31*b1*b3 + 3*a22*a33*b1**2 - 2*a22*b1**2 - 2*a22*b3**2 + 3*a23*a31*b1*b2 - 3*a23*a32*b1**2 + 2*a23*b2*b3 + 2*a31*b1*b3 + 2*a32*b2*b3 - 2*a33*b1**2 - 2*a33*b2**2 + b1**2 + b2**2 + b3**2) + L**2*(-3*a11*a22*b3**2 + 3*a11*a23*b2*b3 + 3*a11*a32*b2*b3 - 3*a11*a33*b2**2 + a11*b2**2 + a11*b3**2 + 3*a12*a21*b3**2 - 3*a12*a23*b1*b3 - 3*a12*a31*b2*b3 + 3*a12*a33*b1*b2 - a12*b1*b2 - 3*a13*a21*b2*b3 + 3*a13*a22*b1*b3 + 3*a13*a31*b2**2 - 3*a13*a32*b1*b2 - a13*b1*b3 - 3*a21*a32*b1*b3 + 3*a21*a33*b1*b2 - a21*b1*b2 + 3*a22*a31*b1*b3 - 3*a22*a33*b1**2 + a22*b1**2 + a22*b3**2 - 3*a23*a31*b1*b2 + 3*a23*a32*b1**2 - a23*b2*b3 - a31*b1*b3 - a32*b2*b3 + a33*b1**2 + a33*b2**2) + L*(a11*a22*b3**2 - a11*a23*b2*b3 - a11*a32*b2*b3 + a11*a33*b2**2 - a12*a21*b3**2 + a12*a23*b1*b3 + a12*a31*b2*b3 - a12*a33*b1*b2 + a13*a21*b2*b3 - a13*a22*b1*b3 - a13*a31*b2**2 + a13*a32*b1*b2 + a21*a32*b1*b3 - a21*a33*b1*b2 - a22*a31*b1*b3 + a22*a33*b1**2 + a23*a31*b1*b2 - a23*a32*b1**2)



def q_AB(A, b, L):
    """

    Denominator of function f.

    """


    a11, a12, a13, a21, a22, a23, a31, a32, a33 = \
    A[0, 0], A[0, 1], A[0, 2], \
    A[1, 0], A[1, 1], A[1, 2], \
    A[2, 0], A[2, 1], A[2, 2]

    b1, b2, b3 = b[0, 0], b[1, 0], b[2, 0]


    return L**3*(-a11*a22*a33 + a11*a22 + a11*a23*a32 + a11*a33 - a11 + a12*a21*a33 - a12*a21 - a12*a23*a31 - a13*a21*a32 + a13*a22*a31 - a13*a31 + a22*a33 - a22 - a23*a32 - a33 + 1) + L**2*(3*a11*a22*a33 - 2*a11*a22 - 3*a11*a23*a32 - 2*a11*a33 + a11 - 3*a12*a21*a33 + 2*a12*a21 + 3*a12*a23*a31 + 3*a13*a21*a32 - 3*a13*a22*a31 + 2*a13*a31 - 2*a22*a33 + a22 + 2*a23*a32 + a33) + L*(-3*a11*a22*a33 + a11*a22 + 3*a11*a23*a32 + a11*a33 + 3*a12*a21*a33 - a12*a21 - 3*a12*a23*a31 - 3*a13*a21*a32 + 3*a13*a22*a31 - a13*a31 + a22*a33 - a23*a32) + a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31



def f_AB(A, b, L):
    """

    """


    return p_AB(A, b, L) / q_AB(A, b, L)



def h_coeffs(A, b):
    """

    Numerator of the derivative of the function f.

    """

    a11, a12, a13, a21, a22, a23, a31, a32, a33 = \
    A[0, 0], A[0, 1], A[0, 2], \
    A[1, 0], A[1, 1], A[1, 2], \
    A[2, 0], A[2, 1], A[2, 2]

    b1, b2, b3 = b[0, 0], b[1, 0], b[2, 0]


    coeff = [a11**2*a22*b2**2 - a11**2*b2**2 - a11*a12*a21*b2**2 - a11*a12*a22*b1*b2 + a11*a12*b1*b2 - a11*a21*a22*b1*b2 + a11*a21*b1*b2 + a11*a22**2*b1**2 - 2*a11*a22*b1**2 - 2*a11*a22*b2**2 + a11*b1**2 + 2*a11*b2**2 + a12**2*a21*b1*b2 + a12*a21**2*b1*b2 - a12*a21*a22*b1**2 + a12*a21*b1**2 + a12*a21*b2**2 + a12*a22*b1*b2 - a12*b1*b2 + a21*a22*b1*b2 - a21*b1*b2 - a22**2*b1**2 + 2*a22*b1**2 + a22*b2**2 - b1**2 - b2**2, -4*a11**2*a22*b2**2 + 2*a11**2*b2**2 + 4*a11*a12*a21*b2**2 + 4*a11*a12*a22*b1*b2 - 2*a11*a12*b1*b2 + 4*a11*a21*a22*b1*b2 - 2*a11*a21*b1*b2 - 4*a11*a22**2*b1**2 + 6*a11*a22*b1**2 + 6*a11*a22*b2**2 - 2*a11*b1**2 - 2*a11*b2**2 - 4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + 4*a12*a21*a22*b1**2 - 4*a12*a21*b1**2 - 4*a12*a21*b2**2 - 2*a12*a22*b1*b2 - 2*a21*a22*b1*b2 + 2*a22**2*b1**2 - 2*a22*b1**2 - 2*a22*b2**2, 6*a11**2*a22*b2**2 - a11**2*b2**2 - 6*a11*a12*a21*b2**2 - 6*a11*a12*a22*b1*b2 + a11*a12*b1*b2 - 6*a11*a21*a22*b1*b2 + a11*a21*b1*b2 + 6*a11*a22**2*b1**2 - 6*a11*a22*b1**2 - 6*a11*a22*b2**2 + a11*b1**2 + 6*a12**2*a21*b1*b2 + 6*a12*a21**2*b1*b2 - 6*a12*a21*a22*b1**2 + 5*a12*a21*b1**2 + 5*a12*a21*b2**2 + a12*a22*b1*b2 + a12*b1*b2 + a21*a22*b1*b2 + a21*b1*b2 - a22**2*b1**2 + a22*b2**2, -4*a11**2*a22*b2**2 + 4*a11*a12*a21*b2**2 + 4*a11*a12*a22*b1*b2 + 4*a11*a21*a22*b1*b2 - 4*a11*a22**2*b1**2 + 2*a11*a22*b1**2 + 2*a11*a22*b2**2 - 4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + 4*a12*a21*a22*b1**2 - 2*a12*a21*b1**2 - 2*a12*a21*b2**2, a11**2*a22*b2**2 - a11*a12*a21*b2**2 - a11*a12*a22*b1*b2 - a11*a21*a22*b1*b2 + a11*a22**2*b1**2 + a12**2*a21*b1*b2 + a12*a21**2*b1*b2 - a12*a21*a22*b1**2]


    return coeff
