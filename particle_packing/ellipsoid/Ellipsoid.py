import numpy as np


class Ellipsoid(object):
    """

    """

    def __init__(self, center, radii, rt_ax, phi):
        """
        """

        center = np.asarray(center, dtype=np.float).flatten()
        if len(center) != 3:
            raise ValueError('incorrect input for center input')
        else:
            self.center = center

        radii = np.asarray(radii, dtype=np.float).flatten()
        if len(radii) == 1:
            self.radii = radii * np.ones(3)
        elif len(radii) == 3:
            self.radii = radii
        else:
            raise ValueError('incorrect input for radii input')


        rt_ax = np.asarray(rt_ax, dtype=np.float).flatten()
        if len(rt_ax) != 3:
            raise ValueError('incorrect input for center rotation axis')

        phi = float(phi)
        self.phi = phi

        v1 = np.cos(phi / 2.) ** 2
        v2 = (np.sin(phi / 2.) * rt_ax[0]) ** 2
        v3 = (np.sin(phi / 2.) * rt_ax[1]) ** 2
        v4 = (np.sin(phi / 2.) * rt_ax[2]) ** 2
        if not np.allclose(v1 + v2 + v3 + v4, 1.0):
            raise ValueError('incorrect input for rotation axis and phi')

        self.rt_ax = rt_ax








    def generate_volume(self, x_ax, y_ax, z_ax):
        """Generate volume for the ellipsoid.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis
        z_ax -- numpy array for z-axis

        Return values:
        vol -- 3D array for the volume

        """

        x_ax = np.asarray(x_ax, dtype=np.float).flatten()
        y_ax = np.asarray(y_ax, dtype=np.float).flatten()
        z_ax = np.asarray(z_ax, dtype=np.float).flatten()

        vol = _generate_ellipsoid_volume(x_ax, y_ax, z_ax,
            self.center, self.radii,
            self.rt_ax[0], self.rt_ax[1], self.rt_ax[2],
            self.phi)

        return vol




    def find_subvolume(self, x_ax, y_ax, z_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        sphere centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis
        z_ax -- numpy array for z-axis

        Return values:
        x_ax_subvol -- x-axis for subvolume
        y_ax_subvol -- y-axis for subvolume
        z_ax_subvol -- z-axis for subvolume

        """

        x_ax_subvol, y_ax_subvol, z_ax_subvol, \
        x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix = \
        _find_sphere_subvolume(x_ax, y_ax, z_ax, self.center, max(self.radii))

        return x_ax_subvol, y_ax_subvol, z_ax_subvol


    def find_subvolume_ix(self, x_ax, y_ax, z_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        sphere centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis
        z_ax -- numpy array for z-axis

        Return values:
        x_ax_subvol_ix -- index values (x_ax) for x-axis for subvolume
        y_ax_subvol_ix -- index values (y_ax) y-axis or subvolume
        z_ax_subvol_ix -- index values (z_ax) z-axis for subvolume

        """

        x_ax_subvol, y_ax_subvol, z_ax_subvol, \
        x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix = \
        _find_sphere_subvolume(x_ax, y_ax, z_ax, self.center, max(self.radii))

        return x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix








    def overlap_potential(self, c):
        """Determine the overlap potential of object self and object c.

        Overlap criterion based on the overlap potential value:
        F(A,B) > 1, A and B are disjoint
        F(A,B) = 0, A and B are externally tangent
        F(A,B) < 1, A and B are overlapping

        Input arguments:
        c -- object to check for overlap with self

        Return values:
        F -- overlap potential value

        """

        if not isinstance(c, Ellipsoid):
            raise ValueError('input is not an ellipsoid')


        F = overlap_potential_py(self.center, self.radii, self.phi, self.rt_ax,
            c.center, c.radii, c.phi, c.rt_ax)

        return F




    def contain_potential(self, c):
        """Determine contain potential for the object.

        Containment criterion based on the contain potential value:
        G(A,B) > 1, A completely inside B
        G(A,B) = 1, A completely inside and tangent to B
        G(A,B) < 1, A at least partially outside B


        Return values:
        G -- contain potential value

        """


        if not isinstance(c, Ellipsoid):
            raise ValueError('input is not an ellipse')

        print 'contain_potential method is not yet implemented.'

        return float('nan')










def _generate_ellipsoid_volume(x, y, z, center, radii, rt_x, rt_y, rt_z, phi):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place an ellipsoid at center having given radii and to quaternion
    rotation given by rotation axis [rt_x, rt_y, rt_z] and rotation angle phi.

    The point x = [xi, yi, zi] is inside an ellipoid if 
    (x - center)^T Q^T A Q (x - center) <= 1.

    center = [xc, yc, zc] is the center point of the ellipsoid

    Q = [ (1-s) * P + (s^2 - 0.5) * I]
    s = cos( phi / 2. )
    p = sin (phi / 2. ) * [ rt_x, rt_y, rt_z ]
    P = p * p^T

    A = [[1/a^2, 0, 0], records major and minor axis radii
         [0, 1/b^2, 0],
         [0, 0, 1/c^2]]

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    z -- extent of the volume along z-axis
    radii -- sphere radius
    center -- center point of the sphere
    rt_x --
    rt_y --
    rt_z --
    phi --

    Return values:
    vol -- generated volume containing ellipse

    """

    # Setup matrices for quadratic form evaluation
    rt_ax = np.asarray([rt_x, rt_y, rt_z])

    B = _characteristic_matrix(radii, rt_ax, phi, exponent=1.)

    B = np.asarray(B)


    # Form square position array for x and y
    X = np.zeros((3, len(y) * len(x) * len(z)))
    X[0] = np.tile(x - center[0], (len(z), len(y), 1)).flatten()
    X[1] = np.tile(y - center[1], (len(z), len(x), 1)).transpose(0, 2, 1).flatten()
    X[2] = np.tile(z - center[2], (len(y), len(x), 1)).transpose(2, 0, 1).flatten()

    # Indicate points inside the ellipse
    vol = (X * (B.dot(X))).sum(0)
    vol = np.reshape(vol, (len(z), len(y), len(x)))
    vol = vol <= 1.
    vol = vol.astype(float)

    return vol





def _find_sphere_subvolume(X, Y, Z, xi, a):
    """Extract sub-volume from the volume with x-, y-, z-axes given by X, Y, Z.
    The sub-volume is just large enough to contain a sphere centered at xi and
    having radius a.

    Keyword arguments:
    X -- extent of the volume along x-axis
    Y -- extent of the volume along y-axis
    Z -- extent of the volume along z-axis
    xi -- sphere center point
    a -- sphere radius


    Return values:
    X_subvol -- extent of the sub-volume along x-axis
    Y_subvol -- extent of the sub-volume along y-axis
    Z_subvol -- extent of the sub-volume along z-axis
    X_subvol_ix -- X_subvol = X[X_subvol_ix.min():X_subvol_ix.max()]
    Y_subvol_ix -- Y_subvol = Y[Y_subvol_ix.min():Y_subvol_ix.max()]
    Z_subvol_ix -- Z_subvol = Z[Z_subvol_ix.min():Y_subvol_ix.max()]

    """

    # Find smallest cube that contains the sphere
    X_subvol_ix = np.nonzero(np.abs(X - xi[0]) <= a)[0]
    Y_subvol_ix = np.nonzero(np.abs(Y - xi[1]) <= a)[0]
    Z_subvol_ix = np.nonzero(np.abs(Z - xi[2]) <= a)[0]

    # Get axis arrays for the sub-volume
    X_subvol = X[X_subvol_ix]
    Y_subvol = Y[Y_subvol_ix]
    Z_subvol = Z[Z_subvol_ix]

    return X_subvol, Y_subvol, Z_subvol, X_subvol_ix, Y_subvol_ix, Z_subvol_ix






def overlap_potential_py(rA, radiiA, phiA, rotaxA, rB, radiiB, phiB, rotaxB):
    """

    Overlap potential function (Python version) provides a distance measure
    for ellipsoids A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of A
    radiiA -- radii of A
    phiA -- rotation angle of A
    rotaxA -- rotation axis for A
    rB -- center of B
    radiiB -- radii of B
    phiB -- rotation angle of B
    rotaxB --

    Return values:
    F -- overlap potential

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




    # What's needed is XA^-1 = Q^T * O^2 * Q
    XA = _characteristic_matrix(radiiA, rotaxA, phiA, exponent=-1.)


    # What's needed is XB^1/2 = Q^T * O^2 * Q
    XB = _characteristic_matrix(radiiB, rotaxB, phiB, exponent=0.5)




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
        if np.isreal(r) and (0. < r < 1.):
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


    coeff = [a11**2*a22**2*a33*b3**2 - a11**2*a22**2*b3**2 - a11**2*a22*a23*a32*b3**2 - a11**2*a22*a23*a33*b2*b3 + a11**2*a22*a23*b2*b3 - a11**2*a22*a32*a33*b2*b3 + a11**2*a22*a32*b2*b3 + a11**2*a22*a33**2*b2**2 - 2*a11**2*a22*a33*b2**2 - 2*a11**2*a22*a33*b3**2 + a11**2*a22*b2**2 + 2*a11**2*a22*b3**2 + a11**2*a23**2*a32*b2*b3 + a11**2*a23*a32**2*b2*b3 - a11**2*a23*a32*a33*b2**2 + a11**2*a23*a32*b2**2 + a11**2*a23*a32*b3**2 + a11**2*a23*a33*b2*b3 - a11**2*a23*b2*b3 + a11**2*a32*a33*b2*b3 - a11**2*a32*b2*b3 - a11**2*a33**2*b2**2 + 2*a11**2*a33*b2**2 + a11**2*a33*b3**2 - a11**2*b2**2 - a11**2*b3**2 - 2*a11*a12*a21*a22*a33*b3**2 + 2*a11*a12*a21*a22*b3**2 + a11*a12*a21*a23*a32*b3**2 + a11*a12*a21*a23*a33*b2*b3 - a11*a12*a21*a23*b2*b3 + a11*a12*a21*a32*a33*b2*b3 - a11*a12*a21*a32*b2*b3 - a11*a12*a21*a33**2*b2**2 + 2*a11*a12*a21*a33*b2**2 + 2*a11*a12*a21*a33*b3**2 - a11*a12*a21*b2**2 - 2*a11*a12*a21*b3**2 + a11*a12*a22*a23*a31*b3**2 + a11*a12*a22*a23*a33*b1*b3 - a11*a12*a22*a23*b1*b3 + a11*a12*a22*a31*a33*b2*b3 - a11*a12*a22*a31*b2*b3 - a11*a12*a22*a33**2*b1*b2 + 2*a11*a12*a22*a33*b1*b2 - a11*a12*a22*b1*b2 - a11*a12*a23**2*a31*b2*b3 - a11*a12*a23**2*a32*b1*b3 - 2*a11*a12*a23*a31*a32*b2*b3 + a11*a12*a23*a31*a33*b2**2 - a11*a12*a23*a31*b2**2 - a11*a12*a23*a31*b3**2 + a11*a12*a23*a32*a33*b1*b2 - a11*a12*a23*a32*b1*b2 - a11*a12*a23*a33*b1*b3 + a11*a12*a23*b1*b3 - a11*a12*a31*a33*b2*b3 + a11*a12*a31*b2*b3 + a11*a12*a33**2*b1*b2 - 2*a11*a12*a33*b1*b2 + a11*a12*b1*b2 + a11*a13*a21*a22*a32*b3**2 + a11*a13*a21*a22*a33*b2*b3 - a11*a13*a21*a22*b2*b3 - 2*a11*a13*a21*a23*a32*b2*b3 - a11*a13*a21*a32**2*b2*b3 + a11*a13*a21*a32*a33*b2**2 - a11*a13*a21*a32*b2**2 - a11*a13*a21*a32*b3**2 - a11*a13*a21*a33*b2*b3 + a11*a13*a21*b2*b3 - a11*a13*a22**2*a31*b3**2 - a11*a13*a22**2*a33*b1*b3 + a11*a13*a22**2*b1*b3 + a11*a13*a22*a23*a31*b2*b3 + a11*a13*a22*a23*a32*b1*b3 + a11*a13*a22*a31*a32*b2*b3 - 2*a11*a13*a22*a31*a33*b2**2 + 2*a11*a13*a22*a31*b2**2 + 2*a11*a13*a22*a31*b3**2 + a11*a13*a22*a32*a33*b1*b2 - a11*a13*a22*a32*b1*b2 + 2*a11*a13*a22*a33*b1*b3 - 2*a11*a13*a22*b1*b3 + a11*a13*a23*a31*a32*b2**2 - a11*a13*a23*a31*b2*b3 - a11*a13*a23*a32**2*b1*b2 - a11*a13*a23*a32*b1*b3 - a11*a13*a31*a32*b2*b3 + 2*a11*a13*a31*a33*b2**2 - 2*a11*a13*a31*b2**2 - a11*a13*a31*b3**2 - a11*a13*a32*a33*b1*b2 + a11*a13*a32*b1*b2 - a11*a13*a33*b1*b3 + a11*a13*b1*b3 + a11*a21*a22*a32*a33*b1*b3 - a11*a21*a22*a32*b1*b3 - a11*a21*a22*a33**2*b1*b2 + 2*a11*a21*a22*a33*b1*b2 - a11*a21*a22*b1*b2 - a11*a21*a23*a32**2*b1*b3 + a11*a21*a23*a32*a33*b1*b2 - a11*a21*a23*a32*b1*b2 - a11*a21*a32*a33*b1*b3 + a11*a21*a32*b1*b3 + a11*a21*a33**2*b1*b2 - 2*a11*a21*a33*b1*b2 + a11*a21*b1*b2 - a11*a22**2*a31*a33*b1*b3 + a11*a22**2*a31*b1*b3 + a11*a22**2*a33**2*b1**2 - 2*a11*a22**2*a33*b1**2 - 2*a11*a22**2*a33*b3**2 + a11*a22**2*b1**2 + 2*a11*a22**2*b3**2 + a11*a22*a23*a31*a32*b1*b3 + a11*a22*a23*a31*a33*b1*b2 - a11*a22*a23*a31*b1*b2 - 2*a11*a22*a23*a32*a33*b1**2 + 2*a11*a22*a23*a32*b1**2 + 2*a11*a22*a23*a32*b3**2 + 2*a11*a22*a23*a33*b2*b3 - 2*a11*a22*a23*b2*b3 + 2*a11*a22*a31*a33*b1*b3 - 2*a11*a22*a31*b1*b3 + 2*a11*a22*a32*a33*b2*b3 - 2*a11*a22*a32*b2*b3 - 2*a11*a22*a33**2*b1**2 - 2*a11*a22*a33**2*b2**2 + 4*a11*a22*a33*b1**2 + 4*a11*a22*a33*b2**2 + 4*a11*a22*a33*b3**2 - 2*a11*a22*b1**2 - 2*a11*a22*b2**2 - 4*a11*a22*b3**2 - a11*a23**2*a31*a32*b1*b2 + a11*a23**2*a32**2*b1**2 - 2*a11*a23**2*a32*b2*b3 - a11*a23*a31*a32*b1*b3 - a11*a23*a31*a33*b1*b2 + a11*a23*a31*b1*b2 - 2*a11*a23*a32**2*b2*b3 + 2*a11*a23*a32*a33*b1**2 + 2*a11*a23*a32*a33*b2**2 - 2*a11*a23*a32*b1**2 - 2*a11*a23*a32*b2**2 - 2*a11*a23*a32*b3**2 - 2*a11*a23*a33*b2*b3 + 2*a11*a23*b2*b3 - a11*a31*a33*b1*b3 + a11*a31*b1*b3 - 2*a11*a32*a33*b2*b3 + 2*a11*a32*b2*b3 + a11*a33**2*b1**2 + 2*a11*a33**2*b2**2 - 2*a11*a33*b1**2 - 4*a11*a33*b2**2 - 2*a11*a33*b3**2 + a11*b1**2 + 2*a11*b2**2 + 2*a11*b3**2 + a12**2*a21**2*a33*b3**2 - a12**2*a21**2*b3**2 - a12**2*a21*a23*a31*b3**2 - a12**2*a21*a23*a33*b1*b3 + a12**2*a21*a23*b1*b3 - a12**2*a21*a31*a33*b2*b3 + a12**2*a21*a31*b2*b3 + a12**2*a21*a33**2*b1*b2 - 2*a12**2*a21*a33*b1*b2 + a12**2*a21*b1*b2 + a12**2*a23**2*a31*b1*b3 + a12**2*a23*a31**2*b2*b3 - a12**2*a23*a31*a33*b1*b2 + a12**2*a23*a31*b1*b2 - a12*a13*a21**2*a32*b3**2 - a12*a13*a21**2*a33*b2*b3 + a12*a13*a21**2*b2*b3 + a12*a13*a21*a22*a31*b3**2 + a12*a13*a21*a22*a33*b1*b3 - a12*a13*a21*a22*b1*b3 + a12*a13*a21*a23*a31*b2*b3 + a12*a13*a21*a23*a32*b1*b3 + a12*a13*a21*a31*a32*b2*b3 + a12*a13*a21*a31*a33*b2**2 - a12*a13*a21*a31*b2**2 - a12*a13*a21*a31*b3**2 - 2*a12*a13*a21*a32*a33*b1*b2 + 2*a12*a13*a21*a32*b1*b2 - a12*a13*a21*a33*b1*b3 + a12*a13*a21*b1*b3 - 2*a12*a13*a22*a23*a31*b1*b3 - a12*a13*a22*a31**2*b2*b3 + a12*a13*a22*a31*a33*b1*b2 - a12*a13*a22*a31*b1*b2 - a12*a13*a23*a31**2*b2**2 + a12*a13*a23*a31*a32*b1*b2 + 2*a12*a13*a23*a31*b1*b3 + a12*a13*a31**2*b2*b3 - a12*a13*a31*a33*b1*b2 + a12*a13*a31*b1*b2 - a12*a21**2*a32*a33*b1*b3 + a12*a21**2*a32*b1*b3 + a12*a21**2*a33**2*b1*b2 - 2*a12*a21**2*a33*b1*b2 + a12*a21**2*b1*b2 + a12*a21*a22*a31*a33*b1*b3 - a12*a21*a22*a31*b1*b3 - a12*a21*a22*a33**2*b1**2 + 2*a12*a21*a22*a33*b1**2 + 2*a12*a21*a22*a33*b3**2 - a12*a21*a22*b1**2 - 2*a12*a21*a22*b3**2 + a12*a21*a23*a31*a32*b1*b3 - 2*a12*a21*a23*a31*a33*b1*b2 + 2*a12*a21*a23*a31*b1*b2 + a12*a21*a23*a32*a33*b1**2 - a12*a21*a23*a32*b1**2 - a12*a21*a23*a32*b3**2 - a12*a21*a23*a33*b2*b3 + a12*a21*a23*b2*b3 - a12*a21*a31*a33*b1*b3 + a12*a21*a31*b1*b3 - a12*a21*a32*a33*b2*b3 + a12*a21*a32*b2*b3 + a12*a21*a33**2*b1**2 + a12*a21*a33**2*b2**2 - 2*a12*a21*a33*b1**2 - 2*a12*a21*a33*b2**2 - 2*a12*a21*a33*b3**2 + a12*a21*b1**2 + a12*a21*b2**2 + 2*a12*a21*b3**2 - a12*a22*a23*a31**2*b1*b3 + a12*a22*a23*a31*a33*b1**2 - a12*a22*a23*a31*b1**2 - a12*a22*a23*a31*b3**2 - a12*a22*a23*a33*b1*b3 + a12*a22*a23*b1*b3 - a12*a22*a31*a33*b2*b3 + a12*a22*a31*b2*b3 + a12*a22*a33**2*b1*b2 - 2*a12*a22*a33*b1*b2 + a12*a22*b1*b2 + a12*a23**2*a31**2*b1*b2 - a12*a23**2*a31*a32*b1**2 + a12*a23**2*a31*b2*b3 + a12*a23**2*a32*b1*b3 + a12*a23*a31**2*b1*b3 + 2*a12*a23*a31*a32*b2*b3 - a12*a23*a31*a33*b1**2 - a12*a23*a31*a33*b2**2 + a12*a23*a31*b1**2 + a12*a23*a31*b2**2 + a12*a23*a31*b3**2 - a12*a23*a32*a33*b1*b2 + a12*a23*a32*b1*b2 + a12*a23*a33*b1*b3 - a12*a23*b1*b3 + a12*a31*a33*b2*b3 - a12*a31*b2*b3 - a12*a33**2*b1*b2 + 2*a12*a33*b1*b2 - a12*b1*b2 + a13**2*a21**2*a32*b2*b3 - a13**2*a21*a22*a31*b2*b3 - a13**2*a21*a22*a32*b1*b3 - a13**2*a21*a31*a32*b2**2 + a13**2*a21*a31*b2*b3 + a13**2*a21*a32**2*b1*b2 + a13**2*a21*a32*b1*b3 + a13**2*a22**2*a31*b1*b3 + a13**2*a22*a31**2*b2**2 - a13**2*a22*a31*a32*b1*b2 - 2*a13**2*a22*a31*b1*b3 - a13**2*a31**2*b2**2 + a13**2*a31*a32*b1*b2 + a13**2*a31*b1*b3 + a13*a21**2*a32**2*b1*b3 - a13*a21**2*a32*a33*b1*b2 + a13*a21**2*a32*b1*b2 - 2*a13*a21*a22*a31*a32*b1*b3 + a13*a21*a22*a31*a33*b1*b2 - a13*a21*a22*a31*b1*b2 + a13*a21*a22*a32*a33*b1**2 - a13*a21*a22*a32*b1**2 - a13*a21*a22*a32*b3**2 - a13*a21*a22*a33*b2*b3 + a13*a21*a22*b2*b3 + a13*a21*a23*a31*a32*b1*b2 - a13*a21*a23*a32**2*b1**2 + 2*a13*a21*a23*a32*b2*b3 + 2*a13*a21*a31*a32*b1*b3 - a13*a21*a31*a33*b1*b2 + a13*a21*a31*b1*b2 + a13*a21*a32**2*b2*b3 - a13*a21*a32*a33*b1**2 - a13*a21*a32*a33*b2**2 + a13*a21*a32*b1**2 + a13*a21*a32*b2**2 + a13*a21*a32*b3**2 + a13*a21*a33*b2*b3 - a13*a21*b2*b3 + a13*a22**2*a31**2*b1*b3 - a13*a22**2*a31*a33*b1**2 + a13*a22**2*a31*b1**2 + a13*a22**2*a31*b3**2 + a13*a22**2*a33*b1*b3 - a13*a22**2*b1*b3 - a13*a22*a23*a31**2*b1*b2 + a13*a22*a23*a31*a32*b1**2 - a13*a22*a23*a31*b2*b3 - a13*a22*a23*a32*b1*b3 - 2*a13*a22*a31**2*b1*b3 - a13*a22*a31*a32*b2*b3 + 2*a13*a22*a31*a33*b1**2 + 2*a13*a22*a31*a33*b2**2 - 2*a13*a22*a31*b1**2 - 2*a13*a22*a31*b2**2 - 2*a13*a22*a31*b3**2 - a13*a22*a32*a33*b1*b2 + a13*a22*a32*b1*b2 - 2*a13*a22*a33*b1*b3 + 2*a13*a22*b1*b3 + a13*a23*a31**2*b1*b2 - a13*a23*a31*a32*b1**2 - a13*a23*a31*a32*b2**2 + a13*a23*a31*b2*b3 + a13*a23*a32**2*b1*b2 + a13*a23*a32*b1*b3 + a13*a31**2*b1*b3 + a13*a31*a32*b2*b3 - a13*a31*a33*b1**2 - 2*a13*a31*a33*b2**2 + a13*a31*b1**2 + 2*a13*a31*b2**2 + a13*a31*b3**2 + a13*a32*a33*b1*b2 - a13*a32*b1*b2 + a13*a33*b1*b3 - a13*b1*b3 - a21*a22*a32*a33*b1*b3 + a21*a22*a32*b1*b3 + a21*a22*a33**2*b1*b2 - 2*a21*a22*a33*b1*b2 + a21*a22*b1*b2 + a21*a23*a32**2*b1*b3 - a21*a23*a32*a33*b1*b2 + a21*a23*a32*b1*b2 + a21*a32*a33*b1*b3 - a21*a32*b1*b3 - a21*a33**2*b1*b2 + 2*a21*a33*b1*b2 - a21*b1*b2 + a22**2*a31*a33*b1*b3 - a22**2*a31*b1*b3 - a22**2*a33**2*b1**2 + 2*a22**2*a33*b1**2 + a22**2*a33*b3**2 - a22**2*b1**2 - a22**2*b3**2 - a22*a23*a31*a32*b1*b3 - a22*a23*a31*a33*b1*b2 + a22*a23*a31*b1*b2 + 2*a22*a23*a32*a33*b1**2 - 2*a22*a23*a32*b1**2 - a22*a23*a32*b3**2 - a22*a23*a33*b2*b3 + a22*a23*b2*b3 - 2*a22*a31*a33*b1*b3 + 2*a22*a31*b1*b3 - a22*a32*a33*b2*b3 + a22*a32*b2*b3 + 2*a22*a33**2*b1**2 + a22*a33**2*b2**2 - 4*a22*a33*b1**2 - 2*a22*a33*b2**2 - 2*a22*a33*b3**2 + 2*a22*b1**2 + a22*b2**2 + 2*a22*b3**2 + a23**2*a31*a32*b1*b2 - a23**2*a32**2*b1**2 + a23**2*a32*b2*b3 + a23*a31*a32*b1*b3 + a23*a31*a33*b1*b2 - a23*a31*b1*b2 + a23*a32**2*b2*b3 - 2*a23*a32*a33*b1**2 - a23*a32*a33*b2**2 + 2*a23*a32*b1**2 + a23*a32*b2**2 + a23*a32*b3**2 + a23*a33*b2*b3 - a23*b2*b3 + a31*a33*b1*b3 - a31*b1*b3 + a32*a33*b2*b3 - a32*b2*b3 - a33**2*b1**2 - a33**2*b2**2 + 2*a33*b1**2 + 2*a33*b2**2 + a33*b3**2 - b1**2 - b2**2 - b3**2, -6*a11**2*a22**2*a33*b3**2 + 4*a11**2*a22**2*b3**2 + 6*a11**2*a22*a23*a32*b3**2 + 6*a11**2*a22*a23*a33*b2*b3 - 4*a11**2*a22*a23*b2*b3 + 6*a11**2*a22*a32*a33*b2*b3 - 4*a11**2*a22*a32*b2*b3 - 6*a11**2*a22*a33**2*b2**2 + 10*a11**2*a22*a33*b2**2 + 10*a11**2*a22*a33*b3**2 - 4*a11**2*a22*b2**2 - 6*a11**2*a22*b3**2 - 6*a11**2*a23**2*a32*b2*b3 - 6*a11**2*a23*a32**2*b2*b3 + 6*a11**2*a23*a32*a33*b2**2 - 6*a11**2*a23*a32*b2**2 - 6*a11**2*a23*a32*b3**2 - 4*a11**2*a23*a33*b2*b3 + 2*a11**2*a23*b2*b3 - 4*a11**2*a32*a33*b2*b3 + 2*a11**2*a32*b2*b3 + 4*a11**2*a33**2*b2**2 - 6*a11**2*a33*b2**2 - 4*a11**2*a33*b3**2 + 2*a11**2*b2**2 + 2*a11**2*b3**2 + 12*a11*a12*a21*a22*a33*b3**2 - 8*a11*a12*a21*a22*b3**2 - 6*a11*a12*a21*a23*a32*b3**2 - 6*a11*a12*a21*a23*a33*b2*b3 + 4*a11*a12*a21*a23*b2*b3 - 6*a11*a12*a21*a32*a33*b2*b3 + 4*a11*a12*a21*a32*b2*b3 + 6*a11*a12*a21*a33**2*b2**2 - 10*a11*a12*a21*a33*b2**2 - 10*a11*a12*a21*a33*b3**2 + 4*a11*a12*a21*b2**2 + 6*a11*a12*a21*b3**2 - 6*a11*a12*a22*a23*a31*b3**2 - 6*a11*a12*a22*a23*a33*b1*b3 + 4*a11*a12*a22*a23*b1*b3 - 6*a11*a12*a22*a31*a33*b2*b3 + 4*a11*a12*a22*a31*b2*b3 + 6*a11*a12*a22*a33**2*b1*b2 - 10*a11*a12*a22*a33*b1*b2 + 4*a11*a12*a22*b1*b2 + 6*a11*a12*a23**2*a31*b2*b3 + 6*a11*a12*a23**2*a32*b1*b3 + 12*a11*a12*a23*a31*a32*b2*b3 - 6*a11*a12*a23*a31*a33*b2**2 + 6*a11*a12*a23*a31*b2**2 + 6*a11*a12*a23*a31*b3**2 - 6*a11*a12*a23*a32*a33*b1*b2 + 6*a11*a12*a23*a32*b1*b2 + 4*a11*a12*a23*a33*b1*b3 - 2*a11*a12*a23*b1*b3 + 4*a11*a12*a31*a33*b2*b3 - 2*a11*a12*a31*b2*b3 - 4*a11*a12*a33**2*b1*b2 + 6*a11*a12*a33*b1*b2 - 2*a11*a12*b1*b2 - 6*a11*a13*a21*a22*a32*b3**2 - 6*a11*a13*a21*a22*a33*b2*b3 + 4*a11*a13*a21*a22*b2*b3 + 12*a11*a13*a21*a23*a32*b2*b3 + 6*a11*a13*a21*a32**2*b2*b3 - 6*a11*a13*a21*a32*a33*b2**2 + 6*a11*a13*a21*a32*b2**2 + 6*a11*a13*a21*a32*b3**2 + 4*a11*a13*a21*a33*b2*b3 - 2*a11*a13*a21*b2*b3 + 6*a11*a13*a22**2*a31*b3**2 + 6*a11*a13*a22**2*a33*b1*b3 - 4*a11*a13*a22**2*b1*b3 - 6*a11*a13*a22*a23*a31*b2*b3 - 6*a11*a13*a22*a23*a32*b1*b3 - 6*a11*a13*a22*a31*a32*b2*b3 + 12*a11*a13*a22*a31*a33*b2**2 - 10*a11*a13*a22*a31*b2**2 - 10*a11*a13*a22*a31*b3**2 - 6*a11*a13*a22*a32*a33*b1*b2 + 4*a11*a13*a22*a32*b1*b2 - 10*a11*a13*a22*a33*b1*b3 + 6*a11*a13*a22*b1*b3 - 6*a11*a13*a23*a31*a32*b2**2 + 4*a11*a13*a23*a31*b2*b3 + 6*a11*a13*a23*a32**2*b1*b2 + 6*a11*a13*a23*a32*b1*b3 + 4*a11*a13*a31*a32*b2*b3 - 8*a11*a13*a31*a33*b2**2 + 6*a11*a13*a31*b2**2 + 4*a11*a13*a31*b3**2 + 4*a11*a13*a32*a33*b1*b2 - 2*a11*a13*a32*b1*b2 + 4*a11*a13*a33*b1*b3 - 2*a11*a13*b1*b3 - 6*a11*a21*a22*a32*a33*b1*b3 + 4*a11*a21*a22*a32*b1*b3 + 6*a11*a21*a22*a33**2*b1*b2 - 10*a11*a21*a22*a33*b1*b2 + 4*a11*a21*a22*b1*b2 + 6*a11*a21*a23*a32**2*b1*b3 - 6*a11*a21*a23*a32*a33*b1*b2 + 6*a11*a21*a23*a32*b1*b2 + 4*a11*a21*a32*a33*b1*b3 - 2*a11*a21*a32*b1*b3 - 4*a11*a21*a33**2*b1*b2 + 6*a11*a21*a33*b1*b2 - 2*a11*a21*b1*b2 + 6*a11*a22**2*a31*a33*b1*b3 - 4*a11*a22**2*a31*b1*b3 - 6*a11*a22**2*a33**2*b1**2 + 10*a11*a22**2*a33*b1**2 + 10*a11*a22**2*a33*b3**2 - 4*a11*a22**2*b1**2 - 6*a11*a22**2*b3**2 - 6*a11*a22*a23*a31*a32*b1*b3 - 6*a11*a22*a23*a31*a33*b1*b2 + 4*a11*a22*a23*a31*b1*b2 + 12*a11*a22*a23*a32*a33*b1**2 - 10*a11*a22*a23*a32*b1**2 - 10*a11*a22*a23*a32*b3**2 - 10*a11*a22*a23*a33*b2*b3 + 6*a11*a22*a23*b2*b3 - 10*a11*a22*a31*a33*b1*b3 + 6*a11*a22*a31*b1*b3 - 10*a11*a22*a32*a33*b2*b3 + 6*a11*a22*a32*b2*b3 + 10*a11*a22*a33**2*b1**2 + 10*a11*a22*a33**2*b2**2 - 16*a11*a22*a33*b1**2 - 16*a11*a22*a33*b2**2 - 16*a11*a22*a33*b3**2 + 6*a11*a22*b1**2 + 6*a11*a22*b2**2 + 8*a11*a22*b3**2 + 6*a11*a23**2*a31*a32*b1*b2 - 6*a11*a23**2*a32**2*b1**2 + 10*a11*a23**2*a32*b2*b3 + 6*a11*a23*a31*a32*b1*b3 + 4*a11*a23*a31*a33*b1*b2 - 2*a11*a23*a31*b1*b2 + 10*a11*a23*a32**2*b2*b3 - 10*a11*a23*a32*a33*b1**2 - 10*a11*a23*a32*a33*b2**2 + 8*a11*a23*a32*b1**2 + 10*a11*a23*a32*b2**2 + 10*a11*a23*a32*b3**2 + 6*a11*a23*a33*b2*b3 - 2*a11*a23*b2*b3 + 4*a11*a31*a33*b1*b3 - 2*a11*a31*b1*b3 + 6*a11*a32*a33*b2*b3 - 2*a11*a32*b2*b3 - 4*a11*a33**2*b1**2 - 6*a11*a33**2*b2**2 + 6*a11*a33*b1**2 + 8*a11*a33*b2**2 + 6*a11*a33*b3**2 - 2*a11*b1**2 - 2*a11*b2**2 - 2*a11*b3**2 - 6*a12**2*a21**2*a33*b3**2 + 4*a12**2*a21**2*b3**2 + 6*a12**2*a21*a23*a31*b3**2 + 6*a12**2*a21*a23*a33*b1*b3 - 4*a12**2*a21*a23*b1*b3 + 6*a12**2*a21*a31*a33*b2*b3 - 4*a12**2*a21*a31*b2*b3 - 6*a12**2*a21*a33**2*b1*b2 + 10*a12**2*a21*a33*b1*b2 - 4*a12**2*a21*b1*b2 - 6*a12**2*a23**2*a31*b1*b3 - 6*a12**2*a23*a31**2*b2*b3 + 6*a12**2*a23*a31*a33*b1*b2 - 6*a12**2*a23*a31*b1*b2 + 6*a12*a13*a21**2*a32*b3**2 + 6*a12*a13*a21**2*a33*b2*b3 - 4*a12*a13*a21**2*b2*b3 - 6*a12*a13*a21*a22*a31*b3**2 - 6*a12*a13*a21*a22*a33*b1*b3 + 4*a12*a13*a21*a22*b1*b3 - 6*a12*a13*a21*a23*a31*b2*b3 - 6*a12*a13*a21*a23*a32*b1*b3 - 6*a12*a13*a21*a31*a32*b2*b3 - 6*a12*a13*a21*a31*a33*b2**2 + 4*a12*a13*a21*a31*b2**2 + 4*a12*a13*a21*a31*b3**2 + 12*a12*a13*a21*a32*a33*b1*b2 - 10*a12*a13*a21*a32*b1*b2 + 6*a12*a13*a21*a33*b1*b3 - 4*a12*a13*a21*b1*b3 + 12*a12*a13*a22*a23*a31*b1*b3 + 6*a12*a13*a22*a31**2*b2*b3 - 6*a12*a13*a22*a31*a33*b1*b2 + 6*a12*a13*a22*a31*b1*b2 + 6*a12*a13*a23*a31**2*b2**2 - 6*a12*a13*a23*a31*a32*b1*b2 - 10*a12*a13*a23*a31*b1*b3 - 4*a12*a13*a31**2*b2*b3 + 4*a12*a13*a31*a33*b1*b2 - 4*a12*a13*a31*b1*b2 + 6*a12*a21**2*a32*a33*b1*b3 - 4*a12*a21**2*a32*b1*b3 - 6*a12*a21**2*a33**2*b1*b2 + 10*a12*a21**2*a33*b1*b2 - 4*a12*a21**2*b1*b2 - 6*a12*a21*a22*a31*a33*b1*b3 + 4*a12*a21*a22*a31*b1*b3 + 6*a12*a21*a22*a33**2*b1**2 - 10*a12*a21*a22*a33*b1**2 - 10*a12*a21*a22*a33*b3**2 + 4*a12*a21*a22*b1**2 + 6*a12*a21*a22*b3**2 - 6*a12*a21*a23*a31*a32*b1*b3 + 12*a12*a21*a23*a31*a33*b1*b2 - 10*a12*a21*a23*a31*b1*b2 - 6*a12*a21*a23*a32*a33*b1**2 + 4*a12*a21*a23*a32*b1**2 + 4*a12*a21*a23*a32*b3**2 + 6*a12*a21*a23*a33*b2*b3 - 4*a12*a21*a23*b2*b3 + 6*a12*a21*a31*a33*b1*b3 - 4*a12*a21*a31*b1*b3 + 6*a12*a21*a32*a33*b2*b3 - 4*a12*a21*a32*b2*b3 - 6*a12*a21*a33**2*b1**2 - 6*a12*a21*a33**2*b2**2 + 10*a12*a21*a33*b1**2 + 10*a12*a21*a33*b2**2 + 8*a12*a21*a33*b3**2 - 4*a12*a21*b1**2 - 4*a12*a21*b2**2 - 4*a12*a21*b3**2 + 6*a12*a22*a23*a31**2*b1*b3 - 6*a12*a22*a23*a31*a33*b1**2 + 6*a12*a22*a23*a31*b1**2 + 6*a12*a22*a23*a31*b3**2 + 4*a12*a22*a23*a33*b1*b3 - 2*a12*a22*a23*b1*b3 + 4*a12*a22*a31*a33*b2*b3 - 2*a12*a22*a31*b2*b3 - 4*a12*a22*a33**2*b1*b2 + 6*a12*a22*a33*b1*b2 - 2*a12*a22*b1*b2 - 6*a12*a23**2*a31**2*b1*b2 + 6*a12*a23**2*a31*a32*b1**2 - 6*a12*a23**2*a31*b2*b3 - 4*a12*a23**2*a32*b1*b3 - 6*a12*a23*a31**2*b1*b3 - 10*a12*a23*a31*a32*b2*b3 + 6*a12*a23*a31*a33*b1**2 + 6*a12*a23*a31*a33*b2**2 - 6*a12*a23*a31*b1**2 - 6*a12*a23*a31*b2**2 - 6*a12*a23*a31*b3**2 + 4*a12*a23*a32*a33*b1*b2 - 4*a12*a23*a32*b1*b2 - 2*a12*a23*a33*b1*b3 - 2*a12*a31*a33*b2*b3 + 2*a12*a33**2*b1*b2 - 2*a12*a33*b1*b2 - 6*a13**2*a21**2*a32*b2*b3 + 6*a13**2*a21*a22*a31*b2*b3 + 6*a13**2*a21*a22*a32*b1*b3 + 6*a13**2*a21*a31*a32*b2**2 - 4*a13**2*a21*a31*b2*b3 - 6*a13**2*a21*a32**2*b1*b2 - 6*a13**2*a21*a32*b1*b3 - 6*a13**2*a22**2*a31*b1*b3 - 6*a13**2*a22*a31**2*b2**2 + 6*a13**2*a22*a31*a32*b1*b2 + 10*a13**2*a22*a31*b1*b3 + 4*a13**2*a31**2*b2**2 - 4*a13**2*a31*a32*b1*b2 - 4*a13**2*a31*b1*b3 - 6*a13*a21**2*a32**2*b1*b3 + 6*a13*a21**2*a32*a33*b1*b2 - 6*a13*a21**2*a32*b1*b2 + 12*a13*a21*a22*a31*a32*b1*b3 - 6*a13*a21*a22*a31*a33*b1*b2 + 6*a13*a21*a22*a31*b1*b2 - 6*a13*a21*a22*a32*a33*b1**2 + 6*a13*a21*a22*a32*b1**2 + 6*a13*a21*a22*a32*b3**2 + 4*a13*a21*a22*a33*b2*b3 - 2*a13*a21*a22*b2*b3 - 6*a13*a21*a23*a31*a32*b1*b2 + 6*a13*a21*a23*a32**2*b1**2 - 10*a13*a21*a23*a32*b2*b3 - 10*a13*a21*a31*a32*b1*b3 + 4*a13*a21*a31*a33*b1*b2 - 4*a13*a21*a31*b1*b2 - 6*a13*a21*a32**2*b2*b3 + 6*a13*a21*a32*a33*b1**2 + 6*a13*a21*a32*a33*b2**2 - 6*a13*a21*a32*b1**2 - 6*a13*a21*a32*b2**2 - 6*a13*a21*a32*b3**2 - 2*a13*a21*a33*b2*b3 - 6*a13*a22**2*a31**2*b1*b3 + 6*a13*a22**2*a31*a33*b1**2 - 6*a13*a22**2*a31*b1**2 - 6*a13*a22**2*a31*b3**2 - 4*a13*a22**2*a33*b1*b3 + 2*a13*a22**2*b1*b3 + 6*a13*a22*a23*a31**2*b1*b2 - 6*a13*a22*a23*a31*a32*b1**2 + 6*a13*a22*a23*a31*b2*b3 + 4*a13*a22*a23*a32*b1*b3 + 10*a13*a22*a31**2*b1*b3 + 6*a13*a22*a31*a32*b2*b3 - 10*a13*a22*a31*a33*b1**2 - 10*a13*a22*a31*a33*b2**2 + 10*a13*a22*a31*b1**2 + 8*a13*a22*a31*b2**2 + 10*a13*a22*a31*b3**2 + 4*a13*a22*a32*a33*b1*b2 - 2*a13*a22*a32*b1*b2 + 6*a13*a22*a33*b1*b3 - 2*a13*a22*b1*b3 - 4*a13*a23*a31**2*b1*b2 + 4*a13*a23*a31*a32*b1**2 + 4*a13*a23*a31*a32*b2**2 - 4*a13*a23*a31*b2*b3 - 4*a13*a23*a32**2*b1*b2 - 4*a13*a23*a32*b1*b3 - 4*a13*a31**2*b1*b3 - 4*a13*a31*a32*b2*b3 + 4*a13*a31*a33*b1**2 + 6*a13*a31*a33*b2**2 - 4*a13*a31*b1**2 - 4*a13*a31*b2**2 - 4*a13*a31*b3**2 - 2*a13*a32*a33*b1*b2 - 2*a13*a33*b1*b3 + 4*a21*a22*a32*a33*b1*b3 - 2*a21*a22*a32*b1*b3 - 4*a21*a22*a33**2*b1*b2 + 6*a21*a22*a33*b1*b2 - 2*a21*a22*b1*b2 - 4*a21*a23*a32**2*b1*b3 + 4*a21*a23*a32*a33*b1*b2 - 4*a21*a23*a32*b1*b2 - 2*a21*a32*a33*b1*b3 + 2*a21*a33**2*b1*b2 - 2*a21*a33*b1*b2 - 4*a22**2*a31*a33*b1*b3 + 2*a22**2*a31*b1*b3 + 4*a22**2*a33**2*b1**2 - 6*a22**2*a33*b1**2 - 4*a22**2*a33*b3**2 + 2*a22**2*b1**2 + 2*a22**2*b3**2 + 4*a22*a23*a31*a32*b1*b3 + 4*a22*a23*a31*a33*b1*b2 - 2*a22*a23*a31*b1*b2 - 8*a22*a23*a32*a33*b1**2 + 6*a22*a23*a32*b1**2 + 4*a22*a23*a32*b3**2 + 4*a22*a23*a33*b2*b3 - 2*a22*a23*b2*b3 + 6*a22*a31*a33*b1*b3 - 2*a22*a31*b1*b3 + 4*a22*a32*a33*b2*b3 - 2*a22*a32*b2*b3 - 6*a22*a33**2*b1**2 - 4*a22*a33**2*b2**2 + 8*a22*a33*b1**2 + 6*a22*a33*b2**2 + 6*a22*a33*b3**2 - 2*a22*b1**2 - 2*a22*b2**2 - 2*a22*b3**2 - 4*a23**2*a31*a32*b1*b2 + 4*a23**2*a32**2*b1**2 - 4*a23**2*a32*b2*b3 - 4*a23*a31*a32*b1*b3 - 2*a23*a31*a33*b1*b2 - 4*a23*a32**2*b2*b3 + 6*a23*a32*a33*b1**2 + 4*a23*a32*a33*b2**2 - 4*a23*a32*b1**2 - 4*a23*a32*b2**2 - 4*a23*a32*b3**2 - 2*a23*a33*b2*b3 - 2*a31*a33*b1*b3 - 2*a32*a33*b2*b3 + 2*a33**2*b1**2 + 2*a33**2*b2**2 - 2*a33*b1**2 - 2*a33*b2**2 - 2*a33*b3**2, 15*a11**2*a22**2*a33*b3**2 - 6*a11**2*a22**2*b3**2 - 15*a11**2*a22*a23*a32*b3**2 - 15*a11**2*a22*a23*a33*b2*b3 + 6*a11**2*a22*a23*b2*b3 - 15*a11**2*a22*a32*a33*b2*b3 + 6*a11**2*a22*a32*b2*b3 + 15*a11**2*a22*a33**2*b2**2 - 20*a11**2*a22*a33*b2**2 - 20*a11**2*a22*a33*b3**2 + 6*a11**2*a22*b2**2 + 6*a11**2*a22*b3**2 + 15*a11**2*a23**2*a32*b2*b3 + 15*a11**2*a23*a32**2*b2*b3 - 15*a11**2*a23*a32*a33*b2**2 + 14*a11**2*a23*a32*b2**2 + 14*a11**2*a23*a32*b3**2 + 6*a11**2*a23*a33*b2*b3 + 6*a11**2*a32*a33*b2*b3 - 6*a11**2*a33**2*b2**2 + 6*a11**2*a33*b2**2 + 6*a11**2*a33*b3**2 - a11**2*b2**2 - a11**2*b3**2 - 30*a11*a12*a21*a22*a33*b3**2 + 12*a11*a12*a21*a22*b3**2 + 15*a11*a12*a21*a23*a32*b3**2 + 15*a11*a12*a21*a23*a33*b2*b3 - 6*a11*a12*a21*a23*b2*b3 + 15*a11*a12*a21*a32*a33*b2*b3 - 6*a11*a12*a21*a32*b2*b3 - 15*a11*a12*a21*a33**2*b2**2 + 20*a11*a12*a21*a33*b2**2 + 20*a11*a12*a21*a33*b3**2 - 6*a11*a12*a21*b2**2 - 6*a11*a12*a21*b3**2 + 15*a11*a12*a22*a23*a31*b3**2 + 15*a11*a12*a22*a23*a33*b1*b3 - 6*a11*a12*a22*a23*b1*b3 + 15*a11*a12*a22*a31*a33*b2*b3 - 6*a11*a12*a22*a31*b2*b3 - 15*a11*a12*a22*a33**2*b1*b2 + 20*a11*a12*a22*a33*b1*b2 - 6*a11*a12*a22*b1*b2 - 15*a11*a12*a23**2*a31*b2*b3 - 15*a11*a12*a23**2*a32*b1*b3 - 30*a11*a12*a23*a31*a32*b2*b3 + 15*a11*a12*a23*a31*a33*b2**2 - 14*a11*a12*a23*a31*b2**2 - 14*a11*a12*a23*a31*b3**2 + 15*a11*a12*a23*a32*a33*b1*b2 - 14*a11*a12*a23*a32*b1*b2 - 6*a11*a12*a23*a33*b1*b3 - 6*a11*a12*a31*a33*b2*b3 + 6*a11*a12*a33**2*b1*b2 - 6*a11*a12*a33*b1*b2 + a11*a12*b1*b2 + 15*a11*a13*a21*a22*a32*b3**2 + 15*a11*a13*a21*a22*a33*b2*b3 - 6*a11*a13*a21*a22*b2*b3 - 30*a11*a13*a21*a23*a32*b2*b3 - 15*a11*a13*a21*a32**2*b2*b3 + 15*a11*a13*a21*a32*a33*b2**2 - 14*a11*a13*a21*a32*b2**2 - 14*a11*a13*a21*a32*b3**2 - 6*a11*a13*a21*a33*b2*b3 - 15*a11*a13*a22**2*a31*b3**2 - 15*a11*a13*a22**2*a33*b1*b3 + 6*a11*a13*a22**2*b1*b3 + 15*a11*a13*a22*a23*a31*b2*b3 + 15*a11*a13*a22*a23*a32*b1*b3 + 15*a11*a13*a22*a31*a32*b2*b3 - 30*a11*a13*a22*a31*a33*b2**2 + 20*a11*a13*a22*a31*b2**2 + 20*a11*a13*a22*a31*b3**2 + 15*a11*a13*a22*a32*a33*b1*b2 - 6*a11*a13*a22*a32*b1*b2 + 20*a11*a13*a22*a33*b1*b3 - 6*a11*a13*a22*b1*b3 + 15*a11*a13*a23*a31*a32*b2**2 - 6*a11*a13*a23*a31*b2*b3 - 15*a11*a13*a23*a32**2*b1*b2 - 14*a11*a13*a23*a32*b1*b3 - 6*a11*a13*a31*a32*b2*b3 + 12*a11*a13*a31*a33*b2**2 - 6*a11*a13*a31*b2**2 - 6*a11*a13*a31*b3**2 - 6*a11*a13*a32*a33*b1*b2 - 6*a11*a13*a33*b1*b3 + a11*a13*b1*b3 + 15*a11*a21*a22*a32*a33*b1*b3 - 6*a11*a21*a22*a32*b1*b3 - 15*a11*a21*a22*a33**2*b1*b2 + 20*a11*a21*a22*a33*b1*b2 - 6*a11*a21*a22*b1*b2 - 15*a11*a21*a23*a32**2*b1*b3 + 15*a11*a21*a23*a32*a33*b1*b2 - 14*a11*a21*a23*a32*b1*b2 - 6*a11*a21*a32*a33*b1*b3 + 6*a11*a21*a33**2*b1*b2 - 6*a11*a21*a33*b1*b2 + a11*a21*b1*b2 - 15*a11*a22**2*a31*a33*b1*b3 + 6*a11*a22**2*a31*b1*b3 + 15*a11*a22**2*a33**2*b1**2 - 20*a11*a22**2*a33*b1**2 - 20*a11*a22**2*a33*b3**2 + 6*a11*a22**2*b1**2 + 6*a11*a22**2*b3**2 + 15*a11*a22*a23*a31*a32*b1*b3 + 15*a11*a22*a23*a31*a33*b1*b2 - 6*a11*a22*a23*a31*b1*b2 - 30*a11*a22*a23*a32*a33*b1**2 + 20*a11*a22*a23*a32*b1**2 + 20*a11*a22*a23*a32*b3**2 + 20*a11*a22*a23*a33*b2*b3 - 6*a11*a22*a23*b2*b3 + 20*a11*a22*a31*a33*b1*b3 - 6*a11*a22*a31*b1*b3 + 20*a11*a22*a32*a33*b2*b3 - 6*a11*a22*a32*b2*b3 - 20*a11*a22*a33**2*b1**2 - 20*a11*a22*a33**2*b2**2 + 24*a11*a22*a33*b1**2 + 24*a11*a22*a33*b2**2 + 24*a11*a22*a33*b3**2 - 6*a11*a22*b1**2 - 6*a11*a22*b2**2 - 4*a11*a22*b3**2 - 15*a11*a23**2*a31*a32*b1*b2 + 15*a11*a23**2*a32**2*b1**2 - 20*a11*a23**2*a32*b2*b3 - 14*a11*a23*a31*a32*b1*b3 - 6*a11*a23*a31*a33*b1*b2 - 20*a11*a23*a32**2*b2*b3 + 20*a11*a23*a32*a33*b1**2 + 20*a11*a23*a32*a33*b2**2 - 12*a11*a23*a32*b1**2 - 18*a11*a23*a32*b2**2 - 18*a11*a23*a32*b3**2 - 6*a11*a23*a33*b2*b3 - 2*a11*a23*b2*b3 - 6*a11*a31*a33*b1*b3 + a11*a31*b1*b3 - 6*a11*a32*a33*b2*b3 - 2*a11*a32*b2*b3 + 6*a11*a33**2*b1**2 + 6*a11*a33**2*b2**2 - 6*a11*a33*b1**2 - 4*a11*a33*b2**2 - 6*a11*a33*b3**2 + a11*b1**2 + 15*a12**2*a21**2*a33*b3**2 - 6*a12**2*a21**2*b3**2 - 15*a12**2*a21*a23*a31*b3**2 - 15*a12**2*a21*a23*a33*b1*b3 + 6*a12**2*a21*a23*b1*b3 - 15*a12**2*a21*a31*a33*b2*b3 + 6*a12**2*a21*a31*b2*b3 + 15*a12**2*a21*a33**2*b1*b2 - 20*a12**2*a21*a33*b1*b2 + 6*a12**2*a21*b1*b2 + 15*a12**2*a23**2*a31*b1*b3 + 15*a12**2*a23*a31**2*b2*b3 - 15*a12**2*a23*a31*a33*b1*b2 + 14*a12**2*a23*a31*b1*b2 - 15*a12*a13*a21**2*a32*b3**2 - 15*a12*a13*a21**2*a33*b2*b3 + 6*a12*a13*a21**2*b2*b3 + 15*a12*a13*a21*a22*a31*b3**2 + 15*a12*a13*a21*a22*a33*b1*b3 - 6*a12*a13*a21*a22*b1*b3 + 15*a12*a13*a21*a23*a31*b2*b3 + 15*a12*a13*a21*a23*a32*b1*b3 + 15*a12*a13*a21*a31*a32*b2*b3 + 15*a12*a13*a21*a31*a33*b2**2 - 6*a12*a13*a21*a31*b2**2 - 6*a12*a13*a21*a31*b3**2 - 30*a12*a13*a21*a32*a33*b1*b2 + 20*a12*a13*a21*a32*b1*b2 - 14*a12*a13*a21*a33*b1*b3 + 6*a12*a13*a21*b1*b3 - 30*a12*a13*a22*a23*a31*b1*b3 - 15*a12*a13*a22*a31**2*b2*b3 + 15*a12*a13*a22*a31*a33*b1*b2 - 14*a12*a13*a22*a31*b1*b2 - 15*a12*a13*a23*a31**2*b2**2 + 15*a12*a13*a23*a31*a32*b1*b2 + 20*a12*a13*a23*a31*b1*b3 + 6*a12*a13*a31**2*b2*b3 - 6*a12*a13*a31*a33*b1*b2 + 6*a12*a13*a31*b1*b2 - 15*a12*a21**2*a32*a33*b1*b3 + 6*a12*a21**2*a32*b1*b3 + 15*a12*a21**2*a33**2*b1*b2 - 20*a12*a21**2*a33*b1*b2 + 6*a12*a21**2*b1*b2 + 15*a12*a21*a22*a31*a33*b1*b3 - 6*a12*a21*a22*a31*b1*b3 - 15*a12*a21*a22*a33**2*b1**2 + 20*a12*a21*a22*a33*b1**2 + 20*a12*a21*a22*a33*b3**2 - 6*a12*a21*a22*b1**2 - 6*a12*a21*a22*b3**2 + 15*a12*a21*a23*a31*a32*b1*b3 - 30*a12*a21*a23*a31*a33*b1*b2 + 20*a12*a21*a23*a31*b1*b2 + 15*a12*a21*a23*a32*a33*b1**2 - 6*a12*a21*a23*a32*b1**2 - 6*a12*a21*a23*a32*b3**2 - 14*a12*a21*a23*a33*b2*b3 + 6*a12*a21*a23*b2*b3 - 14*a12*a21*a31*a33*b1*b3 + 6*a12*a21*a31*b1*b3 - 14*a12*a21*a32*a33*b2*b3 + 6*a12*a21*a32*b2*b3 + 14*a12*a21*a33**2*b1**2 + 14*a12*a21*a33**2*b2**2 - 18*a12*a21*a33*b1**2 - 18*a12*a21*a33*b2**2 - 12*a12*a21*a33*b3**2 + 5*a12*a21*b1**2 + 5*a12*a21*b2**2 + 2*a12*a21*b3**2 - 15*a12*a22*a23*a31**2*b1*b3 + 15*a12*a22*a23*a31*a33*b1**2 - 14*a12*a22*a23*a31*b1**2 - 14*a12*a22*a23*a31*b3**2 - 6*a12*a22*a23*a33*b1*b3 - 6*a12*a22*a31*a33*b2*b3 + 6*a12*a22*a33**2*b1*b2 - 6*a12*a22*a33*b1*b2 + a12*a22*b1*b2 + 15*a12*a23**2*a31**2*b1*b2 - 15*a12*a23**2*a31*a32*b1**2 + 14*a12*a23**2*a31*b2*b3 + 6*a12*a23**2*a32*b1*b3 + 14*a12*a23*a31**2*b1*b3 + 20*a12*a23*a31*a32*b2*b3 - 14*a12*a23*a31*a33*b1**2 - 14*a12*a23*a31*a33*b2**2 + 12*a12*a23*a31*b1**2 + 12*a12*a23*a31*b2**2 + 12*a12*a23*a31*b3**2 - 6*a12*a23*a32*a33*b1*b2 + 6*a12*a23*a32*b1*b2 + 3*a12*a23*b1*b3 + 3*a12*a31*b2*b3 - 2*a12*a33*b1*b2 + a12*b1*b2 + 15*a13**2*a21**2*a32*b2*b3 - 15*a13**2*a21*a22*a31*b2*b3 - 15*a13**2*a21*a22*a32*b1*b3 - 15*a13**2*a21*a31*a32*b2**2 + 6*a13**2*a21*a31*b2*b3 + 15*a13**2*a21*a32**2*b1*b2 + 14*a13**2*a21*a32*b1*b3 + 15*a13**2*a22**2*a31*b1*b3 + 15*a13**2*a22*a31**2*b2**2 - 15*a13**2*a22*a31*a32*b1*b2 - 20*a13**2*a22*a31*b1*b3 - 6*a13**2*a31**2*b2**2 + 6*a13**2*a31*a32*b1*b2 + 6*a13**2*a31*b1*b3 + 15*a13*a21**2*a32**2*b1*b3 - 15*a13*a21**2*a32*a33*b1*b2 + 14*a13*a21**2*a32*b1*b2 - 30*a13*a21*a22*a31*a32*b1*b3 + 15*a13*a21*a22*a31*a33*b1*b2 - 14*a13*a21*a22*a31*b1*b2 + 15*a13*a21*a22*a32*a33*b1**2 - 14*a13*a21*a22*a32*b1**2 - 14*a13*a21*a22*a32*b3**2 - 6*a13*a21*a22*a33*b2*b3 + 15*a13*a21*a23*a31*a32*b1*b2 - 15*a13*a21*a23*a32**2*b1**2 + 20*a13*a21*a23*a32*b2*b3 + 20*a13*a21*a31*a32*b1*b3 - 6*a13*a21*a31*a33*b1*b2 + 6*a13*a21*a31*b1*b2 + 14*a13*a21*a32**2*b2*b3 - 14*a13*a21*a32*a33*b1**2 - 14*a13*a21*a32*a33*b2**2 + 12*a13*a21*a32*b1**2 + 12*a13*a21*a32*b2**2 + 12*a13*a21*a32*b3**2 + 3*a13*a21*b2*b3 + 15*a13*a22**2*a31**2*b1*b3 - 15*a13*a22**2*a31*a33*b1**2 + 14*a13*a22**2*a31*b1**2 + 14*a13*a22**2*a31*b3**2 + 6*a13*a22**2*a33*b1*b3 - 15*a13*a22*a23*a31**2*b1*b2 + 15*a13*a22*a23*a31*a32*b1**2 - 14*a13*a22*a23*a31*b2*b3 - 6*a13*a22*a23*a32*b1*b3 - 20*a13*a22*a31**2*b1*b3 - 14*a13*a22*a31*a32*b2*b3 + 20*a13*a22*a31*a33*b1**2 + 20*a13*a22*a31*a33*b2**2 - 18*a13*a22*a31*b1**2 - 12*a13*a22*a31*b2**2 - 18*a13*a22*a31*b3**2 - 6*a13*a22*a32*a33*b1*b2 - 6*a13*a22*a33*b1*b3 - 2*a13*a22*b1*b3 + 6*a13*a23*a31**2*b1*b2 - 6*a13*a23*a31*a32*b1**2 - 6*a13*a23*a31*a32*b2**2 + 6*a13*a23*a31*b2*b3 + 6*a13*a23*a32**2*b1*b2 + 6*a13*a23*a32*b1*b3 + 6*a13*a31**2*b1*b3 + 6*a13*a31*a32*b2*b3 - 6*a13*a31*a33*b1**2 - 6*a13*a31*a33*b2**2 + 5*a13*a31*b1**2 + 2*a13*a31*b2**2 + 5*a13*a31*b3**2 + 3*a13*a32*b1*b2 + a13*a33*b1*b3 + a13*b1*b3 - 6*a21*a22*a32*a33*b1*b3 + 6*a21*a22*a33**2*b1*b2 - 6*a21*a22*a33*b1*b2 + a21*a22*b1*b2 + 6*a21*a23*a32**2*b1*b3 - 6*a21*a23*a32*a33*b1*b2 + 6*a21*a23*a32*b1*b2 + 3*a21*a32*b1*b3 - 2*a21*a33*b1*b2 + a21*b1*b2 + 6*a22**2*a31*a33*b1*b3 - 6*a22**2*a33**2*b1**2 + 6*a22**2*a33*b1**2 + 6*a22**2*a33*b3**2 - a22**2*b1**2 - a22**2*b3**2 - 6*a22*a23*a31*a32*b1*b3 - 6*a22*a23*a31*a33*b1*b2 + 12*a22*a23*a32*a33*b1**2 - 6*a22*a23*a32*b1**2 - 6*a22*a23*a32*b3**2 - 6*a22*a23*a33*b2*b3 + a22*a23*b2*b3 - 6*a22*a31*a33*b1*b3 - 2*a22*a31*b1*b3 - 6*a22*a32*a33*b2*b3 + a22*a32*b2*b3 + 6*a22*a33**2*b1**2 + 6*a22*a33**2*b2**2 - 4*a22*a33*b1**2 - 6*a22*a33*b2**2 - 6*a22*a33*b3**2 + a22*b2**2 + 6*a23**2*a31*a32*b1*b2 - 6*a23**2*a32**2*b1**2 + 6*a23**2*a32*b2*b3 + 6*a23*a31*a32*b1*b3 + 3*a23*a31*b1*b2 + 6*a23*a32**2*b2*b3 - 6*a23*a32*a33*b1**2 - 6*a23*a32*a33*b2**2 + 2*a23*a32*b1**2 + 5*a23*a32*b2**2 + 5*a23*a32*b3**2 + a23*a33*b2*b3 + a23*b2*b3 + a31*a33*b1*b3 + a31*b1*b3 + a32*a33*b2*b3 + a32*b2*b3 - a33**2*b1**2 - a33**2*b2**2 + a33*b3**2, -20*a11**2*a22**2*a33*b3**2 + 4*a11**2*a22**2*b3**2 + 20*a11**2*a22*a23*a32*b3**2 + 20*a11**2*a22*a23*a33*b2*b3 - 4*a11**2*a22*a23*b2*b3 + 20*a11**2*a22*a32*a33*b2*b3 - 4*a11**2*a22*a32*b2*b3 - 20*a11**2*a22*a33**2*b2**2 + 20*a11**2*a22*a33*b2**2 + 20*a11**2*a22*a33*b3**2 - 4*a11**2*a22*b2**2 - 2*a11**2*a22*b3**2 - 20*a11**2*a23**2*a32*b2*b3 - 20*a11**2*a23*a32**2*b2*b3 + 20*a11**2*a23*a32*a33*b2**2 - 16*a11**2*a23*a32*b2**2 - 16*a11**2*a23*a32*b3**2 - 4*a11**2*a23*a33*b2*b3 - 2*a11**2*a23*b2*b3 - 4*a11**2*a32*a33*b2*b3 - 2*a11**2*a32*b2*b3 + 4*a11**2*a33**2*b2**2 - 2*a11**2*a33*b2**2 - 4*a11**2*a33*b3**2 + 40*a11*a12*a21*a22*a33*b3**2 - 8*a11*a12*a21*a22*b3**2 - 20*a11*a12*a21*a23*a32*b3**2 - 20*a11*a12*a21*a23*a33*b2*b3 + 4*a11*a12*a21*a23*b2*b3 - 20*a11*a12*a21*a32*a33*b2*b3 + 4*a11*a12*a21*a32*b2*b3 + 20*a11*a12*a21*a33**2*b2**2 - 20*a11*a12*a21*a33*b2**2 - 20*a11*a12*a21*a33*b3**2 + 4*a11*a12*a21*b2**2 + 2*a11*a12*a21*b3**2 - 20*a11*a12*a22*a23*a31*b3**2 - 20*a11*a12*a22*a23*a33*b1*b3 + 4*a11*a12*a22*a23*b1*b3 - 20*a11*a12*a22*a31*a33*b2*b3 + 4*a11*a12*a22*a31*b2*b3 + 20*a11*a12*a22*a33**2*b1*b2 - 20*a11*a12*a22*a33*b1*b2 + 4*a11*a12*a22*b1*b2 + 20*a11*a12*a23**2*a31*b2*b3 + 20*a11*a12*a23**2*a32*b1*b3 + 40*a11*a12*a23*a31*a32*b2*b3 - 20*a11*a12*a23*a31*a33*b2**2 + 16*a11*a12*a23*a31*b2**2 + 16*a11*a12*a23*a31*b3**2 - 20*a11*a12*a23*a32*a33*b1*b2 + 16*a11*a12*a23*a32*b1*b2 + 4*a11*a12*a23*a33*b1*b3 + 2*a11*a12*a23*b1*b3 + 4*a11*a12*a31*a33*b2*b3 + 2*a11*a12*a31*b2*b3 - 4*a11*a12*a33**2*b1*b2 + 2*a11*a12*a33*b1*b2 - 20*a11*a13*a21*a22*a32*b3**2 - 20*a11*a13*a21*a22*a33*b2*b3 + 4*a11*a13*a21*a22*b2*b3 + 40*a11*a13*a21*a23*a32*b2*b3 + 20*a11*a13*a21*a32**2*b2*b3 - 20*a11*a13*a21*a32*a33*b2**2 + 16*a11*a13*a21*a32*b2**2 + 16*a11*a13*a21*a32*b3**2 + 4*a11*a13*a21*a33*b2*b3 + 2*a11*a13*a21*b2*b3 + 20*a11*a13*a22**2*a31*b3**2 + 20*a11*a13*a22**2*a33*b1*b3 - 4*a11*a13*a22**2*b1*b3 - 20*a11*a13*a22*a23*a31*b2*b3 - 20*a11*a13*a22*a23*a32*b1*b3 - 20*a11*a13*a22*a31*a32*b2*b3 + 40*a11*a13*a22*a31*a33*b2**2 - 20*a11*a13*a22*a31*b2**2 - 20*a11*a13*a22*a31*b3**2 - 20*a11*a13*a22*a32*a33*b1*b2 + 4*a11*a13*a22*a32*b1*b2 - 20*a11*a13*a22*a33*b1*b3 + 2*a11*a13*a22*b1*b3 - 20*a11*a13*a23*a31*a32*b2**2 + 4*a11*a13*a23*a31*b2*b3 + 20*a11*a13*a23*a32**2*b1*b2 + 16*a11*a13*a23*a32*b1*b3 + 4*a11*a13*a31*a32*b2*b3 - 8*a11*a13*a31*a33*b2**2 + 2*a11*a13*a31*b2**2 + 4*a11*a13*a31*b3**2 + 4*a11*a13*a32*a33*b1*b2 + 2*a11*a13*a32*b1*b2 + 4*a11*a13*a33*b1*b3 - 20*a11*a21*a22*a32*a33*b1*b3 + 4*a11*a21*a22*a32*b1*b3 + 20*a11*a21*a22*a33**2*b1*b2 - 20*a11*a21*a22*a33*b1*b2 + 4*a11*a21*a22*b1*b2 + 20*a11*a21*a23*a32**2*b1*b3 - 20*a11*a21*a23*a32*a33*b1*b2 + 16*a11*a21*a23*a32*b1*b2 + 4*a11*a21*a32*a33*b1*b3 + 2*a11*a21*a32*b1*b3 - 4*a11*a21*a33**2*b1*b2 + 2*a11*a21*a33*b1*b2 + 20*a11*a22**2*a31*a33*b1*b3 - 4*a11*a22**2*a31*b1*b3 - 20*a11*a22**2*a33**2*b1**2 + 20*a11*a22**2*a33*b1**2 + 20*a11*a22**2*a33*b3**2 - 4*a11*a22**2*b1**2 - 2*a11*a22**2*b3**2 - 20*a11*a22*a23*a31*a32*b1*b3 - 20*a11*a22*a23*a31*a33*b1*b2 + 4*a11*a22*a23*a31*b1*b2 + 40*a11*a22*a23*a32*a33*b1**2 - 20*a11*a22*a23*a32*b1**2 - 20*a11*a22*a23*a32*b3**2 - 20*a11*a22*a23*a33*b2*b3 + 2*a11*a22*a23*b2*b3 - 20*a11*a22*a31*a33*b1*b3 + 2*a11*a22*a31*b1*b3 - 20*a11*a22*a32*a33*b2*b3 + 2*a11*a22*a32*b2*b3 + 20*a11*a22*a33**2*b1**2 + 20*a11*a22*a33**2*b2**2 - 16*a11*a22*a33*b1**2 - 16*a11*a22*a33*b2**2 - 16*a11*a22*a33*b3**2 + 2*a11*a22*b1**2 + 2*a11*a22*b2**2 + 20*a11*a23**2*a31*a32*b1*b2 - 20*a11*a23**2*a32**2*b1**2 + 20*a11*a23**2*a32*b2*b3 + 16*a11*a23*a31*a32*b1*b3 + 4*a11*a23*a31*a33*b1*b2 + 2*a11*a23*a31*b1*b2 + 20*a11*a23*a32**2*b2*b3 - 20*a11*a23*a32*a33*b1**2 - 20*a11*a23*a32*a33*b2**2 + 8*a11*a23*a32*b1**2 + 14*a11*a23*a32*b2**2 + 14*a11*a23*a32*b3**2 + 2*a11*a23*a33*b2*b3 + 2*a11*a23*b2*b3 + 4*a11*a31*a33*b1*b3 + 2*a11*a32*a33*b2*b3 + 2*a11*a32*b2*b3 - 4*a11*a33**2*b1**2 - 2*a11*a33**2*b2**2 + 2*a11*a33*b1**2 + 2*a11*a33*b3**2 - 20*a12**2*a21**2*a33*b3**2 + 4*a12**2*a21**2*b3**2 + 20*a12**2*a21*a23*a31*b3**2 + 20*a12**2*a21*a23*a33*b1*b3 - 4*a12**2*a21*a23*b1*b3 + 20*a12**2*a21*a31*a33*b2*b3 - 4*a12**2*a21*a31*b2*b3 - 20*a12**2*a21*a33**2*b1*b2 + 20*a12**2*a21*a33*b1*b2 - 4*a12**2*a21*b1*b2 - 20*a12**2*a23**2*a31*b1*b3 - 20*a12**2*a23*a31**2*b2*b3 + 20*a12**2*a23*a31*a33*b1*b2 - 16*a12**2*a23*a31*b1*b2 + 20*a12*a13*a21**2*a32*b3**2 + 20*a12*a13*a21**2*a33*b2*b3 - 4*a12*a13*a21**2*b2*b3 - 20*a12*a13*a21*a22*a31*b3**2 - 20*a12*a13*a21*a22*a33*b1*b3 + 4*a12*a13*a21*a22*b1*b3 - 20*a12*a13*a21*a23*a31*b2*b3 - 20*a12*a13*a21*a23*a32*b1*b3 - 20*a12*a13*a21*a31*a32*b2*b3 - 20*a12*a13*a21*a31*a33*b2**2 + 4*a12*a13*a21*a31*b2**2 + 4*a12*a13*a21*a31*b3**2 + 40*a12*a13*a21*a32*a33*b1*b2 - 20*a12*a13*a21*a32*b1*b2 + 16*a12*a13*a21*a33*b1*b3 - 4*a12*a13*a21*b1*b3 + 40*a12*a13*a22*a23*a31*b1*b3 + 20*a12*a13*a22*a31**2*b2*b3 - 20*a12*a13*a22*a31*a33*b1*b2 + 16*a12*a13*a22*a31*b1*b2 + 20*a12*a13*a23*a31**2*b2**2 - 20*a12*a13*a23*a31*a32*b1*b2 - 20*a12*a13*a23*a31*b1*b3 - 4*a12*a13*a31**2*b2*b3 + 4*a12*a13*a31*a33*b1*b2 - 4*a12*a13*a31*b1*b2 + 20*a12*a21**2*a32*a33*b1*b3 - 4*a12*a21**2*a32*b1*b3 - 20*a12*a21**2*a33**2*b1*b2 + 20*a12*a21**2*a33*b1*b2 - 4*a12*a21**2*b1*b2 - 20*a12*a21*a22*a31*a33*b1*b3 + 4*a12*a21*a22*a31*b1*b3 + 20*a12*a21*a22*a33**2*b1**2 - 20*a12*a21*a22*a33*b1**2 - 20*a12*a21*a22*a33*b3**2 + 4*a12*a21*a22*b1**2 + 2*a12*a21*a22*b3**2 - 20*a12*a21*a23*a31*a32*b1*b3 + 40*a12*a21*a23*a31*a33*b1*b2 - 20*a12*a21*a23*a31*b1*b2 - 20*a12*a21*a23*a32*a33*b1**2 + 4*a12*a21*a23*a32*b1**2 + 4*a12*a21*a23*a32*b3**2 + 16*a12*a21*a23*a33*b2*b3 - 4*a12*a21*a23*b2*b3 + 16*a12*a21*a31*a33*b1*b3 - 4*a12*a21*a31*b1*b3 + 16*a12*a21*a32*a33*b2*b3 - 4*a12*a21*a32*b2*b3 - 16*a12*a21*a33**2*b1**2 - 16*a12*a21*a33**2*b2**2 + 14*a12*a21*a33*b1**2 + 14*a12*a21*a33*b2**2 + 8*a12*a21*a33*b3**2 - 2*a12*a21*b1**2 - 2*a12*a21*b2**2 + 20*a12*a22*a23*a31**2*b1*b3 - 20*a12*a22*a23*a31*a33*b1**2 + 16*a12*a22*a23*a31*b1**2 + 16*a12*a22*a23*a31*b3**2 + 4*a12*a22*a23*a33*b1*b3 + 2*a12*a22*a23*b1*b3 + 4*a12*a22*a31*a33*b2*b3 + 2*a12*a22*a31*b2*b3 - 4*a12*a22*a33**2*b1*b2 + 2*a12*a22*a33*b1*b2 - 20*a12*a23**2*a31**2*b1*b2 + 20*a12*a23**2*a31*a32*b1**2 - 16*a12*a23**2*a31*b2*b3 - 4*a12*a23**2*a32*b1*b3 - 16*a12*a23*a31**2*b1*b3 - 20*a12*a23*a31*a32*b2*b3 + 16*a12*a23*a31*a33*b1**2 + 16*a12*a23*a31*a33*b2**2 - 10*a12*a23*a31*b1**2 - 10*a12*a23*a31*b2**2 - 10*a12*a23*a31*b3**2 + 4*a12*a23*a32*a33*b1*b2 - 4*a12*a23*a32*b1*b2 + 2*a12*a23*a33*b1*b3 - 2*a12*a23*b1*b3 + 2*a12*a31*a33*b2*b3 - 2*a12*a31*b2*b3 - 2*a12*a33**2*b1*b2 + 2*a12*a33*b1*b2 - 20*a13**2*a21**2*a32*b2*b3 + 20*a13**2*a21*a22*a31*b2*b3 + 20*a13**2*a21*a22*a32*b1*b3 + 20*a13**2*a21*a31*a32*b2**2 - 4*a13**2*a21*a31*b2*b3 - 20*a13**2*a21*a32**2*b1*b2 - 16*a13**2*a21*a32*b1*b3 - 20*a13**2*a22**2*a31*b1*b3 - 20*a13**2*a22*a31**2*b2**2 + 20*a13**2*a22*a31*a32*b1*b2 + 20*a13**2*a22*a31*b1*b3 + 4*a13**2*a31**2*b2**2 - 4*a13**2*a31*a32*b1*b2 - 4*a13**2*a31*b1*b3 - 20*a13*a21**2*a32**2*b1*b3 + 20*a13*a21**2*a32*a33*b1*b2 - 16*a13*a21**2*a32*b1*b2 + 40*a13*a21*a22*a31*a32*b1*b3 - 20*a13*a21*a22*a31*a33*b1*b2 + 16*a13*a21*a22*a31*b1*b2 - 20*a13*a21*a22*a32*a33*b1**2 + 16*a13*a21*a22*a32*b1**2 + 16*a13*a21*a22*a32*b3**2 + 4*a13*a21*a22*a33*b2*b3 + 2*a13*a21*a22*b2*b3 - 20*a13*a21*a23*a31*a32*b1*b2 + 20*a13*a21*a23*a32**2*b1**2 - 20*a13*a21*a23*a32*b2*b3 - 20*a13*a21*a31*a32*b1*b3 + 4*a13*a21*a31*a33*b1*b2 - 4*a13*a21*a31*b1*b2 - 16*a13*a21*a32**2*b2*b3 + 16*a13*a21*a32*a33*b1**2 + 16*a13*a21*a32*a33*b2**2 - 10*a13*a21*a32*b1**2 - 10*a13*a21*a32*b2**2 - 10*a13*a21*a32*b3**2 + 2*a13*a21*a33*b2*b3 - 2*a13*a21*b2*b3 - 20*a13*a22**2*a31**2*b1*b3 + 20*a13*a22**2*a31*a33*b1**2 - 16*a13*a22**2*a31*b1**2 - 16*a13*a22**2*a31*b3**2 - 4*a13*a22**2*a33*b1*b3 - 2*a13*a22**2*b1*b3 + 20*a13*a22*a23*a31**2*b1*b2 - 20*a13*a22*a23*a31*a32*b1**2 + 16*a13*a22*a23*a31*b2*b3 + 4*a13*a22*a23*a32*b1*b3 + 20*a13*a22*a31**2*b1*b3 + 16*a13*a22*a31*a32*b2*b3 - 20*a13*a22*a31*a33*b1**2 - 20*a13*a22*a31*a33*b2**2 + 14*a13*a22*a31*b1**2 + 8*a13*a22*a31*b2**2 + 14*a13*a22*a31*b3**2 + 4*a13*a22*a32*a33*b1*b2 + 2*a13*a22*a32*b1*b2 + 2*a13*a22*a33*b1*b3 + 2*a13*a22*b1*b3 - 4*a13*a23*a31**2*b1*b2 + 4*a13*a23*a31*a32*b1**2 + 4*a13*a23*a31*a32*b2**2 - 4*a13*a23*a31*b2*b3 - 4*a13*a23*a32**2*b1*b2 - 4*a13*a23*a32*b1*b3 - 4*a13*a31**2*b1*b3 - 4*a13*a31*a32*b2*b3 + 4*a13*a31*a33*b1**2 + 2*a13*a31*a33*b2**2 - 2*a13*a31*b1**2 - 2*a13*a31*b3**2 + 2*a13*a32*a33*b1*b2 - 2*a13*a32*b1*b2 + 4*a21*a22*a32*a33*b1*b3 + 2*a21*a22*a32*b1*b3 - 4*a21*a22*a33**2*b1*b2 + 2*a21*a22*a33*b1*b2 - 4*a21*a23*a32**2*b1*b3 + 4*a21*a23*a32*a33*b1*b2 - 4*a21*a23*a32*b1*b2 + 2*a21*a32*a33*b1*b3 - 2*a21*a32*b1*b3 - 2*a21*a33**2*b1*b2 + 2*a21*a33*b1*b2 - 4*a22**2*a31*a33*b1*b3 - 2*a22**2*a31*b1*b3 + 4*a22**2*a33**2*b1**2 - 2*a22**2*a33*b1**2 - 4*a22**2*a33*b3**2 + 4*a22*a23*a31*a32*b1*b3 + 4*a22*a23*a31*a33*b1*b2 + 2*a22*a23*a31*b1*b2 - 8*a22*a23*a32*a33*b1**2 + 2*a22*a23*a32*b1**2 + 4*a22*a23*a32*b3**2 + 4*a22*a23*a33*b2*b3 + 2*a22*a31*a33*b1*b3 + 2*a22*a31*b1*b3 + 4*a22*a32*a33*b2*b3 - 2*a22*a33**2*b1**2 - 4*a22*a33**2*b2**2 + 2*a22*a33*b2**2 + 2*a22*a33*b3**2 - 4*a23**2*a31*a32*b1*b2 + 4*a23**2*a32**2*b1**2 - 4*a23**2*a32*b2*b3 - 4*a23*a31*a32*b1*b3 + 2*a23*a31*a33*b1*b2 - 2*a23*a31*b1*b2 - 4*a23*a32**2*b2*b3 + 2*a23*a32*a33*b1**2 + 4*a23*a32*a33*b2**2 - 2*a23*a32*b2**2 - 2*a23*a32*b3**2, 15*a11**2*a22**2*a33*b3**2 - a11**2*a22**2*b3**2 - 15*a11**2*a22*a23*a32*b3**2 - 15*a11**2*a22*a23*a33*b2*b3 + a11**2*a22*a23*b2*b3 - 15*a11**2*a22*a32*a33*b2*b3 + a11**2*a22*a32*b2*b3 + 15*a11**2*a22*a33**2*b2**2 - 10*a11**2*a22*a33*b2**2 - 10*a11**2*a22*a33*b3**2 + a11**2*a22*b2**2 + 15*a11**2*a23**2*a32*b2*b3 + 15*a11**2*a23*a32**2*b2*b3 - 15*a11**2*a23*a32*a33*b2**2 + 9*a11**2*a23*a32*b2**2 + 9*a11**2*a23*a32*b3**2 + a11**2*a23*a33*b2*b3 + a11**2*a23*b2*b3 + a11**2*a32*a33*b2*b3 + a11**2*a32*b2*b3 - a11**2*a33**2*b2**2 + a11**2*a33*b3**2 - 30*a11*a12*a21*a22*a33*b3**2 + 2*a11*a12*a21*a22*b3**2 + 15*a11*a12*a21*a23*a32*b3**2 + 15*a11*a12*a21*a23*a33*b2*b3 - a11*a12*a21*a23*b2*b3 + 15*a11*a12*a21*a32*a33*b2*b3 - a11*a12*a21*a32*b2*b3 - 15*a11*a12*a21*a33**2*b2**2 + 10*a11*a12*a21*a33*b2**2 + 10*a11*a12*a21*a33*b3**2 - a11*a12*a21*b2**2 + 15*a11*a12*a22*a23*a31*b3**2 + 15*a11*a12*a22*a23*a33*b1*b3 - a11*a12*a22*a23*b1*b3 + 15*a11*a12*a22*a31*a33*b2*b3 - a11*a12*a22*a31*b2*b3 - 15*a11*a12*a22*a33**2*b1*b2 + 10*a11*a12*a22*a33*b1*b2 - a11*a12*a22*b1*b2 - 15*a11*a12*a23**2*a31*b2*b3 - 15*a11*a12*a23**2*a32*b1*b3 - 30*a11*a12*a23*a31*a32*b2*b3 + 15*a11*a12*a23*a31*a33*b2**2 - 9*a11*a12*a23*a31*b2**2 - 9*a11*a12*a23*a31*b3**2 + 15*a11*a12*a23*a32*a33*b1*b2 - 9*a11*a12*a23*a32*b1*b2 - a11*a12*a23*a33*b1*b3 - a11*a12*a23*b1*b3 - a11*a12*a31*a33*b2*b3 - a11*a12*a31*b2*b3 + a11*a12*a33**2*b1*b2 + 15*a11*a13*a21*a22*a32*b3**2 + 15*a11*a13*a21*a22*a33*b2*b3 - a11*a13*a21*a22*b2*b3 - 30*a11*a13*a21*a23*a32*b2*b3 - 15*a11*a13*a21*a32**2*b2*b3 + 15*a11*a13*a21*a32*a33*b2**2 - 9*a11*a13*a21*a32*b2**2 - 9*a11*a13*a21*a32*b3**2 - a11*a13*a21*a33*b2*b3 - a11*a13*a21*b2*b3 - 15*a11*a13*a22**2*a31*b3**2 - 15*a11*a13*a22**2*a33*b1*b3 + a11*a13*a22**2*b1*b3 + 15*a11*a13*a22*a23*a31*b2*b3 + 15*a11*a13*a22*a23*a32*b1*b3 + 15*a11*a13*a22*a31*a32*b2*b3 - 30*a11*a13*a22*a31*a33*b2**2 + 10*a11*a13*a22*a31*b2**2 + 10*a11*a13*a22*a31*b3**2 + 15*a11*a13*a22*a32*a33*b1*b2 - a11*a13*a22*a32*b1*b2 + 10*a11*a13*a22*a33*b1*b3 + 15*a11*a13*a23*a31*a32*b2**2 - a11*a13*a23*a31*b2*b3 - 15*a11*a13*a23*a32**2*b1*b2 - 9*a11*a13*a23*a32*b1*b3 - a11*a13*a31*a32*b2*b3 + 2*a11*a13*a31*a33*b2**2 - a11*a13*a31*b3**2 - a11*a13*a32*a33*b1*b2 - a11*a13*a32*b1*b2 - a11*a13*a33*b1*b3 + 15*a11*a21*a22*a32*a33*b1*b3 - a11*a21*a22*a32*b1*b3 - 15*a11*a21*a22*a33**2*b1*b2 + 10*a11*a21*a22*a33*b1*b2 - a11*a21*a22*b1*b2 - 15*a11*a21*a23*a32**2*b1*b3 + 15*a11*a21*a23*a32*a33*b1*b2 - 9*a11*a21*a23*a32*b1*b2 - a11*a21*a32*a33*b1*b3 - a11*a21*a32*b1*b3 + a11*a21*a33**2*b1*b2 - 15*a11*a22**2*a31*a33*b1*b3 + a11*a22**2*a31*b1*b3 + 15*a11*a22**2*a33**2*b1**2 - 10*a11*a22**2*a33*b1**2 - 10*a11*a22**2*a33*b3**2 + a11*a22**2*b1**2 + 15*a11*a22*a23*a31*a32*b1*b3 + 15*a11*a22*a23*a31*a33*b1*b2 - a11*a22*a23*a31*b1*b2 - 30*a11*a22*a23*a32*a33*b1**2 + 10*a11*a22*a23*a32*b1**2 + 10*a11*a22*a23*a32*b3**2 + 10*a11*a22*a23*a33*b2*b3 + 10*a11*a22*a31*a33*b1*b3 + 10*a11*a22*a32*a33*b2*b3 - 10*a11*a22*a33**2*b1**2 - 10*a11*a22*a33**2*b2**2 + 4*a11*a22*a33*b1**2 + 4*a11*a22*a33*b2**2 + 4*a11*a22*a33*b3**2 - 15*a11*a23**2*a31*a32*b1*b2 + 15*a11*a23**2*a32**2*b1**2 - 10*a11*a23**2*a32*b2*b3 - 9*a11*a23*a31*a32*b1*b3 - a11*a23*a31*a33*b1*b2 - a11*a23*a31*b1*b2 - 10*a11*a23*a32**2*b2*b3 + 10*a11*a23*a32*a33*b1**2 + 10*a11*a23*a32*a33*b2**2 - 2*a11*a23*a32*b1**2 - 4*a11*a23*a32*b2**2 - 4*a11*a23*a32*b3**2 - a11*a31*a33*b1*b3 + a11*a33**2*b1**2 + 15*a12**2*a21**2*a33*b3**2 - a12**2*a21**2*b3**2 - 15*a12**2*a21*a23*a31*b3**2 - 15*a12**2*a21*a23*a33*b1*b3 + a12**2*a21*a23*b1*b3 - 15*a12**2*a21*a31*a33*b2*b3 + a12**2*a21*a31*b2*b3 + 15*a12**2*a21*a33**2*b1*b2 - 10*a12**2*a21*a33*b1*b2 + a12**2*a21*b1*b2 + 15*a12**2*a23**2*a31*b1*b3 + 15*a12**2*a23*a31**2*b2*b3 - 15*a12**2*a23*a31*a33*b1*b2 + 9*a12**2*a23*a31*b1*b2 - 15*a12*a13*a21**2*a32*b3**2 - 15*a12*a13*a21**2*a33*b2*b3 + a12*a13*a21**2*b2*b3 + 15*a12*a13*a21*a22*a31*b3**2 + 15*a12*a13*a21*a22*a33*b1*b3 - a12*a13*a21*a22*b1*b3 + 15*a12*a13*a21*a23*a31*b2*b3 + 15*a12*a13*a21*a23*a32*b1*b3 + 15*a12*a13*a21*a31*a32*b2*b3 + 15*a12*a13*a21*a31*a33*b2**2 - a12*a13*a21*a31*b2**2 - a12*a13*a21*a31*b3**2 - 30*a12*a13*a21*a32*a33*b1*b2 + 10*a12*a13*a21*a32*b1*b2 - 9*a12*a13*a21*a33*b1*b3 + a12*a13*a21*b1*b3 - 30*a12*a13*a22*a23*a31*b1*b3 - 15*a12*a13*a22*a31**2*b2*b3 + 15*a12*a13*a22*a31*a33*b1*b2 - 9*a12*a13*a22*a31*b1*b2 - 15*a12*a13*a23*a31**2*b2**2 + 15*a12*a13*a23*a31*a32*b1*b2 + 10*a12*a13*a23*a31*b1*b3 + a12*a13*a31**2*b2*b3 - a12*a13*a31*a33*b1*b2 + a12*a13*a31*b1*b2 - 15*a12*a21**2*a32*a33*b1*b3 + a12*a21**2*a32*b1*b3 + 15*a12*a21**2*a33**2*b1*b2 - 10*a12*a21**2*a33*b1*b2 + a12*a21**2*b1*b2 + 15*a12*a21*a22*a31*a33*b1*b3 - a12*a21*a22*a31*b1*b3 - 15*a12*a21*a22*a33**2*b1**2 + 10*a12*a21*a22*a33*b1**2 + 10*a12*a21*a22*a33*b3**2 - a12*a21*a22*b1**2 + 15*a12*a21*a23*a31*a32*b1*b3 - 30*a12*a21*a23*a31*a33*b1*b2 + 10*a12*a21*a23*a31*b1*b2 + 15*a12*a21*a23*a32*a33*b1**2 - a12*a21*a23*a32*b1**2 - a12*a21*a23*a32*b3**2 - 9*a12*a21*a23*a33*b2*b3 + a12*a21*a23*b2*b3 - 9*a12*a21*a31*a33*b1*b3 + a12*a21*a31*b1*b3 - 9*a12*a21*a32*a33*b2*b3 + a12*a21*a32*b2*b3 + 9*a12*a21*a33**2*b1**2 + 9*a12*a21*a33**2*b2**2 - 4*a12*a21*a33*b1**2 - 4*a12*a21*a33*b2**2 - 2*a12*a21*a33*b3**2 - 15*a12*a22*a23*a31**2*b1*b3 + 15*a12*a22*a23*a31*a33*b1**2 - 9*a12*a22*a23*a31*b1**2 - 9*a12*a22*a23*a31*b3**2 - a12*a22*a23*a33*b1*b3 - a12*a22*a23*b1*b3 - a12*a22*a31*a33*b2*b3 - a12*a22*a31*b2*b3 + a12*a22*a33**2*b1*b2 + 15*a12*a23**2*a31**2*b1*b2 - 15*a12*a23**2*a31*a32*b1**2 + 9*a12*a23**2*a31*b2*b3 + a12*a23**2*a32*b1*b3 + 9*a12*a23*a31**2*b1*b3 + 10*a12*a23*a31*a32*b2*b3 - 9*a12*a23*a31*a33*b1**2 - 9*a12*a23*a31*a33*b2**2 + 3*a12*a23*a31*b1**2 + 3*a12*a23*a31*b2**2 + 3*a12*a23*a31*b3**2 - a12*a23*a32*a33*b1*b2 + a12*a23*a32*b1*b2 - a12*a23*a33*b1*b3 - a12*a31*a33*b2*b3 + a12*a33**2*b1*b2 + 15*a13**2*a21**2*a32*b2*b3 - 15*a13**2*a21*a22*a31*b2*b3 - 15*a13**2*a21*a22*a32*b1*b3 - 15*a13**2*a21*a31*a32*b2**2 + a13**2*a21*a31*b2*b3 + 15*a13**2*a21*a32**2*b1*b2 + 9*a13**2*a21*a32*b1*b3 + 15*a13**2*a22**2*a31*b1*b3 + 15*a13**2*a22*a31**2*b2**2 - 15*a13**2*a22*a31*a32*b1*b2 - 10*a13**2*a22*a31*b1*b3 - a13**2*a31**2*b2**2 + a13**2*a31*a32*b1*b2 + a13**2*a31*b1*b3 + 15*a13*a21**2*a32**2*b1*b3 - 15*a13*a21**2*a32*a33*b1*b2 + 9*a13*a21**2*a32*b1*b2 - 30*a13*a21*a22*a31*a32*b1*b3 + 15*a13*a21*a22*a31*a33*b1*b2 - 9*a13*a21*a22*a31*b1*b2 + 15*a13*a21*a22*a32*a33*b1**2 - 9*a13*a21*a22*a32*b1**2 - 9*a13*a21*a22*a32*b3**2 - a13*a21*a22*a33*b2*b3 - a13*a21*a22*b2*b3 + 15*a13*a21*a23*a31*a32*b1*b2 - 15*a13*a21*a23*a32**2*b1**2 + 10*a13*a21*a23*a32*b2*b3 + 10*a13*a21*a31*a32*b1*b3 - a13*a21*a31*a33*b1*b2 + a13*a21*a31*b1*b2 + 9*a13*a21*a32**2*b2*b3 - 9*a13*a21*a32*a33*b1**2 - 9*a13*a21*a32*a33*b2**2 + 3*a13*a21*a32*b1**2 + 3*a13*a21*a32*b2**2 + 3*a13*a21*a32*b3**2 - a13*a21*a33*b2*b3 + 15*a13*a22**2*a31**2*b1*b3 - 15*a13*a22**2*a31*a33*b1**2 + 9*a13*a22**2*a31*b1**2 + 9*a13*a22**2*a31*b3**2 + a13*a22**2*a33*b1*b3 + a13*a22**2*b1*b3 - 15*a13*a22*a23*a31**2*b1*b2 + 15*a13*a22*a23*a31*a32*b1**2 - 9*a13*a22*a23*a31*b2*b3 - a13*a22*a23*a32*b1*b3 - 10*a13*a22*a31**2*b1*b3 - 9*a13*a22*a31*a32*b2*b3 + 10*a13*a22*a31*a33*b1**2 + 10*a13*a22*a31*a33*b2**2 - 4*a13*a22*a31*b1**2 - 2*a13*a22*a31*b2**2 - 4*a13*a22*a31*b3**2 - a13*a22*a32*a33*b1*b2 - a13*a22*a32*b1*b2 + a13*a23*a31**2*b1*b2 - a13*a23*a31*a32*b1**2 - a13*a23*a31*a32*b2**2 + a13*a23*a31*b2*b3 + a13*a23*a32**2*b1*b2 + a13*a23*a32*b1*b3 + a13*a31**2*b1*b3 + a13*a31*a32*b2*b3 - a13*a31*a33*b1**2 - a13*a32*a33*b1*b2 - a21*a22*a32*a33*b1*b3 - a21*a22*a32*b1*b3 + a21*a22*a33**2*b1*b2 + a21*a23*a32**2*b1*b3 - a21*a23*a32*a33*b1*b2 + a21*a23*a32*b1*b2 - a21*a32*a33*b1*b3 + a21*a33**2*b1*b2 + a22**2*a31*a33*b1*b3 + a22**2*a31*b1*b3 - a22**2*a33**2*b1**2 + a22**2*a33*b3**2 - a22*a23*a31*a32*b1*b3 - a22*a23*a31*a33*b1*b2 - a22*a23*a31*b1*b2 + 2*a22*a23*a32*a33*b1**2 - a22*a23*a32*b3**2 - a22*a23*a33*b2*b3 - a22*a32*a33*b2*b3 + a22*a33**2*b2**2 + a23**2*a31*a32*b1*b2 - a23**2*a32**2*b1**2 + a23**2*a32*b2*b3 + a23*a31*a32*b1*b3 - a23*a31*a33*b1*b2 + a23*a32**2*b2*b3 - a23*a32*a33*b2**2, -6*a11**2*a22**2*a33*b3**2 + 6*a11**2*a22*a23*a32*b3**2 + 6*a11**2*a22*a23*a33*b2*b3 + 6*a11**2*a22*a32*a33*b2*b3 - 6*a11**2*a22*a33**2*b2**2 + 2*a11**2*a22*a33*b2**2 + 2*a11**2*a22*a33*b3**2 - 6*a11**2*a23**2*a32*b2*b3 - 6*a11**2*a23*a32**2*b2*b3 + 6*a11**2*a23*a32*a33*b2**2 - 2*a11**2*a23*a32*b2**2 - 2*a11**2*a23*a32*b3**2 + 12*a11*a12*a21*a22*a33*b3**2 - 6*a11*a12*a21*a23*a32*b3**2 - 6*a11*a12*a21*a23*a33*b2*b3 - 6*a11*a12*a21*a32*a33*b2*b3 + 6*a11*a12*a21*a33**2*b2**2 - 2*a11*a12*a21*a33*b2**2 - 2*a11*a12*a21*a33*b3**2 - 6*a11*a12*a22*a23*a31*b3**2 - 6*a11*a12*a22*a23*a33*b1*b3 - 6*a11*a12*a22*a31*a33*b2*b3 + 6*a11*a12*a22*a33**2*b1*b2 - 2*a11*a12*a22*a33*b1*b2 + 6*a11*a12*a23**2*a31*b2*b3 + 6*a11*a12*a23**2*a32*b1*b3 + 12*a11*a12*a23*a31*a32*b2*b3 - 6*a11*a12*a23*a31*a33*b2**2 + 2*a11*a12*a23*a31*b2**2 + 2*a11*a12*a23*a31*b3**2 - 6*a11*a12*a23*a32*a33*b1*b2 + 2*a11*a12*a23*a32*b1*b2 - 6*a11*a13*a21*a22*a32*b3**2 - 6*a11*a13*a21*a22*a33*b2*b3 + 12*a11*a13*a21*a23*a32*b2*b3 + 6*a11*a13*a21*a32**2*b2*b3 - 6*a11*a13*a21*a32*a33*b2**2 + 2*a11*a13*a21*a32*b2**2 + 2*a11*a13*a21*a32*b3**2 + 6*a11*a13*a22**2*a31*b3**2 + 6*a11*a13*a22**2*a33*b1*b3 - 6*a11*a13*a22*a23*a31*b2*b3 - 6*a11*a13*a22*a23*a32*b1*b3 - 6*a11*a13*a22*a31*a32*b2*b3 + 12*a11*a13*a22*a31*a33*b2**2 - 2*a11*a13*a22*a31*b2**2 - 2*a11*a13*a22*a31*b3**2 - 6*a11*a13*a22*a32*a33*b1*b2 - 2*a11*a13*a22*a33*b1*b3 - 6*a11*a13*a23*a31*a32*b2**2 + 6*a11*a13*a23*a32**2*b1*b2 + 2*a11*a13*a23*a32*b1*b3 - 6*a11*a21*a22*a32*a33*b1*b3 + 6*a11*a21*a22*a33**2*b1*b2 - 2*a11*a21*a22*a33*b1*b2 + 6*a11*a21*a23*a32**2*b1*b3 - 6*a11*a21*a23*a32*a33*b1*b2 + 2*a11*a21*a23*a32*b1*b2 + 6*a11*a22**2*a31*a33*b1*b3 - 6*a11*a22**2*a33**2*b1**2 + 2*a11*a22**2*a33*b1**2 + 2*a11*a22**2*a33*b3**2 - 6*a11*a22*a23*a31*a32*b1*b3 - 6*a11*a22*a23*a31*a33*b1*b2 + 12*a11*a22*a23*a32*a33*b1**2 - 2*a11*a22*a23*a32*b1**2 - 2*a11*a22*a23*a32*b3**2 - 2*a11*a22*a23*a33*b2*b3 - 2*a11*a22*a31*a33*b1*b3 - 2*a11*a22*a32*a33*b2*b3 + 2*a11*a22*a33**2*b1**2 + 2*a11*a22*a33**2*b2**2 + 6*a11*a23**2*a31*a32*b1*b2 - 6*a11*a23**2*a32**2*b1**2 + 2*a11*a23**2*a32*b2*b3 + 2*a11*a23*a31*a32*b1*b3 + 2*a11*a23*a32**2*b2*b3 - 2*a11*a23*a32*a33*b1**2 - 2*a11*a23*a32*a33*b2**2 - 6*a12**2*a21**2*a33*b3**2 + 6*a12**2*a21*a23*a31*b3**2 + 6*a12**2*a21*a23*a33*b1*b3 + 6*a12**2*a21*a31*a33*b2*b3 - 6*a12**2*a21*a33**2*b1*b2 + 2*a12**2*a21*a33*b1*b2 - 6*a12**2*a23**2*a31*b1*b3 - 6*a12**2*a23*a31**2*b2*b3 + 6*a12**2*a23*a31*a33*b1*b2 - 2*a12**2*a23*a31*b1*b2 + 6*a12*a13*a21**2*a32*b3**2 + 6*a12*a13*a21**2*a33*b2*b3 - 6*a12*a13*a21*a22*a31*b3**2 - 6*a12*a13*a21*a22*a33*b1*b3 - 6*a12*a13*a21*a23*a31*b2*b3 - 6*a12*a13*a21*a23*a32*b1*b3 - 6*a12*a13*a21*a31*a32*b2*b3 - 6*a12*a13*a21*a31*a33*b2**2 + 12*a12*a13*a21*a32*a33*b1*b2 - 2*a12*a13*a21*a32*b1*b2 + 2*a12*a13*a21*a33*b1*b3 + 12*a12*a13*a22*a23*a31*b1*b3 + 6*a12*a13*a22*a31**2*b2*b3 - 6*a12*a13*a22*a31*a33*b1*b2 + 2*a12*a13*a22*a31*b1*b2 + 6*a12*a13*a23*a31**2*b2**2 - 6*a12*a13*a23*a31*a32*b1*b2 - 2*a12*a13*a23*a31*b1*b3 + 6*a12*a21**2*a32*a33*b1*b3 - 6*a12*a21**2*a33**2*b1*b2 + 2*a12*a21**2*a33*b1*b2 - 6*a12*a21*a22*a31*a33*b1*b3 + 6*a12*a21*a22*a33**2*b1**2 - 2*a12*a21*a22*a33*b1**2 - 2*a12*a21*a22*a33*b3**2 - 6*a12*a21*a23*a31*a32*b1*b3 + 12*a12*a21*a23*a31*a33*b1*b2 - 2*a12*a21*a23*a31*b1*b2 - 6*a12*a21*a23*a32*a33*b1**2 + 2*a12*a21*a23*a33*b2*b3 + 2*a12*a21*a31*a33*b1*b3 + 2*a12*a21*a32*a33*b2*b3 - 2*a12*a21*a33**2*b1**2 - 2*a12*a21*a33**2*b2**2 + 6*a12*a22*a23*a31**2*b1*b3 - 6*a12*a22*a23*a31*a33*b1**2 + 2*a12*a22*a23*a31*b1**2 + 2*a12*a22*a23*a31*b3**2 - 6*a12*a23**2*a31**2*b1*b2 + 6*a12*a23**2*a31*a32*b1**2 - 2*a12*a23**2*a31*b2*b3 - 2*a12*a23*a31**2*b1*b3 - 2*a12*a23*a31*a32*b2*b3 + 2*a12*a23*a31*a33*b1**2 + 2*a12*a23*a31*a33*b2**2 - 6*a13**2*a21**2*a32*b2*b3 + 6*a13**2*a21*a22*a31*b2*b3 + 6*a13**2*a21*a22*a32*b1*b3 + 6*a13**2*a21*a31*a32*b2**2 - 6*a13**2*a21*a32**2*b1*b2 - 2*a13**2*a21*a32*b1*b3 - 6*a13**2*a22**2*a31*b1*b3 - 6*a13**2*a22*a31**2*b2**2 + 6*a13**2*a22*a31*a32*b1*b2 + 2*a13**2*a22*a31*b1*b3 - 6*a13*a21**2*a32**2*b1*b3 + 6*a13*a21**2*a32*a33*b1*b2 - 2*a13*a21**2*a32*b1*b2 + 12*a13*a21*a22*a31*a32*b1*b3 - 6*a13*a21*a22*a31*a33*b1*b2 + 2*a13*a21*a22*a31*b1*b2 - 6*a13*a21*a22*a32*a33*b1**2 + 2*a13*a21*a22*a32*b1**2 + 2*a13*a21*a22*a32*b3**2 - 6*a13*a21*a23*a31*a32*b1*b2 + 6*a13*a21*a23*a32**2*b1**2 - 2*a13*a21*a23*a32*b2*b3 - 2*a13*a21*a31*a32*b1*b3 - 2*a13*a21*a32**2*b2*b3 + 2*a13*a21*a32*a33*b1**2 + 2*a13*a21*a32*a33*b2**2 - 6*a13*a22**2*a31**2*b1*b3 + 6*a13*a22**2*a31*a33*b1**2 - 2*a13*a22**2*a31*b1**2 - 2*a13*a22**2*a31*b3**2 + 6*a13*a22*a23*a31**2*b1*b2 - 6*a13*a22*a23*a31*a32*b1**2 + 2*a13*a22*a23*a31*b2*b3 + 2*a13*a22*a31**2*b1*b3 + 2*a13*a22*a31*a32*b2*b3 - 2*a13*a22*a31*a33*b1**2 - 2*a13*a22*a31*a33*b2**2, a11**2*a22**2*a33*b3**2 - a11**2*a22*a23*a32*b3**2 - a11**2*a22*a23*a33*b2*b3 - a11**2*a22*a32*a33*b2*b3 + a11**2*a22*a33**2*b2**2 + a11**2*a23**2*a32*b2*b3 + a11**2*a23*a32**2*b2*b3 - a11**2*a23*a32*a33*b2**2 - 2*a11*a12*a21*a22*a33*b3**2 + a11*a12*a21*a23*a32*b3**2 + a11*a12*a21*a23*a33*b2*b3 + a11*a12*a21*a32*a33*b2*b3 - a11*a12*a21*a33**2*b2**2 + a11*a12*a22*a23*a31*b3**2 + a11*a12*a22*a23*a33*b1*b3 + a11*a12*a22*a31*a33*b2*b3 - a11*a12*a22*a33**2*b1*b2 - a11*a12*a23**2*a31*b2*b3 - a11*a12*a23**2*a32*b1*b3 - 2*a11*a12*a23*a31*a32*b2*b3 + a11*a12*a23*a31*a33*b2**2 + a11*a12*a23*a32*a33*b1*b2 + a11*a13*a21*a22*a32*b3**2 + a11*a13*a21*a22*a33*b2*b3 - 2*a11*a13*a21*a23*a32*b2*b3 - a11*a13*a21*a32**2*b2*b3 + a11*a13*a21*a32*a33*b2**2 - a11*a13*a22**2*a31*b3**2 - a11*a13*a22**2*a33*b1*b3 + a11*a13*a22*a23*a31*b2*b3 + a11*a13*a22*a23*a32*b1*b3 + a11*a13*a22*a31*a32*b2*b3 - 2*a11*a13*a22*a31*a33*b2**2 + a11*a13*a22*a32*a33*b1*b2 + a11*a13*a23*a31*a32*b2**2 - a11*a13*a23*a32**2*b1*b2 + a11*a21*a22*a32*a33*b1*b3 - a11*a21*a22*a33**2*b1*b2 - a11*a21*a23*a32**2*b1*b3 + a11*a21*a23*a32*a33*b1*b2 - a11*a22**2*a31*a33*b1*b3 + a11*a22**2*a33**2*b1**2 + a11*a22*a23*a31*a32*b1*b3 + a11*a22*a23*a31*a33*b1*b2 - 2*a11*a22*a23*a32*a33*b1**2 - a11*a23**2*a31*a32*b1*b2 + a11*a23**2*a32**2*b1**2 + a12**2*a21**2*a33*b3**2 - a12**2*a21*a23*a31*b3**2 - a12**2*a21*a23*a33*b1*b3 - a12**2*a21*a31*a33*b2*b3 + a12**2*a21*a33**2*b1*b2 + a12**2*a23**2*a31*b1*b3 + a12**2*a23*a31**2*b2*b3 - a12**2*a23*a31*a33*b1*b2 - a12*a13*a21**2*a32*b3**2 - a12*a13*a21**2*a33*b2*b3 + a12*a13*a21*a22*a31*b3**2 + a12*a13*a21*a22*a33*b1*b3 + a12*a13*a21*a23*a31*b2*b3 + a12*a13*a21*a23*a32*b1*b3 + a12*a13*a21*a31*a32*b2*b3 + a12*a13*a21*a31*a33*b2**2 - 2*a12*a13*a21*a32*a33*b1*b2 - 2*a12*a13*a22*a23*a31*b1*b3 - a12*a13*a22*a31**2*b2*b3 + a12*a13*a22*a31*a33*b1*b2 - a12*a13*a23*a31**2*b2**2 + a12*a13*a23*a31*a32*b1*b2 - a12*a21**2*a32*a33*b1*b3 + a12*a21**2*a33**2*b1*b2 + a12*a21*a22*a31*a33*b1*b3 - a12*a21*a22*a33**2*b1**2 + a12*a21*a23*a31*a32*b1*b3 - 2*a12*a21*a23*a31*a33*b1*b2 + a12*a21*a23*a32*a33*b1**2 - a12*a22*a23*a31**2*b1*b3 + a12*a22*a23*a31*a33*b1**2 + a12*a23**2*a31**2*b1*b2 - a12*a23**2*a31*a32*b1**2 + a13**2*a21**2*a32*b2*b3 - a13**2*a21*a22*a31*b2*b3 - a13**2*a21*a22*a32*b1*b3 - a13**2*a21*a31*a32*b2**2 + a13**2*a21*a32**2*b1*b2 + a13**2*a22**2*a31*b1*b3 + a13**2*a22*a31**2*b2**2 - a13**2*a22*a31*a32*b1*b2 + a13*a21**2*a32**2*b1*b3 - a13*a21**2*a32*a33*b1*b2 - 2*a13*a21*a22*a31*a32*b1*b3 + a13*a21*a22*a31*a33*b1*b2 + a13*a21*a22*a32*a33*b1**2 + a13*a21*a23*a31*a32*b1*b2 - a13*a21*a23*a32**2*b1**2 + a13*a22**2*a31**2*b1*b3 - a13*a22**2*a31*a33*b1**2 - a13*a22*a23*a31**2*b1*b2 + a13*a22*a23*a31*a32*b1**2]

    return coeff






def _characteristic_matrix(radii, rt_ax, phi, exponent=1.):
    """Find ellipsoid characteristic matrix.

    Keyword arguments:
    radii -- radius values for ellipse
    rt_ax -- rotation axis
    phi -- rotation angle about rotation axis
    exponent -- raise characteristic matrix to this exponential power

    Return values:
    X -- characteristic matrix


    X = Q^T * O^-n*2 * Q, where Q is rotation matrix given by angle phi and
    O is diagonal matrix with radii values along main diagonal


    """

    rt_ax = np.matrix(rt_ax).transpose()

    A = np.matrix( np.diag(radii ** (-2. * exponent)) )

    Q = _quaternion_to_rotation_matrix(rt_ax, phi)

    Q = np.matrix(Q)



    X = Q.T * A * Q



    return X




def _quaternion_to_rotation_matrix(rt_ax, phi):
    """

    Sources:
    Shoemake, K., Animating rotation with quaternion curves, ACM, 19(3), 1985

    """

    w = np.cos(phi / 2.)
    x = np.sin(phi / 2.) * rt_ax[0]
    y = np.sin(phi / 2.) * rt_ax[1]
    z = np.sin(phi / 2.) * rt_ax[2]

    M = np.zeros((3,3))

    M[0, 0] = 1. - 2. * (y ** 2 + z ** 2)
    M[0, 1] = 2. * (x * y + w * z)
    M[0, 2] = 2. * (x * z - w * y)

    M[1, 0] = 2. * (x * y - w * z)
    M[1, 1] = 1. - 2. * (x ** 2 + z ** 2)
    M[1, 2] = 2. * (y * z + w * x)

    M[2, 0] = 2. * (x * z + w * y)
    M[2, 1] = 2. * (y * z - w * x)
    M[2, 2] = 1. - 2. * (x ** 2 + y ** 2)



    return M





def rotation_matrix(alpha, beta, gamma):
    Rx = np.array([[1., 0., 0.],
                   [0., np.cos(alpha), -np.sin(alpha)],
                   [0., np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0., np.sin(beta)],
                   [0., 1., 0.],
                   [-np.sin(beta), 0., np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0.],
                   [np.sin(gamma), np.cos(gamma), 0.],
                   [0., 0., 1.]])
    R = Rz.dot(Ry.dot(Rx))

    return R














def _characteristic_matrix_old(radii, alpha, beta, gamma):
    A = np.diag(1. / radii ** 2)
    Rx = np.array([[1., 0., 0.],
                   [0., np.cos(alpha), -np.sin(alpha)],
                   [0., np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0., np.sin(beta)],
                   [0., 1., 0.],
                   [-np.sin(beta), 0., np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0.],
                   [np.sin(gamma), np.cos(gamma), 0.],
                   [0., 0., 1.]])
    R = Rz.dot(Ry.dot(Rx))
    B = R.transpose().dot(A).dot(R)


    return B










def _generate_ellipsoid_volume_old(x, y, z, center, radii, alpha, beta, gamma):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place an ellipsoid at center having given radii and according to
    alpha, beta, and gamma.

    The point x = [xi, yi, zi] is inside an ellipoid if 
    (x - center)^T R^T A R (x - center) <= 1.

    center = [xc, yc, zc] is the center point of the ellipsoid
    Rx = [[cos(theta), -sin(theta)], is clock-wise rotation of theta radians
         [sin(theta), cos(theta)]]
    Ry =
    Rz =
    A = [[1/a^2, 0, 0], records major and minor axis radii
         [0, 1/b^2, 0],
         [0, 0, 1/c^2]]

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    z -- extent of the volume along z-axis
    radii -- sphere radius
    center -- center point of the sphere
    alpha --
    beta --
    gamma -- 

    Return values:
    vol -- generated volume containing ellipse

    """

    # Setup matrices for quadratic form evaluation
    A = np.diag(1. / radii ** 2)
    Rx = np.array([[1., 0., 0.],
                   [0., np.cos(alpha), -np.sin(alpha)],
                   [0., np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0., np.sin(beta)],
                   [0., 1., 0.],
                   [-np.sin(beta), 0., np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0.],
                   [np.sin(gamma), np.cos(gamma), 0.],
                   [0., 0., 1.]])
    R = Rz.dot(Ry.dot(Rx))
    B = R.transpose().dot(A).dot(R)

    # Form square position array for x and y
    X = np.zeros((3, len(y) * len(x) * len(z)))
    X[0] = np.tile(x - center[0], (len(z), len(y), 1)).flatten()
    X[1] = np.tile(y - center[1], (len(z), len(x), 1)).transpose(0, 2, 1).flatten()
    X[2] = np.tile(z - center[2], (len(y), len(x), 1)).transpose(2, 0, 1).flatten()

    # Indicate points inside the ellipse
    vol = (X * (B.dot(X))).sum(0)
    vol = np.reshape(vol, (len(z), len(y), len(x)))
    vol = vol <= 1.
    vol = vol.astype(float)

    return vol