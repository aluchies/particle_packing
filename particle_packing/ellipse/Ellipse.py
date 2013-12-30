import numpy as np


class Ellipse(object):
    """

    """

    def __init__(self, center, radii, phi):
        """  
        """

        center = np.asarray(center, dtype=np.float).flatten()
        if len(center) != 2:
            raise ValueError('incorrect input for center input')
        else:
            self.center = center

        radii = np.asarray(radii, dtype=np.float).flatten()
        if len(radii) == 1:
            self.radii = radii * np.ones(2)
        elif len(radii) == 2:
            self.radii = radii
        else:
            raise ValueError('incorrect input for radii input')

        phi = float(phi)
        self.phi = phi


    def generate_volume(self, x_ax, y_ax):
        """Generate two-dimensional volume for the ellipse.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis

        Return values:
        vol -- 2D array for the volume

        """

        ndim = len(self.center)

        x_ax = np.asarray(x_ax, dtype=np.float).flatten()
        y_ax = np.asarray(y_ax, dtype=np.float).flatten()

        vol = _generate_ellipse_volume(x_ax, y_ax, self.center, self.radii,
            self.phi)

        return vol



    def find_subvolume(self, x_ax, y_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        sphere centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis

        Return values:
        x_ax_subvol -- x-axis for subvolume
        y_ax_subvol -- y-axis for subvolume

        """

        x_ax_subvol, y_ax_subvol, \
        x_ax_subvol_ix, y_ax_subvol_ix = \
        _find_circle_subvolume(x_ax, y_ax, self.center, max(self.radii))

        return x_ax_subvol, y_ax_subvol


    def find_subvolume_ix(self, x_ax, y_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        sphere centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis

        Return values:
        x_ax_subvol -- x-axis for subvolume
        y_ax_subvol -- y-axis for subvolume

        """

        x_ax_subvol, y_ax_subvol, \
        x_ax_subvol_ix, y_ax_subvol_ix = \
        _find_circle_subvolume(x_ax, y_ax, self.center, max(self.radii))

        return x_ax_subvol_ix, y_ax_subvol_ix


    def overlap_potential(self, c):
        """Determine the overlap potential of object self and object c.

        Overlap criterion based on the overlap potential value:
        F(A,B) > 1, A and B are disjoint
        F(A,B) = 1, A and B are externally tangent
        F(A,B) < 1, A and B are overlapping

        Input arguments:
        c -- object to check for overlap with self

        Return values:
        F -- overlap potential value

        """

        if not isinstance(c, Ellipse):
            raise ValueError('input is not an ellipse')


        F = overlap_potential_py(self.center, self.radii, self.phi,
            c.center, c.radii, c.phi)

        return F



    def square_container_potential(self):
        """Determine if object is contained in the container.

        Containment criterion based on the overlap potential value:
        F(A,B) > 1, object completely inside container
        F(A,B) = 1, object completely inside and tangent to container
        F(A,B) < 1, object at least partially outside container


        Return values:
        F -- overlap potential value

        """

        F = square_container_potential_py(self.center, self.radii, self.phi)

        return F




def _generate_ellipse_volume(x, y, center, radii, phi):
    """Generate the volume having x- and y-axes given by x, and y. In the
    volume, place an ellipse at center having given radii and rotate phi
    radians clockwise.

    The point x = [xi, yi] is inside an ellipse if 
    (x - center)^T R^T A R (x - center) <= 1.

    center = [xc, yc] is the center point of the ellipse
    R = [[cos(phi), -sin(phi)], is clock-wise rotation of phi radians
         [sin(phi), cos(phi)]]
    A = [[1/a^2, 0], records major and minor axis radii
         [0, 1/b^2]]

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    radii -- sphere radius
    center -- center point of the sphere
    phi -- ellipse rotate counter-clockwise in radians

    Return values:
    vol -- generated volume containing ellipse

    """

    # Setup matrices for quadratic form evaluation
    A = np.diag(1. / radii ** 2)
    R = np.array([[np.cos(phi), np.sin(phi)],
        [-np.sin(phi), np.cos(phi)]])
    B = R.transpose().dot(A).dot(R)

    # Form square position array for x and y
    X = np.zeros((2, len(y) * len(x)))
    X[0] = np.tile(x - center[0], (len(y), 1)).flatten()
    X[1] = np.tile(y - center[1], (len(x), 1)).transpose().flatten()

    # Indicate points inside the ellipse
    vol = (X * (B.dot(X))).sum(0)
    vol = np.reshape(vol, (len(y), len(x)))
    vol = vol <= 1.
    vol = vol.astype(float)

    return vol




def _find_circle_subvolume(X, Y, xi, a):
    """Extract sub-volume from the volume with x-, y-axes given by X, Y.
    The sub-volume is just large enough to contain an ellipse centered at xi
    and having maximum radius a.

    Keyword arguments:
    X -- extent of the volume along x-axis
    Y -- extent of the volume along y-axis
    xi -- ellipse center point
    a -- ellipse radius


    Return values:
    X_subvol -- extent of the sub-volume along x-axis
    Y_subvol -- extent of the sub-volume along y-axis
    X_subvol_ix -- X_subvol = X[X_subvol_ix.min():X_subvol_ix.max()]
    Y_subvol_ix -- Y_subvol = Y[Y_subvol_ix.min():Y_subvol_ix.max()]

    """

    # Find smallest cube that contains the sphere
    X_subvol_ix = np.nonzero(np.abs(X - xi[0]) <= a)[0]
    Y_subvol_ix = np.nonzero(np.abs(Y - xi[1]) <= a)[0]

    # Get axis arrays for the sub-volume
    X_subvol = X[X_subvol_ix]
    Y_subvol = Y[Y_subvol_ix]

    return X_subvol, Y_subvol, X_subvol_ix, Y_subvol_ix






















def overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB):
    """

    Overlap potential function (Python version) provides a distance measure
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

    rA = np.asarray(rA).flatten()
    if len(rA) != 2:
        raise ValueError('input error for rA')
    rA = np.matrix(rA).T

    rB = np.asarray(rB).flatten()
    if len(rB) != 2:
        raise ValueError('input error for rB')
    rB = np.matrix(rB).T


    radiiA = np.asarray(radiiA.flatten())
    if len(radiiA) != 2:
        raise ValueError('input error for radiiA')

    radiiB = np.asarray(radiiB.flatten())
    if len(radiiB) != 2:
        raise ValueError('input error for radiiB')

    phiA = float(phiA)
    phiB = float(phiB)







    """

    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and
    O is diagonal matrix with radii values along main diagonal
    
    What's needed is XA^-1 = Q^T * O^2 * Q

    """


    O = np.matrix(np.diag(radiiA ** 2))

    Q = np.matrix([[np.cos(phiA), np.sin(phiA)],
              [-np.sin(phiA), np.cos(phiA)]])

    XA = Q.T * O * Q



    """

    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and
    O is diagonal matrix with radii values along main diagonal

    What's needed is XB^1/2 = Q^T * O^-1 * Q

    """


    O = np.matrix(np.diag(radiiB ** -1))

    Q = np.matrix([[np.cos(phiB), np.sin(phiB)],
              [-np.sin(phiB), np.cos(phiB)]])

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


    a11, a12, a21, a22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    b1, b2 = b[0, 0], b[1, 0]

    return L**3*(a11*b2**2 - a12*b1*b2 - a21*b1*b2 + a22*b1**2 - b1**2 - \
        b2**2) + L**2*(-2*a11*b2**2 + 2*a12*b1*b2 + 2*a21*b1*b2 - \
        2*a22*b1**2 + b1**2 + b2**2) + L*(a11*b2**2 - a12*b1*b2 - a21*b1*b2 + \
        a22*b1**2)

def q_AB(A, b, L):
    """

    Denominator of function f.

    """


    a11, a12, a21, a22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    b1, b2 = b[0, 0], b[1, 0]

    return L**2*(a11*a22 - a11 - a12*a21 - a22 + 1) + L*(-2*a11*a22 + a11 + \
        2*a12*a21 + a22) + a11*a22 - a12*a21

def f_AB(A, b, L):
    """

    """


    return p_AB(A, b, L) / q_AB(A, b, L)

def h_coeffs(A, b):
    """

    Numerator of the derivative of the function f.

    """

    a11, a12, a21, a22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    b1, b2 = b[0, 0], b[1, 0]

    coeff = [

        a11**2*a22*b2**2 - a11**2*b2**2 - a11*a12*a21*b2**2 - \
        a11*a12*a22*b1*b2 + a11*a12*b1*b2 - a11*a21*a22*b1*b2 + \
        a11*a21*b1*b2 + a11*a22**2*b1**2 - 2*a11*a22*b1**2 - \
        2*a11*a22*b2**2 + a11*b1**2 + 2*a11*b2**2 + a12**2*a21*b1*b2 + \
        a12*a21**2*b1*b2 - a12*a21*a22*b1**2 + a12*a21*b1**2 + \
        a12*a21*b2**2 + a12*a22*b1*b2 - a12*b1*b2 + a21*a22*b1*b2 - \
        a21*b1*b2 - a22**2*b1**2 + 2*a22*b1**2 + a22*b2**2 - b1**2 - b2**2,

        -4*a11**2*a22*b2**2 + 2*a11**2*b2**2 + 4*a11*a12*a21*b2**2 + \
        4*a11*a12*a22*b1*b2 - 2*a11*a12*b1*b2 + 4*a11*a21*a22*b1*b2 - \
        2*a11*a21*b1*b2 - 4*a11*a22**2*b1**2 + 6*a11*a22*b1**2 + \
        6*a11*a22*b2**2 - 2*a11*b1**2 - 2*a11*b2**2 - \
        4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + 4*a12*a21*a22*b1**2 - \
        4*a12*a21*b1**2 - 4*a12*a21*b2**2 - 2*a12*a22*b1*b2 - \
        2*a21*a22*b1*b2 + 2*a22**2*b1**2 - 2*a22*b1**2 - 2*a22*b2**2,


        6*a11**2*a22*b2**2 - a11**2*b2**2 - 6*a11*a12*a21*b2**2 - \
        6*a11*a12*a22*b1*b2 + a11*a12*b1*b2 - 6*a11*a21*a22*b1*b2 + \
        a11*a21*b1*b2 + 6*a11*a22**2*b1**2 - 6*a11*a22*b1**2 - \
        6*a11*a22*b2**2 + a11*b1**2 + 6*a12**2*a21*b1*b2 + \
        6*a12*a21**2*b1*b2 - 6*a12*a21*a22*b1**2 + 5*a12*a21*b1**2 + \
        5*a12*a21*b2**2 + a12*a22*b1*b2 + a12*b1*b2 + a21*a22*b1*b2 + \
        a21*b1*b2 - a22**2*b1**2 + a22*b2**2,

        -4*a11**2*a22*b2**2 + 4*a11*a12*a21*b2**2 + 4*a11*a12*a22*b1*b2 +\
        4*a11*a21*a22*b1*b2 - 4*a11*a22**2*b1**2 + 2*a11*a22*b1**2 + \
        2*a11*a22*b2**2 - 4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + \
        4*a12*a21*a22*b1**2 - 2*a12*a21*b1**2 - 2*a12*a21*b2**2,

        a11**2*a22*b2**2 - a11*a12*a21*b2**2 - a11*a12*a22*b1*b2 - \
        a11*a21*a22*b1*b2 + a11*a22**2*b1**2 + a12**2*a21*b1*b2 + \
        a12*a21**2*b1*b2 - a12*a21*a22*b1**2

        ]


    return coeff















def square_container_potential_py(r, radii, phi):
    """

    Container potential function (Python version) determines if an ellipse
    is inside a square container or not.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 1, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of ellipse A
    radiiA -- radii of ellipse A
    phiA -- rotation angle of ellipse A

    Return values:
    F -- overlap potential value


    """


    # Top
    rT = np.array([0.5, 2.]).T
    radiiT = np.array([float("inf"), 1.])
    phiT = 0.
    top = overlap_potential_py(r, radii, phi, rT, radiiT, phiT)

    # Bottom
    rT = np.array([0.5, -1.]).T
    radiiT = np.array([float("inf"), 1.])
    phiT = 0.
    bottom = overlap_potential_py(r, radii, phi, rT, radiiT, phiT)

    # Left
    rL = np.array([-1., 0.5]).T
    radiiL = np.array([1, float("inf")])
    phiL = 0.
    left = overlap_potential_py(r, radii, phi, rL, radiiL, phiL)

    # Right
    rR = np.array([2., 0.5]).T
    radiiR = np.array([1, float("inf")])
    phiR = 0.
    right = overlap_potential_py(r, radii, phi, rR, radiiR, phiR)




    return min(left, right, top, bottom)



