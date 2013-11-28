import numpy as np

class Ellipsoid(object):
    """

    """

    def __init__(self, center, radii, alpha, beta, gamma):
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

        alpha = float(alpha)
        self.alpha = alpha

        beta = float(beta)
        self.beta = beta

        gamma = float(gamma)
        self.gamma = gamma





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
            self.center, self.radii, self.alpha, self.beta, self.gamma)

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









def _generate_ellipsoid_volume(x, y, z, center, radii, alpha, beta, gamma):
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