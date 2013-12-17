import numpy as np
from overlap_potential_py import overlap_potential_py


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

        if not isinstance(c, Ellipse):
            raise ValueError('input is not an ellipse')


        F = overlap_potential_py(self.center, self.radii, self.phi,
            c.center, c.radii, c.phi)

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
    R = np.array([[np.cos(phi), -np.sin(phi)],
        [np.sin(phi), np.cos(phi)]])
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


