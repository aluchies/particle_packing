import numpy as np
from overlap_potential_py import overlap_potential_py

class Circle(object):
    """

    """

    def __init__(self, center, radius):
        """  
        """

        center = np.asarray(center, dtype=np.float).flatten()
        if len(center) != 2:
            raise ValueError('incorrect input for center input')
        else:
            self.center = center

        radius = float(radius)
        self.radius = radius



    def generate_volume(self, x_ax, y_ax):
        """Generate two-dimensional volume for the circle.

        Keyword arguments:
        x_ax -- numpy array for x-axis
        y_ax -- numpy array for y-axis

        Return values:
        vol -- 2D array for the volume

        """

        ndim = len(self.center)

        x_ax = np.asarray(x_ax, dtype=np.float).flatten()
        y_ax = np.asarray(y_ax, dtype=np.float).flatten()

        vol = _generate_circle_volume(x_ax, y_ax, self.center, self.radius)

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
        _find_circle_subvolume(x_ax, y_ax, self.center, self.radius)

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
        _find_circle_subvolume(x_ax, y_ax, self.center, self.radius)

        return x_ax_subvol_ix, y_ax_subvol_ix




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

        if not isinstance(c, Circle):
            raise ValueError('input is not a circle')


        F = overlap_potential_py(self.center, self.radius,
            c.center, c.radius)

        return F







def _generate_circle_volume(x, y, center, radius):
    """Generate the volume having x- and y-axes given by x, and y. In the
    volume, place an circle at center having given radius and rotate phi
    radians clockwise.

    The point x = [xi, yi] is inside an circle if 
    (x - center)^T R^T A R (x - center) <= 1.

    center = [xc, yc] is the center point of the circle
    R = [[cos(phi), -sin(phi)], is clock-wise rotation of phi radians
         [sin(phi), cos(phi)]]
    A = [[1/a^2, 0], records major and minor axis radius
         [0, 1/b^2]]

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    radius -- sphere radius
    center -- center point of the sphere
    phi -- circle rotate counter-clockwise in radians

    Return values:
    vol -- generated volume containing circle

    """




    # Form cubic position array for x, y, z
    X_cube = np.tile(x, (len(y), 1))
    Y_cube = np.tile(y, (len(x), 1)).transpose()


    # Find all points inside sphere inside the cube
    vol = np.sqrt((X_cube - center[0]) ** 2 / radius ** 2 +
        (Y_cube - center[1]) ** 2 / radius ** 2)
    vol = vol <= 1

    return vol.astype(float)







def _find_circle_subvolume(X, Y, xi, a):
    """Extract sub-volume from the volume with x-, y-axes given by X, Y.
    The sub-volume is just large enough to contain an circle centered at xi
    and having maximum radius a.

    Keyword arguments:
    X -- extent of the volume along x-axis
    Y -- extent of the volume along y-axis
    xi -- circle center point
    a -- circle radius


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


