import numpy as np

class Sphere(object):
    """

    """

    def __init__(self, center, radius):
        """  
        """

        center = np.asarray(center, dtype=np.float).flatten()
        if len(center) != 3:
            raise ValueError('incorrect input for center input')
        else:
            self.center = center

        radius = float(radius)
        self.radius = radius



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

        vol = _generate_sphere_volume(x_ax, y_ax, z_ax, self.radius, self.center)

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
        _find_sphere_subvolume(x_ax, y_ax, z_ax, self.center, self.radius)

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
        _find_sphere_subvolume(x_ax, y_ax, z_ax, self.center, self.radius)

        return x_ax_subvol_ix, y_ax_subvol_ix, z_ax_subvol_ix









def _generate_sphere_volume(x, y, z, radius, center):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place a sphere centered at xi and having radius a.

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    z -- extent of the volume along z-axis
    radius -- sphere radius
    center -- center point of the sphere

    Return values:
    subvol -- generated sub-volume containing sphere

    """

    # Form cubic position array for x, y, z
    X_cube = np.tile(x, (len(z), len(y), 1))
    Y_cube = np.tile(y, (len(z), len(x), 1)).transpose(0, 2, 1)
    Z_cube = np.tile(z, (len(y), len(x), 1)).transpose(2, 0, 1)

    # Find all points inside sphere inside the cube
    vol = np.sqrt((X_cube - center[0]) ** 2 / radius ** 2 +
        (Y_cube - center[1]) ** 2 / radius ** 2 +
        (Z_cube - center[2]) ** 2 / radius ** 2)
    vol = vol <= 1

    return vol.astype(float)


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