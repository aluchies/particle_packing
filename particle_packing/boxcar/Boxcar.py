import numpy as np

class Boxcar(object):
    """

    """

    def __init__(self, center, radius):
        """  
        """

        center = float(center)
        self.center = center

        radius = float(radius)
        self.radius = radius



    def generate_volume(self, x_ax):
        """Generate volume for the ellipsoid.

        Keyword arguments:
        x_ax -- numpy array for x-axis


        Return values:
        vol -- 1 array for the volume

        """

        x_ax = np.asarray(x_ax, dtype=np.float).flatten()

        vol = _generate_boxcar_volume(x_ax, self.radius, self.center)

        return vol


    def find_subvolume(self, x_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        boxcar centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis


        Return values:
        x_ax_subvol -- x-axis for subvolume


        """

        x_ax_subvol, x_ax_subvol_ix = \
        _find_boxcar_subvolume(x_ax, self.center, self.radius)

        return x_ax_subvol


    def find_subvolume_ix(self, x_ax):
        """Extract sub-volume from the volume with x-, y-, z-axes given by
        x_ax, y_ax, z_ax. The sub-volume is just large enough to contain a
        boxcar centered at xi and having radius a.

        Keyword arguments:
        x_ax -- numpy array for x-axis


        Return values:
        x_ax_subvol_ix -- index values (x_ax) for x-axis for subvolume


        """

        x_ax_subvol, x_ax_subvol_ix = \
        _find_boxcar_subvolume(x_ax, self.center, self.radius)

        return x_ax_subvol_ix









def _generate_boxcar_volume(x, radius, center):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place a boxcar centered at xi and having radius a.

    Keyword arguments:
    x -- extent of the volume along x-axis
    radius -- boxcar radius
    center -- center point of the boxcar

    Return values:
    subvol -- generated sub-volume containing boxcar

    """

    # Form cubic position array for x, y, z
    X_cube = x.copy()


    # Find all points inside boxcar inside the cube
    vol = np.sqrt((X_cube - center) ** 2 / radius ** 2)
    vol = vol <= 1

    return vol.astype(float)


def _find_boxcar_subvolume(X, Y, Z, xi, a):
    """Extract sub-volume from the volume with x-, y-, z-axes given by X, Y, Z.
    The sub-volume is just large enough to contain a boxcar centered at xi and
    having radius a.

    Keyword arguments:
    X -- extent of the volume along x-axis
    xi -- boxcar center point
    a -- boxcar radius


    Return values:
    X_subvol -- extent of the sub-volume along x-axis
    X_subvol_ix -- X_subvol = X[X_subvol_ix.min():X_subvol_ix.max()]

    """

    # Find smallest cube that contains the boxcar
    X_subvol_ix = np.nonzero(np.abs(X - xi[0]) <= a)[0]

    # Get axis arrays for the sub-volume
    X_subvol = X[X_subvol_ix]

    return X_subvol, X_subvol_ix