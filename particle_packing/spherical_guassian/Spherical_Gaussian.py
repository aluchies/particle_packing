import numpy as np

class Spherical_Gaussian(object):
    """

    """

    def __init__(self, center, radius_eff):
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
        """Generate volume for the spherical gaussian.

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

        vol = _generate_spherical_guassian_volume(x_ax, y_ax, z_ax, \
            self.radius_eff, self.center)

        return vol











def _generate_spherical_gaussian_volume(x, y, z, radius_eff, center):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place a sphere centered at xi and having effective radius radius_eff.

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    z -- extent of the volume along z-axis
    radius_eff -- effective radius
    center -- center point of the sphere

    Return values:
    subvol -- generated sub-volume containing sphere

    """

    # Form cubic position array for x, y, z
    X_cube = np.tile(x, (len(z), len(y), 1))
    Y_cube = np.tile(y, (len(z), len(x), 1)).transpose(0, 2, 1)
    Z_cube = np.tile(z, (len(y), len(x), 1)).transpose(2, 0, 1)

    # Find all points inside sphere inside the cube
    d = (3.0 * np.sqrt(np.pi / 2.0)) ** (1.0 / 3.0) * radius_eff

    vol = np.exp(- ((X_cube - center[0]) ** 2 + (Y_cube - center[1]) ** 2 +
                    (Z_cube - center[2]) ** 2) / d ** 2)

    return vol