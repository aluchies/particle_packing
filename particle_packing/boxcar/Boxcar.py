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

        if not isinstance(c, Boxcar):
            raise ValueError('input is not a boxcar')


        F = overlap_potential_py(self.center, self.radius,
            c.center, c.radius)

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


        if not isinstance(c, Boxcar):
            raise ValueError('input is not a boxcar')


        G = contain_potential_py(self.center, self.radius,
            c.center, c.radius)

        return G









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


def _find_boxcar_subvolume(X, xi, a):
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
    X_subvol_ix = np.nonzero(np.abs(X - xi) <= a)[0]

    # Get axis arrays for the sub-volume
    X_subvol = X[X_subvol_ix]

    return X_subvol, X_subvol_ix










def overlap_potential_py(rA, radiiA, rB, radiiB):
    """

    Overlap potential function (Python version) provides a distance measure
    for boxcars A and B.

    Criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 1, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of boxcar A
    radiiA -- radii of boxcar A
    rB -- center of boxcar B
    radiiB -- radii of boxcar B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """


    """Input argument checking."""

    rA = np.asarray(rA).flatten()
    if len(rA) != 1:
        raise ValueError('input error for rA')
    rA = rA[0]

    rB = np.asarray(rB).flatten()
    if len(rB) != 1:
        raise ValueError('input error for rB')
    rB = rB[0]


    radiiA = np.asarray(radiiA).flatten()
    if len(radiiA) != 1:
        raise ValueError('input error for radiiA')
    radiiA = radiiA[0]

    radiiB = np.asarray(radiiB).flatten()
    if len(radiiB) != 1:
        raise ValueError('input error for radiiB')
    radiiB = radiiB[0]



    return (rA - rB) ** 2 / (radiiA + radiiB) ** 2








def contain_potential_py(rA, radiiA, rB, radiiB):
    """

    Contain potential function (Python version) provides a distance measure
    for boxcars A and B.

    Criterion based on the contain potential value:
    G(A,B) > 1, A entirely inside B
    G(A,B) = 1, A entirely inside and tangent to B
    G(A,B) < 1, A is partly or entirely outside of B

    Keyword arguments:
    rA -- center of boxcar A
    radiiA -- radii of boxcar A
    rB -- center of boxcar B
    radiiB -- radii of boxcar B

    Return values:
    G -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """


    """Input argument checking."""

    rA = np.asarray(rA).flatten()
    if len(rA) != 1:
        raise ValueError('input error for rA')
    rA = rA[0]

    rB = np.asarray(rB).flatten()
    if len(rB) != 1:
        raise ValueError('input error for rB')
    rB = rB[0]


    radiiA = np.asarray(radiiA).flatten()
    if len(radiiA) != 1:
        raise ValueError('input error for radiiA')
    radiiA = radiiA[0]

    radiiB = np.asarray(radiiB).flatten()
    if len(radiiB) != 1:
        raise ValueError('input error for radiiB')
    radiiB = radiiB[0]





    return (radiiB ** 2 - (rA - rB) ** 2) / radiiA ** 2











