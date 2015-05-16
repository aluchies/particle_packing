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



    def generate_volume_spherical_gaussian(self, x_ax, y_ax, z_ax):
        """Generate volume for a spherical Gaussian.

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

        vol = _generate_spherical_gaussian_volume(x_ax, y_ax, z_ax, self.radius, self.center)

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

        if not isinstance(c, Sphere):
            raise ValueError('input is not an sphere')


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


        if not isinstance(c, Sphere):
            raise ValueError('input is not a sphere')


        G = contain_potential_py(self.center, self.radius,
            c.center, c.radius)

        return G







    def container_potential(self, shape):
        """Determine container potential for the object.

        Containment criterion based on the contain potential value:
        G(A,B) > 1, A completely inside the container
        G(A,B) = 1, A completely inside and tangent to container
        G(A,B) < 1, A at least partially outside container

        Keyword arguments:
        shape -- container shape, 'square' or 'circle'


        Return values:
        G -- contain potential value

        """

        # input checking
        if (not shape is 'cube') and (not shape is 'sphere'):
            raise ValueError('shape input is unknown')

        if shape is 'sphere':
            center = [0.5, 0.5, 0.5]
            radius = 0.5
            G = contain_potential_py(self.center, self.radius,
            center, radius)

        elif shape is 'cube':
            G = container_potential_cube_py(self.center, self.radius)

        return G












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

    Notes
    This function is meant to be used in conjuction with 
    _find_sphere_subvolume().

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






def _generate_spherical_gaussian_volume(x, y, z, radius_eff, center):
    """Generate the volume having x-, y-, and z-axes given by x, y, z. In the
    volume, place a spherical Gaussian centered at xi and having effective
    radius radius_eff.

    Keyword arguments:
    x -- extent of the volume along x-axis
    y -- extent of the volume along y-axis
    z -- extent of the volume along z-axis
    radius_eff -- effective radius
    center -- center point of the sphere

    Return values:
    subvol -- generated sub-volume containing sphere

    Notes
    This function is meant to be used in conjuction with 
    _find_sphere_subvolume().

    """

    # Form cubic position array for x, y, z
    X_cube = np.tile(x, (len(z), len(y), 1))
    Y_cube = np.tile(y, (len(z), len(x), 1)).transpose(0, 2, 1)
    Z_cube = np.tile(z, (len(y), len(x), 1)).transpose(2, 0, 1)

    # Find all points inside sphere inside the cube
    sigma = (3.0 * np.sqrt(np.pi / 2.0)) ** (1.0 / 3.0) * radius_eff

    vol = np.exp(- ((X_cube - center[0]) ** 2 + (Y_cube - center[1]) ** 2 +
                    (Z_cube - center[2]) ** 2) / sigma ** 2)

    return vol






def _find_sphere_subvolume(X, Y, Z, center, a):
    """Extract sub-volume from the volume with x-, y-, z-axes given by X, Y, Z.
    The sub-volume is just large enough to contain a sphere centered at xi and
    having radius a.

    Keyword arguments:
    X -- extent of the volume along x-axis
    Y -- extent of the volume along y-axis
    Z -- extent of the volume along z-axis
    center -- sphere center point
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
    X_subvol_ix = np.nonzero(np.abs(X - center[0]) <= a)[0]
    Y_subvol_ix = np.nonzero(np.abs(Y - center[1]) <= a)[0]
    Z_subvol_ix = np.nonzero(np.abs(Z - center[2]) <= a)[0]

    # Get axis arrays for the sub-volume
    X_subvol = X[X_subvol_ix]
    Y_subvol = Y[Y_subvol_ix]
    Z_subvol = Z[Z_subvol_ix]

    return X_subvol, Y_subvol, Z_subvol, X_subvol_ix, Y_subvol_ix, Z_subvol_ix










def overlap_potential_py(rA, radiiA, rB, radiiB):
    """

    Overlap potential function (Python version) provides a distance measure
    for spheres A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 1, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of sphere A
    radiiA -- radii of sphere A
    rB -- center of sphere B
    radiiB -- radii of sphere B

    Return values:
    F -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """


    """Input argument checking."""

    rA = np.asarray(rA).flatten()
    if len(rA) != 3:
        raise ValueError('input error for rA')

    rB = np.asarray(rB).flatten()
    if len(rB) != 3:
        raise ValueError('input error for rB')


    radiiA = np.asarray(radiiA).flatten()
    if len(radiiA) != 1:
        raise ValueError('input error for radiiA')
    radiiA = radiiA[0]

    radiiB = np.asarray(radiiB).flatten()
    if len(radiiB) != 1:
        raise ValueError('input error for radiiB')
    radiiB = radiiB[0]


    rAB = rB - rA

    return (rAB[0] ** 2 + rAB[1] ** 2 + rAB[2] ** 2) / (radiiA + radiiB) ** 2










def contain_potential_py(rA, radiiA, rB, radiiB):
    """

    Contain potential function (Python version) provides a distance measure
    for spheres A and B.

    Criterion based on the contain potential value:
    G(A,B) > 1, A entirely inside B
    G(A,B) = 1, A entirely inside and tangent to B
    G(A,B) < 1, A is partly or entirely outside of B

    Keyword arguments:
    rA -- center of sphere A
    radiiA -- radii of sphere A
    rB -- center of sphere B
    radiiB -- radii of sphere B

    Return values:
    G -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """


    """Input argument checking."""

    rA = np.asarray(rA).flatten()
    if len(rA) != 3:
        raise ValueError('input error for rA')

    rB = np.asarray(rB).flatten()
    if len(rB) != 3:
        raise ValueError('input error for rB')


    radiiA = np.asarray(radiiA).flatten()
    if len(radiiA) != 1:
        raise ValueError('input error for radiiA')
    radiiA = radiiA[0]

    radiiB = np.asarray(radiiB).flatten()
    if len(radiiB) != 1:
        raise ValueError('input error for radiiB')
    radiiB = radiiB[0]


    rAB = rB - rA

    if np.linalg.norm(rAB) < (radiiB - radiiA):
        return 2.
    elif np.linalg.norm(rAB) == (radiiB - radiiA):
        return 1.
    else:
        return 0.
















def container_potential_cube_py(rA, radiiA):
    """

    Container potential function (Python version) provides a distance measure
    for how far A is inside the unit cube

    Overlap criterion based on the overlap potential value:
    G(A,B) > 1, A entirely inside container
    G(A,B) = 1, A entirely inside and tangent to container
    G(A,B) < 1, A is partly or entirely outside of container

    Keyword arguments:
    rA -- center of circle A
    radiiA -- radius of circle A

    Return values:
    G -- overlap potential value

    Sources:
    Donev, A, et. al., Neighbor list collision-driven molecular dynamics
    simulation for nonspherical hard particles. II. Applications to ellipses
    and ellipsoids, J. of Comp. Physics, vol 202, 2004.

    """


    """Input argument checking."""

    rA = np.asarray(rA).flatten()
    if len(rA) != 3:
        raise ValueError('input error for rA')


    radiiA = np.asarray(radiiA).flatten()
    if len(radiiA) != 1:
        raise ValueError('input error for radiiA')
    radiiA = radiiA[0]

    rB = 0.5 * np.ones(3)
    radiiB = 0.5


    rAB = rB - rA


    a = (radiiB ** 2 - rAB[0] ** 2) / radiiA ** 2
    b = (radiiB ** 2 - rAB[1] ** 2) / radiiA ** 2
    c = (radiiB ** 2 - rAB[2] ** 2) / radiiA ** 2

    return max(a, b, c)
