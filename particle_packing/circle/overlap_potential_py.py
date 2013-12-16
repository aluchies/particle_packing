import numpy as np

def overlap_potential_py(rA, radiiA, rB, radiiB):
    """

    Overlap potential function (Python version) provides a distance measure
    for circles A and B.

    Overlap criterion based on the overlap potential value:
    F(A,B) > 1, A and B are disjoint
    F(A,B) = 0, A and B are externally tangent
    F(A,B) < 1, A and B are overlapping

    Keyword arguments:
    rA -- center of circle A
    radiiA -- radii of circle A
    rB -- center of circle B
    radiiB -- radii of circle B

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

    rB = np.asarray(rB).flatten()
    if len(rB) != 2:
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

    return (rAB[0] ** 2 + rAB[1] ** 2) / (radiiA + radiiB) ** 2