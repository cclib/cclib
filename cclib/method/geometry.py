import numpy

def getDistance(coords0, coords1):
    """Distance between atoms 0-1"""
    return numpy.linalg.norm(numpy.subtract(coords0 - coords1))

def getAngle(coords0, coords1, coords2):
    """Angle(in radians) between atoms 0-1-2"""
    return math.acos(((coords1[0]**2 + coords1[1]**2 + coords1[2]**2) - (((coords0[0] + coords2[0]) * coords1[0]) + ((coords0[1] + coords2[1]) * coords1[1]) + ((coords0[2] + coords2[2]) * coords1[2]))) / (getDistance(coords0, coords1) * getDistance(coords1, coords2)))

def getDihedral(coords0, coords1, coords2, coords3):
    """Dihedral angle(in radians) between the planes containing the atoms 0-1 and 2-3"""
    """Praxeolitic formula[https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python]
    1 sqrt, 1 cross product"""

    b0 = -1.0*(coords1 - coords0)
    b1 = coords2 - coords1
    b2 = coords3 - coords2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= numpy.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - numpy.dot(b0, b1)*b1
    w = b2 - numpy.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = numpy.dot(v, w)
    y = numpy.dot(np.cross(b1, v), w)
    return numpy.arctan2(y, x)
