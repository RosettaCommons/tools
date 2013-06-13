########################################################
#Module vector.py
#basic vector functions
#Mike Daily
#created spring 03
#last modified 7/7/03
########################################################

import math

def makeVector(ptA, ptB):
    """vector(ptA, ptB) - return the mathematical vector from A to B
       Given two points, A and B, this function returns the vector going
       from point A to point B, represented as a list of cartesian
       components.  A and be must be lists of coordinates in x, y, and
       z.  That is:

       vector([x1, y1, z1], [x2, y2, z2]) ==> [x2-x1, y2-y1, z2-z1]
    """
    return [ptB[0]-ptA[0], ptB[1]-ptA[1], ptB[2]-ptA[2]]


def vabs(vec):
    """vabs(vec) - return the absolute value of a given vector
       The absolute value of a three-dimensional vector from point A to
       point B is defined as the distance between points A and B.  This
       function evaluates that value for a vector as defined above.
    """
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def unit(vec):
    #returns the unit vector corresponding to vec
    return mult(vec, 1/vabs(vec))

def mult(vec, n):
    #[x,y,z] -> [n*x, n* y, n*z]
    v2 = []
    for item in vec:
        v2.append(n*item)
    return v2

def add(v1, v2):
    #[x1, y1, z1] + [x2, y2, z2] = [x1 + x2, y1 + y2, z1 + z2]
    v3 = []
    for i in range(0, len(v1)):
        v3.append(v1[i] + v2[i])
    return v3

def dot(v1, v2):
    """dot(v1, v2) - return the inner (dot) product of two vectors
       This function evaluates the dot product between two vectors, defined
       as the following:

       dot([x1, y1, z1], [x2, y2, z2]) ==> x1*x2 + y1*y2 + z1*z2

       It can be shown that this is also equal to:

       dot(v1, v2) ==> |v1| * |v2| * cos(theta)

       Where theta is the angle between the vectors.
    """
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def cross(v1, v2):
    """cross(v1, v2) - returns the cross product of two vectors
       This function evaluates the cross product between two cartesian
       vectors in 3d.  For v1 = [x1, y1, z1], v2 = [x2, y2, z2], and
       cartesian unit vectors i, j, and k, this is defined as the
       determinant of the following matrix:

       | i   j   k  |
       | x1  y1  z1 |
       | x2  y2  z2 |

       More simply,

       cross([x1, y1, z1], [x2, y2, z2]) ==>
         [y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2]

       It can also be shown that

--       |cross(v1, v2)| == |v1| * |v2| * sin(theta)

       Where theta is the angle between the two vectors.
    """
    return [v1[1]*v2[2] - v1[2]*v2[1],
            v1[2]*v2[0] - v1[0]*v2[2],
            v1[0]*v2[1] - v1[1]*v2[0]]

def stp(v1, v2, v3):
    """stp(v1, v2, v3) - return the scalar triple product of three vectors
       This function will return the scalar triple product of three vectors,
       defined most simply as dot(v3, cross(v1, v2)).  It can be shown that
       |stp(v1, v2, v3)| is the volume of the parallelpiped with sides
       defined by the vectors.  It follows that the following relationships
       apply:
 
       stp(v1, v2, v3) == stp(v2, v3, v1) == stp(v3, v1, v2)
    """
    return dot(v3, cross(v1, v2))

def angle(v1, v2):
    """angle(v1, v2) - returns the angle (in radians) between two vectors
       If v1 and v2 are two vectors, this function uses the cosine
       relationship given above (in dot) to determine the angle between
       those vectors when one vector is projected on the plane of the other.
       Note that this value can only be the range of the arccos function,
       so theta will be in the range 0 to Pi.
    """
    return math.acos(dot(v1, v2) / vabs(v1) / vabs(v2))

def torsion(A,B,C,D):
    #calculates the B-C torsion angle of any four atoms A,B,C, and D
    AB = makeVector(A,B)
    BC = makeVector(B,C)
    CD = makeVector(C,D)
    phistp = stp(AB,BC,CD)
    norm1 = cross(AB,BC)
    norm2 = cross(BC,CD)
    phi = toDeg(angle(norm1, norm2))
    if phistp >= 0:
        return phi
    else:
        return -phi

def toDeg(ang):
    """toDeg(ang) - converts an angle in radians to an angle in degrees"""
    return ang / math.pi * 180.0

def toRad(ang):
    """toRad(ang) - converts an angle in degrees to an angle in radians"""
    return ang / 180.0 * math.pi

###########################################################################
