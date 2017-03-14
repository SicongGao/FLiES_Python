import ERRCODE
from math import *
from Position import Position

class TreeBoundary:

    face = -1
    io = -1
    distance = -1
    t1 = 1.0e5
    t2 = 1.0e5
    t3 = 1.0e5
    t4 = 1.0e5

    def __init__(self):
        self.face = -1
        self.io = -1
        self.distance = -1
        self.t1 = 1.0e5
        self.t2 = 1.0e5
        self.t3 = 1.0e5
        self.t4 = 1.0e5

    ############################################################
    # This routine calculates the distance between the (x,y,z)
    # and cone external boundary
    ############################################################

    def dealTreeType(self, treeType, phoCoord, vectCoord, tobj):

        # init
        self.__init__()

        if (treeType == 1):
            self.cones(phoCoord, vectCoord, tobj)

        elif (treeType == 2) or (treeType == 4):
            self.cyls(phoCoord, vectCoord, tobj)

        elif (treeType == 3):
            self.elpss(phoCoord, vectCoord, tobj)

        elif (treeType == 5):
            self.helpss(phoCoord, vectCoord, tobj)

    def cones(self, phoCoord, vectCoord, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z), phoCoord
        # photon direction (ux, uy, uz), vectCoord
        # cone apex(cx,cy,cz),crownCoord
        # cone height h and radius bottom circle r
        # rp=(r/h)
        # face: the face of cylinder 1=side;2=bottom circle
        # d = quadratic judgements
        # t1,t2,t3 = distance

        # we have to solve following eq.
        # a*t**2+b*t+c=0
        # for cone boundary

        a = b = c = 0.0

        # set initial value
        t = 0.0

        # crown related parameters
        crownCoord = Position()
        crownCoord.setPosition(tobj[1], tobj[2], tobj[3])

        h = tobj[4]
        r = tobj[5]
        rp = r / h

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = vectCoord.x ** 2 + vectCoord.y ** 2 - (rp * vectCoord.z) ** 2
        a = copysign(max(abs(a), 0.0001), a)

        b = 2.0 * (vectCoord.x * (phoCoord.x - crownCoord.x) + vectCoord.y * (phoCoord.y - crownCoord.y) -
                   vectCoord.z * (phoCoord.z - crownCoord.z) * rp ** 2)
        c = (phoCoord.x - crownCoord.x) ** 2 + (phoCoord.y - crownCoord.y) ** 2 - \
            (rp * (phoCoord.z - crownCoord.z)) ** 2

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            # if t do not meet the following conditions,
            # penalties are added in the distance.

            if (((crownCoord.z - h - 1.0e-4) <= (phoCoord.z + self.t1 * vectCoord.z)) and
                    ((phoCoord.z + self.t1 * vectCoord.z) <= (crownCoord.z + 1.0e-4))):

                if (self.t1 < 0):
                    self.t1 = 1.0e5
            else:
                self.t1 = 1.0e5

            if (((crownCoord.z - h - 1.0e-4) <= (phoCoord.z + self.t2 * vectCoord.z)) and
                    ((phoCoord.z + self.t2 * vectCoord.z) <= (crownCoord.z + 1.0e-4))):

                if (self.t2 < 0):
                    self.t2 = 1.0e5
            else:
                self.t2 = 1.0e5

        # if ray - line cross the bottom circle of the cone crown.
        vectCoord.z = copysign(max(abs(vectCoord.z), 1.0e-4), vectCoord.z)
        self.t3 = (crownCoord.z - h - phoCoord.z) / vectCoord.z

        circle = (phoCoord.x + self.t3 * vectCoord.x - crownCoord.x) ** 2 + \
                 (phoCoord.y + self.t3 * vectCoord.y - crownCoord.y) ** 2

        if ((self.t3 < 0) or (circle > r ** 2)):
            self.t3 = 1.0e5

        # determination of the sb
        t = self.t1
        self.face = 1

        if (t >= self.t2):
            t = self.t2

        if (t >= self.t3):
            t = self.t3
            self.face = 2

        # there is no distance, if t = 1.0e5
        self.distance = t
        if (self.distance > 0.9e5):
            self.face = -1

        # inout check
        if (((crownCoord.z - h) < phoCoord.z) and (phoCoord.z < crownCoord.z)):
            if ((phoCoord.x - crownCoord.x) ** 2 + (phoCoord.y - crownCoord.y) ** 2 <
                        (rp * (crownCoord.z - phoCoord.z)) ** 2):
                self.io = 0
            else:
                self.io = 1
        else:
            self.io = 1

        # prevent numerical error
        if (self.distance > 1.0e5):
            self.io = 1

        return ERRCODE.SUCCESS

    def cyls(self, phoCoord, vectCoord, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z), phoCoord, phoCoord
        # photon direction (ux, uy, uz), vectCoord
        # cone apex(cx,cy,cz),crownCoord
        # cone height h and radius bottom circle r
        # face: the face of cylinder 1=side;2=bottom circle; 3=upper circle
        # d = quadratic judgements
        # t1,t2,t3 = distance

        # we have to solve following eq.
        # a*t**2+b*t+c=0
        # for cone boundary

        a = b = c = 0.0

        # set initial value
        t = 0.0

        # crown related parameters
        crownCoord = Position()
        crownCoord.setPosition(tobj[1], tobj[2], tobj[3])
        h = tobj[4]
        r = tobj[5]

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = vectCoord.x ** 2 + vectCoord.y ** 2
        a = copysign(max(abs(a), 0.0001), a)

        b = 2.0 * (vectCoord.x * (phoCoord.x - crownCoord.x) + vectCoord.y * (phoCoord.y - crownCoord.y))
        c = (phoCoord.x - crownCoord.x) ** 2 + (phoCoord.y - crownCoord.y) ** 2 - r ** 2

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            # if t do not meet the following conditions,
            # penalties are added in the distance.

            if (((crownCoord.z - h - 1.0e-4) <= (phoCoord.z + self.t1 * vectCoord.z)) and
                    ((phoCoord.z + self.t1 * vectCoord.z) <= (crownCoord.z + 1.0e-4))):

                if (self.t1 < 0):
                    self.t1 = 1.0e5
            else:
                self.t1 = 1.0e5

            if (((crownCoord.z - h - 1.0e-4) <= (phoCoord.z + self.t2 * vectCoord.z)) and
                    ((phoCoord.z + self.t2 * vectCoord.z) <= (crownCoord.z + 1.0e-4))):

                if (self.t2 < 0):
                    self.t2 = 1.0e5
            else:
                self.t2 = 1.0e5

        # if ray - line cross the bottom circle of the cylinder crown.
        vectCoord.z = copysign(max(abs(vectCoord.z), 1.0e-4), vectCoord.z)
        self.t3 = (crownCoord.z - h - phoCoord.z) / vectCoord.z

        circle = (phoCoord.x + self.t3 * vectCoord.x - crownCoord.x) ** 2 + \
                 (phoCoord.y + self.t3 * vectCoord.y - crownCoord.y) ** 2

        if ((self.t3 < 0) or (circle > r ** 2)):
            self.t3 = 1.0e5

        # if ray - line cross the upper circle of the cylinder crown.
        vectCoord.z = copysign(max(abs(vectCoord.z), 1.0e-4), vectCoord.z)
        self.t4 = (crownCoord.z - phoCoord.z) / vectCoord.z

        circle = (phoCoord.x + self.t4 * vectCoord.x - crownCoord.x) ** 2 + \
                 (phoCoord.y + self.t4 * vectCoord.y - crownCoord.y) ** 2

        if ((self.t4 < 0) or (circle > r ** 2)):
            self.t4 = 1.0e5

        # determination of the s
        t = self.t1
        self.face = 1

        if (t >= self.t2):
            t = self.t2

        if (t >= self.t3):
            t = self.t3
            self.face = 2

        if (t >= self.t4):
            t = self.t4
            self.face = 3

        # there is no distance, t = 1.0e5
        self.distance = t
        if (self.distance > 0.9e5):
            self.face = -1

        # inout check
        if ((((crownCoord.z - h) < phoCoord.z) and (phoCoord.z < crownCoord.z)) and
            ((phoCoord.x - crownCoord.x) ** 2 + (phoCoord.y - crownCoord.y) ** 2 < r ** 2)):
            self.io = 0
        else:
            self.io = 1

        # prevent numerical error
        if (self.distance > 1.0e5):
            self.io = 1

        return ERRCODE.SUCCESS

    def elpss(self, phoCoord, vectCoord, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z), phoCoord
        # photon direction (ux, uy, uz), vectCoord
        # cone apex(cx,cy,cz),crownCoord
        # radius for z-axis : r1
        # radius for x - y plane: r2
        # d = quadratic judgements
        # t1,t2,t3 = distance

        # we have to solve following eq.
        # a*t**2+b*t+c=0
        # for rotational elliptical sphere external boundary

        a = b = c = 0.0

        # set initial value
        t = 0.0

        # crown related parameters
        crownCoord = Position()
        crownCoord.setPosition(tobj[1], tobj[2], tobj[3])
        r1 = tobj[4]
        r2 = tobj[5]
        r1Square = r1 ** 2
        r2Square = r2 ** 2

        self.distance = 0.0
        elps = 0.0

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = r1Square * vectCoord.x ** 2
        a += r1Square * vectCoord.y ** 2
        a += r2Square * vectCoord.z ** 2
        a = max(a, 0.0001)

        b = r1Square * vectCoord.x * (phoCoord.x - crownCoord.x)
        b += r1Square * vectCoord.y * (phoCoord.y - crownCoord.y)
        b += r2Square * vectCoord.z * (phoCoord.z - crownCoord.z)
        b *= 2.0

        c = r1Square * (phoCoord.x - crownCoord.x) ** 2
        c += r1Square * (phoCoord.y - crownCoord.y) ** 2
        c += r2Square * (phoCoord.z - crownCoord.z) ** 2
        c -= r1Square * r2Square

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            if (self.t1 < 0):
                self.t1 = 1.0e5
            if (self.t2 < 0):
                self.t2 = 1.0e5

        self.distance = min(self.t1, self.t2)

        elps = (phoCoord.x - crownCoord.x) ** 2 / r2Square
        elps += (phoCoord.y - crownCoord.y) ** 2 / r2Square
        elps += (phoCoord.z - crownCoord.z) ** 2 / r1Square

        if ((elps < 1.0) and (self.distance < 1.0e5)):
            self.io = 0
        else:
            self.io = 1

        return ERRCODE.SUCCESS

    def helpss(self, phoCoord, vectCoord, tobj):

        # semi-spheroid
        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z), phoCoord
        # photon direction (ux, uy, uz), vectCoord
        # cone apex(cx,cy,cz),crownCoord
        # radius for z-axis : r1
        # radius for x - y plane: r2
        # d = quadratic judgements
        # t1,t2,t3 = distance

        # we have to solve following eq.
        # a*t**2+b*t+c=0
        # for rotational elliptical sphere external boundary
        # face: the face of half elpsd 1=upper elpsd;2=bottom circle

        a = b = c = 0.0

        # set initial value
        t = 0.0
        z1 = z2 = 0.0
        r1 = r2 = 0.0

        # crown related parameters
        crownCoord = Position()
        crownCoord.setPosition(tobj[1], tobj[2], tobj[3])
        r1 = tobj[4]
        r2 = tobj[5]
        r1Square = r1 ** 2
        r2Square = r2 ** 2

        self.distance = 0.0
        elps = 0.0

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = r1Square * vectCoord.x ** 2
        a += r1Square * vectCoord.y ** 2
        a += r2Square * vectCoord.z ** 2
        a = max(a, 0.0001)

        b = r1Square * vectCoord.x * (phoCoord.x - crownCoord.x)
        b += r1Square * vectCoord.y * (phoCoord.y - crownCoord.y)
        b += r2Square * vectCoord.z * (phoCoord.z - crownCoord.z)
        b *= 2.0

        c = r1Square * (phoCoord.x - crownCoord.x) ** 2
        c += r1Square * (phoCoord.y - crownCoord.y) ** 2
        c += r2Square * (phoCoord.z - crownCoord.z) ** 2
        c -= r1Square * r2Square

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            if (self.t1 < 0):
                self.t1 = 1.0e5
            if (self.t2 < 0):
                self.t2 = 1.0e5

            z1 = phoCoord.z + self.t1 * vectCoord.z
            z2 = phoCoord.z + self.t2 * vectCoord.z

            if (self.t1 < 0) or (z1 < crownCoord.z):
                self.t1 = 1.0e5
            if (self.t2 < 0) or (z2 < crownCoord.z):
                self.t2 = 1.0e5

        # if ray-line cross the bottom circle of the cone crown
        vectCoord.z = copysign(max(max(abs(vectCoord.z), 1.0e-4)), vectCoord.z)
        self.t3 = (crownCoord.z - phoCoord.z) / vectCoord.z

        circle = (phoCoord.x + self.t3 * vectCoord.x - crownCoord.x) ** 2 + \
                 (phoCoord.y + self.t3 * vectCoord.y - crownCoord.y) ** 2

        if (self.t3 < 0) or (circle > r2 * r2):
            self.t3 = 1.0e5

        # determination of the s
        t = self.t1
        self.face = 1

        if (t > self.t2):
            t = self.t2

        if (t > self.t3):
            t = self.t3
            self.face = 2

        self.distance = t
        if (self.distance > 0.9e5):
            self.face = -1

        # in out check
        elps = (phoCoord.x - crownCoord.x) ** 2 / r2Square
        elps += (phoCoord.y - crownCoord.y) ** 2 / r2Square
        elps += (phoCoord.z - crownCoord.z) ** 2 / r1Square

        if ((elps < 1.0) and (self.distance < 1.0e5)):
            self.io = 0
        else:
            self.io = 1

        return ERRCODE.SUCCESS
