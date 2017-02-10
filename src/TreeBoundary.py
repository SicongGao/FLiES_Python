import ERRCODE
from math import *


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

    def dealTreeType(self, treeType, x, y, z, ux, uy, uz, tobj):
        if (treeType == 1):
            self.cones(x, y, z, ux, uy, uz, tobj)

        elif (treeType == 2) or (treeType == 4):
            self.cyls(x, y, z, ux, uy, uz, tobj)

        elif (treeType == 3):
            self.elpss(x, y, z, ux, uy, uz, tobj)

        elif (treeType == 5):
            self.helpss(x, y, z, ux, uy, uz, tobj)

    def cones(self, x, y, z, ux, uy, uz, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z)
        # photon direction (ux, uy, uz)
        # cone apex(cx,cy,cz)
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
        cx = tobj[1]
        cy = tobj[2]
        cz = tobj[3]
        h = tobj[4]
        r = tobj[5]
        rp = r / h

        # init
        self.__init__()

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = ux ** 2 + uy ** 2 - (rp * uz) ** 2
        a = copysign(max(abs(a), 0.0001), a)

        b = 2.0 * (ux * (x - cx) + uy * (y - cy) - uz * (z - cz) * rp ** 2)
        c = (x - cx) ** 2 + (y - cy) ** 2 - (rp * (z - cz)) ** 2

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            # if t do not meet the following conditions,
            # penalties are added in the distance.

            if (((cz - h - 1.0e-4) <= (z + self.t1 * uz)) and
                    ((z + self.t1 * uz) <= (cz + 1.0e-4))):

                if (self.t1 < 0):
                    self.t1 = 1.0e5
            else:
                self.t1 = 1.0e5

            if (((cz - h - 1.0e-4) <= (z + self.t2 * uz)) and
                    ((z + self.t2 * uz) <= (cz + 1.0e-4))):

                if (self.t2 < 0):
                    self.t2 = 1.0e5
            else:
                self.t2 = 1.0e5

        # if ray - line cross the bottom circle of the cone crown.
        uz = copysign(max(abs(uz), 1.0e-4), uz)
        self.t3 = (cz - h - z) / uz

        circle = (x + self.t3 * ux - cx) ** 2 + (y + self.t3 * uy - cy) ** 2

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
        if (((cz - h) < z) and (z < cz)):
            if ((x - cx) ** 2 + (y - cy) ** 2 < (rp * (cz - z)) ** 2):
                self.io = 0
            else:
                self.io = 1
        else:
            self.io = 1

        # prevent numerical error
        if (self.distance > 1.0e5):
            self.io = 1

        return ERRCODE.SUCCESS

    def cyls(self, x, y, z, ux, uy, uz, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z)
        # photon direction (ux, uy, uz)
        # cone apex(cx,cy,cz)
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
        cx = tobj[1]
        cy = tobj[2]
        cz = tobj[3]
        h = tobj[4]
        r = tobj[5]

        # init
        self.__init__()

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = ux ** 2 + uy ** 2
        a = copysign(max(abs(a), 0.0001), a)

        b = 2.0 * (ux * (x - cx) + uy * (y - cy))
        c = (x - cx) ** 2 + (y - cy) ** 2 - r ** 2

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            # if t do not meet the following conditions,
            # penalties are added in the distance.

            if (((cz - h - 1.0e-4) <= (z + self.t1 * uz)) and
                    ((z + self.t1 * uz) <= (cz + 1.0e-4))):

                if (self.t1 < 0):
                    self.t1 = 1.0e5
            else:
                self.t1 = 1.0e5

            if (((cz - h - 1.0e-4) <= (z + self.t2 * uz)) and
                    ((z + self.t2 * uz) <= (cz + 1.0e-4))):

                if (self.t2 < 0):
                    self.t2 = 1.0e5
            else:
                self.t2 = 1.0e5

        # if ray - line cross the bottom circle of the cylinder crown.
        uz = copysign(max(abs(uz), 1.0e-4), uz)
        self.t3 = (cz - h - z) / uz

        circle = (x + self.t3 * ux - cx) ** 2 + (y + self.t3 * uy - cy) ** 2

        if ((self.t3 < 0) or (circle > r ** 2)):
            self.t3 = 1.0e5

        # if ray - line cross the upper circle of the cylinder crown.
        uz = copysign(max(abs(uz), 1.0e-4), uz)
        self.t4 = (cz - z) / uz

        circle = (x + self.t4 * ux - cx) ** 2 + (y + self.t4 * uy - cy) ** 2

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
        if ((((cz - h) < z) and (z < cz)) and
            ((x - cx) ** 2 + (y - cy) ** 2 < r ** 2)):
            self.io = 0
        else:
            self.io = 1

        # prevent numerical error
        if (self.distance > 1.0e5):
            self.io = 1

        return ERRCODE.SUCCESS

    def elpss(self, x, y, z, ux, uy, uz, tobj):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z)
        # photon direction (ux, uy, uz)
        # cone apex(cx,cy,cz)
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
        cx = tobj[1]
        cy = tobj[2]
        cz = tobj[3]
        r1 = tobj[4]
        r2 = tobj[5]
        r1Square = r1 ** 2
        r2Square = r2 ** 2

        self.distance = 0.0
        elps = 0.0

        # init
        self.__init__()

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = r1Square * ux ** 2
        a += r1Square * uy ** 2
        a += r2Square * uz ** 2
        a = max(a, 0.0001)

        b = r1Square * ux * (x - cx)
        b += r1Square * uy * (y - cy)
        b += r2Square * uz * (z - cz)
        b *= 2.0

        c = r1Square * (x - cx) ** 2
        c += r1Square * (y - cy) ** 2
        c += r2Square * (z - cz) ** 2
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

        elps = (x - cx) ** 2 / r2Square
        elps += (y - cy) ** 2 / r2Square
        elps += (z - cz) ** 2 / r1Square

        if ((elps < 1.0) and (self.distance < 1.0e5)):
            self.io = 0
        else:
            self.io = 1

        return ERRCODE.SUCCESS

    def helpss(self, x, y, z, ux, uy, uz, tobj):

        # semi-spheroid
        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z)
        # photon direction (ux, uy, uz)
        # cone apex(cx,cy,cz)
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
        cx = tobj[1]
        cy = tobj[2]
        cz = tobj[3]
        r1 = tobj[4]
        r2 = tobj[5]
        r1Square = r1 ** 2
        r2Square = r2 ** 2

        self.distance = 0.0
        elps = 0.0

        # init
        self.__init__()

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = r1Square * ux ** 2
        a += r1Square * uy ** 2
        a += r2Square * uz ** 2
        a = max(a, 0.0001)

        b = r1Square * ux * (x - cx)
        b += r1Square * uy * (y - cy)
        b += r2Square * uz * (z - cz)
        b *= 2.0

        c = r1Square * (x - cx) ** 2
        c += r1Square * (y - cy) ** 2
        c += r2Square * (z - cz) ** 2
        c -= r1Square * r2Square

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            self.t1 = (-b - sqrt(d)) / (2.0 * a)
            self.t2 = (-b + sqrt(d)) / (2.0 * a)

            if (self.t1 < 0):
                self.t1 = 1.0e5
            if (self.t2 < 0):
                self.t2 = 1.0e5

            z1 = z + self.t1 * uz
            z2 = z + self.t2 * uz

            if (self.t1 < 0) or (z1 < cz):
                self.t1 = 1.0e5
            if (self.t2 < 0) or (z2 < cz):
                self.t2 = 1.0e5

        # if ray-line cross the bottom circle of the cone crown
        uz = copysign(max(max(abs(uz), 1.0e-4)), uz)
        self.t3 = (cz - z) / uz

        circle = (x + self.t3 * ux - cx) ** 2 + (y + self.t3 * uy - cy) ** 2

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
        elps = (x - cx) ** 2 / r2Square
        elps += (y - cy) ** 2 / r2Square
        elps += (z - cz) ** 2 / r1Square

        if ((elps < 1.0) and (self.distance < 1.0e5)):
            self.io = 0
        else:
            self.io = 1

        return ERRCODE.SUCCESS
