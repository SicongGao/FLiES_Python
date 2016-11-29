import ERRCODE
from math import *
import common as comm


class TreeBoundary:

    ############################################################
    # This routine calculates the distance between the (x,y,z)
    # and cone external boudary
    ############################################################

    def cones(self, s, x, y, z, ux, uy, uz, tobj, face, io):

        # Parameter Definition
        # inout : the position inside object, then =0, outside=1
        # photon position (x,y,z)
        # photon direction (ux, uy, uz)
        # cone apex(cx,cy,cz)
        # cone height h and radius bottom circle r
        # rp=(r/h)
        # face: the face of cylinder 1=side;2=bottom circle
        # d = quadratic judgements
        # t1,t2,t3 = solution

        # we have to solve following eq.
        # a*t**2+b*t+c=0
        # for cone boundary

        a = b = c = 0.0

        # set initial value
        t = 0.0
        t1 = t2 = t3 = 1.0e5

        # crown related parameters
        cx = tobj[1]
        cy = tobj[2]
        cz = tobj[3]
        h = tobj[4]
        r = tobj[5]
        rp = r / h

        # calculate the quadratic parameters
        # a*t^2+b*t+c=0
        a = ux ** 2 * uy ** 2 - (rp * uz) ** 2
        a = copysign(max(abs(a), 0.0001), a)

        b = 2.0 * (ux * (x - cx) + uy * (y - cy) - uz * (z - cz) * rp ** 2)
        c = (x - cx) ** 2 + (y - cy) ** 2 - (rp * (z - cz)) ** 2

        d = b ** 2 - 4.0 * a * c

        if (d >= 0):
            t1 = (-b - sqrt(d)) / (2.0 * a)
            t2 = (-b + sqrt(d)) / (2.0 * a)

            # if t do not meet the following conditions,
            # penalties are added in the solution.

            if (((cz - h - 1.0e-4) <= (z + t1 * uz))
                and ((z + t1 * uz) <= (cz + 1.0e-4))):

                if (t1 < 0):
                    t1 = 1.0e5
            else:
                t1 = 1.0e5

            if (((cz - h - 1.0e-4) <= (z + t2 * uz))
                and ((z + t2 * uz) <= (cz + 1.0e-4))):

                if (t2 < 0):
                    t2 = 1.0e5
            else:
                t2 = 1.0e5

        # if ray - line cross the bottom circle of the cone crown.
        uz = copysign(max(abs(uz), 1.0e-4), uz)
        t3 = (cz - h - z) / uz

        circle = (x + t3 * ux - cx) ** 2 +(y + t3 * uy - cy) ** 2

        if ((t3 < 0) or (circle > r ** 2)):
            t3 = 1.0e5

        # determination of the sb
        t = t1
        face = 1

        if (t >= t2):
            t = t2

        if (t >= t3):
            t = t3
            face = 2

        # there is no solution, t = 1.0e5
        s = t
        if (s > 0.9e5):
            face = -1

        # inout check
        if (((cz - h) < z) and (z < cz)):
            if ( (x - cx) ** 2 + (y - cy) ** 2 < (rp * (cz - z)) ** 2):
                io = 0
            else:
                io = 1
        else:
            io = 1

        # prevent numerical error
        if (s > 1.0e5):
            io = 1

        return ERRCODE.SUCCESS

    def cyls(self, s, x, y, z, ux, uy, uz, tobj, face, io):
        return ERRCODE.SUCCESS

    def elpss(self, s, x, y, z, ux, uy, uz, tobj, face, io):
        return ERRCODE.SUCCESS

    def helpss(self, s, x, y, z, ux, uy, uz, tobj, face, io):
        return ERRCODE.SUCCESS
