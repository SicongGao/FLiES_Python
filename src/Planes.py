import ERRCODE
from math import *


# ***********************************************************
# This routine calculates the distance between the (x,y,x)
# and cube
# ***********************************************************


class Planes:

    face = -1
    x = 0.0
    y = 0.0
    z = 0.0
    distance = 0.0

    def deal_X_Y_Plane(self, x0, y0, z0, ux, uy, uz, x1, y1, z1, intv):

        distanceBottom = (z1 - z0) / uz
        distanceUpper = (z1 - z0 + intv[3]) / uz

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 5
        else:
            self.distance = distanceUpper
            self.face = 6

        self.x = x0 + self.distance * ux
        if((self.x >= x1) and (self.x <= x1 + intv[1])):
            self.y = y0 + self.distance * uy
            if ((self.y >= y1) and (self.y <= y1 + intv[2])):
                print("In x - y plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    def deal_Y_Z_Plane(self, x0, y0, z0, ux, uy, uz, x1, y1, z1, intv):

        distanceBottom = (x1 - x0) / ux
        distanceUpper = (x1 - x0 + intv[1]) / ux

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 1
        else:
            self.distance = distanceUpper
            self.face = 2

        self.y = y0 + self.distance * uy
        if((self.y >= y1) and (self.y <= y1 + intv[2])):
            self.z = z0 + self.distance * uz
            if ((self.z >= z1) and (self.z <= z1 + intv[3])):
                print("In y - z plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    def deal_X_Z_Plane(self, x0, y0, z0, ux, uy, uz, x1, y1, z1, intv):

        distanceBottom = (y1 - y0) / uy
        distanceUpper = (y1 - y0 + intv[2]) / uy

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 3
        else:
            self.distance = distanceUpper
            self.face = 4

        self.z = z0 + self.distance * uz
        if ((self.z >= z1) and (self.z <= z1 + intv[3])):
            self.x = x0 + self.distance * ux
            if ((self.x >= x1) and (self.x <= x1 + intv[1])):
                print("In x - z plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    ###############################################################
    # face 1:x = x1, 2:x = x1 + intv, 3:y = y1,
    #      4:y = y1 + intv, 5:z = z1(ground), 6:z = z1 + intv(sky)
    # (x1, y1, z1) are minimum position of cubic apex
    # intv(1):x0, intv(2):y0, intv(3):z0
    ###############################################################
    def calPlanes(self, d, x0, y0, z0, ux, uy, uz, x1, y1, z1, intv):

        errCode = 0
        MIN_VALUE = 1.0e-6
        self.distance = d

        self.x = self.y = self.z = 0.0

        # if not parallel
        # check manual P.51-52
        if (abs(ux) >= MIN_VALUE):
            errCode = self.deal_Y_Z_Plane(x0, y0, z0, ux, uy, uz, x1, y1, z1, intv)
            if (errCode == 0):
                return ERRCODE.SUCCESS

        if (abs(uy) >= MIN_VALUE):
            errCode = self.deal_X_Z_Plane(x0, y0, z0, ux, uy, uz, x1, y1, z1, intv)
            if (errCode == 0):
                return ERRCODE.SUCCESS

        errCode = self.deal_X_Y_Plane(x0, y0, z0, ux, uy, uz, x1, y1, z1, intv)
        if (errCode == 0):
            return ERRCODE.SUCCESS

        # if cannot find the intersection due to the limitation of numerical calculation
        # move the photon position slightly then back to the top of this code

        self.x = x0 + 0.1 * copysign(1.0, (x1 + intv[1]) - x0)
        self.y = y0 + 0.1 * copysign(1.0, (y1 + intv[2]) - y0)
        self.z = z0 + 0.1 * copysign(1.0, (z1 + intv[3]) - z0)

        print("can't find cube intersection")
        return ERRCODE.CANNOT_FIND
