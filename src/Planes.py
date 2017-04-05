import ERRCODE
from math import *
from Position import Position

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

    def deal_X_Y_Plane(self, phoCoord, vecCoord, objCoord, intv):

        distanceBottom = (objCoord.z - phoCoord.z) / vecCoord.x
        distanceUpper = (objCoord.z - phoCoord.z + intv[3]) / vecCoord.z

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 5
        else:
            self.distance = distanceUpper
            self.face = 6

        self.x = phoCoord.x + self.distance * vecCoord.x
        if((self.x >= objCoord.x) and (self.x <= objCoord.x + intv[1])):
            self.y = phoCoord.y + self.distance * vecCoord.y
            if ((self.y >= objCoord.y) and (self.y <= objCoord.y + intv[2])):
                #print("In x - y plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    def deal_Y_Z_Plane(self, phoCoord, vecCoord, objCoord, intv):

        distanceBottom = (objCoord.x - phoCoord.x) / vecCoord.x
        distanceUpper = (objCoord.x - phoCoord.x + intv[1]) / vecCoord.x

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 1
        else:
            self.distance = distanceUpper
            self.face = 2

        self.y =  phoCoord.y + self.distance * vecCoord.y
        if((self.y >= objCoord.y) and (self.y <= objCoord.y + intv[2])):
            self.z =phoCoord.z + self.distance * vecCoord.z
            if ((self.z >= objCoord.z) and (self.z <= objCoord.z + intv[3])):
                #print("In y - z plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    def deal_X_Z_Plane(self, phoCoord, vecCoord, objCoord, intv):

        distanceBottom = (objCoord.y - phoCoord.y) / vecCoord.y
        distanceUpper = (objCoord.y - phoCoord.y + intv[2]) / vecCoord.y

        if (distanceBottom > distanceUpper):
            self.distance = distanceBottom
            self.face = 3
        else:
            self.distance = distanceUpper
            self.face = 4

        self.z = phoCoord.z + self.distance * vecCoord.z
        if ((self.z >= objCoord.z) and (self.z <= objCoord.z + intv[3])):
            self.x = phoCoord.x + self.distance * vecCoord.x
            if ((self.x >= objCoord.x) and (self.x <= objCoord.x + intv[1])):
                #print("In x - z plane")
                return ERRCODE.SUCCESS

        return ERRCODE.CANNOT_FIND

    ###############################################################
    # face 1:x = x1, 2:x = x1 + intv, 3:y = y1,
    #      4:y = y1 + intv, 5:z = z1(ground), 6:z = z1 + intv(sky)
    # (x1, y1, z1) are minimum position of cubic apex
    # intv(1):x0, intv(2):y0, intv(3):z0
    ###############################################################
    def calPlanes(self, phoCoord, vecCoord, objCoord, intv):

        errCode = 0
        MIN_VALUE = 1.0e-6

        self.x = self.y = self.z = 0.0

        # if not parallel
        # check manual P.51-52
        if (abs(vecCoord.x) >= MIN_VALUE):
            errCode = self.deal_Y_Z_Plane(phoCoord, vecCoord, objCoord, intv)
            if (errCode == 0):
                return ERRCODE.SUCCESS

        if (abs(vecCoord.y) >= MIN_VALUE):
            errCode = self.deal_X_Z_Plane(phoCoord, vecCoord, objCoord, intv)
            if (errCode == 0):
                return ERRCODE.SUCCESS

        errCode = self.deal_X_Y_Plane(phoCoord, vecCoord, objCoord, intv)
        if (errCode == ERRCODE.SUCCESS):
            return ERRCODE.SUCCESS

        # if cannot find the intersection due to the limitation of numerical calculation
        # move the photon position slightly then back to the top of this code

        self.x = phoCoord.x + 0.1 * copysign(1.0, (objCoord.x + intv[1]) - phoCoord.x)
        self.y = phoCoord.y + 0.1 * copysign(1.0, (objCoord.y + intv[2]) - phoCoord.y)
        self.z = phoCoord.z + 0.1 * copysign(1.0, (objCoord.z + intv[3]) - phoCoord.z)

        print("can't find cube intersection")
        return ERRCODE.CANNOT_FIND
