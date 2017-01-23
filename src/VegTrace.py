import ERRCODE
import common as comm
from math import *
from Planes import Planes
from TreeBoundary import TreeBoundary
from G_Function import G_Function

# #####################################################
# simulate the optical thickness in the canopy
# along the photon trajectory
#
# Written by H. Kobayashi
# Last modified 08/04/02
# #####################################################


class VegTrace:

    sFlag = 1
    tau = 0.0

    def trace(self, x, y, z, uxr, uyr, uzr):

        # marginal value
        mgn = 1.0e-2
        iNObj = 1

        # weight minimum limit
        conv = 1.0e-8

        # initial distance from object
        distanceObj = 1.0e5
        distancePho = 1.0e5

        # sflg: stem flag = 0 stem collision, 1 = no stem
        self.sFlag = 1

        # optical thickness
        self.tau = tauc12 = tauc = taub = 0.0

        x0 = y0 = z0 = 0.0

        x0 = x
        y0 = y
        z0 = max(z, mgn)

        oFace = pFace = 0
        rb12 = 1.0
        io = io2 = io12 = 0
        distance2 = distance12 = 0.0

        intv = [50.0 / comm.RES] * 4

        # z boundary for trance flag = 1 (zlim = 0.0), flag = 2 (zlim = zmax)
        zlim = [0.0, 0.0, comm.Z_MAX]
        flag = int(1.5 + copysign(0.5, uzr))

        tobj = [0.0] * 6
        tobjb = [0.0] * 6
        tobj12 = [0.0] * 6

        planes = Planes()
        treeBoundary = TreeBoundary()
        gFunction = G_Function()
        gFunction.igtbl()

        # do while photon reaches the terminal point
        while(1):
            x1 = trunc(x0 / intv[1])
            y1 = trunc(y0 / intv[2])
            z1 = trunc(z0 / intv[3])

            iVOX = comm.IX_MAX * comm.IY_MAX * int(z1)
            iVOX += int(y1) * comm.IY_MAX
            iVOX += int(x1) + 1

            x1 *= intv[1]
            y1 *= intv[2]
            z1 *= intv[3]

            io = 1
            iNobj = -1
            index = -1

            # check the photon intersection with big-voxel walls
            errCode = planes.calPlanes(x0, y0, z0, uxr, uyr, uzr, x1, y1, z1, intv)
            distancePho = planes.distance
            if (errCode == ERRCODE.CANNOT_FIND):
                # update the x0, y0, z0
                x0 = planes.x
                y0 = planes.y
                z0 = planes.z

            # check the photon intersection with objects
            if (comm.N_DIVS[iVOX] != 0):

                distanceObj = 1.0e5
                for idiv in range(1, comm.N_DIVS[iVOX]):

                    index = comm.DIVS[iVOX][idiv]

                    for l in range(1, 6):
                        tobj[l] = comm.OBJ[index][l-1]

                    treeBoundary.dealTreeType(comm.S_OBJ[index], x0, y0, z0, uxr, uyr, uzr, tobj)

                    tempDistance = treeBoundary.distance
                    tempIO = treeBoundary.io

                    if (tempDistance < distanceObj):
                        distanceObj = tempDistance
                        iNobj = index
                        io = tempIO

                    if (io == 0):
                        break

                # if stem collision, return
                if (comm.S_OBJ[iNobj] == 4) and (distanceObj < 1.0e5):
                    self.sFlag = 0
                    return ERRCODE.FAILURE

            # increment of the optical path
            if (io == 0):
                # check branch optical thickness
                tobjb[1] = tobj[1]
                tobjb[2] = tobj[2]
                tobjb[3] = tobj[3] - tobj[4] * (1.0 - comm.RB) * min(1, abs(comm.S_OBJ[iNObj] - 5))
                tobjb[4] = tobj[4] * comm.RB
                tobjb[5] = tobj[5] * comm.RB

                tobj12[1] = tobj[1]
                tobj12[2] = tobj[2]
                tobj12[3] = tobj[3] - tobj[4] * (1.0 - comm.RB) * min(1, abs(comm.S_OBJ[iNObj] - 5))
                tobj12[4] = tobj[4] * rb12
                tobj12[5] = tobj[5] * rb12

                treeBoundary.dealTreeType(comm.S_OBJ[index], x0, y0, z0, uxr, uyr, uzr, tobj)
                io2 = treeBoundary.io
                distance2 = treeBoundary.distance

                treeBoundary.dealTreeType(comm.S_OBJ[index], x0, y0, z0, uxr, uyr, uzr, tobj12)
                io12 = treeBoundary.io
                distance12 = treeBoundary.distance

                # if outside of branch (io12 = 1) go into the branch
                if (io12 == 1):
                    xb = x0 + (distance12 + mgn) * uxr
                    yb = y0 + (distance12 + mgn) * uyr
                    zb = z0 + (distance12 + mgn) * uzr

                    treeBoundary.dealTreeType(comm.S_OBJ[index], xb, yb, zb, uxr, uyr, uzr, tobj12)
                    distance12 = treeBoundary.distance

                # if outside of branch (io2 = 1) go into the branch
                if (io2 == 1):
                    xb = x0 + (distance2 + mgn) * uxr
                    yb = y0 + (distance2 + mgn) * uyr
                    zb = z0 + (distance2 + mgn) * uzr

                    treeBoundary.dealTreeType(comm.S_OBJ[index], xb, yb, zb, uxr, uyr, uzr, tobj12)
                    distance2 = treeBoundary.distance

                ################################
                # calculation of optical path
                ################################
                # Don't know what's the meaning
                ################################
                th = acos(uzr)
                ith = int(radians(th))
                rio = 1.0 - float(io2)
                rio12 = 1.0 - float(io12)
                cf = comm.S_BAR[comm.I_OBJ[index]]
                cf12 = comm.S_BAR[comm.I_OBJ[index]]

                taub = rio * (distance2 + mgn) * comm.BAD[comm.I_OBJ[index]] * gFunction.GT_BLB[ith] * comm.BP2
                taub += rio * (distance2 + mgn) * comm.U[comm.I_OBJ[index]] * gFunction.GT_BLC[ith] * 4.0 * cf * (1.0 - comm.BP2)

                tauc12 = ((distance12 + mgn) - (distance12 + mgn) * rio) * rio12
                tauc12 = comm.U[comm.I_OBJ[index]] * gFunction.GT_BLC[ith] * 4.0 * cf12 * (1.0 - comm.BP1) \
                         + (comm.BAD[comm.I_OBJ[index]] * gFunction.GT_BLB[ith]) * comm.BP1

                tauc = (distanceObj + mgn) - (distance12 + mgn) * rio12
                tauc = comm.U[comm.I_OBJ[index]] * gFunction.GT_BLC[ith] * 4.0 * cf * (1.0 - comm.BP1) \
                         + (comm.BAD[comm.I_OBJ[index]] * gFunction.GT_BLB[ith]) * comm.BP1

                self.tau += tauc + tauc12 + taub

            # refresh photon position
            d = distanceObj * (1.0 - float(io)) + min(distanceObj, distancePho) * float(io)
            x0 += (d + mgn) * uxr
            y0 += (d + mgn) * uyr
            z0 += (d + mgn) * uzr

            # check upper or bottom boundary condition
            if (copysign(1.0, uzr) * (z0 - zlim[flag]) >= 0.0):
                return ERRCODE.SUCCESS

            x0 -= (trunc(x0 / comm.X_MAX) - 0.5 + copysign(0.5, x0)) * comm.X_MAX
            y0 -= (trunc(y0 / comm.Y_MAX) - 0.5 + copysign(0.5, y0)) * comm.Y_MAX

        return ERRCODE.SUCCESS
