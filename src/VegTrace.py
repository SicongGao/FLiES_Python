import ERRCODE
import common as comm
from math import *
from Planes import Planes
from TreeBoundary import TreeBoundary
from G_Function import G_Function
from Position import Position

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

    def trace(self, coord, vectCoord):

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

        # x0 = y0 = z0 = 0.0
        #
        # x0 = x
        # y0 = y
        # z0 = max(z, mgn)
        objCoord = Position()
        phoCoord = Position()
        phoCoord.setPosition(coord.x, coord.y, max(coord.z, mgn))

        oFace = pFace = 0
        rb12 = 1.0
        io = io2 = io12 = 0
        distance2 = distance12 = 0.0

        intv = [50.0 / comm.RES] * 4

        # z boundary for trance flag = 1 (zlim = 0.0), flag = 2 (zlim = zmax)
        zlim = [0.0, 0.0, comm.Z_MAX]
        flag = int(1.5 + copysign(0.5, vectCoord.z))

        tobj = [0.0] * 6
        tobjb = [0.0] * 6
        tobj12 = [0.0] * 6

        planes = Planes()
        treeBoundary = TreeBoundary()
        gFunction = G_Function()
        gFunction.igtbl()

        print("Vegetation trace start...")

        # do while photon reaches the terminal point
        while(1):
            # x1 = trunc(x0 / intv[1])
            # y1 = trunc(y0 / intv[2])
            # z1 = trunc(z0 / intv[3])
            objCoord.setPosition(trunc(phoCoord.x / intv[1]), trunc(phoCoord.y / intv[2]), trunc(phoCoord.z / intv[3]))

            iVOX = comm.IX_MAX * comm.IY_MAX * int(objCoord.z)
            iVOX += int(objCoord.y) * comm.IY_MAX
            iVOX += int(objCoord.x) + 1

            objCoord.x *= intv[1]
            objCoord.y *= intv[2]
            objCoord.z *= intv[3]

            io = 1
            iNobj = -1
            index = -1

            # check the photon intersection with big-voxel walls
            errCode = planes.calPlanes(phoCoord, vectCoord, objCoord, intv)
            distancePho = planes.distance
            if (errCode == ERRCODE.CANNOT_FIND):
                # update the x0, y0, z0
                phoCoord.x = planes.x
                phoCoord.y = planes.y
                phoCoord.z = planes.z

            # check the photon intersection with objects
            if (comm.N_DIVS[iVOX] != 0):

                distanceObj = 1.0e5
                for idiv in range(1, comm.N_DIVS[iVOX]):

                    index = comm.DIVS[iVOX][idiv]

                    for l in range(1, 6):
                        tobj[l] = comm.OBJ[index][l-1]

                    treeBoundary.dealTreeType(comm.OBJ_Shape[index], phoCoord, vectCoord, tobj)

                    tempDistance = treeBoundary.distance
                    tempIO = treeBoundary.io

                    if (tempDistance < distanceObj):
                        distanceObj = tempDistance
                        iNobj = index
                        io = tempIO

                    if (io == 0):
                        break

                # if stem collision, return
                if (comm.OBJ_Shape[iNobj] == 4) and (distanceObj < 1.0e5):
                    self.sFlag = 0
                    return ERRCODE.FAILURE

            # increment of the optical path
            if (io == 0):
                # check branch optical thickness
                tobjb[1] = tobj[1]
                tobjb[2] = tobj[2]
                tobjb[3] = tobj[3] - tobj[4] * (1.0 - comm.RB) * min(1, abs(comm.OBJ_Shape[iNObj] - 5))
                tobjb[4] = tobj[4] * comm.RB
                tobjb[5] = tobj[5] * comm.RB

                tobj12[1] = tobj[1]
                tobj12[2] = tobj[2]
                tobj12[3] = tobj[3] - tobj[4] * (1.0 - comm.RB) * min(1, abs(comm.OBJ_Shape[iNObj] - 5))
                tobj12[4] = tobj[4] * rb12
                tobj12[5] = tobj[5] * rb12

                treeBoundary.dealTreeType(comm.OBJ_Shape[index], phoCoord, vectCoord, tobj)
                io2 = treeBoundary.io
                distance2 = treeBoundary.distance

                treeBoundary.dealTreeType(comm.OBJ_Shape[index], phoCoord, vectCoord, tobj12)
                io12 = treeBoundary.io
                distance12 = treeBoundary.distance

                # if outside of branch (io12 = 1) go into the branch
                if (io12 == 1):
                    xb = phoCoord.x + (distance12 + mgn) * vectCoord.x
                    yb = phoCoord.y + (distance12 + mgn) * vectCoord.y
                    zb = phoCoord.z + (distance12 + mgn) * vectCoord.z
                    branchCoord = Position()
                    branchCoord.setPosition(xb, yb, zb)

                    treeBoundary.dealTreeType(comm.OBJ_Shape[index], branchCoord, vectCoord, tobj12)
                    distance12 = treeBoundary.distance

                # if outside of branch (io2 = 1) go into the branch
                if (io2 == 1):
                    xb = phoCoord.x + (distance2 + mgn) * vectCoord.x
                    yb = phoCoord.y + (distance2 + mgn) * vectCoord.y
                    zb = phoCoord.z + (distance2 + mgn) * vectCoord.z
                    branchCoord = Position()
                    branchCoord.setPosition(xb, yb, zb)

                    treeBoundary.dealTreeType(comm.OBJ_Shape[index], branchCoord, vectCoord, tobj12)
                    distance2 = treeBoundary.distance

                ################################
                # calculation of optical path
                ################################
                # Don't know what's the meaning
                ################################
                th = acos(vectCoord.z)
                ith = int(radians(th))
                rio = 1.0 - float(io2)
                rio12 = 1.0 - float(io12)
                cf = comm.S_BAR[comm.OBJ_Group[index]]
                cf12 = comm.S_BAR[comm.OBJ_Group[index]]

                taub = rio * (distance2 + mgn) * comm.BAD[comm.OBJ_Group[index]] * gFunction.GT_BLB[ith] * comm.BP2
                taub += rio * (distance2 + mgn) * comm.U[comm.OBJ_Group[index]] * gFunction.GT_BLC[ith] * 4.0 * cf * (1.0 - comm.BP2)

                tauc12 = ((distance12 + mgn) - (distance12 + mgn) * rio) * rio12
                tauc12 = comm.U[comm.OBJ_Group[index]] * gFunction.GT_BLC[ith] * 4.0 * cf12 * (1.0 - comm.BP1) \
                         + (comm.BAD[comm.OBJ_Group[index]] * gFunction.GT_BLB[ith]) * comm.BP1

                tauc = (distanceObj + mgn) - (distance12 + mgn) * rio12
                tauc = comm.U[comm.OBJ_Group[index]] * gFunction.GT_BLC[ith] * 4.0 * cf * (1.0 - comm.BP1) \
                         + (comm.BAD[comm.OBJ_Group[index]] * gFunction.GT_BLB[ith]) * comm.BP1

                self.tau += tauc + tauc12 + taub

            # refresh photon position
            d = distanceObj * (1.0 - float(io)) + min(distanceObj, distancePho) * float(io)
            phoCoord.x += (d + mgn) * vectCoord.x
            phoCoord.y += (d + mgn) * vectCoord.y
            phoCoord.z += (d + mgn) * vectCoord.z

            # check upper or bottom boundary condition
            if (copysign(1.0, vectCoord.z) * (phoCoord.z - zlim[flag]) >= 0.0):
                return ERRCODE.SUCCESS

            phoCoord.x -= (trunc(phoCoord.x / comm.X_MAX) - 0.5 + copysign(0.5, phoCoord.x)) * comm.X_MAX
            phoCoord.y -= (trunc(phoCoord.y / comm.Y_MAX) - 0.5 + copysign(0.5, phoCoord.y)) * comm.Y_MAX

        print("Vegetation trace finish.")
        return ERRCODE.SUCCESS
