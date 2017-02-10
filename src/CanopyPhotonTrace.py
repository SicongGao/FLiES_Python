import ERRCODE
import common as comm
from math import *
from Planes import Planes
from TreeBoundary import TreeBoundary
from G_Function import G_Function
from MonteCarlo import MonteCarlo

# ##################################################
# Simulator for forest light environment
# start MC simulation
#
# by H. Kobayashi
# modified 08/03/26
# ##################################################


class CanopyPhotonTrace:

    sFlag = 1
    tau = 0.0

    def __init__(self):
        return

    def trace(self, x, y, uxr, uyr, uzr):

        iVOX = 0
        distancePho = 0.0

        MIN_VALUE = 1.0e-8  # wieght minimum limit
        distanceObj = 1.0e5    # initial distance from object
        intv = [50.0 / comm.RES] * 4
        tObj = [0.0] * 6

        # marginal value
        mgn = 1.0e-2

        x0 = y0 = z0 = 0.0

        x0 = x
        y0 = y
        z0 = comm.Z_MAX - mgn

        planes = Planes()
        treeBoundary = TreeBoundary()
        gFunction = G_Function()
        gFunction.igtbl()
        mcSimulation = MonteCarlo()

        # do wile photon exit from canopy space
        while (1):
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
            distanceObj = 1.0e5
            
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
                        tObj[l] = comm.OBJ[index][l - 1]

                    treeBoundary.dealTreeType(comm.T_OBJ[index], x0, y0, z0, uxr, uyr, uzr, tObj)

                    tempDistance = treeBoundary.distance
                    tempIO = treeBoundary.io

                    if (tempDistance < distanceObj):
                        distanceObj = tempDistance
                        iNobj = index
                        io = tempIO

                    if (io == 0):
                        break

            # canopy interaction
            if ((distanceObj < distancePho) or (io == 0)):

                # stem interaction
                if (comm.T_OBJ[iNobj] == 4):
                    # if the photon is in the stem object by mistake, exit from stem
                    # in other case, photon go to the stem surface

                    return

                # canopy interaction [Monte Carlo in canopy media]
                else:
                    return

            # big-voxel wall interaction
            else:
                # Monte Carlo in forest floor
                if():
                    return

                # sky (exit from canopy space)
                elif():
                    return

                # refresh the x, y position using the boundary condition
                else:
                    return







        return ERRCODE.SUCCESS