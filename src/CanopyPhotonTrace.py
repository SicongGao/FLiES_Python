import ERRCODE
import common as comm
from math import *
from Planes import Planes
from TreeBoundary import TreeBoundary
from G_Function import G_Function
from MonteCarlo import MonteCarlo
from Position import Position
from iparam import Parameters

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
    weight = 0
    cNscat = 0
    cIchi = 0
    cIkd = 0

    def __init__(self):
        return

    def save(self, w, nscat, ichi, ikd):
        self.weight = w
        self.cNscat = nscat
        self.cIchi = ichi
        self.cIkd = ikd

        return ERRCODE.SUCCESS

    def trace(self, phoCoord, vectCoord, w, wq, nscat, ichi, ikd, iParameter):

        iVOX = 0
        distancePho = 0.0

        MIN_VALUE = 1.0e-8  # weight minimum limit
        distanceObj = 1.0e5    # initial distance from object
        intv = [50.0 / comm.RES] * 4
        tObj = [0.0] * 6
        face = 0

        # marginal value
        mgn = 1.0e-2

        #x = x0
        #y = y0
        phoCoord.z = comm.Z_MAX - mgn
        objCoord = Position()

        planes = Planes()
        treeBoundary = TreeBoundary()
        gFunction = G_Function()
        gFunction.igtbl()
        mcSimulation = MonteCarlo()

        # do wile photon exit from canopy space
        while (1):
            # objCoord.x = trunc(phoCoord.x / intv[1])
            # y1 = trunc(phoCoord.y / intv[2])
            # z1 = trunc(phoCoord.z / intv[3])
            objCoord.setPosition(trunc(phoCoord.x / intv[1]),
                                 trunc(phoCoord.y / intv[2]),
                                 trunc(phoCoord.z / intv[3]))

            iVOX = comm.IX_MAX * comm.IY_MAX * int(objCoord.z)
            iVOX += int(objCoord.y) * comm.IY_MAX
            iVOX += int(objCoord.x) + 1

            objCoord.x *= intv[1]
            objCoord.y *= intv[2]
            objCoord.z *= intv[3]

            io = 1
            iNobj = -1
            index = -1
            distanceObj = 1.0e5
            
            # check the photon intersection with big-voxel walls
            errCode = planes.calPlanes(phoCoord, vectCoord, objCoord, intv)
            distancePho = planes.distance
            if (errCode == ERRCODE.CANNOT_FIND):
                # update the x, y, z
                phoCoord.x = planes.x
                phoCoord.y = planes.y
                phoCoord.z = planes.z

            # check the photon intersection with objects
            if (comm.N_DIVS[iVOX] != 0):

                distanceObj = 1.0e5
                for idiv in range(1, comm.N_DIVS[iVOX]):

                    index = comm.DIVS[iVOX][idiv]

                    tObj[1:6] = comm.OBJ[index][0:5]

                    treeBoundary.dealTreeType(comm.T_OBJ[index], phoCoord, vectCoord, tObj)

                    tempDistance = treeBoundary.distance
                    tempIO = treeBoundary.io
                    face = treeBoundary.face

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
                    phoCoord.x += vectCoord.x * distanceObj
                    phoCoord.y += vectCoord.y * distanceObj
                    phoCoord.z += vectCoord.z * distanceObj

                    tObj[1:6] = comm.OBJ[iNobj][0:5]

                    # !!!! modified the input parameters here !!!!
                    mcSimulation.stem(w, wq, phoCoord, vectCoord, nscat,
                                      tObj, face, ichi, ikd, iParameter)
                    ichi = mcSimulation.cIchi
                    ikd = mcSimulation.cIkd
                    nscat = mcSimulation.cNscat
                    w = mcSimulation.weight

                    # !!!! modified the weight here !!!!
                    if (w < MIN_VALUE):
                        return ERRCODE.CANNOT_FIND

                    phoCoord.x += mgn * vectCoord.x
                    phoCoord.y += mgn * vectCoord.y
                    phoCoord.z += mgn * vectCoord.z
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



        self.save(w, nscat)

        return ERRCODE.SUCCESS

# a = [1,2,3,4,5,6,7,8]
# b = [0]*10
#
# print(a)
# print(b)
# b[1:5] = a[1:5]
# print(b)
#
# class coo:
#     x = 0
#     y = 0
#     z = 0
#
# def change(a):
#     print(id(a))
#     a+= "a"
#     print(a)
#
# a = Position()
# b = Position()
# a.setPosition(1,2,3)
# b.setPosition(4,5,6)
# print(b.x)
# b = a
# del a
# print(b.x)
