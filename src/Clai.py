from TreeBoundary import TreeBoundary
import ERRCODE
import common as comm
from math import *
from VegTrace import VegTrace
from Position import Position

# *****************************************
# LAI calculation
# spn is the sampling grid unit (m)
# *****************************************


class CLAI:

    LAI = 0.0
    pLAI = []
    crownCover = 0.0

    def calLAI(self, spn):

        self.LAI = 0.0
        self.crownCover = 0.0

        sFlag = 1
        tau = 0.0
        tLai = 0.0
        intv = 50.0 / comm.RES

        # ux = 0.0
        # uy = 0.0
        # uz = 1.0
        vectCoord = Position()
        vectCoord.setPosition(0.0, 0.0, 1.0)
        # x = y = z = 0.0
        phoCoord = Position()
        phoCoord.setPosition(0.0, 0.0, 1.0)

        m = 1
        its = 0
        fc = io1 = io2 = 0
        tobj = [0] * 6
        tobjb = [0] * 6

        treeBoundary = TreeBoundary()

        for k in range(1, int(comm.Z_MAX / spn) + 1):
            phoCoord.z = (float(k) - 0.5) * spn
            if (fmod(k, int(1.0 / spn)) == 0):
                m += 1

            for i in range(1, int(comm.X_MAX / spn) + 1):
                phoCoord.x = (float(i) - 0.5) * spn

                for j in range(1, int(comm.Y_MAX / spn) + 1):
                    phoCoord.y = (float(j) - 0.5) * spn

                    x1 = trunc(phoCoord.x / intv)
                    y1 = trunc(phoCoord.y / intv)
                    z1 = trunc(phoCoord.z / intv)

                    ivox = comm.IX_MAX * comm.IY_MAX * z1
                    ivox += int(y1) * comm.IY_MAX
                    ivox += int(x1)
                    # print(ivox, x1, y1)
                    if (comm.N_DIVS[ivox] != 0):

                        for idiv in range(comm.N_DIVS[ivox]):

                            # selected object number
                            index = comm.DIVS[ivox, idiv]

                            tobj[1:6] = comm.OBJ[index][0:5]

                            # define the branch dominant region
                            tobjb[1] = tobj[1]
                            tobjb[2] = tobj[2]
                            tobjb[3] = tobj[3] - tobj[4] * comm.RB
                            tobjb[4] = tobj[4] * comm.RB
                            tobjb[5] = tobj[5] * comm.RB

                            treeBoundary.dealTreeType(comm.OBJ_Shape[i], phoCoord, vectCoord, tobj)
                            io1 = treeBoundary.io
                            treeBoundary.dealTreeType(comm.OBJ_Shape[i], phoCoord, vectCoord, tobjb)
                            io2 = treeBoundary.io

                            if (io2 == 0):
                                tLai += comm.U[comm.OBJ_Group[index]] * (1.0 - comm.BP2)
                                its += 1
                            elif (io1 == 0):
                                tLai += comm.U[comm.OBJ_Group[index]] * (1.0 - comm.BP1)
                                its += 1

                        if (its != 0):
                            self.pLAI[m] += tLai / float(its)

                        its = 0
                        tLai = 0.0

        multiSpnXYMAX = spn * spn * spn / (comm.X_MAX * comm.Y_MAX)

        # for i in range(1, 100):
        #     self.pLAI[i] *= multiSpnXYMAX
        #     self.LAI += self.pLAI[i]
        self.LAI = sum(self.pLAI) * multiSpnXYMAX

        # crown cover calculation
        print("crown cover calculation ...")

        phoCoord.z = 0.01

        # dummy clumping factor and lea area density
        for i in range(1, 5):
            comm.S_BAR[i] = 0.25
            comm.U[i] = 1.0

        vegTrace = VegTrace()
        for i in range(1, int(comm.X_MAX / spn)):
            phoCoord.x = (float(i) - 0.5) * spn

            for j in range(1, int(comm.Y_MAX / spn)):
                phoCoord.y = (float(j) - 0.5) * spn
                vegTrace.trace(phoCoord, vectCoord)
                if (vegTrace.tau != 0.0):
                    self.crownCover += 1.0

        self.crownCover /= float(int(comm.X_MAX / spn) * int(comm.Y_MAX / spn))

        return ERRCODE.SUCCESS
