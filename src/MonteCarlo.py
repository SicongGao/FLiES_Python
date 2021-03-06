import ERRCODE
from Random import Random
from math import *
import common as comm
from iparam import Parameters
from Position import Position
import RussianRoulette
import TransformCoordinate
from VegRadiation import VegRadiation
import G_Function
from TreeBoundary import TreeBoundary
from Planes import Planes
import logging


class MonteCarlo:

    # nscat = 0
    # cNscat = 0
    # cIchi = 0
    # cIkd = 0
    # weight = -1

    def __init__(self):

        self.cNscat = 0
        self.cIchi = 0
        self.cIkd = 0
        self.weight = -1

    def save(self, nscat, w):
        self.cNscat = nscat
        self.weight = w
        string = "weight = " + str(self.weight) + ", nscat = " + str(self.cNscat)
        logging.debug(string)
        return ERRCODE.SUCCESS

    # recollision loop for shoot clumping effect
    def recollisionShootClump(self):




        return ERRCODE.SUCCESS

    # ******************************************
    # scattering in the surface boundary
    # ******************************************
    def scatterReflectSurface(self, vectCoord, mode):

        MIN_VALUE = 0.0174524
        # lambertian reflection
        if (mode == 1):
            randMethod = Random()
            th = 0.5 * acos(1.0 - 2.0 * randMethod.getRandom())
            ph = 2.0 * pi * randMethod.getRandom()
            vectCoord.setPosition(sin(th) * cos(ph),
                                  sin(th) * sin(ph),
                                  cos(th))
            if (abs(vectCoord.z) < MIN_VALUE):
                vectCoord.z = copysign(MIN_VALUE, vectCoord.z)

            self.weight = RussianRoulette.roulette(self.weight)
            self.cNscat += 1

        return ERRCODE.SUCCESS
    # ****************************************************
    # Scattering direction
    # Hideki Kobayashi
    # We use the method of Frontier Technical Report No7
    # p90-91, based on Rejection method
    # ****************************************************
    def scatterDirection(self, leaf_reflectance, leaf_transmittance, vectCoord, m):

        fglm = (0, 1.0, 4.0 / pi, 4.0 / pi)
        phl = thl = 0
        # step 1: determination of the leaf normal vector
        cRandom = Random()
        rand = 1
        cosa = rand - 1
        ref = - 1
        vectLCoord = Position()
        vectObjCoord = Position()

        while (rand > abs(cosa)):

            rand = cRandom.getRandom()
            phl = 2.0 * pi * rand
            # make sure it could go to the inner loop
            ref = rand - 1
            while(rand > ref):
                rand = cRandom.getRandom()
                thl = 0.5 * pi * rand
                ref = G_Function.fgl(thl, m) / fglm[m]
                rand = cRandom.getRandom()

            vectLCoord.setPosition(sin(thl) * cos(phl),
                                   sin(thl) * sin(phl),
                                   cos(thl))

            # step 2: Adjustment to follow the |omega*omegaL|=cosa
            cosa = vectCoord.x * vectLCoord.x + vectCoord.y * vectLCoord.y + vectCoord.z * vectLCoord.z
            rand = cRandom.getRandom()

        # step 3: Determination of the scattering direction on the leaf
        # Reflection or transmission
        b = 2.0 * pi * cRandom.getRandom()
        a = sqrt(cRandom.getRandom())
        a = acos(a)

        # reflection
        if (cRandom.getRandom() <= (leaf_reflectance / (leaf_reflectance + leaf_transmittance))):
            if (cosa >= 0.0):
                b += pi
                if (b > (2.0 * pi)):
                    b -= 2.0 * pi
                a = pi - a
        # transmission
        elif (cosa < 0.0):
            b += pi
            if (b > (2.0 * pi)):
                b -= 2.0 * pi
            a = pi - a

        # For tho and pho, coordinate transformation
        TransformCoordinate.transformCoordinate(vectLCoord, a, b, vectCoord)

        return ERRCODE.SUCCESS

    def stem(self, w, wq, phoCoord, vectCoord, nscat, tObj, face, trunkRef, ichi, ikd, para):
        MGN = 1.0e-2
        MIN_VALUE = 1.0e-8
        FD = 0.0
        CB = 3
        UZM = 0.0174524
        a = 0
        fd = 0.0

        randomMethod = Random()
        vegRadiation = VegRadiation()

        logging.debug("Monte Carlo stem simulation start...")
        string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
        logging.debug(string)
        # reflectance at the side of stem
        if (face == 1):
            # stem normal vector
            rx = phoCoord.x - tObj[1]
            ry = phoCoord.y - tObj[2]
            rr = sqrt(rx ** 2 + ry ** 2)

            tempCoord = Position()
            tempCoord.setPosition(rx / rr, ry / rr, 0)

            th = 0.5 * acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()

            a = acos(tempCoord.x)
            CB = 6

            TransformCoordinate.transformCoordinate(tempCoord, th, ph, vectCoord)

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

        # reflectance at the bottom of stem
        elif (face == 2):
            th = 0.5 * pi + acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()
            vectCoord.setPosition(sin(th) * cos(ph),
                                  sin(th) * sin(ph),
                                  cos(th))
            phoCoord.z -= MGN

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

            a = 1.0

        # reflectance at the top of stem
        else:
            th = 0.5 * acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()
            vectCoord.setPosition(sin(th) * cos(ph),
                                  sin(th) * sin(ph),
                                  cos(th))
            phoCoord.z += MGN

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

            a = 1.0

        # fpar samping (leave or branch)
        # every tree specie may have the same reflectance (trunkRef)
        para.B_FPR += w * wq * (1.0 - trunkRef)
        # ipara.AP_NP range (0,100)
        iz = int(phoCoord.z) + 1
        if (iz > 99):
            logging.ERROR("Z overhight: MonteCarlo.py")
            self.save(nscat, w)
            return ERRCODE.OUT_OF_RANGE
        else:
            para.AP_NP[iz] += w * wq * (1.0 - trunkRef)

        w *= trunkRef
        nscat += 1

        # Russian Roulette
        w = RussianRoulette.roulette(w)

        # call vegrad()
        vegRadiation.simulate(phoCoord, vectCoord, w, 1.0, 0.0, CB, a, fd, ichi, ikd, para)
        if (CB == 6):
            a = vegRadiation.save_a

        logging.debug("Monte Carlo stem simulation finish.")
        self.save(nscat, w)
        return ERRCODE.SUCCESS

    def canopy(self, w, wq, phoCoord, vectCoord, nscat, tObj, inobj, trunkRef, ichi, ikd, leaf_reflectance, leaf_transmittance, para):

        MIN_VALUE = 1.0e-8
        mgn = 1.0e-2
        UZ_MIN = 0.0174524

        temp = comm.OBJ_Group[inobj]
        cf = comm.S_BAR[temp]
        branchAreaDensity = comm.BAD[temp]
        la = comm.U[temp]
        obj_shape = comm.OBJ_Shape[inobj]

        cf12 = comm.S_BAR[comm.OBJ_Group[inobj]]
        rb12 = 1.0
        fd = 0.0

        # define the second canopy area
        tObj12 = (0, tObj[1],
                     tObj[2],
                     tObj[3] - tObj[4] * (1.0 - rb12) * min(1, abs(obj_shape - 5)),
                     tObj[4] * rb12,
                     tObj[5] * rb12)

        # define the branch dominant region
        tObjb = (0, tObj[1],
                    tObj[2],
                    tObj[3] - tObj[4] * (1.0 - comm.RB) * min(1, abs(obj_shape - 5)),
                    tObj[4] * comm.RB,
                    tObj[5] * comm.RB)

        treeBounday = TreeBoundary()
        randMethod = Random()
        vegRadiant = VegRadiation()

        logging.debug("Monte Carlo canopy simulation start...")
        string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
        logging.debug(string)
        # check status:leaf dominant region, branch dominant or irregulary outside
        while (1):
            io1 = 0
            io2 = 0
            io12 = 0

            treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObj)
            distance1 = treeBounday.distance
            io1 = treeBounday.io

            treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObj12)
            distance12 = treeBounday.distance
            io12 = treeBounday.io

            treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObjb)
            distance2 = treeBounday.distance
            io2 = treeBounday.io

            string = "phoCoord = " + str(phoCoord.x) + ", " + str(phoCoord.y) + ", " + str(phoCoord.z)
            logging.debug(string)
            string = "vectCoord = " + str(vectCoord.x) + ", " + str(vectCoord.y) + ", " + str(vectCoord.z)
            logging.debug(string)
            string = "io1, io2, io12 = " + str(io1) + ", " + str(io2) + ", " + str(io12)
            logging.debug(string)
            string = "s1, s2, s12 = " + str(distance1) + ", " + str(distance2) + ", " + str(distance12)
            logging.debug(string)

            if (io1 == 1):
                self.save(nscat, w)
                logging.debug("Monte Carlo canopy simulation finish. (io1 == 1)")
                return ERRCODE.OUTSIDE

            # if photon inside branch dominant region
            if (io2 == 0):
                logging.debug("MC Canopy: io2 == 0")
                while (1):
                    rand = randMethod.getRandom()

                    # branch scattering or leaf scattering
                    th = acos(vectCoord.z)
                    ith = int(degrees(th))
                    bp = 0.5 + copysign(0.5, comm.BP2 - rand)

                    sgm = comm.GT_BLB[ith] * branchAreaDensity * bp
                    sgm += 4.0 * cf * comm.GT_BLC[ith] * la * (1.0 - bp)
                    sgm = sgm * bp + (sgm / comm.FE) * (1.0 - bp)
                    sgm = max(1.0e-5, sgm)
                    ref = trunkRef * bp + leaf_reflectance * (1.0 - bp)
                    tr = leaf_transmittance * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance < distance12):
                        phoCoord.movePositionDistance(distance, vectCoord, comm.X_MAX, comm.Y_MAX)
                        # re-collision loop for shoot clumping effect
                        while (1):
                            # collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar sampling (leave or branch)
                            para.B_FPR += w * wq * (1.0 - ssa) * bp
                            para.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            par = w * wq * (1.0 - ssa)
                            para.AP[ix, iy, iz] += (1.0 - bp) * par
                            para.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            para.AP_D[ix, iy, iz] += (1.0 - bp) * min(nscat, 1) * par
                            para.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            w = RussianRoulette.roulette(w)

                            if (w < MIN_VALUE):
                                self.save(nscat, w)
                                logging.debug("Monte Carlo canopy simulation finish. (low w)[1]")
                                return ERRCODE.LOW_WEIGHT

                            psh = (1.0 - 4.0 * cf) * (1.0 - bp)
                            if (randMethod.getRandom() > psh):
                                break

                        # radiance sampling
                        cb = int((1.0 - bp) + 2.0 * bp)
                        vegRadiant.simulate(phoCoord, vectCoord, w, ref, tr, cb, 1.0, fd, ichi, ikd, para)

                        # new photon direction after scattering
                        if (randMethod.getRandom() > fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_B)
                            if (abs(vectCoord.z) < UZ_MIN):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        # check status
                        treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObjb)
                        io2 = treeBounday.io
                        distance2 = treeBounday.distance

                        # for unexpected case(photon is outside of canopy)
                        if (io2 == 1):
                            break
                    else:
                        phoCoord.movePositionDistance(distance2, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break

            # if photon inside the canopy dominant region
            elif (io12 == 0):
                logging.debug("MC Canopy: io12 == 0")
                while (1):
                    rand = randMethod.getRandom()
                    # branch scattering or leaf scattering
                    th = acos(vectCoord.z)
                    ith = int(degrees(th))
                    bp = 0.5 + copysign(0.5, comm.BP1 - rand)

                    sgm = comm.GT_BLB[ith] * branchAreaDensity * bp
                    sgm += 4.0 * cf * comm.GT_BLC[ith] * la * (1.0 - bp)
                    sgm = sgm * bp + (sgm / comm.FE) * (1.0 - bp)
                    sgm = max(1.0e-5, sgm)

                    ref = trunkRef * bp + leaf_reflectance * (1.0 - bp)
                    tr = leaf_transmittance * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance > distance2):
                        phoCoord.movePositionDistance(distance2, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    elif (distance > distance12):
                        phoCoord.movePositionDistance(distance12, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    else:
                        phoCoord.movePositionDistance(distance, vectCoord, comm.X_MAX, comm.Y_MAX)

                        while (1):
                            #  collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar samping (leave or branch)
                            para.B_FPR += w * wq * (1.0 - ssa) * bp
                            para.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            ix = min(ix, comm.SIZE - 1)
                            iy = min(iy, comm.SIZE - 1)

                            par = w * wq * (1.0 - ssa)
                            para.AP[ix, iy, iz] += par * (1.0 - bp)
                            para.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            para.AP_D[ix, iy, iz] += par * (1.0 - bp) * min(nscat, 1)
                            para.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            psh = (1.0 - 4.0 * cf) * (1.0 - bp)
                            if (randMethod.getRandom() >= psh):
                                break

                        w = RussianRoulette.roulette(w)

                        if (w < MIN_VALUE):
                            self.save(nscat, w)
                            logging.debug("Monte Carlo canopy simulation finish. (low w)[2]")
                            return ERRCODE.LOW_WEIGHT

                        cb = int((1.0 - bp) + 2.0 * bp)

                        vegRadiant.simulate(phoCoord, vectCoord, w, leaf_reflectance, leaf_transmittance, cb, 1.0, fd, ichi, ikd, para)

                        rand = randMethod.getRandom()
                        if (rand >= fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_C)
                            if (abs(vectCoord.z) < UZ_MIN):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObj12)
                        io12 = treeBounday.io
                        distance12 = treeBounday.distance
                        treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObjb)
                        io2 = treeBounday.io
                        distance2 = treeBounday.distance

                        # for unexpected case (photon is outside of canopy)
                        if (io12 == 1):
                            break

            # if photon inside the canopy dominant region
            elif (io1 == 0):
                logging.debug("MC Canopy: io1 == 0")
                while (1):
                    rand = randMethod.getRandom()
                    # branch scattering or leaf scattering
                    th = acos(vectCoord.z)
                    ith = int(degrees(th))

                    bp = 0.5 + copysign(0.5, float(comm.BP1 - rand))

                    sgm = comm.GT_BLB[ith] * branchAreaDensity * bp
                    sgm += 4.0 * cf * comm.GT_BLC[ith] * la * (1.0 - bp)
                    sgm = sgm * bp + (sgm / comm.FE) * (1.0 - bp)
                    sgm = max(1.0e-5, sgm)

                    ref = trunkRef * bp + leaf_reflectance * (1.0 - bp)
                    tr = leaf_transmittance * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance > distance12):
                        phoCoord.movePositionDistance(distance12, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    elif (distance > distance1):
                        phoCoord.x += (distance1 + mgn) * vectCoord.x
                        phoCoord.y += (distance1 + mgn) * vectCoord.y
                        phoCoord.z += (distance1 + mgn) * vectCoord.z
                        self.save(nscat, w)

                        logging.debug("Monte Carlo canopy simulation finish. distance > distance1")

                        return ERRCODE.OUTSIDE
                    else:
                        phoCoord.movePositionDistance(distance, vectCoord, comm.X_MAX, comm.Y_MAX)

                        while (1):

                            #  collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar samping (leave or branch)
                            para.B_FPR += w * wq * (1.0 - ssa) * bp
                            para.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            ix = min(ix, comm.SIZE - 1)
                            iy = min(iy, comm.SIZE - 1)

                            par = w * wq * (1.0 - ssa)
                            para.AP[ix, iy, iz] += par * (1.0 - bp)
                            para.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            para.AP_D[ix, iy, iz] += par * (1.0 - bp) * min(nscat, 1)
                            para.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            psh = (1.0 - 4.0 * cf) * (1.0 - bp)
                            if (randMethod.getRandom() >= psh):
                                break

                        w = RussianRoulette.roulette(w)
                        self.save(nscat, w)

                        if (w < MIN_VALUE):
                            self.save(nscat, w)

                            logging.debug("Monte Carlo canopy simulation finish. (low w)[3]")

                            return ERRCODE.LOW_WEIGHT

                        cb = int((1.0 - bp) + 2.0 * bp)
                        vegRadiant.simulate(phoCoord, vectCoord, w, leaf_reflectance, leaf_transmittance, cb, 1.0, fd, ichi, ikd, para)

                        if (randMethod.getRandom() >= fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_C)
                            if (vectCoord.z < UZ_MIN):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        # check status
                        treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObj)
                        distance1 = treeBounday.distance
                        io1 = treeBounday.io
                        treeBounday.dealTreeType(obj_shape, phoCoord, vectCoord, tObj12)
                        distance12 = treeBounday.distance
                        io12 = treeBounday.io

                        if (io1 == 1):
                            self.save(nscat, w)

                            logging.debug("Monte Carlo canopy simulation finish. (io1 == 1)")

                            return ERRCODE.OUTSIDE

        self.save(nscat, w)
        #print("Monte Carlo canopy simulation finish.")

        logging.debug("Monte Carlo canopy simulation finish.")
        string = "weight = " + str(self.weight) + ", nscat = " + str(self.cNscat)
        logging.debug(string)

        return ERRCODE.SUCCESS

    # ##################################################################################
    # This subroutine calculates photon trajectory in the forest floor
    # rand     : random numbers
    # distance : traveling distance
    # sgm      : extinction
    # w        : weight of the photon
    # z (m)    : Height of the photon position from the upper grass boundary (-1.0<=z<0)
    # ##################################################################################S
    def floor(self, w, wq, phoCoord, vectCoord, nscat, floor_reflectance, floor_transmittance, sor, ichi, ikd, para):

        # tentative fd
        fd = 0.0

        MIN_VALUE = 1.0e-8
        mgn = 1.0e-2
        MIN_Z = 0.0174524

        # clumping factor
        cf = comm.S_BAR[1]

        # mode = 1 Lambertian
        mode = 1

        # upper and bottom of the forest floor layer
        zUpper = 0.0
        zBottom = -1.0

        intv = (0, comm.X_MAX, comm.Y_MAX, abs(zUpper - zBottom))

        randMethod = Random()
        plane = Planes()
        vegRadiant = VegRadiation()

        logging.debug("Monte Carlo floor simulation start...")
        string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
        logging.debug(string)
        string = "vectCoord = " + str(vectCoord.x) + ", " + str(vectCoord.y) + ", " + str(vectCoord.z)
        logging.debug(string)
        # Monte Carlo loop
        while (1):
            # string = "Every while: vectCoord = " + str(vectCoord.x) + ", " + str(vectCoord.y) + ", " + str(vectCoord.z)
            # logging.debug(string)
            rand = randMethod.getRandom()
            th = acos(vectCoord.z)
            ith = int(degrees(th))
            sgm = comm.GT_BLF[ith] * comm.G_LAI
            sgm = max(1.0e-5, sgm)
            distance = -log(rand) / sgm

            # check the intersection
            objCoord = Position()
            objCoord.setPosition(0.0, 0.0, -1.0)
            plane.calPlanes(phoCoord, vectCoord, objCoord, intv)

            distancePlane = plane.distance

            if (distance < distancePlane):
                # scattering
                phoCoord.moveVector(distance, vectCoord)

                # recollision
                while (1):
                    para.F_FPR += w * wq * (1.0 - floor_reflectance - floor_transmittance)

                    ix = int(phoCoord.x * comm.RES) + 1
                    iy = int(phoCoord.y * comm.RES) + 1

                    ix = min(ix, comm.SIZE - 1)
                    iy = min(iy, comm.SIZE - 1)

                    para.AP_F[ix, iy] += w * wq * (1.0 - floor_transmittance - floor_reflectance)
                    para.AP_FD[ix, iy] += w * wq * (1.0 - floor_transmittance - floor_reflectance) * min(float(nscat), 1.0)

                    w *= (floor_transmittance + floor_reflectance)
                    nscat += 1

                    w = RussianRoulette.roulette(w)
                    if (w < MIN_VALUE):
                        self.save(nscat, w)

                        logging.debug("Monte Carlo floor simulation finish.")

                        return ERRCODE.LOW_WEIGHT

                    psh = 1.0 - 4.0 * cf
                    if (randMethod.getRandom() >= psh):
                        break

                vegRadiant.simulate(phoCoord, vectCoord, w, floor_reflectance, floor_transmittance, 4, 1.0, fd, ichi, ikd, para)

                # new direction
                self.scatterDirection(floor_reflectance, floor_transmittance, vectCoord, comm.M_F)
                # string = "scatter: vectCoord = " + str(vectCoord.x) + ", " + str(vectCoord.y) + ", " + str(vectCoord.z)
                # logging.debug(string)
                if (abs(vectCoord.z) < MIN_Z):
                    vectCoord.z = copysign(MIN_Z, vectCoord.z)

            else:
                logging.debug("Monte Carlo floor: else")
                # string = "before: phoCoord = " + str(phoCoord.x) + ", " + str(phoCoord.y) + ", " + str(phoCoord.z)
                # logging.debug(string)
                # string = "vectCoord = " + str(vectCoord.x) + ", " + str(vectCoord.y) + ", " + str(vectCoord.z)
                # logging.debug(string)

                phoCoord.x += (distancePlane + mgn) * vectCoord.x
                phoCoord.y += (distancePlane + mgn) * vectCoord.y
                phoCoord.z += distancePlane * vectCoord.z
                phoCoord.x -= (int(phoCoord.x / comm.X_MAX) - 0.5 + copysign(0.5, phoCoord.x)) * comm.X_MAX
                phoCoord.y -= (int(phoCoord.y / comm.Y_MAX) - 0.5 + copysign(0.5, phoCoord.y)) * comm.Y_MAX

                # string = "after: phoCoord = " + str(phoCoord.x) + ", " + str(phoCoord.y) + ", " + str(phoCoord.z)
                # logging.debug(string)
                # logging.debug("sp = " + str(distancePlane))

                if (abs(phoCoord.z) < MIN_VALUE):
                    phoCoord.z = MIN_VALUE

                if (phoCoord.z <= zBottom):
                    # surface boundary reflectance mode 1: Lambertian, 2: RPV model, 3: DSM model
                    # fpar sampling
                    para.S_FPR += w * wq * (1.0 - sor)

                    ix = int(phoCoord.x * comm.RES) + 1
                    iy = int(phoCoord.y * comm.RES) + 1
                    ix = min(ix, comm.SIZE - 1)
                    iy = min(iy, comm.SIZE - 1)

                    para.AP_S[ix, iy] += w * wq * (1.0 - floor_reflectance - floor_reflectance)

                    # surface downward flux
                    para.SF_DIR[ix, iy] += w * wq * (1.0 - min(nscat, 1))
                    para.SF_DIF[ix, iy] += w * wq * min(nscat, 1)
                    w *= sor

                    if (w < MIN_VALUE):
                        self.save(nscat, w)
                        logging.debug("Monte Carlo floor simulation finish. (low w)")
                        return ERRCODE.LOW_WEIGHT

                    vegRadiant.simulate(phoCoord, vectCoord, w, floor_reflectance, floor_transmittance, 5, 1.0, fd, ichi, ikd, para)

                    self.save(nscat, w)
                    self.scatterReflectSurface(vectCoord, mode)
                    nscat = self.cNscat
                    w = self.weight

                    phoCoord.z = zBottom + mgn

                elif (phoCoord.z >= zUpper):
                    self.save(nscat, w)
                    logging.debug("Monte Carlo floor simulation finish.(higher)")
                    return ERRCODE.OUT_OF_RANGE

        self.save(nscat, w)

        logging.debug("Monte Carlo floor simulation finish.")

        return ERRCODE.SUCCESS


# ##############################
# FOR TEST
# Monte Carlo in Calculating PI
# ##############################

def test():
    R = 0.5
    inside = 0
    allNum = 10000
    for i in range(allNum):
        import Random
        rX = Random.uniform(0, 1)
        rY = Random.uniform(0, 1)
        rX = (rX - R) * (rX - R)
        rY = (rY - R) * (rY - R)
        rNEW = math.sqrt(rX + rY)

        if (rNEW <= R):
            inside += 1

    PI = 4 * inside / allNum
    print("PI = ", PI)
    print("ERR = ", abs((math.pi - PI) / math.pi))
