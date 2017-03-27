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
    def scatterDirection(self, lr, lt, vectCoord, m):

        fglm = (0, 1.0, 4.0 / pi, 4.0 / pi)
        phl = thl = 0
        # step 1: determination of the leaf normal vector
        cRandom = Random()
        rnd = 1
        cosa = rnd - 1
        ref = rnd - 1
        vectLCoord = Position()
        vectObjCoord = Position()

        while (rnd < abs(cosa)):

            while(rnd > ref):
                phl = 2.0 * pi * cRandom.getRandom()
                thl = 0.5 * pi * cRandom.getRandom()
                # TODO check fgl
                ref = G_Function.fgl(thl, m) / fglm[m]
                rnd = cRandom > ref

            vectLCoord.setPosition(sin(thl) * cos(phl),
                                   sin(thl) * sin(phl),
                                   cos(thl))

            # step 2: Adjustment to follow the |omega*omegaL|=cosa
            cosa = vectCoord.x * vectLCoord.x + vectCoord.y * vectLCoord.y + vectCoord.z * vectLCoord.z
            rnd = cRandom.getRandom()

        # step 3: Determination of the scattering direction on the leaf
        # Reflection or transmission
        b = 2.0 * pi * cRandom.getRandom()
        a = sqrt(cRandom.getRandom())
        a = acos(a)

        # reflection
        if (cRandom.getRandom() <= (lr / (lr + lt))):
            if (cosa >= 0.0):
                b += pi
                if (cosa > 0.0):
                    b -= 2.0 * pi
                    a = pi - a
        # transmission
        elif (cosa < 0.0):
            b += pi
            if (cosa > 0.0):
                b -= 2.0 * pi
                a = pi - a

        # For tho and pho, coordinate transformation
        TransformCoordinate.transformCoordinate(vectLCoord, a, b, vectCoord)

        return ERRCODE.SUCCESS

    def stem(self, w, wq, phoCoord, vectCoord, nscat, tObj, face, truncRef, ichi, ikd):
        MGN = 1.0e-2
        MIN_VALUE = 1.0e-8
        FD = 0.0
        CB = 3
        UZM = 0.0174524
        a = 0
        fd = 0.0

        randomMethod = Random()
        vegRadiation = VegRadiation()

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
            cb = 6

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
        # every tree specie may have the same reflectance (truncRef)
        comm.B_FPR += w * wq * (1.0 - truncRef)
        # ipara.AP_NP range (0,100)
        iz = int(phoCoord.z) + 1
        if (iz > 99):
            print("Z overhight: MonteCarlo.py")
            self.save(nscat, w)
            return ERRCODE.OUT_OF_RANGE
        else:
            comm.AP_NP[iz] += w * wq * (1.0 - truncRef)

        w *= truncRef
        nscat += 1

        # Russian Roulette
        w = RussianRoulette.roulette(w)

        # call vegrad()
        vegRadiation.simulate(phoCoord, vectCoord, CB, 1.0, 0.0, a, w, fd, ichi, ikd)
        a = vegRadiation.save_a

        self.save(nscat, w)
        return ERRCODE.SUCCESS

    def canopy(self, w, wq, phoCoord, vectCoord, nscat, tObj, inobj, truncRef, ichi, ikd, lr, lt):

        MIN_VALUE = 1.0e-8
        mgn = 1.0e-2
        UZ_MIN = 0.0174524

        temp = comm.I_OBJ[inobj]
        cf = comm.S_BAR[temp]
        branchAreaDensity = comm.BAD[temp]
        la = comm.U[temp]
        tsobj = comm.T_OBJ[inobj]

        cf12 = comm.S_BAR[comm.I_OBJ[inobj]]
        rb12 = 1.0
        fd = 0.0

        # define the second canopy area
        tObj12 = {0, tObj[1],
                     tObj[2],
                     tObj[3] - tObj[4] * (1.0 - rb12) * min(1, abs(tsobj - 5)),
                     tObj[4] * rb12,
                     tObj[5] * rb12}

        # define the branch dominant region
        tObjb = {0, tObj[1],
                    tObj[2],
                    tObj[3] - tObj[4] * (1.0 - comm.RB) * min(1, abs(tsobj - 5)),
                    tObj[4] * comm.RB,
                    tObj[5] * comm.RB}

        treeBounday = TreeBoundary()
        randMethod = Random()
        vegRadiant = VegRadiation()

        # check status:leaf dominant region, branch dominant or irregulary outside
        while (1):
            io1 = 0
            io2 = 0
            io12 = 0

            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObj)
            distance1 = treeBounday.distance
            io1 = treeBounday.io

            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObj12)
            distance12 = treeBounday.distance
            io12 = treeBounday.io

            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObjb)
            distance2 = treeBounday.distance
            io2 = treeBounday.io

            if (io1 == 1):
                self.save(nscat, w)
                return ERRCODE.OUTSIDE

            # if photon inside branch dominant region
            if (io2 == 0):
                while (1):
                    rand = randMethod.getRandom()

                    # branch scattering or leaf scattering
                    th = acos(vectCoord.z)
                    ith = int(degrees(th))
                    bp = 0.5 + copysign(0.5, comm.BP2 - rand)

                    sgm = comm.GT_BLB(ith) * branchAreaDensity * bp
                    sgm += 4.0 * cf * comm.GT_BLC[ith] * la * (1.0 - bp)
                    sgm = sgm * bp + (sgm / comm.FE) * (1.0 - bp)
                    sgm = max(1.0e-5, sgm)
                    ref = truncRef * bp + lr * (1.0 - bp)
                    tr = lt * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance < distance12):
                        phoCoord.movePosition(distance, vectCoord, comm.X_MAX, comm.Y_MAX)
                        # re-collision loop for shoot clumping effect
                        while (1):
                            # collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar samping (leave or branch)
                            comm.B_FPR += w * wq * (1.0 - ssa) * bp
                            comm.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            par = w * wq * (1.0 - ssa)
                            comm.AP[ix, iy, iz] += (1.0 - bp) * par
                            comm.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            comm.AP_D[ix, iy, iz] += (1.0 - bp) * min(nscat, 1) * par
                            comm.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            w = RussianRoulette.roulette(w)

                            if (w < MIN_VALUE):
                                self.save(nscat, w)
                                return ERRCODE.LOW_WEIGHT

                            psh = (1.0 - 4.0 * cf)
                            if (randMethod.getRandom() > psh):
                                break

                        # radiance sampling
                        cb = int((1.0 - bp) + 2.0 * bp)
                        vegRadiant.simulate(phoCoord, vectCoord, w, ref, tr, cb, 1.0, fd, ichi, ikd)

                        # new photon direction after scattering
                        if (randMethod.getRandom() > fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_B)
                            if (abs(vectCoord.z) < UZ_MIN):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        # check status
                        if (tsobj == 1):
                            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObjb)
                            io2 = treeBounday.io
                            distance2 = treeBounday.distance

                        # for unexpected case(photon is outside of canopy)
                        if (io2 == 1):
                            break
                    else:
                        phoCoord.movePosition(distance, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break

            # if photon inside the canopy dominant region
            elif (io12 == 0):
                while (1):
                    rand = randMethod.getRandom()
                    # branch scattering or leaf scattering
                    th = acos(vectCoord.z)
                    ith = int(degrees(th))
                    bp = 0.5 + copysign(0.5, comm.BP1 - rand)

                    sgm = comm.GT_BLB[ith] * branchAreaDensity * bp
                    sgm += 4.0 * cf * comm.GT_BLC[ith] * la * (1.0 - bp)
                    sgm= sgm * bp + (sgm / comm.FE) * (1.0 - bp)
                    sgm = max(1.0e-5, sgm)

                    ref = truncRef * bp + lr * (1.0 - bp)
                    tr = lt * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance > distance2):
                        phoCoord.movePosition(distance2, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    elif (distance > distance12):
                        phoCoord.movePosition(distance12, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    else:
                        phoCoord.movePosition(distance, vectCoord, comm.X_MAX, comm.Y_MAX)

                        while (1):
                            #  collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar samping (leave or branch)
                            comm.B_FPR += w * wq * (1.0 - ssa) * bp
                            comm.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            ix = min(ix, comm.SIZE)
                            iy = min(iy, comm.SIZE)

                            par = w * wq * (1.0 - ssa)
                            comm.AP[ix, iy, iz] += par * (1.0 - bp)
                            comm.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            comm.AP_D[ix, iy, iz] += par * (1.0 - bp) * min(nscat, 1)
                            comm.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            psh = (1.0 - 4.0 * cf) * (1.0 - bp)
                            if (randMethod.getRandom() >= psh):
                                break

                        w = RussianRoulette.roulette(w)

                        if (w < MIN_VALUE):
                            self.save(nscat, w)
                            return ERRCODE.LOW_WEIGHT

                        cb = int((1.0 - bp) + 2.0 * bp)
                        vegRadiant.simulate(phoCoord, vectCoord, w, lr, lt, cb, 1.0, fd, ichi, ikd)

                        if (randMethod.getRandom() >= fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_C)
                            if (abs(vectCoord.z < UZ_MIN)):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        if (tsobj == 1):
                            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObj12)
                            io12 = treeBounday.io
                            distance12 = treeBounday.distance
                            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObjb)
                            io2 = treeBounday.io
                            distance2 = treeBounday.distance

                        # for unexpected case (photon is outside of canopy)
                        if (io12 == 1):
                            break

            # if photon inside the canopy dominant region
            elif (io1 == 0):
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

                    ref = truncRef * bp + lr * (1.0 - bp)
                    tr = lt * (1.0 - bp)

                    rand = randMethod.getRandom()
                    distance = -log(rand) / sgm
                    distance = min(distance, 0.9e5)

                    if (distance > distance12):
                        phoCoord.movePosition(distance12, vectCoord, comm.X_MAX, comm.Y_MAX)
                        break
                    elif (distance > distance1):
                        mgn = 1.0e-2
                        phoCoord.x += (distance1 + mgn) * vectCoord.x
                        phoCoord.y += (distance1 + mgn) * vectCoord.y
                        phoCoord.z += (distance1 + mgn) * vectCoord.z
                        self.save(nscat, w)
                        return ERRCODE.OUTSIDE
                    else:
                        phoCoord.movePosition(distance, vectCoord, comm.X_MAX, comm.Y_MAX)

                        while (1):

                            #  collision forcing parameters
                            ssa = 1.0 - (1.0 - ref - tr) * comm.FE
                            ssa *= (1.0 - bp)
                            ssa += (ref + tr) * bp
                            fd = (1.0 - comm.FE) / ssa
                            fd *= (1.0 - bp)

                            # fpar samping (leave or branch)
                            comm.B_FPR += w * wq * (1.0 - ssa) * bp
                            comm.C_FPR += w * wq * (1.0 - ssa) * (1.0 - bp)

                            ix = int(phoCoord.x * comm.RES) + 1
                            iy = int(phoCoord.y * comm.RES) + 1
                            iz = int(phoCoord.z) + 1

                            ix = min(ix, comm.SIZE)
                            iy = min(iy, comm.SIZE)

                            par = w * wq * (1.0 - ssa)
                            comm.AP[ix, iy, iz] += par * (1.0 - bp)
                            comm.AP_B[ix, iy, iz] += (1.0 - bp) * (1.0 - min(nscat, 1))
                            comm.AP_D[ix, iy, iz] += par * (1.0 - bp) * min(nscat, 1)
                            comm.AP_NP[iz] += par * bp

                            w *= ssa
                            nscat += 1

                            psh = (1.0 - 4.0 * cf) * (1.0 - bp)
                            if (randMethod.getRandom() >= psh):
                                break

                        w = RussianRoulette.roulette(w)
                        self.save(nscat, w)

                        if (w < MIN_VALUE):
                            self.save(nscat, w)
                            return ERRCODE.LOW_WEIGHT

                        cb = int((1.0 - bp) + 2.0 * bp)
                        vegRadiant.simulate(phoCoord, vectCoord, w, lr, lt, cb, 1.0, fd, ichi, ikd)

                        if (randMethod.getRandom() >= fd):
                            self.scatterDirection(ref, tr, vectCoord, comm.M_C)
                            if (vectCoord.z < UZ_MIN):
                                vectCoord.z = copysign(UZ_MIN, vectCoord.z)

                        # check status
                        if (tsobj == 1):
                            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObj)
                            distance1 = treeBounday.distance
                            io1 = treeBounday.io
                            treeBounday.dealTreeType(tsobj, phoCoord, vectCoord, tObj12)
                            distance12 = treeBounday.distance
                            io12 = treeBounday.io
                        if (io1 == 1):
                            self.save(nscat, w)
                            return ERRCODE.OUTSIDE

        self.save(nscat, w)
        return ERRCODE.SUCCESS

    # ##################################################################################
    # This subroutine calculates photon trajectory in the forest floor
    # rand     : random numbers
    # distance : traveling distance
    # sgm      : extinction
    # w        : weight of the photon
    # z (m)    : Height of the photon position from the upper grass boundary (-1.0<=z<0)
    # ##################################################################################S
    def floor(self, w, wq, phoCoord, vectCoord, nscat, ulr, ult, sor, ichi, ikd):

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

        # Monte Carlo loop
        while (1):

            rand = randMethod.getRandom()
            th = acos(vectCoord.z)
            ith = degrees(th)
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
                    comm.F_FPR += w * wq * (1.0 - ulr - ult)

                    ix = int(phoCoord.x * comm.RES) + 1
                    iy = int(phoCoord.y * comm.RES) + 1

                    ix = min(ix, comm.SIZE)
                    iy = min(iy, comm.SIZE)

                    comm.AP_F[ix, iy] += w * wq * (1.0 - ult - ulr)
                    comm.AP_FD[ix, iy] += w * wq * (1.0 - ult - ulr) * min(float(nscat), 1.0)

                    w *= (ult + ulr)
                    nscat += 1

                    w = RussianRoulette.roulette(w)
                    if (w < MIN_VALUE):
                        self.save(nscat, w)
                        return ERRCODE.LOW_WEIGHT

                    psh = 1.0 - 4.0 * cf
                    if (randMethod.getRandom() >= psh):
                        break

                vegRadiant.simulate(phoCoord, vectCoord, w, ulr, ult, 4, 1.0, fd, ichi, ikd)

                # new direction
                self.scatterDirection(ulr, ult, vectCoord, comm.M_F)
                if (abs(vectCoord.z) < MIN_Z):
                    vectCoord.z = copysign(MIN_Z, vectCoord.z)

            else:
                phoCoord.x += (distancePlane + mgn) * vectCoord.z
                phoCoord.y += (distancePlane + mgn) * vectCoord.y
                phoCoord.z += distancePlane * vectCoord.z
                phoCoord.x -= (trunc(phoCoord.x / comm.X_MAX) - 0.5 + copysign(0.5, phoCoord.x)) * comm.X_MAX
                phoCoord.y -= (trunc(phoCoord.y / comm.Y_MAX) - 0.5 + copysign(0.5, phoCoord.y)) * comm.Y_MAX

                if (phoCoord.z <= zBottom):
                    # surface boundary reflectance mode 1: Lambertian, 2: RPV model, 3: DSM model
                    # fpar sampling
                    comm.S_FPR += w * wq * (1.0 - sor)

                    ix = int(phoCoord.x * comm.RES) + 1
                    iy = int(phoCoord.y * comm.RES) + 1
                    ix = min(ix, comm.SIZE)
                    iy = min(iy, comm.SIZE)

                    comm.AP_S += w * wq * (1.0 - ulr - ulr)

                    # surface downward flux
                    comm.SF_DIR[ix, iy] += w * wq * (1.0 - min(nscat, 1))
                    comm.SF_DIF[ix, iy] += w * wq * min(nscat, 1)
                    w *= sor

                    if (w < MIN_VALUE):
                        self.save(nscat, w)
                        return ERRCODE.LOW_WEIGHT

                    vegRadiant.simulate(phoCoord, vectCoord, w, ulr, ult, 5, 1.0, fd, ichi, ikd)

                    self.save(nscat, w)
                    self.scatterReflectSurface(vectCoord, mode)
                    nscat = self.cNscat
                    w = self.weight

                    phoCoord.z = zBottom + mgn
                elif (phoCoord.z >= zUpper):
                    self.save(nscat, w)
                    return ERRCODE.OUT_OF_RANGE

        self.save(nscat, w)
        return ERRCODE.SUCCESS


# ##############################
# FOR TEST
# Monte Carlo in Calculating PI
# ##############################
#
# def test():
#     R = 0.5
#     inside = 0
#     allNum = 10000
#     for i in range(allNum):
#         rX = random.uniform(0, 1)
#         rY = random.uniform(0, 1)
#         rX = (rX - R) * (rX - R)
#         rY = (rY - R) * (rY - R)
#         rNEW = math.sqrt(rX + rY)
#
#         if (rNEW <= R):
#             inside += 1
#
#     PI = 4 * inside / allNum
#     print("PI = ", PI)
#     print("ERR = ", abs((math.pi - PI) / math.pi))
# a = [1,2,3,4,5,6,7]
# print(a)
# a.remove(3)
# print(a)
# print(a[1 + 9])
