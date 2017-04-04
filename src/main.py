import common as comm
from iparam import Parameters
import numpy as np
import ERRCODE
from G_Function import G_Function
from ipf import ipf
from idivspace import idivspace
from Clai import CLAI
import datetime
from Random import Random
from math import *
from Position import Position
from CanopyPhotonTrace import CanopyPhotonTrace
from MonteCarlo_1D import MonteCarlo_1D
import input_parameters

# only in main
knmix = 10
knang = 2000
knzext = 200
fpv = [0.0] * 101
fpc = fpf = 0.0
ip = iz = num = 0

plai = [0.0] * 100
wkd = [0.0] * comm.K_NKD
ang = [0.0] * knang
re = [0.0] * knmix
Qext_ref = [0.0] * knmix
Qext = [0.0] * knmix
G = [0.0] * knmix
Qabs = [0.0] * knmix
omg = [0.0] * knmix
phs = np.zeros(knmix * knang, dtype=float).reshape(knmix, knang)
ext = np.zeros(knmix * knzext, dtype=float).reshape(knmix, knzext)

MIN_VALUE = 1.0e-8
MIN_UZ = 0.0174524

# photon data
nscat = nscata = ichi = ikd = 0
w = 0.0
PhotonCoord = Position()
VectorCoord = Position()

# function parameters
# real fsin, fcos, facos, r_acos,frnd

Nid = 0
Nprod = 1  # comm.T_SIN[-1] = 5
# print(comm.T_SIN)

rand = Random()


def simulateATM(para):
    global PhotonCoord, VectorCoord
    global nscat, nscata
    global ichi
    global rand

    mc1D = MonteCarlo_1D()
    canopyTrace = CanopyPhotonTrace()

    #print("simulate without atmosphere")

    for iwl in range(1, para.nwl):
        print("This is the ", iwl, " of ", para.nwl, " wavelength.")
        para.preparAtm(iwl)
        mc1D.optics(para.ext, para.omg, para.phs, para.ang, para.nmix + para.cflg)
        wq = para.wq

        if (para.npl[iwl] == 0):
            print("NPL ERROR")
            return ERRCODE.FAILURE

        # loop for single wavelength
        for ip in range(1, para.npl[iwl] + 1):
            if (fmod(ip, para.nPhoton) == 0):
                print(str(ip) + " of " + str(para.npl[iwl]))

            w = 1.0
            PhotonCoord.setPosition(comm.X_MAX * float(rand.getRandom()),
                                    comm.Y_MAX * float(rand.getRandom()),
                                    comm.Z_GRD[comm.N_Z])
            VectorCoord.setPosition(para.sin_q0 * para.cos_f0,
                                    para.sin_q0 * para.sin_f0,
                                    para.cos_q0)

            ftau = -log(max(1.0e-35, float(rand.getRandom())))
            chi = 1.0
            ichi = 1
            iz = comm.N_Z

            # determination of the k-term
            randNum = rand.getRandom()
            ikd = 1
            while (randNum > para.wkd[ikd]):
                ikd += 1

            # monte carlo core loop
            while (1):
                mc1D.trace(PhotonCoord, VectorCoord, w, wq, ftau, chi, ikd, iz, nscat, ichi)
                # load
                w = mc1D.weight
                nscat = mc1D.nscat
                iz = mc1D.iz

                nscata = nscat
                if ((w <= 0.0) or (iz > comm.N_Z)):
                    break

                comm.scmpf[1, iwl] += w
                comm.scmpf[2, iwl] += w * (1.0 - min(float(nscat), 1.0))
                comm.scmpf[3, iwl] += w * min(float(nscat), 1.0)

                comm.scmpp[1, iwl] += w * wq
                comm.scmpp[2, iwl] += w * wq * (1.0 - min(float(nscat), 1.0))
                comm.scmpp[3, iwl] += w * wq * min(float(nscat), 1.0)

                # surface reflection
                if (para.surfaceType == 1):
                    # lambertian
                    w *= para.alb[1]
                    if (w < MIN_VALUE):
                        break
                    iz = 1

                    th = 0.5 * acos(0.01 * rand.getRandom())
                    ph = 2 * pi * rand.getRandom()

                    VectorCoord.setPosition(sin(th) * cos(ph),
                                            sin(th) * sin(ph),
                                            cos(th))
                    if (abs(VectorCoord.z) < 0.0174524):
                        VectorCoord.z = copysign(0.0174524, VectorCoord.z)

                    para.rflx += w
                    para.rbflx += w * (1.0 - min(float(nscata), 1.0))
                    para.rdflx += w * min(float(nscata), 1.0)

                else:
                    # 3-D surface
                    tLR = [0.0] * (para.nts + 1)
                    tLT = [0.0] * (para.nts + 1)
                    tSTR = [0.0] * (para.nts + 1)
                    for i in range(para.nts):
                        tLR[i] = para.lr[i, 1]
                        tLT[i] = para.lt[i, 1]
                        tSTR[i] = para.truncRef[i, 1]

                    # call the canopy radiation transfer module
                    canopyTrace.trace(PhotonCoord, VectorCoord, w, para.wq, nscat, ichi, ikd, tSTR, para.soilRef[1],
                                      tLR, tLT, para.ulr[1], para.ult[1])

                    w = canopyTrace.weight
                    nscat = canopyTrace.cNscat

                    if (w < MIN_VALUE):
                        break

                    iz = 1
                    para.rflx += w
                    para.rbflx += w * (1 - min(float(nscata), 1.0))
                    para.rdflx += w * min(float(nscata), 1.0)

        # summary of spectral flx
        para.tflx += comm.scmpf[1, iwl]
        para.bflx += comm.scmpf[2, iwl]
        para.dflx += comm.scmpf[3, iwl]

        para.tpfd += comm.scmpp[1, iwl]
        para.bpfd += comm.scmpp[2, iwl]
        para.dpfd += comm.scmpp[3, iwl]

    return ERRCODE.SUCCESS


def simulateNoATM(para):

    global PhotonCoord, VectorCoord
    global nscat, nscata
    global ichi
    global rand

    #global th, ph

    canopyTrace = CanopyPhotonTrace()

    for iPhoton in range(para.nPhotonProcess):

        print("Current photon: ", iPhoton + 1)

        w = 1.0     # initial weight of photon
        ikd = 0     # initialization of CDK (correlated k-dist)

        #  selection of beam or diffuse flux
        if (rand.getRandom() > para.diffuse):
            # beam
            VectorCoord.setPosition(para.sin_q0 * para.cos_f0,
                                    para.sin_q0 * para.sin_f0,
                                    para.cos_q0)

            if (abs(VectorCoord.z) < MIN_UZ):
                VectorCoord.z = copysign(MIN_UZ, VectorCoord.z)

            nscat = 0
            nscata = nscat

        else:
            # diffuse
            th = 0.5 * pi + 0.5 * acos(1.0 - 2.0 * rand.getRandom())
            ph = 2.0 * pi * rand.getRandom()

            VectorCoord.setPosition(sin(th) * cos(ph),
                                    sin(th) * sin(ph),
                                    cos(th))

            if (abs(VectorCoord.z) < MIN_UZ):
                VectorCoord.z = copysign(MIN_UZ, VectorCoord.z)

            nscat = 1
            nscata = nscat

        # initial position (x, y)
        PhotonCoord.setPosition(comm.X_MAX * rand.getRandom(),
                                comm.Y_MAX * rand.getRandom(),
                                0)
        print("Initial potion x = ", PhotonCoord.x, ", y = ", PhotonCoord.y)

        para.tflx += w
        para.tpfd += w * para.wq

        para.bflx += w * (1.0 - min(float(nscat), 1.0))
        para.bpfd += w * para.wq * (1.0 - min(float(nscat), 1.0))
        para.dflx += w * (1.0 - min(float(nscat), 1.0))
        para.dpfd += w * para.wq * (1.0 - min(float(nscat), 1.0))

        # surface reflection
        if (para.surfaceType == 1):
            # lambertian
            w *= para.alb[1]
            th = 0.5 * acos(1.0 - 2.0 * rand.getRandom())
            ph = 2 * pi * rand.getRandom()

            VectorCoord.setPosition(sin(th) * cos(ph),
                                    sin(th) * sin(ph),
                                    cos(th))

            if (abs(VectorCoord.z) < MIN_VALUE):
                VectorCoord.z = copysign(MIN_VALUE, VectorCoord.z)
        else:
            # 3D surface
            tLR = [0.0] * (para.nts + 1)
            tLT = [0.0] * (para.nts + 1)
            tSTR = [0.0] * (para.nts + 1)
            for i in range(para.nts):
                tLR[i] = para.lr[i, 1]
                tLT[i] = para.lt[i, 1]
                tSTR[i] = para.truncRef[i, 1]

            # call the canopy radiation transfer module
            canopyTrace.trace(PhotonCoord, VectorCoord, w, para.wq, nscat, ichi, ikd, tSTR, para.soilRef[1],
                              tLR, tLT, para.ulr[1], para.ult[1])

            w = canopyTrace.weight
            nscat = canopyTrace.cNscat

            para.rflx += w
            para.rbflx += w * (1 - min(float(nscata), 1.0))
            para.rdflx += w * min(float(nscata), 1.0)

    return ERRCODE.SUCCESS


def main(**args):
    if (Nid == 0):
        print("*********************************************")
        print("3D canopy-atmosphere radiative transfer model")
        print("Forest Light Environmental Simulator (FLiES) ")
        print("                         by Hideki Kobayashi ")
        print("        Special thanks to Hironobu Iwabuchi  ")
        print("                      Since Ver2.00 2008/5/1 ")
        print("*********************************************")

    # ---- Preparation of the initial condition -------
    # set by iparam
    print("iparam")
    print("Parameters initialization...")
    para = Parameters()
    errCode = para.readParameters(**args)

    if (errCode):
        return errCode

    # Initialize math function
    print("imath")

    # number the initial file read divide (5 for standard input)
    num = 5

    # Read required parameters

    # igtbl
    if (para.surfaceType == 2):
        print("G-function")     # G-function LUT
        gFunction = G_Function()
        errCode = gFunction.igtbl()

        if errCode:
            return errCode

        print("idivspace")      # Initialize space division
        errCode = idivspace()

        print("ipf")            # LUT for phase function
        ipf()

    # fish eye image at the forest floor
    # call local estimation for fish eye image
    if (para.nPhoton == -4):
        tgx = int(input("Input the x, y position of view\n"))
        tgy = int(input())

        # to be continue, not important

        return ERRCODE.SUCCESS

    #  LAI calculation
    if (para.nPhoton == -5):
        # spn = float(input("input sampling grid (m), suggested value is 0.1 (m)\n"))
        # debug
        spn = 0.1
        cLAI = CLAI()
        cLAI.calLAI(spn)

        print("LAI = ", cLAI.LAI)
        print("Crown cover = ", cLAI.crownCover)
        print("pLAI is")
        print(cLAI.pLAI)
        print("Finish print LAI.")

    # ---- start simulation --------

    print("start simulation...")

    # #################################
    # Without atmosphere
    # #################################
    if (para.AtmMode == 2):
        print("without atmosphere simulation.")
        errCode = simulateNoATM(para)
        if (errCode != ERRCODE.SUCCESS):
            return errCode

    # #################################
    # With atmosphere
    # #################################
    else:
        print("with atmosphere simulation.")
        errCode = simulateATM(para)
        if (errCode != ERRCODE.SUCCESS):
            return errCode

    print("end simulation.")

    # ---- end simulation --------
    # Output
    errCode = para.writeData()
    print("Start to write results...")
    return errCode


# ##################################################################
# # errCode = main()
# np.set_printoptions(threshold=1000)
# para = Parameters()
# gFunction = G_Function()
# errCode = gFunction.igtbl(para)
# para.nts = 2
# para.process202()
# errCode = idivspace.idivspace()
# print(comm.Z_MAX)
# ERRCODE.printMessage(errCode)
# import math
# time.sleep(2)

# a = np.zeros((6 + 1) * (6 + 1), dtype=float).reshape((6 + 1), (6 + 1))
#
# print(a[1][1])
# print(a[1, 1])

# #################################################################
# TEST
# #################################################################
START_TIME = datetime.datetime.now()

err = main(**input_parameters.Atmosphere_Mode_Args)
ERRCODE.printMessage(err)

END_TIME = datetime.datetime.now()
print("TIME USED (HOURS, MINUTES, SECONDS): ", (END_TIME - START_TIME))



