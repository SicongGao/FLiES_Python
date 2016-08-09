# Parameters Initialization

import common
import math
import numpy as np
import ERRCODE

class Parameters:
    __ix = [0] * 1025
    wl0 = 0.0  # Default target wavelength(micron)
    wls = 0.0  # Default start point of wavelength integration
    nwl = 0  # Default sampling number of the wavelength
    taur = 0.0  # default AOT(tau at 550 nm)
    ctaur = 0.0  # default COT(tau at 550 nm)
    d = 0.0  # default scale height(m)
    atmType = 0  # default atmospheric type
    aerosolType = 0  # default aerosol type
    cloudType = 0
    nkd = 0  # of subbands for the k-distribution
    th0 = 0.0  # solar zenith angle in degree
    ph0 = 0.0  # solar azimuth angle in degree
    tm = 0.0  # atmopsheric transmittance
    hfov = 0.0  # detector half widt of the field of view
    cflg = 0  # cflg = 0 -> cloud - free, cflg = 1 -> cloud
    tflx = 0.0  # Total downward flux at top of canopy
    bflx = 0.0  # Beam downward flux at TOC
    dflx = 0.0  # Diffuse downward flux at TOC
    obflx = 0.0  # Beam downawr flux at TOC(observed)
    odflx = 0.0  # Diffuse downward flux at TOC(observed)
    rflx = 0.0  # total reflcted irradiaice
    rbflx = 0.0  # beam reflected irradiance
    rdflx = 0.0  # diffuse reflected irradiance
    tpfd = 0.0  # Total downward PFD at top of canopy
    bpfd = 0.0  # Beam downward PFD at top of canopy
    dpfd = 0.0  # Diffuse downward PDF at top of canopy
    obpfd = 0.0  # Beam downawr PFD at TOC(observed)
    odpfd = 0.0  # Diffuse downward PFD at TOC(observed)
    cbnz = 0  # Number of cloud bottom layer
    ctnz = 0  # Number of cloud top layer
    RF = 0.0  # incident irradiance and photon flux density at TOA
    RQ = 0.0  # incident irradiance and photon flux density at TOC
    fpc = 0.0  # fpar for canopy floor
    fpf = 0.0  # fpar for forest floor
    span = [0] * 101  # Wavelength span
    Nprocess = 1  # number of processor
    nPhotonProcess = 1

    # read parameters
    nPhoton = nmix = amode = imode = cmode = surfaceType = 0
    dif = phi = th = ph = tgx = tgy = 0.0
    wq = sinf0 = cosf0 = cosq0 = sinq0 = 0

    npl = [0] * 200
    alb = [0.0] * 100
    spcf = [0.0] * 400
    spcq = [0.0] * 400
    ulr = [0.0] * 100
    ult = [0.0] * 100
    sor = [0.0] * 100
    tlr = [0.0] * 5
    tlt = [0.0] * 5
    tstr = [0.0] * 5

    lr = np.zeros(5 * 100, dtype=float).reshape(5, 100)
    lt = np.zeros(5 * 100, dtype=float).reshape(5, 100)
    str = np.zeros(5 * 100, dtype=float).reshape(5, 100)

    fname = []
    rfname = ""

    nts = 0  # group of tree species
    bound = 1

    # for math
    T_SIN = [0.0] * common.ANGLE_SHIFT * 2
    T_COS = [1.0] * common.ANGLE_SHIFT * 2
    T_ACOS = [0.0] * common.ACOOS_SHIFT * 2
    T_EXP = [0.0] * common.ACOOS_SHIFT
    DLT = np.zeros((6 + 1) * (6 + 1), dtype=float).reshape((6 + 1), (6 + 1))

    def __init__(self):
        self.initParameters()

    def initParameters(self):

        common.N_Z = 12  # of layers
        common.X_MAX = common.Y_MAX = 30.0  # X domain size (m), Y domain size (m)
        common.RES = common.SIZE / common.X_MAX  # inverse of the spatial resolution
        common.WRR = 0.1  # ideal weight (used for Russian roulette in atm.

        # branch portion in region 1(outer canopy) and 2(internal canopy)
        # 1.0 = 100 % branch, 0.0: 100 % leaf
        common.bp1 = 0.0
        common.bp2 = 1.0

        # leaf angle distribution of canopy, branches, and floor
        common.M_C = common.M_B = common.M_F = 1

        # #   ---  local parameters ---
        self.wl0 = 0.55  # Default target wavelength(micron)
        self.wls = 0.3  # Default start point of wavelength integration
        self.nwl = 75  # Default sampling number of the wavelength
        self.taur = 0.3  # default AOT(tau at 550 nm)
        self.ctaur = 0.00  # default COT(tau at 550 nm)
        self.d = 8000.0  # default scale height(m)
        self.atmType = 1  # default atmospheric and aerosol type
        self.aerosolType = 2
        self.nkd = 3  # of subbands for the k-distribution
        self.th0 = 10.0  # solar zenith angle in degree
        self.ph0 = 0.0  # solar azimuth angle in degree
        self.tm = 1.0  # atmopsheric transmittance
        self.hfov = 1.0  # detector half widt of the field of view
        self.hfov = math.radians(self.hfov)
        self.hfov = math.cos(self.hfov)
        self.cflg = 0  # cflg = 0 -> cloud - free, cflg = 1 -> cloud

        self.span[1: 20] = [0.02] * 20  # Wavelength span (0.3-0.7 micron)
        self.span[21: 100] = [0.1] * 80  # Wavelength span (0.7- micron)

        #     Irradiance
        self.tflx = 0.0  # Total downward flux at top of canopy
        self.bflx = 0.0  # Beam downward flux at TOC
        self.dflx = 0.0  # Diffuse downward flux at TOC
        self.obflx = 0.0  # Beam downawr flux at TOC(observed)
        self.odflx = 0.0  # Diffuse downward flux at TOC(observed)

        #     reflected radiance from canopy
        self.rflx = 0.0  # total reflcted irradiaice
        self.rbflx = 0.0  # beam reflected irradiance
        self.rdflx = 0.0  # diffuse reflected irradiance

        #     Photon flux density
        self.tpfd = 0.0  # Total downward PFD at top of canopy
        self.bpfd = 0.0  # Beam downward PFD at top of canopy
        self.dpfd = 0.0  # Diffuse downward PDF at top of canopy
        self.obpfd = 0.0  # Beam downawr PFD at TOC(observed)
        self.odpfd = 0.0  # Diffuse downward PFD at TOC(observed)

        #     Number of cloud bottom and top layer
        self.cbnz = 0
        self.ctnz = 1

        #     incident irradiance and photon flux density at TOA or TOC
        self.RF = 0.0
        self.RQ = 0.0

        #     initialization of random number generator
        # flg = frndi(ix(Nid + 1))
        # self.flg = 0

        self.initMathparameters()

        return ERRCODE.SUCCESS

    def initMathparameters(self):
        for i in range(0, 62832 * 2):
            self.T_SIN[i] = math.sin(float(i - 62832) * 0.0001)
            self.T_COS[i] = math.cos(float(i - 62832) * 0.0001)

        for i in range(0, common.ACOOS_SHIFT * 2):
            self.T_ACOS[i] = math.acos(float(i - common.ACOOS_SHIFT) * 0.0001)

        for i in range(1, 6):
            common.DLT[i, i] = 1.0

        return ERRCODE.SUCCESS

    def getInputLR_LT_ULR_ULT(self, i):
        for j in range(self.nts):
            self.lr[j, i] = float(input())
        for j in range(self.nts):
            self.lt[j, i] = float(input())
        self.ulr[i] = float(input())
        self.ult[i] = float(input())
        for j in range(self.nts):
            self.str[j] = float(input())
        self.sor[i] = float(input())
        return ERRCODE.SUCCESS

    def process100(self):

        if (self.cloudType >= 2):
            common.N_ANG_C = 1
            common.UX_RC[0] = 0.0
            common.UY_RC[0] = 0.0
            common.UZ_RC[0] = 1.0

        self.process201()

        return ERRCODE.SUCCESS

    def process202(self):
        umax = 0.0

        print("u: leaf area density 1,2,3...# tree species")
        for i in range(self.nts):
            common.U[i] = float(input())

        if (not((self.nPhoton == -4) or (self.nPhoton == -5))):
            common.G_LAI = float(input("gLAI: forest floor LAI\n"))

        print("BAD: branch area density 1,2,3... # of tree species")
        for i in range(self.nts):
            common.BAD[i] = float(input())

        for i in range(self.nts):
            if ((common.U[i] < 0.0) or (common.U[i] > 8.0)):
                print(str(i) + "th leaf area density " + str(common.U[i]) + " should be set in the range (0.0-8.0)")
                print("EXIT")
                return ERRCODE.INPUT_ERROR

        if ((common.G_LAI < 0.0) or (common.G_LAI > 8.0)):
            print(str(common.G_LAI) + " should be set in the range (0.0-8.0)")
            print("EXIT")
            return ERRCODE.INPUT_ERROR

        if (self.nPhoton == -5):
            print("sbar: Spherical ave. shoot silhouette to total needle area ratio")
            print("1,2,3... # of tree species (0.0-0.25)")
            print("For broadleaves, please input 0.25")
            for i in range(self.nts):
                common.S_BAR = float(input())

        umax = common.U[0]
        for i in range(1, self.nts):
            if (umax < common.U[i]):
                umax = common.U[i]

        if (umax > 1.0):
            common.FE = 1.0
        elif (umax <= 0.01):
            common.FE = 0.01
        else:
            common.FE = umax

        # canopy object parameters
        # object id initialization
        common.I_OBJ = [1] * 6000

        # input obj_nt
        result = []
        with open("../data/crowndata.txt", "r") as file:
            common.N_OBJ = int(file.readline())

            if (common.N_OBJ == 0):
                common.N_OBJ = common.S_OBJ = 1
                common.OBJ[1, 1] = 0.01
                common.OBJ[1, 2] = 0.01
                common.OBJ[1, 3] = 0.01
                common.OBJ[1, 4] = 1E-5
                common.OBJ[1, 5] = 1E-5
            else:
                result = []
                result = file.readlines()
                result = np.loadtxt(result)
        file.close()

        obj_nt = common.N_OBJ
        for i in range(len(result)):
            common.S_OBJ[i] = result[i][0]
            common.OBJ[i][0:5] = result[i][1: 6]
            common.I_OBJ[i] = result[i][6]
            if ( result[i][0] != 4):
                if ((result[i][4] < 0.01) or (result[i][5] < 0.01)):
                    print(str(i + 1) + "th canopy neglected!")
                    obj_nt -= 1

        common.N_OBJ = obj_nt

        # check the obj id range
        for i in range(obj_nt):
            if ((common.I_OBJ[i] <= 0) or (common.I_OBJ[i] > self.nts)):
                print("species id should be the range betweem 0-" + str(self.nts))
            return ERRCODE.OUT_OF_RANGE

        print("Total Object is " + str(obj_nt))

        # change from height to radius
        for i in range(obj_nt):
            if (common.S_OBJ[i] == 3):
                common.OBJ[i, 3] /= 2

        # in case periodic boudary
        # add objects that are partially in the outside from the simulation scence
        a = [1.0, 1.0, 0.0, 0.0]
        b = [0.0, 0.0, 1.0, 1.0]
        # preparation of the i-th grid for space divided method
        c = [0.0, common.X_MAX, 0.0, common.Y_MAX]
        cc = [0.0, 1.0, 0.0, 1.0]

        if (self.bound == 1):
            # set distance calculation parameters
            kobj = common.N_OBJ

            # definition of rectangular of the i-th object
            for j in (kobj):
                xr = common.OBJ[j, 0]
                yr = common.OBJ[j, 1]

                # check the intersection on the x - y plane
                for k in range(4):
                    d = abs(xr * a[k] + yr * b[k] - c[k])

                    if (d <= common.OBJ[j, 4]):
                        dd = math.sqrt(common.OBJ[j, 4] * common.OBJ[j, 4] - d * d)
                        p1 = b[k] * xr + a[k] * yr - dd
                        p2 = b[k] * xr + a[k] * yr + dd
                        min = 0.0
                        max = common.X_MAX

                        if (((min < p1) and (max > p1)) or
                                ((min < p2) and (max > p2))):
                            common.N_OBJ += 1

                            n = common.N_OBJ
                            common.OBJ[n, 0] = common.OBJ[j, 0] + a[k] * common.X_MAX * (1.0 - 2.0 * cc[k])
                            common.OBJ[n, 1] = common.OBJ[j, 1] + a[k] * common.Y_MAX * (1.0 - 2.0 * cc[k])
                            common.OBJ[n, 2] = common.OBJ[j, 2]
                            common.OBJ[n, 3] = common.OBJ[j, 3]
                            common.OBJ[n, 4] = common.OBJ[j, 4]

                            common.S_OBJ[n] = common.S_OBJ[j]
                            common.I_OBJ[n] = common.I_OBJ[j]

                for k in range(2):
                    for l in range(2):
                        rx = xr - c[k]
                        ry = yr - c[l + 2]
                        rr = math.sqrt(rx * rx + ry * ry)

                        if (rr <= common.OBJ[j, 4]):
                            common.N_OBJ += 1

                            n = common.N_OBJ
                            common.OBJ[n, 0] = common.OBJ[j, 0] + common.X_MAX - 2.0 * cc[k]
                            common.OBJ[n, 1] = common.OBJ[j, 1] + common.Y_MAX - 2.0 * cc[l + 2]
                            common.OBJ[n, 2] = common.OBJ[j, 2]
                            common.OBJ[n, 3] = common.OBJ[j, 3]
                            common.OBJ[n, 4] = common.OBJ[j, 4]

                            common.S_OBJ[n] = common.S_OBJ[j]
                            common.I_OBJ[n] = common.I_OBJ[j]

            # determination of the epsi(epsiron) that is used to perform the
            # Russian Roulette
            # epsi is setted to be c * (leaf_r + leaf_t)
            if (self.nPhoton == -4):
                avelr = 0.0
                avelt = 0.0

                for i in range(self.nwl):
                    for j in range(self.nts):
                        avelr += self.lr[j, i]
                        avelt += self.lt[j, i]

                avelr /= float(self.nwl * self.nts)
                avelt /= float(self.nwl * self.nts)

        return ERRCODE.SUCCESS

    def process201(self):
        ispc = 0
        # get the number of forest species
        self.nts = int(input("nts: # of group of tree species\n"))
        if ((self.nPhoton == -4) or (self.nPhoton == -5)):
            return self.process202()

        common.U[5] = float(self.nts)

        # currently max nts should be less than 5
        if ((self.nts <= 0) or (self.nts >= 5)):
            print("Error # of tree species")
            print("tree species should be 1-5")
            return ERRCODE.INPUT_ERROR

        # read leaf reflectance/transmittance
        ramda = self.wls
        if(self.nPhoton == -4):
            self.nwl = 1

        for i in range(self.nwl + 1):
            if (self.nwl == 1):
                print("Input surface optical properties")
            else:
                print("Input surface optical properties in")
                if (i <= 20):
                    print(str(ramda) + " and " + str(ramda + self.span[1]))
                else:
                    print(str(ramda) + " and " + str(ramda + 5 * self.span[1]))

            print("lr1 lr2.. lt1 lt2.. ulr ult str1 str2.. sor")
            print("(lr, lt:leaf refl. & leaf transm.")
            print("ulr,ult:understory leaf refl. & transm.")
            print("stmr: stem refl., soilr: soil refl.)")

            if (self.nwl == 1):
                self.getInputLR_LT_ULR_ULT(i)

            else:
                if (i == 1):
                    print("Input PAR average values:")
                    self.getInputLR_LT_ULR_ULT(i)
                    ispc = i

                elif (i == 21):
                    print("Input NIR average values")
                    self.getInputLR_LT_ULR_ULT(i)
                    ispc = i

                else:
                    for j in range(self.nts):
                        self.lr[j, i] = self.lr[j, ispc]
                        self.lt[j, i] = self.lt[j, ispc]
                        self.str[j, i] = self.str[j, ispc]

                    self.ulr[i] = self.ult[ispc]
                    self.ult[i] = self.ult[ispc]
                    self.sor[i] = self.sor[ispc]

            # check the parameter range
            for j in range(self.nts):
                if (self.lr[j,i] + self.lt[j,i] > 0.99):
                    print("canopy leaf reflectance+transmittance is too large, exit!")
                    return ERRCODE.OUT_OF_RANGE

                if (self.lr[j, i] + self.lt[j, i] < 0.0001):
                    self.lr[j,i] = 0.0001
                    self.lt[j,i] = 0.0001

            if (self.ulr[i] + self.ult[i] > 0.99):
                print("floor leaf reflectance+transmittance is too large, exit!")
                return ERRCODE.OUT_OF_RANGE

            if (self.ulr[i] + self.ult[i] < 0.0001):
                self.ulr[i] = 0.0001
                self.ult[i] = 0.0001

            for j in range(self.nts):
                if (self.str[j, i] > 1.00):
                    print("stem reflectance is too large, exit!")
                    return ERRCODE.OUT_OF_RANGE

                if (self.str[j, i] < 0.0001):
                    self.str[j, i] = 0.0001

            if (self.sor[i] > 1.0):
                print("soil reflectance is too large, exit!")
                return ERRCODE.OUT_OF_RANGE

            if (self.nwl <= 20):
                ramda += self.span[1]
            else:
                ramda += 5 * self.span[i]

        self.process202()

        return ERRCODE.SUCCESS

    def readVegParameters(self):

        if ((self.nPhoton == -4) or (self.nPhoton == -5)):
            self.process201()

        # input output mode
        print("cmode: calculation mode")
        self.cloudType = int(input("1: BRF only 2: BRF Nadir Image 3: 3D APAR\n"))

        if ((self.cloudType <= 0) or (self.cloudType >= 4)):
            print("Bad mode selection exit")
            return ERRCODE.INPUT_ERROR

        # read the condition file
        if ((self.bound != 1) and (self.bound != 2)):
            print("boundary condition should be 1 or 2.")
            print("  1: Periodic  2: Non-periodic")
            return ERRCODE.INPUT_ERROR

        if (self.cloudType != 1):
            self.process100()

        print("nth, angt: # of angle anfinished process 201 200d zenith angle(max 18)for BRF")
        print("eg. 5 10. 20. 45. 50. 70.")
        for i in range(common.N_TH):
            common.ANG_T[i] = float(input())

        print("nph, angp:# of angle and azimuth angle(max 36)for BRF")
        print("eg. 3 0. 90. 180.")
        for i in range(common.N_TH):
            common.ANG_P[i] = float(input())

        if (common.N_TH * common.N_PH > 700):
            print("The number of sampling angles are too huge! It should be theta*phi<348")
            return ERRCODE.INPUT_ERROR

        if (common.N_TH > 18):
            print("The number of sampling theta are too huge! It should be theta*phi<100")
            return ERRCODE.INPUT_ERROR

        if (common.N_PH > 36):
            print("The number of sampling phi are too huge! It should be theta*phi<100")
            return ERRCODE.INPUT_ERROR

        k = 0

        for i in range(common.N_PH):
            for j in range(common.N_TH):
                if (common.ANG_T[j] > 80.0):
                    print("Zenith angle should be less than 80")
                    print(str(common.ANG_T[j]) + " is ignored !")
                else:
                    common.UX_RC[k] = math.sin(math.radians(common.ANG_T[j])) * math.cos(math.radians(common.ANG_P[i]))
                    common.UY_RC[k] = math.sin(math.radians(common.ANG_T[j])) * math.sin(math.radians(common.ANG_P[i]))
                    common.UZ_RC[k] = math.cos(math.radians(common.ANG_T[j]))

                    k += 1

        common.N_ANG_C = k

        self.process100()

        return ERRCODE.SUCCESS

    def readGeoParameters(self):

        f0 = q0 = 0.0
        ntha = npha = 0
        th = [0] * 18
        ph = [0] * 36

        # solar zenith angle
        self.th0 = float(input("the0: Solar zenith angle(degree)\n"))
        q0 = math.pi - math.radians(self.th0)

        if (self.ph0 <= math.pi):
            f0 = math.radians(self.ph0) + math.pi
        else:
            f0 = math.radians(self.ph0) - math.pi

        sin_q0 = math.sin(q0)
        cos_q0 = math.cos(q0)
        sin_fo = math.sin(f0)
        cos_fo = math.cos(f0)

        if (self.amode == 2):
            return

        # radiance smapling angle
        print("Radiance at the bottom atmospheric boundary")
        print("ntha, th: # of angle and zenith angle(degree,max 18)")
        print("ex. 5 100. 120. 140. 160. 170.")
        ntha = int(input())
        for i in range(ntha):
            th[i] = int(input())

        print("npha, ph: # of angle and azimuth angle(degree,max 36)")
        print("ex. 3 0. 90. 180.")
        npha = int(input())
        for i in range(npha):
            ph[i] = int(input())

        k = 0
        for i in range(1,npha):
            for j in range(1, ntha):
                if (th[j] < 91):
                    print("Zenith angle should be greater than 91")
                    print(str(th[j]) + " is ignored !")
                else:
                    k = k + 1
                    common.UX_RTAB[k] = math.sin(math.radians(th[j])) * math.cos(math.radians(ph[i]))
                    common.UY_RTAB[k] = math.sin(math.radians(th[j])) * math.sin(math.radians(ph[i]))
                    common.UX_RTAB[k] = math.cos(math.radians(th[j]))

        common.N_RDC = k
        return ERRCODE.SUCCESS

    def readAtmParameters(self):
        imode = nspc = 0
        ctop = cbot = 0.0
        parF = 531.2593
        parQ = 2423.93994
        swF = 1351.81531
        swQ = 10213.1367

        spcf = [0.0] * 400
        spcq = [0.0] * 400
        npl = [0] * 200
        wl0d = []
        spcdf = []
        spcdq = []

        ch = ["", "hi", "lo", "lo"]

        if (self.amode == 2):
            # "Only monochro wavelength calculation!"
            imode = 1
            self.RF = 1000.0
            self.RQ = 1000.0
            self.wq = 1.0
            self.nwl = 1

        else:
            imode = int(input("imode: Integration mode 1:Monochro 2:PAR 3:SW\n"))
            if ((imode < 0) or (imode > 4)):
                print("Mode ERR: number should less than 4 and larger than 0.")
                return ERRCODE.INPUT_ERROR

            elif (imode == 1):
                self.wl0 = float(input("wl0:wavelength (micron ex 0.55)\n"))
                self.wls = 0.2
                self.nwl = int(round(4.0 - 0.3) / self.span[1])
                npl[1] = self.nPhotonProcess

            elif (imode == 2):
                self.wls = 0.4
                self.nwl = int(round(0.7 - 0.4) / self.span[1])
                self.RF = parF
                self.RQ = parQ

            elif (imode == 3):
                self.wls = 0.3
                self.nwl = int(round(0.7 - 0.3) / self.span[1])
                self.nwl += int((4.0 - 0.7) / (5.0 * self.span[1]))
                self.RF = swF
                self.RQ = swQ

            # read solar radiation
            contentFile = np.loadtxt("..\data\solar_rad")
            for i in range(len(contentFile)):
                wl0d.append(contentFile[i][0])
                spcdf.append(contentFile[i][1])
                spcdq.append(contentFile[i][2])
            nspc = len(contentFile)

            # search a spectral irradiance data in monochromatic calculation
            j = 1
            wlsd = self.wls + 0.0005

            if (imode == 1):
                for i in range(nspc):
                    if ((wlsd > self.wl0 - 0.0025) and (wlsd < self.wl0 + 0.0025)):
                        self.RF = (wl0d[i + 1] * spcdf[i] + wl0d[i] * spcdf[i + 1]) / (wl0d[i] + wl0d[i + 1])
                        self.RQ = (wl0d[i + 1] * spcdq[i] + wl0d[i] * spcdq[i + 1]) / (wl0d[i] + wl0d[i + 1])
                        spcf[j] = self.RF
                        spcq[j] = self.RQ
                        self.wq = 1.0
                        self.nwl = 1
                        break
                    wlsd += 0.001

            elif (imode == 2):
                for i in range(nspc):
                    if (abs(wl0d[i] - wlsd) < 1.0E-4):
                        for k in range(20):
                            spcf[j] += spcdf[i + k]
                            spcq[j] += spcdq[i + k]

                        spcf[j] /= self.span[1] * 1000
                        spcq[j] /= self.span[1] * 1000
                        j += 1

                        if (j> self.nwl):
                            break
                        wlsd += self.span[1]

            elif (imode == 3):
                for i in range(nspc):
                    if (abs(wl0d[i] - wlsd) < 1.0E-4):
                        if (wlsd < 0.7):
                            for k in range(20):
                                spcf[j] += spcdf[i + k]
                                spcq[j] += spcdq[i + k]

                            spcf[j] /= self.span[1] * 1000
                            spcq[j] /= self.span[1] * 1000

                        else:
                            for k in range(100):
                                spcf[j] += spcdf[i + k]
                                spcq[j] += spcdq[i + k]

                            spcf[j] /= self.span[1] * 5000
                            spcq[j] /= self.span[1] * 5000

                        j += 1
                        if (j > self.nwl):
                            break

                        if (wlsd < 0.7):
                            wlsd += self.span[1]
                        else:
                            wlsd += 5.0 * self.span[1]

            # make input photon weight for each wavelength
            dum = 0.0

            for i in range(self.nwl + 1):
                if (i < 20):
                    dum += spcf[i]
                else:
                    dum += spcf[i] * 5.0

            for i in range(self.nwl):
                if (i < 20):
                    npl[i] = int(float(self.nPhotonProcess) * spcf[i] / dum)
                else:
                    npl[i] = int(float(self.nPhotonProcess) * 5.0 * spcf[i] / dum)

            dum = 0.0

            for i in range(self.nwl):
                dum += npl[i]

            self.nPhotonProcess = int(dum)
            self.nPhoton = self.nPhotonProcess * self.Nprocess
            print("Actural number of photon is: " + str(self.nPhoton))

            # with/without atmospheric
            if (self.amode == 2):
                self.atmType = 0
                self.aerosolType = 0
                self.taur = 0.0
                ctype = 0

            else:
                # read z profile
                common.Z_GRD = np.loadtxt("..\data\zgrd")
                print("atmType: Atmospheric profile")
                print(" 1: Tropical")
                print(" 2: Mid latitude summer")
                print(" 3: Mid latitude winter")
                print(" 4: High latitude summer")
                print(" 5: High latitude winter")
                print(" 6: US standard atm.")
                self.atmType = int(input())

                if (self.atmType == 1):
                    self.rfname = "Data/gas_TR_" + ch[imode]
                elif (self.atmType == 2):
                    self.rfname = "Data/gas_TMS_" + ch[imode]
                elif (self.atmType == 3):
                    self.rfname = "Data/gas_MW_" + ch[imode]
                elif (self.atmType == 4):
                    self.rfname = "Data/gas_HS_" + ch[imode]
                elif (self.atmType == 5):
                    self.rfname = "Data/gas_HW_" + ch[imode]
                elif (self.atmType == 6):
                    self.rfname = "Data/gas_US_" + ch[imode]
                else:
                    print("Input error!")
                    return ERRCODE.INPUT_ERROR

                # read the aerosol data
                print("aerosolType: aerosol type")
                print(" 1:  Continental clean")
                print(" 2:  Continental average")
                print(" 3:  Continental polluted")
                print(" 4:  Urban")
                print(" 5:  Desert")
                print(" 6:  Maritime clean")
                print(" 7:  Maritime polluted")
                print(" 8:  Maritime Tropical")
                print(" 9:  Arctic")
                print("10:  Antactic")
                print("11:  Smoke")
                self.aerosolType = int(input())

                #AOT
                taur = int(input("tauref: AOT at 0.550 micron\n"))
                print(" - this version uses a extinction with RH=70% -")

                #current version uses a mixed extinction coef. by Iwabushi san
                self.nmix = 2

                if (self.aerosolType == 1):
                    self.fname.append("../data/opt_type1_rh0.70_" + ch[imode])
                    self.d = 8000.0     # scale height (m)
                elif (self.aerosolType == 2):
                    self.fname.append("../data/opt_type2_rh0.70_" + ch[imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 3):
                    self.fname.append("../data/opt_type3_rh0.70_" + ch[imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 4):
                    self.fname.append("../data/opt_type4_rh0.70_" + ch[imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 5):
                    self.fname.append("../data/opt_type5_rh0.70_" + ch[imode])
                    self.d = 2000.0  # scale height (m)
                elif ((self.aerosolType == 6) or (self.aerosolType == 8)):
                    self.fname.append("../data/opt_type6_rh0.70_" + ch[imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 7):
                    self.fname.append("../data/opt_type7_rh0.70_" + ch[imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 8):
                    self.fname.append("../data/opt_type8_rh0.70_" + ch[imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 9):
                    self.fname.append("../data/opt_type9_rh0.70_" + ch[imode])
                    self.d = 99000.0  # scale height (m)
                elif (self.aerosolType == 10):
                    self.fname.append("../data/opt_type10_rh0.70_" + ch[imode])
                    self.d = 99000.0  # scale height (m)
                elif (self.aerosolType == 11):
                    print("smoke aerosol under construction !!")
                    self.d = 8000.0  # scale height (m)
                    return ERRCODE.INPUT_ERROR
                else:
                    print("Input error !!")
                    return ERRCODE.INPUT_ERROR

                # read cloud data
                print("ctype: cloud type")
                print(" 0:  Cloud-free")
                print(" 1:  Stratus continental")
                print(" 2:  Stratus maritime")
                print(" 3:  Cumulus continental clean")
                print(" 4:  Culumus continental pulluted")
                print(" 5:  Culumus maritime")
                print(" 6:  Fog")
                print(" 7:  Cirrus 1 (-25degC)")
                print(" 8:  Cirrus 2 (-50 degC)")
                print(" 9:  Cirrus 3 (-50 degC + small particles)")

                self.cloudType = int(input())

                if (self.cloudType != 0):
                    self.ctaur = float(input("ctauref:COT at 0.55 micron\n" ))
                    ctop = float(input("cloud top height (m)\n"))
                    cbot = float(input("cloud bottom height (m)\n"))

                    self.cflg = 1

                    if (ctop < cbot):
                        print("cloud top should be greater than cloud bottom")
                        return ERRCODE.INPUT_ERROR

                    if (self.cloudType == 1):
                        self.fname.append("../data/opt_type101_" + ch[imode])
                    elif (self.cloudType == 2):
                        self.fname.append("../data/opt_type102_" + ch[imode])
                    elif (self.cloudType == 3):
                        self.fname.append("../data/opt_type103_" + ch[imode])
                    elif (self.cloudType == 4):
                        self.fname.append("../data/opt_type104_" + ch[imode])
                    elif (self.cloudType == 5):
                        self.fname.append("../data/opt_type105_" + ch[imode])
                    elif (self.cloudType == 6):
                        self.fname.append("../data/opt_type106_" + ch[imode])
                    elif (self.cloudType == 7):
                        self.fname.append("../data/opt_type107_" + ch[imode])
                    elif (self.cloudType == 8):
                        self.fname.append("../data/opt_type108_" + ch[imode])
                    elif (self.cloudType == 9):
                        self.fname.append("../data/opt_type109_" + ch[imode])
                    else:
                        print("Input error !!")
                        return ERRCODE.INPUT_ERROR

                    for i in range(common.N_Z):
                        if (common.Z_GRD[i] < cbot):
                            self.cbnz = i
                        if (common.Z_GRD[i] > ctop):
                            self.ctnz = i
                            break

                    print(common.Z_GRD[self.cbnz], common.Z_GRD[self.ctnz])
                    print("clouds are located between " + str(self.ctnz) + " and " + str(self.cbnz))
        return ERRCODE.SUCCESS

    def readParameters(self):
        self.nPhoton = int(input("np: Input number of photon\n"))

        # fish eye simulation mode
        if (self.nPhoton == -4):
            if (self.Nprocess > 1):
                print("fish eye mode cannot work with multi-processors")
                print("please execute under single processor")
                return ERRCODE.INPUT_ERROR

            print("fish eye mode - selected")
            self.readVegParameters()
            self.surfaceType = 2
            return ERRCODE.SUCCESS

        # LAI calculation mode
        if (self.nPhoton == -5):
            if (self.Nprocess > 1):
                print("LAI mode cannot work with multi-processors")
                print("please execute under single processor")
                return ERRCODE.INPUT_ERROR

            print("LAI calculation mode - selected")
            self.readVegParameters()
            self.surfaceType = 2
            return ERRCODE.SUCCESS

        self.nPhotonProcess = int(self.nPhoton / self.Nprocess)

        # Atmospheric mode
        self.amode = int(input("amode: 1:atmospheric mode, 2: without atmospheric\n"))

        if (self.amode == 2):
            self.dif = int(input("dif: Input frac. of diffuse radiation (0-1)\n"))
            if ((self.dif < 0.0) or (self.dif > 1.0)):
                print("fraction of diffuse should be 0.0-1.0... exit ")
                return ERRCODE.INPUT_ERROR

        # Get geometrical parameters
        self.readGeoParameters()

        # Get atmospheric parameters
        self.readAtmParameters()

        # vegetation parameters read / initialize
        self.surfaceType = int(input("stye: Surface mode 1: Lambertian, 2: 3-D Vegetation\n"))
        if (self.surfaceType == 1):
            for iwl in range(1, self.nwl):
                if (self.imode == 1):
                    print("Input albedo")
                else:
                    if (iwl <= 20):
                        ab1 = self.wls + self.span[iwl] * float(iwl -1)
                        ab2 = self.wls + self.span[iwl] * float(iwl)
                    else:
                        ab1 = 0.7 + self.span[iwl] * float(iwl - 21)
                        ab2 = 0.7 + self.span[iwl] * float(iwl - 20)
                    print("Input albedo: " + str(ab1) + " and " + str(ab2))
                self.alb[iwl] = float(input())

        elif (self.surfaceType == 2):
            self.readVegParameters()

        else:
            print("Bad surface type selection exit")
            return ERRCODE.INPUT_ERROR

        return ERRCODE.SUCCESS
