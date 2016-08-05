# Parameters Initialization

import common
import math
import numpy as np


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
    nPhoton = nmix = amode = imode = cmode = stype = 0
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
        self.flg = 0

    def readVegParameters(self):






        return

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
        return 0

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
                return

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
                    return -1

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
                    return -1
                else:
                    print("Input error !!")
                    return -1

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
                        return -1

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
                        return -1

                    for i in range(common.N_Z):
                        if (common.Z_GRD[i] < cbot):
                            self.cbnz = i
                        if (common.Z_GRD[i] > ctop):
                            self.ctnz = i
                            break

                    print(common.Z_GRD[self.cbnz], common.Z_GRD[self.ctnz])
                    print("clouds are located between " + str(self.ctnz) + " and " + str(self.cbnz))
        return 0

    def readParameters(self):
        self.nPhoton = int(input("np: Input number of photon\n"))

        # fish eye simulation mode
        if (self.nPhoton == -4):
            if (self.Nprocess > 1):
                print("fish eye mode cannot work with multi-processors")
                print("please execute under single processor")
                return

            print("fish eye mode - selected")
            self.readVegParameters()
            self.stype = 2
            return

        # LAI calculation mode
        if (self.nPhoton == -5):
            if (self.Nprocess > 1):
                print("LAI mode cannot work with multi-processors")
                print("please execute under single processor")
                return

            print("LAI calculation mode - selected")
            self.readVegParameters()
            self.stype = 2
            return

        self.nPhotonProcess = int(self.nPhoton / self.Nprocess)

        # Atmospheric mode
        self.amode = int(input("amode: 1:atmospheric mode, 2: without atmospheric\n"))

        if (self.amode == 2):
            self.dif = int(input("dif: Input frac. of diffuse radiation (0-1)\n"))
            if ((self.dif < 0.0) or (self.dif > 1.0)):
                print("fraction of diffuse should be 0.0-1.0... exit ")
                return

        # Get geometrical parameters
        self.readGeoParameters()

        # Get atmospheric parameters
        self.readAtmParameters()

        # vegetation parameters read / initialize
        self.stype = int(input("stye: Surface mode 1: Lambertian, 2: 3-D Vegetation\n"))
        if (self.stype == 1):
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

        elif (self.stype == 2):
            self.readVegParameters()

        else:
            print("Bad surface type selection exit")

        return 0
