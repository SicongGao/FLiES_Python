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
    rtype = 0  # default atmospheric and aerosol type
    atype = 0  # default aerosol type
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
    cbnz = 0  # Number of cloud bottom and top layer
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

    fname = ""
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
        self.rtype = 1  # default atmospheric and aerosol type
        self.atype = 2
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
        parF = 531.2593
        parQ = 2423.93994
        swF = 1351.81531
        swQ = 10213.1367

        spcf = [0.0] * 400
        spcq = [0.0] * 400
        npl = [0] * 200

        ch = ["hi", "lo", "lo"]

        if (self.amode == 2):
            # "Only monochro wavelength calculation!"
            imode = 1
            self.RF = 1000.0
            self.RQ = 1000.0
            self.wq = 1.0
            self.nwl = 1

        else:
            imode = int(input("imode: Integration mode 1:Monochro 2:PAR 3:SW"))
            if ((imode < 0) or (imode > 4)):
                print("Mode ERR: number should less than 4 and larger than 0.")
                return

            elif (imode == 1):
                self.wl0 = float(input("wl0:wavelength (micron ex 0.55)"))
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
