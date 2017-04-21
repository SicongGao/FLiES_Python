import common as comm
from math import *
import numpy as np
import ERRCODE
import config.config as config
import logging
from Position import Position

# Parameters Initialization
class Parameters:
    ix = [0] * 1025
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
    nPhoton = nmix = AtmMode = imode = cmode = surfaceType = 0
    diffuse = phi = th = ph = tgx = tgy = 0.0
    wq = sin_f0 = cos_f0 = cos_q0 = sin_q0 = 0

    npl = [0] * 200
    alb = [0.0] * 100
    spcf = [0.0] * 400
    spcq = [0.0] * 400
    ulr = [0.0] * 100 # reflectance of forest floor vegetaion
    ult = [0.0] * 100 # tansmittance of forest floor vegetaion
    soilRef = [0.0] * 100 # soil reflectance
    tlr = [0.0] * 5
    tlt = [0.0] * 5
    tstr = [0.0] * 5

    lr = np.zeros(5 * 100, dtype=float).reshape(5, 100) # reflectance of crown foliage
    lt = np.zeros(5 * 100, dtype=float).reshape(5, 100) # transmittance of crown foliage
    truncRef = np.zeros(5 * 100, dtype=float).reshape(5, 100) # trunk reflectance
    ext = np.zeros(10 * 200, dtype=float).reshape(10, 200)
    wkd = [0.0] * 200

    fname = [""]
    rfname = ""

    nts = 0  # group of tree species
    bound = 1

    # for math
    T_SIN = [0.0] * comm.ANGLE_SHIFT * 2
    T_COS = [1.0] * comm.ANGLE_SHIFT * 2
    T_ACOS = [0.0] * comm.ACOOS_SHIFT * 2
    T_EXP = [0.0] * comm.ACOOS_SHIFT
    DLT = np.zeros((6 + 1) * (6 + 1), dtype=float).reshape((6 + 1), (6 + 1))

    # for output
    #AP_NP = [0.0] * 100

    def __init__(self):
        self.initParameters()

    def initParameters(self):
        comm.N_Z = 12  # of layers
        comm.X_MAX = comm.Y_MAX = 30.0  # X domain size (m), Y domain size (m)
        comm.RES = (comm.SIZE - 1) / comm.X_MAX  # inverse of the spatial resolution
        comm.WRR = 0.1  # ideal weight (used for Russian roulette in atm.

        # branch portion in region 1(outer canopy) and 2(internal canopy)
        # 1.0 = 100 % branch, 0.0: 100 % leaf
        comm.BP1 = 0.0
        comm.BP2 = 1.0

        # leaf angle distribution of canopy, branches, and floor
        comm.M_C = comm.M_B = comm.M_F = 1

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
        self.hfov = radians(self.hfov)
        self.hfov = cos(self.hfov)
        self.cflg = 0  # cflg = 0 -> cloud - free, cflg = 1 -> cloud

        self.span[1: 20] = [0.02] * 20  # Wavelength span (0.3-0.7 micron)
        self.span[21: 100] = [0.1] * 80  # Wavelength span (0.7- micron)

        # Irradiance
        self.tflx = 0.0  # Total downward flux at top of canopy
        self.bflx = 0.0  # Beam downward flux at TOC
        self.dflx = 0.0  # Diffuse downward flux at TOC
        self.obflx = 0.0  # Beam downawr flux at TOC(observed)
        self.odflx = 0.0  # Diffuse downward flux at TOC(observed)

        # reflected radiance from canopy
        self.rflx = 0.0  # total reflcted irradiaice
        self.rbflx = 0.0  # beam reflected irradiance
        self.rdflx = 0.0  # diffuse reflected irradiance

        # Photon flux density
        self.tpfd = 0.0  # Total downward PFD at top of canopy
        self.bpfd = 0.0  # Beam downward PFD at top of canopy
        self.dpfd = 0.0  # Diffuse downward PDF at top of canopy
        self.obpfd = 0.0  # Beam downawr PFD at TOC(observed)
        self.odpfd = 0.0  # Diffuse downward PFD at TOC(observed)

        # Number of cloud bottom and top layer
        self.cbnz = 0
        self.ctnz = 1

        # incident irradiance and photon flux density at TOA or TOC
        self.RF = 0.0
        self.RQ = 0.0

        self.ix = [0, 715327539,474399417,182768403,
       780769431,576730197,937392842,445151859,688527083,160158012,
       750394380,704034435,445074075,668471240,896236437,567474752,
       686366724,313996988,336207690,785978227,778020274,918287307,
       165417294,345438969,907890635,347416159,564662778,423338472,
       322334997,537865632,861550217,847868806,131114898,991933631,
       711366790,789845275,154507895,445727944,529958587,952987849,
       941223627,914187765,387129643,731071531,106928367,428804814,
       339530062,987977910,349373608,666589081,249856479,762473708,
       492770269,984295237,762602078,309875385,181659605,554070627,
       544579014,145075585,518001246,574235904,470872688,375689646,
       919288736,793184089,576772630,583573585,630078959,718095433,
       922435718,260494938,899714994,845035588,550636380,614489316,
       129748380,342385542,752870810,561946344,901663726,773463362,
       291476356,142718061,857835626,599125462,822165143,148115470,
       102699640,719609755,114280930,315919730,632002210,662364387,
       475737383,521125644,473855143,713211613,581004834,124188662,
       385785588,212828838,395410802,861018067,166374165,861383384,
       600531530,232779699,631097769,686852794,999512106,415913715,
       881422436,984974586,539663088,857359588,472863933,448952791,
       946881961,785362821,238347978,227637797,822765690,738837653,
       262806318,959824693,686128383,326913756,243790771,282925300,
       347129324,483847659,666342264,911506146,217384530,323057584,
       898291599,182956265,233336405,245978702,227818927,221598428,
       509776547,447233226,214040942,414571726,517627480,507070145,
       293723551,952738082,125537664,628467971,608408802,524136102,
       279781483,918029707,548632100,640354341,235301552,938506978,
       780259168,659470605,956230694,325113028,469281736,218069674,
       669464838,897083532,428805136,686481469,360384586,358491134,
       472625645,958973413,731986004,450382757,828050124,528788635,
       516020488,711633133,549766916,820648300,744078099,756772887,
       181891301,327751624,988430935,213427722,866142696,964985638,
       496034786,735881912,572488659,743238997,879294425,713844186,
       767544317,516990455,924970513,495397278,293978777,835805356,
       387324693,959867930,279599267,563193517,901016831,583683180,
       467890018,888909125,290367062,825984930,897427982,667439448,
       193090536,669374716,796016824,842227637,893353706,100791099,
       874773830,225275495,700911164,505312681,384159177,495752993,
       185093925,754601478,990426230,540048843,117183163,284735697,
       654715317,720082843,589553463,999204993,292192117,956943678,
       520778593,787604343,475951343,238244029,738028538,450868237,
       882337337,852777892,714244800,845837140,319648818,163957206,
       158674488,640909075,952790546,799179255,480213829,353940284,
       619044995,555360287,428642246,562991225,297469578,458329311,
       272641456,355631524,310239267,815292841,100535538,449940916,
       726618528,110184604,638394987,845519244,222809445,274570728,
       777646374,989143276,722952055,386900100,173563846,259192337,
       804837858,150796596,944107592,241958057,536792668,172954717,
       245518931,765096962,292466963,598520195,342293274,742080068,
       732690459,332958489,487763014,331451785,453715831,304306945,
       670086896,387259837,228740735,573510956,726406580,433203637,
       671245718,170436887,221798990,328530055,528306588,384406611,
       432767805,562542760,704605746,837753498,786313664,629190021,
       191473640,666614669,755747532,379664975,702568554,695219945,
       570526945,565430802,113455511,729851019,708613234,846314841,
       461649617,449048089,973078024,396699738,389404076,724621301,
       591522955,566275215,346788763,514627617,805478692,874203109,
       984912198,859577828,624179011,987852382,776539105,645075726,
       205693244,235354499,950277531,601152193,119875329,678681015,
       703088581,926062923,501421090,737043684,835391223,339952000,
       736254954,626308470,726111215,589620465,246772944,575693416,
       244674022,664175254,804692804,492974600,502179351,105707521,
       912971061,712505763,772011214,409435805,514390751,434718710,
       568864083,693147993,108412769,717491453,163571384,748691016,
       193417558,103828531,244246532,738387846,294243296,629307234,
       450758266,322059829,351363754,160164791,367319893,587539178,
       447192484,850313907,566630554,149750910,831939166,588153833,
       151727745,870573383,992379736,351560574,619998413,357348138,
       398057472,510176438,798292303,778476464,269448219,235280376,
       546121954,946146821,171925198,713047999,751387333,145542333,
       540686538,213339531,646850329,254061029,567989629,101491674,
       511445713,976020085,765062844,854291462,126419931,320127592,
       726073664,501840615,784555584,416186872,298616638,204913941,
       530023872,333112770,711398977,421232482,105471433,698953098,
       718673020,211273633,639423400,481902280,834338939,337258231,
       174255883,958625638,597575950,439795354,988194686,725370764,
       288606071,652712941,172317773,325139743,940627425,492783117,
       544820171,488717451,606341838,582640069,792202669,592277354,
       302321077,258376838,831743204,792493635,310322348,249334776,
       769693785,506666097,202308101,126671762,852585899,243781933,
       212450419,658941084,781093174,142874742,976827108,774326443,
       954241406,423692041,437779057,549959793,687359732,584795171,
       687710458,396150502,210256204,749684453,794749748,148324950,
       497469338,535762405,520967233,329461720,307514692,695794153,
       976011556,567257386,908122807,379152056,132028931,968366134,
       158820642,132730855,394504806,284807178,713499200,565936505,
       857418435,116150068,922622507,299500587,244074843,634191107,
       434319168,498511296,211463896,254241408,509565162,174777144,
       619098263,461522963,979462164,194666354,479420781,654949367,
       546923801,163301688,194472229,701938289,729811966,704713141,
       172904969,337296640,573704075,450154206,821816778,352260041,
       542208772,740407979,513786503,360789814,481501424,877278476,
       618555653,294747538,852938073,941717422,425875252,316396425,
       585554182,501336869,220635233,326374846,482671993,951440107,
       742191863,967128026,140806645,884695804,254039706,525141957,
       268142265,567372882,102905902,645042306,319375286,843252885,
       412742513,558692395,284106625,157423258,815817373,698219943,
       309580288,735815125,137278393,790358489,945752859,237983600,
       396490848,782421189,500067272,753798800,435267785,335600680,
       821875894,309346842,522548636,936234825,363221210,816270077,
       582875460,850904583,545167812,463222974,339507237,308870281,
       487556564,668586951,768946576,840601950,380444666,970218682,
       130818789,979514467,574695205,740414148,875068712,155558241,
       950725245,243294508,298463980,521267265,632429379,645401185,
       729721897,136967290,199017487,444539350,997838729,839200335,
       587381196,880762827,409252047,128474914,763457381,833866065,
       703891903,492486625,938723272,988592243,272785383,440529581,
       120519592,180093396,280580070,786109495,777856391,936727976,
       123705718,169591884,630828422,871188521,159661253,927312511,
       555372625,239872492,571394586,245277479,566895931,616924226,
       641473895,565217298,578284591,411284571,396620988,226639439,
       978504025,722013819,707072252,627414184,728702336,629999244,
       236095204,545612818,233551371,128490893,929005736,541735309,
       146455038,388578063,701168870,497181886,605335366,593862646,
       649999511,700344359,596646457,941726326,131553557,196462468,
       728614306,552047276,187689284,466783955,250581520,345926004,
       675570994,721695333,538164564,320628093,120164535,848130375,
       442103126,742325812,733032494,804702085,865477508,751452833,
       299490810,944428598,186884447,739405047,928227144,685114991,
       352093368,644984370,794827264,487401399,580778777,538013046,
       577440607,531804564,556137967,210654774,909569001,431705918,
       649663913,347817739,942342323,686363077,104334047,311862446,
       307242770,404658630,444671663,637728619,189879047,492842099,
       351035720,497138702,637451654,536757424,474827873,203027128,
       282233773,932952767,410118830,789557313,872889310,839059573,
       472144806,824421149,713342022,561833637,732224237,326279896,
       976913100,614424353,758227449,574671709,377494746,615525883,
       927542752,151950633,727296161,154128218,822941052,471023160,
       481720185,882544672,442547461,615864646,615197741,573103207,
       793041825,521291270,878148531,450254708,830705833,349478349,
       832974821,895334196,196964927,245983436,348982274,768813431,
       911096143,891245651,180867790,567402976,478034266,628175771,
       945884341,229315733,569857144,186785695,165582142,475116720,
       707223582,686165130,367357766,253716620,739102226,569993185,
       792424112,621876281,606350421,228971847,890372592,244926230,
       444357630,404276227,150045567,875459617,892705684,314217399,
       397682526,665889990,201763533,671684849,485985970,890571665,
       426160129,234326247,737729150,718820112,939559102,128263398,
       713078951,537823629,430655351,549691599,975228619,907891225,
       505661314,774230527,423441469,467359691,535485547,314569331,
       873291748,917821514,231576071,766682898,216533385,431382417,
       276676692,243621818,451869961,978292721,966844195,831742131,
       747027510,444085198,995860546,728259825,370109578,439736077,
       814043951,270593723,233908681,640119379,542656674,618111264,
       529876700,421904590,486525496,405207249,814850115,681845700,
       548767498,400876072,362372827,531800755,383228987,229728204,
       819131296,700170713,172289831,906871557,850126636,864350444,
       627107125,450550422,675152140,868928599,211951635,166159857,
       337720108,891153866,361899819,642024445,439627662,143827704,
       199016615,240285606,443759444,359789028,815392673,564878320,
       209565602,388550195,928746581,489589545,230407714,563131290,
       518506333,278271935,680408090,512241283,480297219,282204416,
       855593729,850283330,630477321,671157473,363739091,857839488,
       432136413,724229806,108647099,781859052,169139583,780270111,
       773479938,837574487,465207675,747239726,194121423,436023390,
       145433472,474004194,421773108,983829176,436824458,275210985,
       940699148,166465293,877754408,439134913,850040054,698777306,
       410009986,570765233,848371505,693898797,757403367,623156607,
       516752839,263076202,703176987,697732532,349000808,699375867,
       956577932,462465894,851222318,581388550,957187867,852737927,
       399635469,145038761,117024413,636334460,604822957,573563528,
       414576151,678724038,665047293,898052990,110377544,318905645,
       496241262,498133292,660918354,947286063,582360690,637050825,
       704483383,474406659,929325616,747815275,391085982,936901086,
       708458739,204045931,692734181,129042436,696603751,869574207,
       853373664,836296576,153963825,433753517,322276310,642164617,
       700613063,695282870,907520276,246662250,926170587,616911351,
       882269155,618787074,538490641,105029377,852397340,770021498,
       665877866,276560163,350916147,205514039,416253310,498361334,
       854810738,282842473,952407902,452399832,966281896,295668850,
       871440219,612405884,906151652,634743267,182620781,631412714,
       345198562]

        # initialization of random number generator
        # self.flg = frndi(ix(Nid + 1)) move to the Random.py _init_()
        self.flg = 0

        self.initMathparameters()

        return ERRCODE.SUCCESS

    def initMathparameters(self):
        for i in range(0, 62832 * 2):
            self.T_SIN[i] = sin(float(i - 62832) * 0.0001)
            self.T_COS[i] = cos(float(i - 62832) * 0.0001)

        for i in range(0, comm.ACOOS_SHIFT * 2):
            self.T_ACOS[i] = acos(float(i - comm.ACOOS_SHIFT) * 0.0001)

        for i in range(1, 6):
            comm.DLT[i, i] = 1.0

        return ERRCODE.SUCCESS

    def getInputLR_LT_ULR_ULT(self, i, **args):
        if (config.INPUT_MODE):
            tempList = list(args.get("optical_parameters"))
            #print(args)
            for j in range(1, self.nts + 1):
                self.lr[j,i] = tempList.pop(0)

            for j in range(1, self.nts + 1):
                self.lt[j, i] = tempList.pop(0)

            self.ulr[i] = tempList.pop(0)
            self.ult[i] = tempList.pop(0)

            for j in range(1, self.nts + 1):
                self.truncRef[j] = tempList.pop(0)

            self.soilRef[i] = tempList.pop(0)
        else:
            for j in range(1, self.nts + 1):
                self.lr[j, i] = float(input())

            for j in range(1, self.nts + 1):
                self.lt[j, i] = float(input())

            self.ulr[i] = float(input())
            self.ult[i] = float(input())

            for j in range(1, self.nts + 1):
                self.truncRef[j] = float(input())
            self.soilRef[i] = float(input())

        return ERRCODE.SUCCESS

    def process100(self, **args):

        if (self.cmode >= 2):
            comm.N_ANG_C = 1
            comm.URC_coord[1].setPosition(0.0, 0.0, 1.0)

        self.process201(**args)

        return ERRCODE.SUCCESS

    def process202(self, **args):
        umax = 0.0

        logging.info("u: leaf area density 1,2,3...# tree species")
        if (config.INPUT_MODE):
            comm.U = list(args.get("leaf_area_density"))
            logging.info(comm.U)
            comm.U.insert(0, -1)
        else:
            for i in range(1, self.nts + 1):
                comm.U[i] = float(input())

        if (not((self.nPhoton == -4) or (self.nPhoton == -5))):
            logging.info("gLAI: forest floor LAI")
            if (config.INPUT_MODE):
                comm.G_LAI = args.get("forest_floor_LAI")
                logging.info(comm.G_LAI)
            else:
                comm.G_LAI = float(input())

        logging.info("BAD: branch area density 1,2,3... # of tree species")
        if (config.INPUT_MODE):
            comm.BAD = list(args.get("branch_area_density"))
            logging.info(comm.BAD)
            comm.BAD.insert(0, -1)
        else:
            for i in range(1, self.nts + 1):
                comm.BAD[i + 1] = float(input())

        for i in range(1, self.nts + 1):
            if ((comm.U[i] < 0.0) or (comm.U[i] > 8.0)):
                logging.ERROR(str(i) + "th leaf area density " + str(comm.U[i]) + " should be set in the range (0.0-8.0)")
                logging.ERROR("EXIT")
                return ERRCODE.INPUT_ERROR

        if ((comm.G_LAI < 0.0) or (comm.G_LAI > 8.0)):
            logging.ERROR(str(comm.G_LAI) + " should be set in the range (0.0-8.0)")
            logging.ERROR("EXIT")
            return ERRCODE.INPUT_ERROR

        if (self.nPhoton != -5):
            logging.info("sbar: Spherical ave. shoot silhouette to total needle area ratio")
            logging.info("1,2,3... # of tree species (0.0-0.25)")
            logging.info("For broadleaves, please input 0.25")
            if (config.INPUT_MODE):
                comm.S_BAR = list(args.get("sbar"))
                logging.info(comm.S_BAR)
                comm.S_BAR.insert(0, -1)
            else:
                for i in range(1, self.nts + 1):
                    comm.S_BAR[i] = float(input())

        umax = comm.U[1]
        for i in range(1, self.nts):
            if (umax < comm.U[i]):
                umax = comm.U[i]

        if (umax > 1.0):
            comm.FE = 1.0
        elif (umax <= 0.01):
            comm.FE = 0.01
        else:
            comm.FE = umax

        # canopy object parameters
        # object id initialization
        comm.OBJ_Group = [1] * 6000

        # input obj_nt
        result = ''
        with open("../data/crowndata.txt", "r") as file:
            comm.N_OBJ = int(file.readline())

            if (comm.N_OBJ == 0):
                comm.N_OBJ = comm.OBJ_Shape = 1
                comm.OBJ[0, 0] = 0.01
                comm.OBJ[0, 1] = 0.01
                comm.OBJ[0, 2] = 0.01
                comm.OBJ[0, 3] = 1E-5
                comm.OBJ[0, 4] = 1E-5
            else:
                result = file.readlines()
                result = np.loadtxt(result)
        file.close()

        obj_nt = comm.N_OBJ
        for i in range(len(result)):
            comm.OBJ_Shape[i] = int(result[i][0])
            comm.OBJ[i][0:5] = result[i][1: 6]
            comm.OBJ_Group[i] = int(result[i][6])
            if (result[i][0] != 4):
                if ((result[i][4] < 0.01) or (result[i][5] < 0.01)):
                    logging.info(str(i + 1) + "th canopy neglected!")
                    obj_nt -= 1

        comm.N_OBJ = obj_nt

        # check the obj id range
        for i in range(obj_nt):
            if ((comm.OBJ_Group[i] <= 0) or (comm.OBJ_Group[i] > self.nts)):
                logging.info("species id should be the range betweem 0-" + str(self.nts))
                return ERRCODE.OUT_OF_RANGE

        logging.info("Total Object is " + str(obj_nt))

        # change from height to radius
        for i in range(obj_nt):
            if (comm.OBJ_Shape[i] == 3):
                comm.OBJ[i, 3] /= 2

        # in case periodic boundary
        # add objects that are partially in the outside from the simulation scene
        a = [1.0, 1.0, 0.0, 0.0]
        b = [0.0, 0.0, 1.0, 1.0]
        # preparation of the i-th grid for space divided method
        c = [0.0, comm.X_MAX, 0.0, comm.Y_MAX]
        cc = [0.0, 1.0, 0.0, 1.0]
        iMin = 0.0
        iMax = comm.X_MAX

        if (self.bound == 1):
            # set distance calculation parameters
            kobj = comm.N_OBJ

            # definition of rectangular of the i-th object
            for j in range(kobj):
                xr = comm.OBJ[j, 0]
                yr = comm.OBJ[j, 1]

                # check the intersection on the x - y plane
                for k in range(4):
                    d = abs(xr * a[k] + yr * b[k] - c[k])

                    if (d <= comm.OBJ[j, 4]):   # partially in the out of range
                        dd = sqrt(comm.OBJ[j, 4] ** 2 - d ** 2)
                        p1 = b[k] * xr + a[k] * yr - dd
                        p2 = b[k] * xr + a[k] * yr + dd

                        if (((iMin < p1) and (iMax > p1)) or
                                ((iMin < p2) and (iMax > p2))):

                            n = comm.N_OBJ
                            comm.OBJ[n, 0] = comm.OBJ[j, 0] + a[k] * comm.X_MAX * (1.0 - 2.0 * cc[k])
                            comm.OBJ[n, 1] = comm.OBJ[j, 1] + b[k] * comm.Y_MAX * (1.0 - 2.0 * cc[k])
                            comm.OBJ[n, 2] = comm.OBJ[j, 2]
                            comm.OBJ[n, 3] = comm.OBJ[j, 3]
                            comm.OBJ[n, 4] = comm.OBJ[j, 4]

                            comm.OBJ_Shape[n] = comm.OBJ_Shape[j]
                            comm.OBJ_Group[n] = comm.OBJ_Group[j]

                            comm.N_OBJ += 1

                for k in range(2):
                    for l in range(2):
                        rx = xr - c[k]
                        ry = yr - c[l + 2]
                        rr = sqrt(rx * rx + ry * ry)

                        if (rr <= comm.OBJ[j, 4]):

                            n = comm.N_OBJ
                            comm.OBJ[n, 0] = comm.OBJ[j, 0] + comm.X_MAX - 2.0 * c[k]
                            comm.OBJ[n, 1] = comm.OBJ[j, 1] + comm.Y_MAX - 2.0 * c[l + 2]
                            comm.OBJ[n, 2] = comm.OBJ[j, 2]
                            comm.OBJ[n, 3] = comm.OBJ[j, 3]
                            comm.OBJ[n, 4] = comm.OBJ[j, 4]

                            comm.OBJ_Shape[n] = comm.OBJ_Shape[j]
                            comm.OBJ_Group[n] = comm.OBJ_Group[j]

                            comm.N_OBJ += 1

            # determination of the epsi(epsiron) that is used to perform the
            # Russian Roulette
            # epsi is setted to be c * (leaf_r + leaf_t)
            if (self.nPhoton != -4):
                avelr = 0.0
                avelt = 0.0

                for i in range(1, self.nwl + 1):
                    for j in range(1, self.nts + 1):
                        avelr += self.lr[j, i]
                        avelt += self.lt[j, i]

                avelr /= float(self.nwl * self.nts)
                avelt /= float(self.nwl * self.nts)
                comm.EPSI = 0.5 * 0.1 * (avelr + avelt)

        return ERRCODE.SUCCESS

    def process201(self, **args):
        ispc = 0
        # get the number of forest species
        logging.info("nts: # of group of tree species")
        if (config.INPUT_MODE):
            self.nts = args.get("tree_species")
            logging.info(self.nts)
        else:
            self.nts = int(input())

        comm.N_TS = self.nts

        if ((self.nPhoton == -4) or (self.nPhoton == -5)):
            return self.process202()

        #comm.U[5] = float(self.nts)

        # currently max nts should be less than 5
        if ((self.nts <= 0) or (self.nts >= 5)):
            logging.ERROR("Error # of tree species")
            logging.ERROR("tree species should be 1-5")
            return ERRCODE.INPUT_ERROR

        # read leaf reflectance/transmittance
        ramda = self.wls
        if(self.nPhoton == -4):
            self.nwl = 1

        for i in range(1, self.nwl + 1):
            if (self.nwl == 1):
                logging.info("Input surface optical properties")
            else:
                logging.info("Input surface optical properties in")
                if (i <= 20):
                    logging.info(str(ramda) + " and " + str(ramda + self.span[1]))
                else:
                    logging.info(str(ramda) + " and " + str(ramda + 5 * self.span[1]))

            if (i == 1):
                logging.info("lr1 lr2.. lt1 lt2.. ulr ult str1 str2.. sor")
                logging.info("(lr, lt:leaf refl. & leaf transm.")
                logging.info("ulr,ult:understory leaf refl. & transm.")
                logging.info("stmr: stem refl., soilr: soil refl.)")

            if (self.nwl == 1):
                self.getInputLR_LT_ULR_ULT(i, **args)

            else:
                if (i == 1):
                    logging.info("Input PAR average values:")
                    self.getInputLR_LT_ULR_ULT(i, **args)
                    ispc = i

                elif (i == 21):
                    logging.info("Input NIR average values")
                    self.getInputLR_LT_ULR_ULT(i, **args)
                    ispc = i

                else:
                    for j in range(1, self.nts + 1):
                        self.lr[j, i] = self.lr[j, ispc]
                        self.lt[j, i] = self.lt[j, ispc]
                        self.truncRef[j, i] = self.truncRef[j, ispc]

                    self.ulr[i] = self.ult[ispc]
                    self.ult[i] = self.ult[ispc]
                    self.soilRef[i] = self.soilRef[ispc]

            # check the parameter range
            for j in range(1, self.nts + 1):
                if (self.lr[j,i] + self.lt[j,i] > 0.99):
                    logging.ERROR("canopy leaf reflectance+transmittance is too large, exit!")
                    return ERRCODE.OUT_OF_RANGE

                if (self.lr[j, i] + self.lt[j, i] < 0.0001):
                    self.lr[j,i] = 0.0001
                    self.lt[j,i] = 0.0001

            if (self.ulr[i] + self.ult[i] > 0.99):
                logging.ERROR("floor leaf reflectance+transmittance is too large, exit!")
                return ERRCODE.OUT_OF_RANGE

            if (self.ulr[i] + self.ult[i] < 0.0001):
                self.ulr[i] = 0.0001
                self.ult[i] = 0.0001

            for j in range(1, self.nts + 1):
                if (self.truncRef[j, i] > 1.00):
                    logging.ERROR("stem reflectance is too large, exit!")
                    return ERRCODE.OUT_OF_RANGE

                if (self.truncRef[j, i] < 0.0001):
                    self.truncRef[j, i] = 0.0001

            if (self.soilRef[i] > 1.0):
                logging.ERROR("soil reflectance is too large, exit!")
                return ERRCODE.OUT_OF_RANGE

            if (self.nwl <= 20):
                ramda += self.span[1]
            else:
                ramda += 5 * self.span[i]

        self.process202(**args)

        return ERRCODE.SUCCESS

    def readVegParameters(self, **args):

        if ((self.nPhoton == -4) or (self.nPhoton == -5)):
            return self.process201()

        # input output mode
        logging.info("cmode: calculation mode")
        logging.info("1: BRF only 2: BRF Nadir Image 3: 3D APAR")
        if (config.INPUT_MODE):
            self.cmode = args.get("calculation_mode")
            logging.info(self.cmode)
        else:
            self.cmode = int(input())

        if ((self.cmode <= 0) or (self.cmode >= 4)):
            logging.ERROR("Bad mode selection exit")
            return ERRCODE.INPUT_ERROR

        # read the condition file
        if ((self.bound != 1) and (self.bound != 2)):
            logging.info("boundary condition should be 1 or 2.")
            logging.info("  1: Periodic  2: Non-periodic")
            return ERRCODE.INPUT_ERROR

        if (self.cmode != 1):
            self.process100()

        logging.info("nth, angt: # of angle anfinished process 201 200d zenith angle(max 18)for BRF")
        logging.info("eg. 5 10. 20. 45. 50. 70.")
        if (config.INPUT_MODE):
            comm.ANG_T = list(args.get("BRF_zenith_angles"))
            comm.N_TH = len(comm.ANG_T)
            logging.info(str(comm.N_TH) + ": " + str(comm.ANG_T))
            comm.ANG_T.insert(0, -1)
        else:
            self.N_TH = int(input())
            for i in range(1, comm.N_TH + 1):
                comm.ANG_T[i] = float(input())

        logging.info("nph, angp:# of angle and azimuth angle(max 36)for BRF")
        logging.info("eg. 3 0. 90. 180.")
        if (config.INPUT_MODE):
            comm.ANG_P = list(args.get("BRF_azimuth_angles"))
            comm.N_PH = len(comm.ANG_P)
            logging.info(str(comm.N_PH) + ": " + str(comm.ANG_P))
            comm.ANG_P.insert(0, -1)
        else:
            self.N_PH = int(input())
            for i in range(1, comm.N_PH + 1):
                comm.ANG_P[i] = float(input())

        if (comm.N_TH * comm.N_PH > 700):
            logging.ERROR("The number of sampling angles are too huge! It should be theta*phi<348")
            return ERRCODE.INPUT_ERROR

        if (comm.N_TH > 18):
            logging.ERROR("The number of sampling theta are too huge! It should be theta*phi<100")
            return ERRCODE.INPUT_ERROR

        if (comm.N_PH > 36):
            logging.ERROR("The number of sampling phi are too huge! It should be theta*phi<100")
            return ERRCODE.INPUT_ERROR

        k = 0

        temp = Position()
        temp.setPosition(0, 0, 0)
        comm.URC_coord.append(temp)
        for i in range(1, comm.N_PH + 1):
            for j in range(1, comm.N_TH + 1):
                if (comm.ANG_T[j] > 80.0):
                    logging.warning("Zenith angle should be less than 80")
                    logging.warning(str(comm.ANG_T[j]) + " is ignored !")
                else:
                    k += 1
                    temp = Position()
                    temp.setPosition(sin(radians(comm.ANG_T[j])) * cos(radians(comm.ANG_P[i])),
                                     sin(radians(comm.ANG_T[j])) * sin(radians(comm.ANG_P[i])),
                                     cos(radians(comm.ANG_T[j])))
                    comm.URC_coord.append(temp)

        comm.N_ANG_C = k

        self.process100(**args)

        return ERRCODE.SUCCESS

    def readGeoParameters(self, **args):

        f0 = q0 = 0.0
        ntha = npha = 0
        th = [0] * 18
        ph = [0] * 36

        # solar zenith angle

        #self.th0 = float(input("the0: Solar zenith angle(degree)\n"))
        #comm.Z_MIN = float(input(" \n"))

        logging.info("the0: Solar zenith angle(degree)")
        if (config.INPUT_MODE):
            self.th0 = args.get("solar_angle")
            logging.info(self.th0)
        else:
            self.th0 = float(input())

        logging.info("Solar elevation")
        if (config.INPUT_MODE):
            comm.Z_MIN = args.get("solar_elevation")
            logging.info(comm.Z_MIN)
        else:
            comm.Z_MIN = float(input())

        q0 = pi - radians(self.th0)

        if (self.ph0 <= pi):
            f0 = radians(self.ph0) + pi
        else:
            f0 = radians(self.ph0) - pi

        self.sin_q0 = sin(q0)
        self.cos_q0 = cos(q0)
        self.sin_f0 = sin(f0)
        self.cos_f0 = cos(f0)
        #self.sin_f0 =

        if (self.AtmMode == 2):
            return

        # radiance smapling angle
        logging.info("Radiance at the bottom atmospheric boundary")
        logging.info("ntha, th: # of angle and zenith angle(degree,max 18)")
        logging.info("ex. 5 100. 120. 140. 160. 170.")
        if (config.INPUT_MODE):
            th = list(args.get("zenith_angle"))
            ntha = len(th)
            logging.info(str(ntha) + ": " + str(th))
            th.insert(0, -1)
        else:
            ntha = int(input())
            for i in range(1, ntha + 1):
                th[i] = int(input())

        logging.info("npha, ph: # of angle and azimuth angle(degree,max 36)")
        logging.info("ex. 3 0. 90. 180.")
        if (config.INPUT_MODE):
            ph = list(args.get("azimuth_angle"))
            npha = len(ph)
            logging.info(str(npha) + ": " + str(ph))
            ph.insert(0, -1)
        else:
            npha = int(input())
            for i in range(1, npha + 1):
                ph[i] = int(input())

        k = 0
        for i in range(1, npha + 1):
            for j in range(1, ntha + 1):
                if (th[j] < 91):
                    logging.warning("Zenith angle should be greater than 91")
                    logging.warning(str(th[j]) + " is ignored !")
                else:
                    k = k + 1
                    comm.UX_RTAB[k] = sin(radians(th[j])) * cos(radians(ph[i]))
                    comm.UY_RTAB[k] = sin(radians(th[j])) * sin(radians(ph[i]))
                    comm.UZ_RTAB[k] = cos(radians(th[j]))

        comm.N_RDC = k
        return ERRCODE.SUCCESS

    def readAtmParameters(self, **args):
        nspc = 0
        ctop = cbot = 0.0
        parF = 531.2593
        parQ = 2423.93994
        swF = 1351.81531
        swQ = 10213.1367

        wl0d = [0.0]
        spcdf = [0.0]
        spcdq = [0.0]

        ch = ["", "hi", "lo", "lo"]

        if (self.AtmMode == 2):
            # "Only monochro wavelength calculation!"
            self.imode = 1
            self.RF = 1000.0
            self.RQ = 1000.0
            self.wq = 1.0
            self.nwl = 1

        else:
            logging.info("imode: Integration mode 1:Monochro 2:PAR 3:Snow")
            if (config.INPUT_MODE):
                self.imode = args.get("integration_mode")
                logging.info(self.imode)
            else:
                self.imode = int(input())

            if ((self.imode < 0) or (self.imode > 4)):
                logging.ERROR("Mode ERR: number should less than 4 and larger than 0.")
                return ERRCODE.INPUT_ERROR

            elif (self.imode == 1):
                logging.info("wl0:wavelength (micron ex 0.55)")
                self.wl0 = float(input())
                self.wls = 0.2
                self.nwl = int(round((4.0 - 0.3) / self.span[1]))
                self.npl[1] = self.nPhotonProcess

            elif (self.imode == 2):
                self.wls = 0.4
                self.nwl = int(round((0.7 - 0.4) / self.span[1]))
                self.RF = parF
                self.RQ = parQ

            elif (self.imode == 3):
                self.wls = 0.3
                self.nwl = int(round((0.7 - 0.3) / self.span[1]))
                self.nwl += int((4.0 - 0.7) / (5.0 * self.span[1]))
                self.RF = swF
                self.RQ = swQ

            # read solar radiation
            contentFile = np.loadtxt("../data/solar_rad")
            nspc = len(contentFile)
            for i in range(nspc):
                wl0d.append(contentFile[i][0])
                spcdf.append(contentFile[i][1])
                spcdq.append(contentFile[i][2])

            # search a spectral irradiance data in monochromatic calculation
            j = 1
            wlsd = self.wls + 0.0005

            if (self.imode == 1):
                for i in range(1, nspc + 1):
                    if ((wlsd > self.wl0 - 0.0025) and (wlsd < self.wl0 + 0.0025)):
                        self.RF = (wl0d[i + 1] * spcdf[i] + wl0d[i] * spcdf[i + 1]) / (wl0d[i] + wl0d[i + 1])
                        self.RQ = (wl0d[i + 1] * spcdq[i] + wl0d[i] * spcdq[i + 1]) / (wl0d[i] + wl0d[i + 1])
                        self.spcf[j] = self.RF
                        self.spcq[j] = self.RQ
                        self.wq = 1.0
                        self.nwl = 1
                        break
                    wlsd += 0.001

            elif (self.imode == 2):
                for i in range(1, nspc + 1):
                    if (abs(wl0d[i] - wlsd) < 1.0E-4):
                        for k in range(20):
                            self.spcf[j] += spcdf[i + k]
                            self.spcq[j] += spcdq[i + k]

                        self.spcf[j] /= self.span[1] * 1000.0
                        self.spcq[j] /= self.span[1] * 1000.0
                        j += 1

                        if (j> self.nwl):
                            break
                        wlsd += self.span[1]

            elif (self.imode == 3):
                for i in range(1, nspc + 1):
                    if (abs(wl0d[i] - wlsd) < 1.0E-4):
                        if (wlsd < 0.7):
                            for k in range(20):
                                self.spcf[j] += spcdf[i + k]
                                self.spcq[j] += spcdq[i + k]

                            self.spcf[j] /= self.span[1] * 1000
                            self.spcq[j] /= self.span[1] * 1000

                        else:
                            for k in range(100):
                                self.spcf[j] += spcdf[i + k]
                                self.spcq[j] += spcdq[i + k]

                            self.spcf[j] /= self.span[1] * 5000
                            self.spcq[j] /= self.span[1] * 5000

                        j += 1
                        if (j > self.nwl):
                            break

                        if (wlsd < 0.7):
                            wlsd += self.span[1]
                        else:
                            wlsd += 5.0 * self.span[1]

            # make input photon weight for each wavelength
            sum = 0.0

            for i in range(1, self.nwl + 1):
                if (i <= 20):
                    sum += self.spcf[i]
                else:
                    sum += self.spcf[i] * 5.0

            for i in range(1, self.nwl + 1):
                if (i <= 20):
                    self.npl[i] = int(round(float(self.nPhotonProcess) * self.spcf[i] / sum))
                else:
                    self.npl[i] = int(round(float(self.nPhotonProcess) * 5.0 * self.spcf[i] / sum))

            sum = 0.0

            for i in range(1, self.nwl + 1):
                sum += self.npl[i]

            self.nPhotonProcess = int(sum)
            self.nPhoton = self.nPhotonProcess * self.Nprocess
            logging.info("Actural number of photon is: " + str(self.nPhoton))

            # with/without atmospheric
            if (self.AtmMode == 2):
                self.atmType = 0
                self.aerosolType = 0
                self.taur = 0.0
                ctype = 0

            else:
                # read z profile
                comm.Z_GRD = np.loadtxt("../data/zgrd")

                # rescaling of zgrd according to the elevation
                for i in range(comm.N_Z + 1):
                    comm.Z_GRD_BACK[i] = comm.Z_GRD[i]
                    comm.Z_GRD[i] = (comm.Z_GRD[comm.N_Z] - comm.Z_MIN) * (comm.Z_GRD[i] / comm.Z_GRD[comm.N_Z]) + comm.Z_MIN
                    comm.K_LAYER[i] = i

                # mapping the iz to the actual height level corrected above
                for i in range(comm.N_Z + 1):
                    for j in range(comm.N_Z + 1):
                        if (comm.Z_GRD[i] < comm.Z_GRD_BACK[j]):
                            comm.K_LAYER[i] = j - 1
                            break

                # make the middle of the layer height of the original zgrd (zgrd_back)
                for i in range(1, comm.N_Z + 1):
                    comm.Z_GRD_M[i] = 0.5 * (comm.Z_GRD_BACK[i] + comm.Z_GRD_BACK[i - 1])

                logging.info("atmType: Atmospheric profile")
                logging.info(" 1: Tropical")
                logging.info(" 2: Mid latitude summer")
                logging.info(" 3: Mid latitude winter")
                logging.info(" 4: High latitude summer")
                logging.info(" 5: High latitude winter")
                logging.info(" 6: US standard atm.")

                if (config.INPUT_MODE):
                    self.atmType = args.get("atmosphere_type")
                    logging.info(self.atmType)
                else:
                    self.atmType = int(input())

                if (self.atmType == 1):
                    self.rfname = "../data/gas_TR_" + ch[self.imode]
                elif (self.atmType == 2):
                    self.rfname = "../data/gas_MS_" + ch[self.imode]
                elif (self.atmType == 3):
                    self.rfname = "../data/gas_MW_" + ch[self.imode]
                elif (self.atmType == 4):
                    self.rfname = "../data/gas_HS_" + ch[self.imode]
                elif (self.atmType == 5):
                    self.rfname = "../data/gas_HW_" + ch[self.imode]
                elif (self.atmType == 6):
                    self.rfname = "../data/gas_US_" + ch[self.imode]
                else:
                    logging.ERROR("Input error!")
                    return ERRCODE.INPUT_ERROR

                # read the aerosol data
                logging.info("aerosolType: aerosol type")
                logging.info(" 1:  Continental clean")
                logging.info(" 2:  Continental average")
                logging.info(" 3:  Continental polluted")
                logging.info(" 4:  Urban")
                logging.info(" 5:  Desert")
                logging.info(" 6:  Maritime clean")
                logging.info(" 7:  Maritime polluted")
                logging.info(" 8:  Maritime Tropical")
                logging.info(" 9:  Arctic")
                logging.info("10:  Antactic")
                logging.info("11:  Smoke")

                if (config.INPUT_MODE):
                    self.aerosolType = args.get("aerosol_type")
                    logging.info(self.aerosolType)
                else:
                    self.aerosolType = int(input())

                #AOT
                logging.info("tauref: AOT at 0.550 micron")
                logging.info(" - this version uses a extinction with RH=70% -")
                if (config.INPUT_MODE):
                    self.taur = args.get("AOT")
                    logging.info(self.taur)
                else:
                    self.taur = float(input())

                #current version uses a mixed extinction coef. by Iwabushi san
                self.nmix = 2

                if (self.aerosolType == 1):
                    self.fname.append("../data/opt_type1_rh0.70_" + ch[self.imode])
                    self.d = 8000.0     # scale height (m)
                elif (self.aerosolType == 2):
                    self.fname.append("../data/opt_type2_rh0.70_" + ch[self.imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 3):
                    self.fname.append("../data/opt_type3_rh0.70_" + ch[self.imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 4):
                    self.fname.append("../data/opt_type4_rh0.70_" + ch[self.imode])
                    self.d = 8000.0  # scale height (m)
                elif (self.aerosolType == 5):
                    self.fname.append("../data/opt_type5_rh0.70_" + ch[self.imode])
                    self.d = 2000.0  # scale height (m)
                elif ((self.aerosolType == 6) or (self.aerosolType == 8)):
                    self.fname.append("../data/opt_type6_rh0.70_" + ch[self.imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 7):
                    self.fname.append("../data/opt_type7_rh0.70_" + ch[self.imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 8):
                    self.fname.append("../data/opt_type8_rh0.70_" + ch[self.imode])
                    self.d = 1000.0  # scale height (m)
                elif (self.aerosolType == 9):
                    self.fname.append("../data/opt_type9_rh0.70_" + ch[self.imode])
                    self.d = 99000.0  # scale height (m)
                elif (self.aerosolType == 10):
                    self.fname.append("../data/opt_type10_rh0.70_" + ch[self.imode])
                    self.d = 99000.0  # scale height (m)
                elif (self.aerosolType == 11):
                    logging.ERROR("smoke aerosol under construction !!")
                    self.d = 8000.0  # scale height (m)
                    return ERRCODE.INPUT_ERROR
                else:
                    logging.ERROR("Input error !!")
                    return ERRCODE.INPUT_ERROR

                # read cloud data
                logging.info("ctype: cloud type")
                logging.info(" 0:  Cloud-free")
                logging.info(" 1:  Stratus continental")
                logging.info(" 2:  Stratus maritime")
                logging.info(" 3:  Cumulus continental clean")
                logging.info(" 4:  Culumus continental pulluted")
                logging.info(" 5:  Culumus maritime")
                logging.info(" 6:  Fog")
                logging.info(" 7:  Cirrus 1 (-25degC)")
                logging.info(" 8:  Cirrus 2 (-50 degC)")
                logging.info(" 9:  Cirrus 3 (-50 degC + small particles)")

                if (config.INPUT_MODE):
                    self.cloudType = args.get("cloud_type")
                    logging.info(self.cloudType)
                else:
                    self.cloudType = int(input())
                comm.CloudType = self.cloudType

                if (self.cloudType != 0):
                    logging.info("ctauref:COT at 0.55 micron")
                    self.ctaur = float(input())
                    logging.info("cloud top height (m)")
                    ctop = float(input())
                    logging.info("cloud bottom height (m)")
                    cbot = float(input())

                    self.cflg = 1

                    if (ctop < cbot):
                        logging.info("cloud top should be greater than cloud bottom")
                        return ERRCODE.INPUT_ERROR

                    if (self.cloudType == 1):
                        self.fname.append("../data/opt_type101_" + ch[self.imode])
                    elif (self.cloudType == 2):
                        self.fname.append("../data/opt_type102_" + ch[self.imode])
                    elif (self.cloudType == 3):
                        self.fname.append("../data/opt_type103_" + ch[self.imode])
                    elif (self.cloudType == 4):
                        self.fname.append("../data/opt_type104_" + ch[self.imode])
                    elif (self.cloudType == 5):
                        self.fname.append("../data/opt_type105_" + ch[self.imode])
                    elif (self.cloudType == 6):
                        self.fname.append("../data/opt_type106_" + ch[self.imode])
                    elif (self.cloudType == 7):
                        self.fname.append("../data/opt_type107_" + ch[self.imode])
                    elif (self.cloudType == 8):
                        self.fname.append("../data/opt_type108_" + ch[self.imode])
                    elif (self.cloudType == 9):
                        self.fname.append("../data/opt_type109_" + ch[self.imode])
                    else:
                        logging.ERROR("Input error !!")
                        return ERRCODE.INPUT_ERROR

                    for i in range(comm.N_Z + 1):
                        if (comm.Z_GRD[i] < cbot):
                            self.cbnz = i
                        if (comm.Z_GRD[i] > ctop):
                            self.ctnz = i
                            break

                    logging.info(comm.Z_GRD[self.cbnz], comm.Z_GRD[self.ctnz])
                    logging.info("clouds are located between " + str(self.ctnz) + " and " + str(self.cbnz))
        return ERRCODE.SUCCESS

    def readParameters(self, **args):

        logging.info("np: Input number of photon")
        if (config.INPUT_MODE):
            self.nPhoton = args.get("number_photon")
            logging.info(self.nPhoton)
        else:
            self.nPhoton = int(input())

        # fish eye simulation mode
        if (self.nPhoton == -4):
            if (self.Nprocess > 1):
                logging.ERROR("fish eye mode cannot work with multi-processors")
                logging.ERROR("please execute under single processor")
                return ERRCODE.INPUT_ERROR

            logging.info("fish eye mode - selected")
            self.readVegParameters(**args)
            self.surfaceType = 2
            return ERRCODE.SUCCESS

        # LAI calculation mode
        if (self.nPhoton == -5):
            if (self.Nprocess > 1):
                logging.ERROR("LAI mode cannot work with multi-processors")
                logging.ERROR("please execute under single processor")
                return ERRCODE.INPUT_ERROR

            logging.info("LAI calculation mode - selected")
            self.readVegParameters()
            self.surfaceType = 2
            return ERRCODE.SUCCESS

        self.nPhotonProcess = int(self.nPhoton / self.Nprocess)

        # Atmospheric mode
        logging.info("AtmMode: 1:atmospheric mode, 2: without atmospheric")
        if (config.INPUT_MODE):
            self.AtmMode = args.get("atm_type")
            logging.info(self.AtmMode)
        else:
            self.AtmMode = int(input())

        if (self.AtmMode == 2):
            logging.info("diffuse: Input frac. of diffuse radiation (0-1)")
            self.diffuse = int(input())
            if ((self.diffuse < 0.0) or (self.diffuse > 1.0)):
                logging.ERROR("fraction of diffuse should be 0.0-1.0... exit ")
                return ERRCODE.INPUT_ERROR

        # Get geometrical parameters
        self.readGeoParameters(**args)

        # Get atmospheric parameters
        self.readAtmParameters(**args)

        # vegetation parameters read / initialize
        logging.info("stye: Surface mode 1: Lambertian, 2: 3-D Vegetation")
        if (config.INPUT_MODE):
            self.surfaceType = args.get("surface_type")
            logging.info(self.surfaceType)
        else:
            self.surfaceType = int(input())

        if (self.surfaceType == 1):
            for iwl in range(1, self.nwl):
                if (self.imode == 1):
                    logging.info("Input albedo")
                else:
                    if (iwl <= 20):
                        ab1 = self.wls + self.span[iwl] * float(iwl -1)
                        ab2 = self.wls + self.span[iwl] * float(iwl)
                    else:
                        ab1 = 0.7 + self.span[iwl] * float(iwl - 21)
                        ab2 = 0.7 + self.span[iwl] * float(iwl - 20)
                    logging.info("Input albedo: " + str(ab1) + " and " + str(ab2))
                self.alb[iwl] = float(input())

        elif (self.surfaceType == 2):
            self.readVegParameters(**args)

        else:
            logging.ERROR("Bad surface type selection exit")
            return ERRCODE.INPUT_ERROR

        return ERRCODE.SUCCESS

    def preparAtm(self, iWL):

        if (self.imode != 1):
            # weight of photon for PPFD
            self.wq = self.span[iWL] * self.spcq[iWL] / self.RQ
            self.wq *= float(self.nPhotonProcess) / float(self.npl[iWL])

            if(iWL <= 20):
                self.wl0 = self.wls + float(iWL - 1) * self.span[iWL] + self.span[iWL] / 2.0
            else:
                self.wl0 = 0.7 + float(iWL - 21) * self.span[iWL] + self.span[iWL] / 2.0

        # get atomospheric parameters
        ext_back = np.zeros(10 * 200, dtype=float).reshape(10, 200)
        absg1d_back = np.zeros(200 * 4, dtype=float).reshape(200, 4)
        zmed = 0
        with open(self.rfname, "r") as file:
            for i in range(740):
                file.readline() # read: ##---

                result = file.readline()
                result = result.split()
                wl1 = float(result[1])
                wl2 = float(result[2])
                file.readline()

                for iz in range(1, comm.N_Z + 1):
                    result = file.readline()
                    result = result.split()

                    # #print("i=", i)
                    # print(iz)
                    # print(self.ext[1, iz])
                    self.ext[1, iz] = float(result[1])

                    for j in range(self.nkd):
                        comm.ABS_G1D[iz, j + 1] = float(result[2 + j])
                    for j in range(self.nkd):
                        self.wkd[j + 1] = float(result[5 + j])

                    # the ext(1,iz) and absg1d values are adjusted according to the
                    # actual height level
                    # caution!! this is very rough mapping, I didn't do a extrapolation
                    # I just used the values from closest layers.

                    for iz in range(1, comm.N_Z + 1):
                        ext_back[1, iz] = self.ext[1, iz]
                        nl = comm.K_LAYER[iz]
                        zmed = 0.5 * (comm.Z_GRD[iz] + comm.Z_GRD[iz - 1])
                        # print("iz = "+ str(iz))
                        # print("nl = " + str(nl))
                        self.ext[1, iz] = (comm.Z_GRD_M[nl + 1] - zmed) * self.ext[1, nl]
                        self.ext[1, iz] += (zmed - comm.Z_GRD_M[nl]) * self.ext[1, nl + 1]
                        self.ext[1, iz] /= (comm.Z_GRD_M[nl + 1] - comm.Z_GRD_M[nl])

                        if (iz == comm.N_Z):
                            self.ext[1, iz] = self.ext[1, comm.K_LAYER[iz]]

                        for j in range(1, 4):
                            absg1d_back[iz, j] = comm.ABS_G1D[iz, j]
                            comm.ABS_G1D[iz, j] = (comm.Z_GRD_M[nl + 1] - zmed) * comm.ABS_G1D[nl, j]
                            comm.ABS_G1D[iz, j] += (zmed - comm.Z_GRD_M[nl]) * comm.ABS_G1D[nl + 1, j]
                            comm.ABS_G1D[iz, j] /= (comm.Z_GRD_M[nl + 1] - comm.Z_GRD_M[nl])

                            if (iz == comm.N_Z):
                                comm.ABS_G1D[iz, j] = comm.ABS_G1D[comm.K_LAYER[iz], j]

                if ((wl2 > self.wl0) and (self.wl0 >= wl1)):
                    self.wkd[self.nkd] = 1.000001
                    break
            file.close()

        # read optic file
        self.Qext_ref = [0.0] * 10
        self.Qext = [0.0] * 10
        self.Qabs = [0.0] * 10
        self.omg = [0.0] * 10
        G = [0.0] * 10
        self.ang = [0]
        #dphs = [0]
        self.phs = np.zeros(10 * 201, dtype=float).reshape(10, 201)

        for i in range(2, self.nmix + 1):
            with open(self.fname[i - 1], "r") as file:

                if (self.imode == 1):

                    for j in range(740):
                        file.readline()

                        result = file.readline().split()
                        wl1 = float(result[0])
                        dum = float(result[1])

                        for ii in range(4):
                            file.readline()

                        result = file.readline().split()
                        self.Qext_ref[i - 1] = float(result[0])
                        self.Qext[i - 1] = float(result[1])
                        self.Qabs[i - 1] = float(result[2])
                        self.omg[i - 1] = float(result[3])
                        G[i - 1] = float(result[4])

                        file.readline()

                        result = file.readline().split()
                        nAng = float(result[0])

                        file.readline()

                        for k in range(nAng):
                            result = file.readline().split()
                            self.ang.append(float(result[0]))
                            self.phs[i, k + 1] = float(result[1])

                        if ((wl1 + 0.005) > self.wl0):
                            file.readline()

                            result = file.readline()
                            result = result.split()
                            wl12 = float(result[0])
                            dum = float(result[1])

                            for ii in range(4):
                                file.readline()

                            result = file.readline().split()
                            self.Qext_ref[i - 1] = float(result[0])
                            self.Qext = float(result[1])
                            self.Qabs = float(result[2])
                            self.omg = float(result[3])
                            dG = float(result[4])

                            file.readline()

                            result = file.readline().split()
                            nAng = float(result[0])

                            file.readline()

                            for k in range(nAng):
                                result = file.readline().split()
                                self.ang.append(float(result[0]))
                                dphs.append(float(result[1]))

                            self.Qext[i - 1] = ((self.wl0 - wl1) * self.Qext + (wl2 - self.wl0) * self.Qext[i + 1]) \
                                          * (1.0 / 0.005)
                            self.Qabs[i - 1] = ((self.wl0 - wl1) * self.Qabs + (wl2 - self.wl0) * self.Qabs[i + 1]) \
                                          * (1.0 / 0.005)
                            self.omg[i] = ((self.wl0 - wl1) * self.omg + (wl2 - self.wl0) * self.omg[i + 1]) \
                                          * (1.0 / 0.005)
                            G[i - 1] = ((self.wl0 - wl1) * dG + (wl2 - self.wl0) * G[i + 1]) \
                                          * (1.0 / 0.005)

                            for k in range(nAng):
                                self.phs[i, k + 1] = ((self.wl0 - wl1) * dphs + (wl2 - self.wl0) * self.phs[i + 1, k + 1]) \
                                          * (1.0 / 0.005)

                            break

                else:
                    for j in range(53):
                        file.readline()

                        result = file.readline().split()
                        wl1 = float(result[0])
                        wl2 = float(result[1])

                        # skip reading RvDry, ReDry, rhoDry, RvWet, ReWet, rhoWet
                        # because of not using
                        for ii in range(4):
                            file.readline()

                        result = file.readline().split()
                        self.Qext_ref[i - 1] = float(result[0])
                        self.Qext[i - 1] = float(result[1])
                        self.Qabs[i - 1] = float(result[2])
                        self.omg[i] = float(result[3])
                        G[i - 1] = float(result[4])

                        file.readline()
                        result = file.readline().split()
                        nAng = int(result[0])
                        comm.N_ANG = nAng

                        file.readline()

                        for k in range(nAng):
                            result = file.readline().split()
                            self.ang.append(float(result[0]))
                            self.phs[i, k + 1] = float(result[1])

                        if ((self.wl0 > wl1) and (self.wl0 < wl2)):
                            break

                file.close()

        # make rayleigh scattering albedo
        self.omg[1] = 1.0

        fac = 3.0 / (16.0 * pi)

        #print("nang = ", nAng)
        for i in range(nAng):
            self.phs[1, i + 1] = fac * (1.0 + cos(radians(self.ang[i + 1])) ** 2)

        # get the ratio of extinction corf. of each aerosol
        rat = (0, 1.0, 2.0)
        zmd = [0]

        # Calcualte the extinction coef. from total tau
        for i in range(1, comm.N_Z + 1):
            temp = (0.5 * (exp(- comm.Z_GRD[i] / self.d) + exp(- comm.Z_GRD[i - 1] / self.d)))
            temp = -self.d * log(temp)
            zmd.append(temp)

        # Calculate the scale factor (sfc)
        sfc = [0.0] * 10
        zsum = 0.0
        for i in range(1, comm.N_Z + 1):
            zsum += exp(- zmd[i] / self.d)

        for imix in range(2, self.nmix + 1):
            sfc[imix - 1] = rat[imix - 1] * self.taur / zsum

        for imix in range(2, self.nmix + 1):
            for i in range(1, comm.N_Z + 1):
                self.ext [imix, i] = sfc[imix - 1] * exp(- zmd[i] / self.d) /\
                                         (comm.Z_GRD[i] - comm.Z_GRD[i - 1])
                # Convert to ext from ext_ref
                self.ext[imix, i] *= self.Qext[imix - 1] / self.Qext_ref[imix - 1]

        if (self.cflg == 1):
            for i in range(self.cbnz + 1, self.ctnz + 1):
                self.ext[self.nmix + self.cflg, i] = self.ctaur /\
                                            (comm.Z_GRD[self.ctnz] - comm.Z_GRD[self.cbnz])
                self.ext[self.nmix + self.cflg, i] *= self.Qext(self.nmix) / \
                                             self.Qext_ref[self.nmix]

        return ERRCODE.SUCCESS

    def printFishEye():
        logging.warning("Haven't finish finish eye function.")
        return ERRCODE.SUCCESS

    def writeData(self):

        app = [0.0] * int(comm.Z_MAX + 2)
        # back to 300
        comm.SIZE -= 1
        # summary of the radiative flux

        # summary calculation of fulx
        pixelNP = float(self.nPhoton) / float(comm.SIZE ** 2)
        tm = self.tflx / float(self.nPhoton)

        for i in range(1, self.nwl + 1):
            comm.scmpf[1, i] *= 100.0 / self.tflx
            comm.scmpp[1, i] *= 100.0 / self.tpfd
            comm.scmpf[2, i] *= 100.0 / self.bflx
            comm.scmpp[2, i] *= 100.0 / self.bpfd
            comm.scmpf[3, i] *= 100.0 / self.dflx
            comm.scmpp[3, i] *= 100.0 / self.dpfd

        self.tflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)
        self.bflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)
        self.dflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)

        self.rflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)
        self.rbflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)
        self.rdflx *= self.RF * abs(self.cos_q0) / float(self.nPhoton)

        self.tpfd *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        self.bpfd *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        self.dpfd *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)

        comm.C_FPR *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        comm.B_FPR *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        comm.F_FPR *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        comm.S_FPR *= self.RQ * abs(self.cos_q0) / float(self.nPhoton)
        comm.T_FPR = comm.C_FPR + comm.B_FPR + comm.F_FPR

        comm.T_FPR /= self.tpfd
        comm.C_FPR /= self.tpfd
        comm.B_FPR /= self.tpfd
        comm.F_FPR /= self.tpfd
        comm.S_FPR /= self.tpfd

        # summary of 3D far etc
        th = acos(self.cos_q0)
        ith = int(degrees(th))

        for k in range(int(comm.Z_MAX) + 1, 0, -1):
            for i in range(1, comm.SIZE + 1):
                for j in range(1, comm.SIZE + 1):
                    comm.AP[i, j, k] *= self.RQ * abs(self.cos_q0) / pixelNP
                    comm.AP_D[i, j, k] *= self.RQ * abs(self.cos_q0) / pixelNP
                    comm.AP_B[i, j, k] *= abs(self.cos_q0) / pixelNP
                    comm.AP_B[i, j, k] *= (1.0 / comm.GT_BLC[ith])
                    comm.AP_B[i, j, k] = min(comm.AP_B[i, j, k], 1.0)
                    app[k] += comm.AP[i, j, k]

            app[k] /= (self.tflx * float(comm.SIZE ** 2))
            comm.AP_NP[k] *= self.RQ * abs(self.cos_q0) / float(self.nPhoton) / self.tflx

        # summary of forest floor far etc
        for i in range(1, comm.SIZE + 1):
            for j in range(1, comm.SIZE + 1):
                comm.AP_F[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP
                comm.AP_FD[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP

        # summary of forest floor and surface downward flux
        ffSum = 0.0
        sfSum = 0.0

        for i in range(1, comm.SIZE + 1):
            for j in range(1, comm.SIZE + 1):
                comm.SF_DIR[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP
                comm.SF_DIF[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP
                comm.FF_DIR[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP
                comm.FF_DIF[i, j] *= self.RQ * abs(self.cos_q0) / pixelNP

                ffSum += comm.FF_DIR[i, j] + comm.FF_DIF[i, j]
                sfSum += comm.SF_DIR[i, j] + comm.SF_DIF[i, j]

        ffSum /= float(comm.SIZE ** 2)
        sfSum /= float(comm.SIZE ** 2)
        # writting
        f = open(config.OUTPUT_PATH + "flxsum.txt", "w")
        f.write("-- summary of radiative quantity --" + "\n")
        f.write("")
        f.write("-- simulation mode --" + "\n")

        if (self.AtmMode == 1):
            f.write(" Atmospheric module:      Yes" + "\n")
        else:
            f.write(" Atmospheric module:      No" + "\n")

        if (self.AtmMode == 1):
            if (self.imode == 1):
                f.write("Spectral domain: single wavelength, " + str(self.wl0) + "\n")
            elif (self.imode == 2):
                f.write("Spectral domain:       PAR" + "\n")
            elif (self.imode == 3):
                f.write("Spectral domain:       Shortwave total" + "\n")

        if (self.surfaceType == 1):
            f.write("Surface type:          lambertian" + "\n")
        else:
            f.write("Surface type:          3-D canopy" + "\n")


        f.write("\n" + format("Solar zenith angle", '22') + format("Elv(m)", '12') + format("atm. transmittance", '20')+ "\n")
        f.write(format(degrees(pi - acos(self.cos_q0)), '12.2f') + " " +
                format(comm.Z_MIN, '17.5f') + "  " +
                format(tm, '12.5f') + "\n")

        f.write("\n-- Downward Irradiance at TOC --")

        if (self.AtmMode == 2):
            f.write("\nwarning ! irradiance & PPFD are normalized by 1000.*cos(th) !!"+ "\n")

        f.write("\nEnergy unit (W/m2)"+ "\n")
        f.write(format("Total", '12') + format("Beam", '12') + format("Diffuse", '12') + format("Diff.Ratio", '12') + "\n")
        f.write(format(self.tflx, '8.5f') + format(self.bflx, '12.5f') + format(self.dflx, '12.5f') + format(
            self.dflx / self.tflx, '12.5f') + "\n")

        f.write("\nPhoton unit (W/m2)"+ "\n")
        f.write(format("Total", '12') + format("Beam", '12') + format("Diffuse", '12') + format("Diff.Ratio", '12') + "\n")
        f.write(format(self.tpfd, '8.5f') + format(self.bpfd, '12.5f') + format(self.dpfd, '12.5f') + format(
            self.dpfd / self.tpfd, '12.5f') + "\n")

        f.write("\nAlbedo (A)"+ "\n")
        f.write(format("Actual", '12') + format("black", '12') + format("white", '12') + "\n")
        f.write(format("", '12') + format("(beam)", '12') + format("(diffuse)", '12') + "\n")
        f.write(format(self.rflx / self.tflx, '8.5f') + format(self.rbflx / self.bflx, '11.5f') + format(
            self.rdflx / self.dflx, '11.5f') + "\n")

        if (self.surfaceType == 2):
            f.write("\n-- fraction of downward flux (T) -- "+ "\n")
            f.write(format("Forest Floor", '15') + format("Surface Floor", '15') + "\n")
            f.write(format(ffSum / self.tpfd, '12.5f') + format(sfSum / self.tpfd, '16.5f') + "\n")

            f.write("\n-- fraction of absorbed radiation (far) --"+ "\n")
            f.write(format("total", '13') + format("leaf", '13') + format("non-leaf", '15') + format("floor", '13') + "\n")
            f.write(format(comm.T_FPR, '10.7f') + format(comm.C_FPR, '12.7f') + format(comm.B_FPR, '14.7f') + format(
                comm.F_FPR, '13.7f') + "\n")

            f.write("\nVertical profile of far" + "\n")
            f.write(format("H(m)", '12') + format("total", '12') + format("leaf", '12') + format("non-leaf", '8') + "\n")
            string = ""
            for k in range(int(comm.Z_MAX) + 1, 0, -1):
                string += format(k, '5') + format(app[k] + comm.AP_NP[k], '12.6f') + format(app[k], '14.6f') + \
                          format(comm.AP_NP[k], '13.6f') + "\n"
            f.write(string + "\n")

        f.write("\n-- spctrl fraction. energy & photon flux(%) --" + "\n")
        f.write(format("wl(um)", '10') + format("Etot", '10') + format("Ebeam", '10') + format("Ediffuse", '14')
                + format("Ptot", '12') + format("Pbeam", '10') + format("Pdiffuse", '10') + "\n")
        for i in range(1, self.nwl + 1):
            if (i <= 20):
                temp = self.wls + self.span[i] * float(i) - self.span[i] / 2.0
            else:
                temp = 0.7 + self.wls + self.span[i] * float(i - 20) - self.span[i] / 2.0

            f.write(format(temp, '5.3f'))
            string = ""
            for j in range(3):
                string += format(comm.scmpf[j + 1, i], '12.5f')
            for j in range(3):
                string += format(comm.scmpp[j + 1, i], '12.5f')

            f.write(string + "\n")

        f.write("\n-- Downward radiance at the top of canopy --\n")
        f.write(format("idrc", '10') + format("theta", '10') + format("phi", '9') + format("radiance_F", '13')
                + format("stdev_F", '10') + format("radiance_Q", '13') + format("radiance_Q", '13') + "\n")

        pixelNP_A = float(self.nPhoton) / float(comm.K_NXR * comm.K_NYR)

        for k in range(1, comm.N_RDC + 1):
            sumF = sum2F = 0.0
            pF_Max = 0.0
            sumQ = sum2Q = 0.0
            pQ_Max = 0.0

            for j in range(1, comm.K_NYR + 1):
                for i in range(1, comm.K_NXR + 1):
                    pF = comm.PROC_F[i, j, k] / pixelNP_A
                    pQ = comm.PROC_Q[i, j, k] / pixelNP_A

                    sumF += pF
                    sum2F += pF ** 2
                    pF_Max = max(pF_Max, pF)

                    sumQ += pQ
                    sum2Q += pQ ** 2
                    pQ_Max = max(pQ_Max, pQ)

            aveF = sumF / float(comm.K_NXR * comm.K_NYR)
            # print("sum2F / float(comm.K_NXR * comm.K_NYR) - aveF ** 2 = ", sum2F / float(comm.K_NXR * comm.K_NYR) - aveF ** 2)
            sdvF = sqrt(sum2F / float(comm.K_NXR * comm.K_NYR) - aveF ** 2)
            aveQ = sumQ / float(comm.K_NXR * comm.K_NYR)
            sdvQ = sqrt(sum2Q / float(comm.K_NXR * comm.K_NYR) - aveQ ** 2)

            uxc = copysign(max(1.0e-5, abs(comm.UX_RTAB[k])), comm.UX_RTAB[k])
            phi = sqrt(comm.UX_RTAB[k] ** 2 + comm.UY_RTAB[k])
            phi = max(1.0e-5, phi)
            phi = degrees(acos(uxc / phi))

            f.write(format(k, "4") + format(degrees(acos(comm.UZ_RTAB[k])), '11.3f') + format(phi, '9.3f') +
                    format(aveF, '12.5f') + format(sdvF, '12.5f') + format(aveQ, '12.5f') + format(sdvQ, '12.5f') + "\n")

        f.close()

        if (self.surfaceType == 1):
            return ERRCODE.SUCCESS

        if (self.cmode == 1):
            # BRDF result output
            for i in range(1, comm.N_ANG_C + 1):
                comm.BRF_C[1, i] /= comm.BRF[1, i]
                comm.BRF_S[1, i] /= comm.BRF[1, i]
                comm.BRF_F[1, i] /= comm.BRF[1, i]

                if (self.AtmMode == 1):
                    comm.BRF[1, i] *= pi / (float(self.nPhoton) * tm)
                else:
                    comm.BRF[1, i] *= pi / float(self.nPhoton)

                comm.BRF_C[2, i] /= comm.BRF[2, i]
                comm.BRF_S[2, i] /= comm.BRF[2, i]
                comm.BRF_F[2, i] /= comm.BRF[2, i]
                comm.BRF[2, i] *= pi / float(self.nPhoton)

            f = open(config.OUTPUT_PATH + "brfsum.txt", 'w')

            if (self.AtmMode == 1):
                f.write(format("Theta", '12') + format("Phi", '6') + format("BRF(TOC)", '15') + format("BRF(TOA)",'15') +
                        format("rover", '13') + format("rbark", '13') + format("rfloor", '12') + "\n")
            else:
                f.write(format("Theta", '12') + format("Phi", '6') + format("BRF(TOC)", '15') +
                        format("rover", '13') + format("rbark", '13') + format("rfloor", '12') + "\n")

            k = 1

            for j in range(1, comm.N_PH + 1):
                for i in range(1, comm.N_TH + 1):
                    if (self.AtmMode == 1):
                        f.write(format(comm.ANG_T[i], '8.2f') + format(comm.ANG_P[j], '8.2f') +
                                format(comm.BRF[1, k], '10.6f') + format(comm.BRF[2, k], '14.6f') +
                                format(comm.BRF_C[1, k], '14.6f') + format(comm.BRF_S[1, k], '14.6f') +
                                format(comm.BRF_F[1, k], '14.6f') + "\n")
                    else:
                        f.write(format(comm.ANG_T[i], '8.2f') + format(comm.ANG_P[j], '8.2f') +
                                format(comm.BRF[1, k], '10.6f') +
                                format(comm.BRF_C[1, k], '14.6f') + format(comm.BRF_S[1, k], '14.6f') +
                                format(comm.BRF_F[1, k], '14.6f') + "\n")

                    k += 1

            f.close()

        # Nadir view image
        if (self.cmode != 1):
            f = open(config.OUTPUT_PATH + "TOC_IMG.txt")

            for j in range(comm.SIZE, 0, -1):
                string = ""
                if (self.atmType == 1):
                    for i in range(1, comm.SIZE + 1):
                        string += format(pi * comm.REFL[1, i, j] / (pixelNP * tm), '10.5f')
                else:
                    for i in range(1, comm.SIZE + 1):
                        string += format(pi * comm.REFL[1, i, j] / pixelNP, '10.5f')
                f.write(s)

            f.close()

            f = open(config.OUTPUT_PATH + "nTOC_IMG.txt", 'w')
            for j in range(comm.SIZE, 0, -1):
                string = ""
                for i in range(1, comm.SIZE + 1):
                    string += format(pi * comm.I_REFL[1, i, j] / pixelNP, '10.5f')
            f.close()

        if (self.cmode == 3):
            f = open(config.OUTPUT_PATH + "apar.txt", 'w')
            for k in range(int(comm.Z_MAX) + 1, 0, -1):
                for j in range(comm.SIZE, 0, -1):
                    string = ""
                    for i in range(1, comm.SIZE + 1):
                        string += format(comm.AP[i, j, k], '12.5f')
                    f.write(string)
            f.close()

            f = open(config.OUTPUT_PATH + "aparb.txt", 'w')
            for k in range(int(comm.Z_MAX) + 1, 0, -1):
                for j in range(comm.SIZE, 0, -1):
                    string = ""
                    for i in range(1, comm.SIZE + 1):
                        string += format(comm.AP_B[i, j, k], '12.5f')
                    f.write(string)
            f.close()

            f = open(config.OUTPUT_PATH + "apard.txt", 'w')
            for k in range(int(comm.Z_MAX) + 1, 0, -1):
                for j in range(comm.SIZE, 0, -1):
                    string = ""
                    for i in range(1, comm.SIZE + 1):
                        string += format(comm.AP_D[i, j, k], '12.5f')
                    f.write(string)
            f.close()

            f = open(config.OUTPUT_PATH + "aparf.txt", 'w')
            for j in range(comm.SIZE, 0, -1):
                string = ""
                for i in range(1, comm.SIZE + 1):
                    string += format(comm.AP_F[i, j, k], '12.5f')
                f.write(string)
            f.close()

            f = open(config.OUTPUT_PATH + "aparfd.txt", 'w')
            for j in range(comm.SIZE, 0, -1):
                string = ""
                for i in range(1, comm.SIZE + 1):
                    string += format(comm.AP_FD[i, j, k], '12.5f')
                f.write(string)
            f.close()

        # print fish eye data
        if (self.nPhoton == -4):
            self.printFishEye()

        return ERRCODE.SUCCESS
