import numpy as np
import math

# DEFINE CONST
OBJ_NUM_MAX = 6000
ANGLE_SHIFT = 65000
ACOOS_SHIFT = 10000

# space definition
SIZE = 300
N_DICS = [0] * 720
DIVS = np.zeros(720 * 300, dtype=int).reshape(720, 300)
M_DIV = 0
IX_MAX = IY_MAX = IZ_MAX = 0
X_MAX = Y_MAX = Z_MAX = 0
RES = EPSI = 0

# atmospheric common
K_NZ = 200
K_NKD = 3
K_NCHI = 10
K_NXR = K_NYR = 300
K_NRDC = 700
N_LUT = 2000
N_LUT1 = N_LUT + 1
N_ANG = 0
Z_GRD = [0.0] * K_NZ
CHI_GRD = [0.0] * 100
TRU_LUT = np.zeros(6 * K_NCHI * K_NZ, dtype=float).reshape(6, K_NCHI, K_NZ)
ABS_G1D = np.zeros(K_NZ * K_NKD, dtype=float).reshape(K_NZ, K_NKD)
EXT_T1D = np.zeros(K_NZ * K_NCHI, dtype=float).reshape(K_NZ, K_NCHI)
ABS_T1D = [0.0] * K_NZ
PF_LUT = np.zeros(N_LUT1 * K_NZ, dtype=float).reshape(N_LUT1, K_NZ)
SA_LUT = np.zeros(N_LUT1 * K_NZ, dtype=float).reshape(N_LUT1, K_NZ)
UX_RTAB = [0.0] * K_NRDC
UY_RTAB = [0.0] * K_NRDC
UZ_RTAB = [0.0] * K_NRDC
PROC_F = np.zeros(K_NXR * K_NYR * K_NRDC, dtype=float).reshape(K_NXR, K_NYR, K_NRDC)
PROC_Q = np.zeros(K_NXR * K_NYR * K_NRDC, dtype=float).reshape(K_NXR, K_NYR, K_NRDC)
WRR = FS_ANG = 0.0
N_Z = N_CHI = N_RDC = 0

# canopy common
N_ANG_C = N_TH = N_PH = N_OBJ = 0
I_OBJ = [0] * OBJ_NUM_MAX
S_OBJ = [0] * OBJ_NUM_MAX
DIF_TYPE = DIR_FLAG = BOUND = 0
REFL = np.zeros(2 * SIZE * SIZE, dtype=float).reshape(2, SIZE, SIZE)
I_REFL = np.zeros(2 * SIZE * SIZE, dtype=int).reshape(2, SIZE, SIZE)
M_C = M_B = M_F = N_TS = 0
OBJ = np.zeros(OBJ_NUM_MAX * 5, dtype=float).reshape(OBJ_NUM_MAX, 5)
S_BAR = [0.0] * 5
ANG_T = [0.0] * 100
ANG_P = [0.0] * 100
BRF = np.zeros(2 * 700, dtype=float).reshape(2, 700)
BRF_C = np.zeros(2 * 700, dtype=float).reshape(2, 700)
BRF_S = np.zeros(2 * 700, dtype=float).reshape(2, 700)
BRF_F = np.zeros(2 * 700, dtype=float).reshape(2, 700)
F_EYE = np.zeros(90 * 360, dtype=float).reshape(90, 360)
RF_EYE = np.zeros(90 * 360, dtype=float).reshape(90, 360)
AP = np.zeros(SIZE * SIZE * 100, dtype=float).reshape(SIZE, SIZE, 100)
AP_D = np.zeros(SIZE * SIZE * 100, dtype=float).reshape(SIZE, SIZE, 100)
AP_B = np.zeros(SIZE * SIZE * 100, dtype=float).reshape(SIZE, SIZE, 100)
AP_F = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
AP_S = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
AP_FD = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
AP_NP = [0.0] * 100
T_FPR = C_FPR = B_FPR = F_FPR = S_FPR = 0.0
FF_DIR = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
FF_DIF = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
SF_DIR = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
SF_DIF = np.zeros(SIZE * SIZE, dtype=float).reshape(SIZE, SIZE)
G_LAI = 0.0
U = [0.0] * 6
BAD = [0.0] * 5
IR = FE = BP1 = BP2 = 0.0
UX_RC = [0.0] * K_NRDC
UY_RC = [0.0] * K_NRDC
UZ_RC = [0.0] * K_NRDC
RB = 0.5

# scattering angle parameters
G_MRC = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)
G_MTC = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)
G_MRB = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)
G_MTB = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)
G_MRF = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)
G_MTF = np.zeros(19 * 19 * 37, dtype=float).reshape(19, 19, 37)

# BACK TRACE
N_B = B_G = 0

# MATH
T_SIN = [0.0] * ANGLE_SHIFT * 2
T_COS = [1.0] * ANGLE_SHIFT * 2
T_ACOS = [0.0] * ACOOS_SHIFT * 2
T_EXP = [0.0] * ACOOS_SHIFT
DLT = np.zeros((6 + 1) * (6 + 1), dtype=float).reshape((6 + 1), (6 + 1))

for i in range(0, 62832 * 2):
    T_SIN[i] = math.sin(float(i - 62832) * 0.0001)
    T_COS[i] = math.cos(float(i - 62832) * 0.0001)

for i in range(0, ACOOS_SHIFT * 2):
    T_ACOS[i] = math.acos(float(i - ACOOS_SHIFT) * 0.0001)

for i in range(1,6):
    DLT[i, i] = 1.0

def getSin(n):
    return T_SIN[int((n + 62382) * 1E4)]

def getCos(n):
    return T_COS[int((n + 62382) * 1E4)]

def getACos(n):
    return T_SIN[int((n + ACOOS_SHIFT) * 1E4)]
