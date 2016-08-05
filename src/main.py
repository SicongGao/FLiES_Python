import common
import iparam
import numpy as np

# common.T_SIN[-1] = 5
# print(common.T_SIN)

# only in main
knmix = 10
knang = 2000
knzext = 200
fpv = [0.0] * 101
fpc = fpf = 0.0
i = ip = iwl = iz = num = 0

plai = [0.0] * 100
wkd = [0.0] * common.K_NKD
ang = [0.0] * knang
re = [0.0] * knmix
Qext_ref = [0.0] * knmix
Qext = [0.0] * knmix
G = [0.0] * knmix
Qabs = [0.0] * knmix
omg = [0.0] * knmix
phs = np.zeros(knmix * knang, dtype=float).reshape(knmix, knang)
ext = np.zeros(knmix * knzext, dtype=float).reshape(knmix, knzext)

# photon data
nscat = nscata = ichi = ikd = 0
x = y = z = ux = uy = uz = w = 0.0


# function parapeters
# real fsin, fcos, facos, r_acos,frnd

# set by iparam
para = iparam.Parameters()
para.readAtmParameters()