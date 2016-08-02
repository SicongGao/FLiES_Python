import common
import iparam
import numpy as np
# common.T_SIN[-1] = 5
# print(common.T_SIN)

knmix = 10
knang = 2000
knzext = 200
ext = np.zeros(knmix * knzext, dtype=float).reshape(knmix, knzext)
fpv = [0.0] * 101
fpc = fpf = 0.0




a = [0] * 5
b = [0] *5

for i in range(0,4):
    a[i]= 2

for i in range(0, 4):
    b[i] = 2

print(a,b)

