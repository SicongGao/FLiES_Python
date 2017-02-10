import common as comm
import math
import ERRCODE


def idivspace():
    intv = 50.0 / comm.RES

    comm.IX_MAX = 6
    comm.IY_MAX = 6
    comm.IZ_MAX = 20

    ix = iy = iz = 1
    xr = yr = zu = zb = 0.0
    a = [1.0, 1.0, 0.0, 0.0]
    b = [0.0, 0.0, 1.0, 1.0]

    # initialize the div array
    # has done in common.py

    idiv = comm.IX_MAX * comm.IY_MAX * comm.IZ_MAX

    print("idiv = ", idiv)

    # start the idiv loop
    for i in range(idiv):
        n = 0

        # preparation of the i-th grid for space divided method
        divX = (float(ix) - 1.0) * intv
        divY = (float(iy) - 1.0) * intv
        divZ = (float(iz) - 1.0) * intv

        c = [divX, divX + intv, divY, divY + intv, divZ, divZ + intv]

        # definition of rectangular of the i-th object
        for j in range(comm.N_OBJ):
            flag = 0

            if (comm.T_OBJ[j] == 1):                            # cone
                xr = comm.OBJ[j][0]
                yr = comm.OBJ[j][1]
                zu = comm.OBJ[j][2]
                zb = comm.OBJ[j][2] - comm.OBJ[j][3]
            elif (comm.T_OBJ[j] == 2) or (comm.T_OBJ[j] == 4):  # cylinder, 4 for trunk
                xr = comm.OBJ[j][0]
                yr = comm.OBJ[j][1]
                zu = comm.OBJ[j][2]
                zb = comm.OBJ[j][2] - comm.OBJ[j][3]            # zb is always 0 for trunk
            elif (comm.T_OBJ[j] == 3):                          # ellipsoid
                xr = comm.OBJ[j][0]
                yr = comm.OBJ[j][1]
                zu = comm.OBJ[j][2] + comm.OBJ[j][3]            # zb, zu is bottom and upper coordinate of objects
                zb = comm.OBJ[j][2] - comm.OBJ[j][3]            # if T_OBJ[3] == 3, T_OBJ[3] =/ 2 in iparam.py, line 280
            elif (comm.T_OBJ[j] == 5):                          # half ellipsoid
                xr = comm.OBJ[j][0]
                yr = comm.OBJ[j][1]
                zu = comm.OBJ[j][2] + comm.OBJ[j][3]
                zb = comm.OBJ[j][2]

            # check the intersection on the x-y plane
            for k in range(4):
                d = abs(xr * a[k] + yr * b[k] - c[k])

                if (d <= comm.OBJ[j][4]):   # (xr, yr) near the boundary
                    dd = math.sqrt(comm.OBJ[j][4] * comm.OBJ[j][4] - d * d)
                    p1 = b[k] * xr + a[k] * yr - dd
                    p2 = b[k] * xr + a[k] * yr + dd
                    iMin = b[k] * divX + a[k] * divY
                    iMax = b[k] * divX + a[k] * divY + intv

                    if ((iMin <= p1) and (max >= p1)):
                        flag = 1
                    if ((iMax <= p2) and (max >= p2)):
                        flag = 1

            for k in range(2):
                for l in range(2):
                    rx = xr - c[k]
                    ry = yr - c[l + 2]
                    rr = math.sqrt(rx * rx + ry * ry)
                    if (rr <= comm.OBJ[j][4]):
                        flag = 1

            if (flag == 0):
                if (((xr >= c[0]) and (xr <= c[1])) and
                        ((yr >= c[2]) and (yr <= c[3]))):
                    flag = 1

            # check the intersection for z - axis
            # s: object[j] inter the voxel
            if (flag == 1):
                if ((zu > c[4]) and (zu <= c[5])):
                    flag = 2
                elif ((zb >= c[4]) and (zb < c[5])):
                    flag = 2
                elif ((zu >= c[5]) and (zu <= c[4])):
                    flag = 2

            # input data number for the ndivs & divs
            # s: current voxel contains obj, then record
            if (flag == 2):
                comm.N_DIVS[i] += 1
                comm.DIVS[i][n] = j
                n += 1
                comm.M_DIV = i

        # count
        ix += 1
        if (ix > comm.IX_MAX):
            ix = 1
            iy += 1

            if (iy > comm.IY_MAX):
                iy = 1
                iz += 1

    # determination of zmax at the boundary of the big voxel
    comm.Z_MAX = intv * (1.0 + float((comm.M_DIV - 1) / (comm.IX_MAX * comm.IY_MAX)))

    print("N_DIVS")
    print(comm.N_DIVS)
    print("DIVS")
    print(comm.DIVS)

    return ERRCODE.SUCCESS
