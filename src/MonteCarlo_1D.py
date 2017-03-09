import ERRCODE
import common as comm
from Position import Position
from math import *
from VegTrace import VegTrace

class MonteCarlo_1D:

    weight = -1
    tau = -1

    def __init__(self):
        self.weight = -1
        self.tau = -1

    def save(self):

        return ERRCODE.SUCCESS

    # Horizontal shift
    def arthshift(self, x, ux, path, xmax):

        x += path * ux

        if ((x < 0.0) or (x > xmax)):
            x = x % xmax
            if (x < 0.0):
                x = x + xmax

        return x

    def optics(self, phoCoord, vectCoord):
        return ERRCODE.SUCCESS

    # Traces a trajectory in plane-parallel vertically-inhomogeneous atmosphere
    def escape(self, phoCoord, vectCoord, iz, ichi, ikd, resultCoord):

        resultCoord.setPosition(phoCoord.x, phoCoord.y, phoCoord.z)

        # TAU integration
        tau = 0.0
        if (ikd == 0):
            return ERRCODE.SUCCESS

        if (vectCoord.z >= 0.0):
            for izr in range(iz, comm.N_Z):
                tau += (comm.Z_GRD(izr) - resultCoord.z) * (comm.ABS_G1D[izr, ikd] + comm.EXT_T1D[izr, ichi])
                resultCoord.z = comm.Z_GRD[izr]

            tau /= vectCoord.z
        else:
            for izr in range(iz, 0, -1):
                tau += (resultCoord.z - comm.Z_GRD(izr - 1)) * (comm.ABS_G1D[izr, ikd] + comm.EXT_T1D[izr, ichi])
                resultCoord.z = comm.Z_GRD[izr - 1]

            tau /= (-1.0 * vectCoord.z)

        # escape location
        fpath = (resultCoord.z - phoCoord.z) / abs(vectCoord.z)
        resultCoord.x = self.arthshift(resultCoord.x, vectCoord.x, fpath, comm.X_MAX)
        resultCoord.y = self.arthshift(resultCoord.y, vectCoord.y, fpath, comm.Y_MAX)

        return tau

    def trace(self, phoCoord, vectCoord):
        return ERRCODE.SUCCESS


