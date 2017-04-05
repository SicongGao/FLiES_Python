from math import *
import ERRCODE


def transformCoordinate(oriVecCoord, angleA, angleB, destVecCoord):

    MIN_VALUE = 1.0e-15

    if (oriVecCoord.x == 0 and oriVecCoord.y == 0):
        print("coord got wrong, check TransformCoordinate.py, line. 8")
        exit()

    sinA = sin(angleA)
    sinB = sin(angleB)
    cosA = cos(angleA)
    cosB = cos(angleB)
    sinT = sqrt(oriVecCoord.x ** 2 + oriVecCoord.y ** 2)
    sinP = oriVecCoord.y / sinT
    cosP = oriVecCoord.x / sinT

    if (oriVecCoord.x ** 2 + oriVecCoord.y ** 2 > MIN_VALUE):
        destVecCoord.x = cosA * oriVecCoord.x + sinA * (cosB * oriVecCoord.z * cosP - sinB * sinP)
        destVecCoord.y = cosA * oriVecCoord.y + sinA * (cosB * oriVecCoord.z * sinP + sinB * cosP)
        destVecCoord.z = cosA * oriVecCoord.z - sinA * cosB * sinT

    else:
        destVecCoord.x = sinA * cosB * copysign(1.0, oriVecCoord.z)
        destVecCoord.x = sinA * sinB * copysign(1.0, oriVecCoord.z)
        destVecCoord.z = cosA * copysign(1.0, oriVecCoord.z)

    # convert to unit vector
    c = destVecCoord.x ** 2 + destVecCoord.y ** 2 + destVecCoord.z ** 2
    c = sqrt(c)
    destVecCoord.x /= c
    destVecCoord.y /= c
    destVecCoord.z /= c

    return ERRCODE.SUCCESS