from math import *

class Position:
    x = 0.0
    y = 0.0
    z = 0.0

    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

    def setPosition(self, x1, y1, z1):
        self.x = x1
        self.y = y1
        self.z = z1

    def clearPotion(self):
        self.__init__()

    def movePositon(self, distance, vectCoord, x_max, y_max):
        mgn = 1.0e-2
        self.x += (distance + mgn) * vectCoord.x
        self.y += (distance + mgn) * vectCoord.y
        self.z += (distance + mgn) * vectCoord.z

        self.x -= (trunc(self.x / x_max) - 0.5 + copysign(0.5, self.x)) * x_max
        self.y -= (trunc(self.y / y_max) - 0.5 + copysign(0.5, self.y)) * y_max

    def moveVector(self, distance, vectCoord):
        self.x += distance * vectCoord.x
        self.y += distance * vectCoord.y
        self.z += distance * vectCoord.z