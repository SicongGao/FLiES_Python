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