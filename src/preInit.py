import pandas as pd
import numpy as np
import csv

X_MAX = Y_MAX = 30
outputArray = []
treeObjFilePath = "C:/Users/12442063/Documents/software/fusion/Example_UTM/TEST2/output/"
treeObjFileName = "Example_Out.csv"
outputPath = "../data/"


def makeTreeTable():
    numObj = 0

    with open(treeObjFilePath + treeObjFileName) as file:
        fileContent = csv.reader(file)
        next(fileContent)
        for rowContent in fileContent:
            #rowContent = int(rowContent)
            if ((float(rowContent[1]) <= X_MAX) and (float(rowContent[2]) <= Y_MAX)):
                # [objType, x, y, z, high, radius, species]
                # crown length was assumed to be 50% of tree height, and the maximum projected crown radius
                # A satellite-based method for monitoring seasonality in the overstoryleaf area index of Siberian
                # larch forest.pdf
                # Hideki

                height = format(float(rowContent[3]) / 2.0, "3.5f")
                rowContent[3] = format(float(rowContent[3]), "3.5f")
                rowContent[5] = format(float(rowContent[5]), "3.5f")
                tempArr = [3, rowContent[1], rowContent[2], rowContent[3], height, rowContent[5], 1]
                outputArray.append(tempArr)

                tempArr = [4, rowContent[1], rowContent[2], height, height, 0.1, 1]
                outputArray.append(tempArr)

                numObj += 2

    file = open(outputPath + "crowndata.txt", "w")
    file.write(str(numObj) + "\n")
    for i in range(numObj):
        string = str(outputArray[i][0]) + " "

        for j in range(1, 6):
            string += format(float(outputArray[i][j]), "12.5f")

        string += "   " + str(outputArray[i][6]) + "\n"

        file.write(string)
    file.close()

    return outputArray

import pyglet
from pyglet.gl import *
from pyglet.window import mouse
import math

class DrawTrees:

    def __init__(self, treeData):
        self.WIDTH = 800
        self.HEIGHT = 600
        self.RED = (255, 0, 0)
        self.GREEN = (0, 255, 0)
        self.BLUE = (0, 0, 255)
        self.DARKBLUE = (0, 0, 128)
        self.WHITE = (255, 255, 255)
        self.BLACK = (0, 0, 0)
        self.PINK = (255 / 255.0, 200 / 255.0, 200 / 255.0)
        self.DARKGREEN = (156 / 255.0, 184 / 255.0, 28 / 255.0)
        self.BROWN = (172 / 255.0, 103 / 255.0, 13 / 255.0)

        self.window = pyglet.window.Window()
        self.window.set_caption("Research Area")
        self.window.set_size(self.WIDTH, self.HEIGHT)

        self.on_draw = self.window.event(self.on_draw)
        self.on_mouse_drag = self.window.event(self.on_mouse_drag)
        self.on_mouse_scroll = self.window.event(self.on_mouse_scroll)

        self.sensitivity = 0.1
        self.camX = 0.0
        self.camZ = 5.0
        self.camVecX = 0.0
        self.camVecZ = -1.0
        self.camAngle = 0
        self.deltaAngle = 0.0
        self.deltaMove = 0
        self.scrollScale = 1

        self.treeData = treeData
        self.processTreeData()
        self.mouseFlag = 0


        glClearColor(0.0, 0.0, 0.0, 0.0)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0, 0, self.WIDTH, self.HEIGHT)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, self.WIDTH / self.HEIGHT, 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def processTreeData(self):
        i = 0
        while (i < len(self.treeData)):
            if (self.treeData[i][0] == 4):
                self.treeData.pop(i)
            else:
                for j in range(7):
                    self.treeData[i][j] = float(self.treeData[i][j]) * 10
                i += 1


    def drawCrown(self, height, radius):
        return True

    def drawTrunk(self, height):
        # glRotatef(24, 0.0, 0.0, 0.0)
        # glLineWidth(12.5)
        #
        # glPushMatrix()
        # glColor3f(1.0, 0.0, 0.0)
        # glBegin(GL_LINES)
        # glVertex3f(0.0, 0.0, 0.0)
        # glVertex3f(150, 0, 0)
        # glEnd()
        # glPopMatrix()
        #
        #
        # glPushMatrix()
        # glColor3f(1.0, 10.0, 0.0)
        # glBegin(GL_LINES)
        # glVertex3f(0.0, 0.0, 0.0)
        # glVertex3f(0, 150, 0)
        # glEnd()
        # glPopMatrix()
        #
        # # glPushMatrix()
        # # glColor3f(0.0, 10.0, 19.0)
        # glBegin(GL_LINES)
        # glVertex3f(0.0, 0.0, 0.0)
        # glVertex3f(0, 0, 150)
        # glEnd()

        # glPopMatrix()
        glColor3f(1, 0, 0)
        # cylinder = gluNewQuadric()
        # gluQuadricDrawStyle(cylinder, GLU_FILL)
        # gluQuadricNormals(cylinder, GL_SMOOTH)
        # gluQuadricOrientation(cylinder, GLU_INSIDE)
        # # gluCylinder(cylinder, base radius, top radius, length, slice, slice)
        # gluCylinder(cylinder, 20.0, 20.0, height, 100, 100)
        #
        # gluDeleteQuadric(cylinder)

        #self.drawGround()

    def drawSingleTree(self, posX, posY, heightTrunk, heightCrown, radius):

        glPushMatrix()

        glTranslatef(posX, posY, 0)

        self.drawTrunk(heightTrunk * 10)

        # move to crown
        # glTranslatef(0, 0, heightTrunk)
        # self.drawCrown(heightCrown, radius)

        glPopMatrix()

    def drawGround(self):
        area_width = 300
        area_long = 300
        area_start = 0
        ground_vertices = [(area_start, area_start, 0),
                           (area_start + area_width, area_start, 0),
                           (area_start + area_width, area_start + area_long, 0),
                           (area_start, area_start + area_long, 0)]

        glPushMatrix()
        glBegin(GL_QUADS)

        # for vertex in ground_vertices:
        #     glColor3f(self.DARKGREEN[0], self.DARKGREEN[1], self.DARKGREEN[2])
        #     glVertex3f(vertex[0], vertex[1], vertex[2])
        glVertex3f(-100.0, 0.0, -100.0)
        glVertex3f(-100.0, 0.0, 100.0)
        glVertex3f(100.0, 0.0, 100.0)
        glVertex3f(100.0, 0.0, -100.0)

        glEnd()
        glPopMatrix()

    def drawWholeArea(self):

        self.drawGround()


        glTranslatef(100, 100, 0)
        #glTranslatef(0, -10, -10)
        glLineWidth(10)


        glPushMatrix()
        glColor3f(0.0, 1.0, 0.0) # Y: green
        glBegin(GL_LINES)
        glVertex3f(0, -100, 0)
        glVertex3f(0, 100, 0)
        glEnd()
        glPopMatrix()

        glPushMatrix()
        glColor3f(0.0, 0.0, 1.0) # X: BLUE
        glBegin(GL_LINES)
        glVertex3f(-100, 0, 0)
        glVertex3f(100, 0, 0)
        glEnd()
        glPopMatrix()

        glPushMatrix()
        glColor3f(1.0, 0.0, 0.0)  # z: red
        glBegin(GL_LINES)
        glVertex3f(0, 0, -100)
        glVertex3f(0, 0, 100)
        glEnd()
        glPopMatrix()

    def doRotateAndZoom(self):
        # transX = 100
        # transY = 100
        # transZ = 0
        #
        # glTranslatef(0, 0, -self.scrollScale)
        # glTranslatef(transX, transY, transZ)
        # #glRotatef(self.rotateY, 1.0, 0.0, 0.0)
        # glRotatef(self.rotateX, 0.0, 1.0, 0.0)
        # glTranslatef(-transX, -transY, -transZ)
        return True

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):

        if (buttons & mouse.LEFT):
            # # print(x, y, dx, dy, buttons, modifiers)
            # self.rotateX += dx * self.sensitivity
            # self.rotateY += dy * self.sensitivity

            self.deltaAngle = dx * 0.001

            self.camVecX = math.sin(self.camAngle + self.deltaAngle)
            self.camVecZ = -math.cos(self.camAngle + self.deltaAngle)
            self.mouseFlag = 1


    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        self.scrollScale -= self.scrollScale * scroll_y * self.sensitivity

    def on_draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        gluLookAt(0, 0, 2, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)


        self.drawWholeArea()

        # glRotatef(9, 0.0, 0.0, 1.0)

    def startDrawing(self):
        pyglet.app.run()

treeData = makeTreeTable()
drawTrees = DrawTrees(treeData)
drawTrees.startDrawing()
