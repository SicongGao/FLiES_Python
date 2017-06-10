import pandas as pd
import numpy as np
import csv

X_MAX = Y_MAX = 30
outputArray = []
treeObjFilePath = "C:/Users/12442063/Documents/software/fusion/Example_UTM/TEST2/output/"
treeObjFileName = "Example_Out.csv"
outputPath = "../data/"
treeObjFilePath = "../result/"

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
#
# import pyglet
# from pyglet.gl import *
# from pyglet.window import mouse

# treeData = makeTreeTable()
# drawTrees = DrawTrees(treeData)
# drawTrees.startDrawing()

# author: Somsubhra Bairi (201101056)

# Draw sphere with QUAD_STRIP
# Controls: UP/DOWN - scale up/down
#           LEFT/RIGHT - rotate left/right
#           F1 - Toggle surface as SMOOTH or FLAT

# Python imports
from math import *

# OpenGL imports for python
try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *
except:
    print("OpenGL wrapper for python not found")

# Last time when sphere was re-displayed
last_time = 0

# The sphere class
class DrawTrees:

    # Constructor for the sphere class
    def __init__(self, treeData):
        #
        # # Radius of sphere
        # self.radius = radius
        #
        # # Number of latitudes in sphere
        # self.lats = 100
        #
        # # Number of longitudes in sphere
        # self.longs = 100

        self.user_theta = 0
        self.user_height = 0

        # Direction of light
        self.direction = [0.0, 2.0, -1.0, 1.0]

        # Intensity of light
        self.intensity = [0.7, 0.7, 0.7, 1.0]

        # Intensity of ambient light
        self.ambient_intensity = [0.3, 0.3, 0.3, 1.0]

        # The surface type(Flat or Smooth)
        self.surface = GL_FLAT

        self.WIDTH = 800
        self.HEIGHT = 600
        self.RED = (255, 0, 0)
        self.GREEN = (0, 255, 0)
        self.BLUE = (0, 0, 255/255.0)
        self.DARKBLUE = (0, 0, 128)
        self.WHITE = (255, 255, 255)
        self.BLACK = (0, 0, 0)
        self.PINK = (255 / 255.0, 200 / 255.0, 200 / 255.0)
        self.DARKGREEN = (156 / 255.0, 184 / 255.0, 28 / 255.0)
        self.BROWN = (172 / 255.0, 103 / 255.0, 13 / 255.0)

        self.treeData = treeData
        self.processTreeData()
    # Initialize
    def init(self):

        # Set background color to black
        glClearColor(0.0, 0.0, 0.0, 0.0)

        self.compute_location()

        # Set OpenGL parameters
        glEnable(GL_DEPTH_TEST)

        # Enable lighting
        glEnable(GL_LIGHTING)

        # Set light model
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, self.ambient_intensity)

        # Enable light number 0
        glEnable(GL_LIGHT0)

        # Set position and intensity of light
        glLightfv(GL_LIGHT0, GL_POSITION, self.direction)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, self.intensity)

        # Setup the material
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)

    def processTreeData(self):
        i = 0
        scale = 0.05
        while(i < len(self.treeData)):
            if (self.treeData[i][0] == 4):
                self.treeData.pop(i)
            else:
                for j in range(7):
                    self.treeData[i][j] = float(self.treeData[i][j]) * scale
                i += 1

    # Compute location
    def compute_location(self):
        x = 2 * cos(self.user_theta)
        y = 2 * sin(self.user_theta)
        z = self.user_height
        d = sqrt(x * x + y * y + z * z)

        # Set matrix mode
        glMatrixMode(GL_PROJECTION)

        # Reset matrix
        glLoadIdentity()
        glFrustum(-d * 0.5, d * 0.5, -d * 0.5, d * 0.5, d - 1.1, d + 1.1)

        # Set camera
        gluLookAt(x, y, z, 0, 0, 0, 0, 0, 1)

    # Display the sphere
    def display(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Set color to white
        glColor3f(1.0, 1.0, 1.0)

        # Set shade model
        glShadeModel(self.surface)

        self.draw()
        glutSwapBuffers()

    def drawCrown(self, height, radius):
        glColor3f(self.DARKGREEN[0], self.DARKGREEN[1], self.DARKGREEN[2])
        cylinder = gluNewQuadric()
        #glutSolidSphere(radius, 20, 20)
        gluSphere(cylinder, radius, 20, 20)

    def drawTrunk(self, height):
        glRotatef(90, 0, 1, 0)
        glColor3f(self.BROWN[0], self.BROWN[1], self.BROWN[2])
        cylinder = gluNewQuadric()
        gluQuadricDrawStyle(cylinder, GLU_FILL)
        gluQuadricNormals(cylinder, GL_SMOOTH)
        gluQuadricOrientation(cylinder, GLU_INSIDE)
        #glEnable(GL_LIGHTING)
        gluCylinder(cylinder, 0.03, 0.03, height, 40, 20)
        gluDeleteQuadric(cylinder)

    def drawSingleTree(self, posX, posY, heightTrunk, heightCrown, radius):
        glPushMatrix()

        glTranslatef(0, posX, posY)
        self.drawTrunk(heightTrunk)
        glTranslatef(0, 0, heightTrunk + radius)
        self.drawCrown(heightCrown, radius)

        glPopMatrix()

    def drawGround(self):
        area = 100
        glPushMatrix()
        glBegin(GL_QUADS)
        glColor3f(self.PINK[0], self.PINK[1], self.PINK[2])
        glVertex3f(0, 0, 0)
        glVertex3f(0, area, 0)
        glVertex3f(0, area, area)
        glVertex3f(0, 0, area)
        glEnd()
        glPopMatrix()

    # Draw the sphere
    def draw(self):
        self.drawCoordinate()
        self.drawGround()

        for tree in self.treeData:
            self.drawSingleTree(tree[1], tree[2], tree[4], tree[3], tree[5])

    def drawCoordinate(self):
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

    # Keyboard controller for sphere
    def special(self, key, x, y):

        # Scale the sphere up or down
        if key == GLUT_KEY_UP:
            self.user_height += 0.3
        if key == GLUT_KEY_DOWN:
            self.user_height -= 0.3

        # Rotate the cube
        if key == GLUT_KEY_LEFT:
            self.user_theta += 0.3
        if key == GLUT_KEY_RIGHT:
            self.user_theta -= 0.3

        # Toggle the surface
        if key == GLUT_KEY_F1:
            if self.surface == GL_FLAT:
                self.surface = GL_SMOOTH
            else:
                self.surface = GL_FLAT

        self.compute_location()
        glutPostRedisplay()

    # The idle callback
    def idle(self):
        global last_time
        time = glutGet(GLUT_ELAPSED_TIME)

        if last_time == 0 or time >= last_time + 40:
            last_time = time
            glutPostRedisplay()

    # The visibility callback
    def visible(self, vis):
        if vis == GLUT_VISIBLE:
            glutIdleFunc(self.idle)
        else:
            glutIdleFunc(None)

    def openGL(self):
        # Initialize the OpenGL pipeline
        glutInit(sys.argv)

        # Set OpenGL display mode
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)

        # Set the Window size and position
        glutInitWindowSize(800, 800)
        glutInitWindowPosition(50, 100)

        # Create the window with given title
        glutCreateWindow('Research Area')

        self.init()

        # Set the callback function for display
        glutDisplayFunc(self.display)

        # Set the callback function for the visibility
        glutVisibilityFunc(self.visible)

        # Set the callback for special function
        glutSpecialFunc(self.special)

        # Run the OpenGL main loop
        glutMainLoop()

treeData = makeTreeTable()
drawTree = DrawTrees(treeData)
drawTree.openGL()
