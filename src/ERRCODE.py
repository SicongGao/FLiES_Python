# Error code
SUCCESS = 0
FAILURE = 1
OUT_OF_RANGE = 2
INPUT_ERROR = 3
CALCULATE_ERROR = 4
CANNOT_FIND = 5
LOW_WEIGHT = 6

# Error messages
ERR_MESSAGE = ["Run Success!",
               "Run Failure!",
               "Input number out of range!",
               "Input Error!",
               "Calculate Error!",
               "Cannot Find",
               "Low weight, return."]


def printMessage(itr):
    if (itr >= len(ERR_MESSAGE)):
        print("Return unknown error message number! Please check!")
    else:
        print(ERR_MESSAGE[itr])
        print("Error Code: " + str(itr))

import common as comm
import numpy as np

gmrx = np.zeros(2 * 2 * 2, dtype=float).reshape(2, 2, 2)

for i in range(10,1, -1):
    print(i)