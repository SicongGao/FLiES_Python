# Error code
SUCCESS = 0
FAILURE = 1
OUT_OF_RANGE = 2
INPUT_ERROR = 3
CALCULATE_ERROR = 4

# Error messages
ERR_MESSAGE = ["Run Success!",
               "Run Failure!",
               "Input number out of range!",
               "Input Error!",
               "Calculate Error!"]


def printMessage(itr):
    if (itr >= len(ERR_MESSAGE)):
        print("Return unknown error message number! Please check!")
    else:
        print(ERR_MESSAGE[itr])
        print("Error Code is " + str(itr))