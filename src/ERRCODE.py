import logging

# Error code
SUCCESS = 0
FAILURE = 1
OUT_OF_RANGE = 2
INPUT_ERROR = 3
CALCULATE_ERROR = 4
CANNOT_FIND = 5
LOW_WEIGHT = 6
OUTSIDE = 7


# Error messages
ERR_MESSAGE = ["Run Success!",
               "Run Failure!",
               "Input number out of range!",
               "Input Error!",
               "Calculate Error!",
               "Cannot Find",
               "Low weight.",
               "Outside."]


def printMessage(itr):
    if (itr >= len(ERR_MESSAGE)):
        logging.ERROR("Return unknown error message number! Please check!")
    else:
        logging.info(ERR_MESSAGE[itr])
        logging.info("Error Code: " + str(itr))
