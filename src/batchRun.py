import shutil
import os
import subprocess
import time
import re
import csv

INPUT_DIR = "/Users/Gauss/Dropbox/Share_with_MAC/PhD_Studies/Model_work/Experiment_1_Alice_spring_seasonal" \
      "/materials/input_files/"
FLIES_DIR = "/Users/Gauss/Dropbox/Share_with_MAC/FLiES/FLiES_v2.46/"
OUTPUT_DIR = "/Users/Gauss/Dropbox/Share_with_MAC/PhD_Studies/Model_work/Experiment_1_Alice_spring_seasonal" \
      "/materials/"


def moveFile(fileName, inputFileName):

    print("Move file: " + fileName + ".txt")
    newName = fileName + "_" + inputFileName
    newName = newName.replace("init_", "")
    shutil.move(FLIES_DIR + fileName + ".txt", FLIES_DIR + "/results/" + newName)

def batchRun():

    if not os.path.exists(FLIES_DIR + "/results"):
        os.mkdir(FLIES_DIR + "/results")

    for file in os.listdir(INPUT_DIR):
        print("Processing inputfile: ", file)
        command = FLIES_DIR + "flies" + " < " + INPUT_DIR + file
        print(command)

        process = subprocess.Popen(command, cwd=FLIES_DIR, shell=True)
        process.wait()
        print("\n###Finish: " + file)
        time.sleep(2)

        moveFile("flxsum", file)
        moveFile("brfsum", file)

    return True


def getBandName(fileName):
    BAND_NAME = ["RED", "BLUE", "GREEN", "NIR"]
    for name in BAND_NAME:
        if name in fileName:
            return name

def getDate(fileName):
    # print(fileName)
    date = re.search('\d{4}_\d{2}', fileName)
    return date.group(0)

def writeCSV(bandName, bandsReflectance):
    bandReflectance = list(bandsReflectance.get(bandName))
    bandDate = list(bandsReflectance.get(bandName + "_DATE"))

    with open(OUTPUT_DIR + bandName + ".csv", "w") as file:
        wr = csv.writer(file, dialect="excel")
        wr.writerow(["DATE", bandName])
        for i in range(len(bandDate)):
            wr.writerow([str(bandDate[i]), str(bandReflectance[i])])

        file.close()



def extractData():
    bandsReflectance = {"RED": [],
                        "RED_DATE": [],
                        "GREEN": [],
                        "GREEN_DATE": [],
                        "BLUE": [],
                        "BLUE_DATE": [],
                        "NIR": [],
                        "NIR_DATE": [],
                        }

    for fileName in os.listdir(FLIES_DIR + "results/"):
        print(FLIES_DIR + "results")
        print(fileName)
        if (fileName[0] == '.'):
            continue
        date = getDate(fileName)
        band = getBandName(fileName)

        #print(date)
        #print(band)
        with open(FLIES_DIR + "results/" + fileName) as file:
            flag = False
            string = ""

            for line in file:
                if line.startswith(" # idrc theta phi radiance_F stdev_F radiance_Q stdev_Q"):
                    flag = True
                    string = next(file)
                    break

            if flag:
                reflectance = float(string.split()[3])
                #print(reflectance)
                bandsReflectance[band].append(reflectance)
                bandsReflectance[band + "_DATE"].append(date)

        file.close()
    print(bandsReflectance)

    # write to csv
    writeCSV("RED", bandsReflectance)
    writeCSV("GREEN", bandsReflectance)
    writeCSV("BLUE", bandsReflectance)
    writeCSV("NIR", bandsReflectance)


#batchRun()
extractData()
# print(getDate("init_RED_2012_04.txt"))