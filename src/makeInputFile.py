from config import config
import csv

input_parameters = {

    "number_photon": 1000000,
    "atm_type": 1,
    "diffuse": 0.1,
    "solar_angle": 45.0,
    "solar_elevation": 0.0,
    #"zenith_angle": [100.0, 120.0, 140.0, 150.0, 160.0],
    "zenith_angle": [140.0],
    "azimuth_angle": [180],
    "integration_mode": 3,  # short wavelength
    "single_wavelength": 0.82,
    "atmosphere_type": 2,  # Mid latitude summer
    "aerosol_type": 2,  # Continental average
    "AOT": 0.3,
    "cloud_type": 0,  # cloud free
    "surface_type": 2,  # 3-D canopy
    "calculation_mode": 2,  # NIR
    "BRF_zenith_angles": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
    "BRF_azimuth_angles": [0.0, 180],
    "tree_species": 2,
    "optical_parameters": [0.1, 0.1,
                           0.05, 0.05,
                           0.15, 0.1,  # floor reflectance, transmittance
                           0.1, 0.0001,  # trunk reflectance
                           0.05],  # soil reflectance
    "optical_parameters_NIR": [0.1, 0.1,
                               0.05, 0.05,
                               0.15, 0.1,  # floor reflectance, transmittance
                               0.1, 0.0001,  # trunk reflectance
                               0.05],  # soil reflectance
    "leaf_area_density": [0.5, 0.5],
    "forest_floor_LAI": 0.12666,
    "branch_area_density": [0.5, 0.1],
    "sbar": [0.25, 0.25],
}

LAI_DATA = {
    "DATE": [],
    "LAI_total": [],
    "LAI_canopy": [],
    "LAI_understory": [],
}

BANDS_SETTING = {
    "RED": 0.65,
    "GREEN": 0.55,
    "BLUE": 0.47,
    "NIR": 0.80,
}

SUMMER = [1, 2, 3, 9, 10, 11, 12]
# SUMMER = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

filePath = "C:/Users/12442063/Dropbox/Share with MAC/PhD Studies/Model work/" \
           "Experiment_1_Alice_spring_seasonal/materials/"
fileName_totalLAI = "total_LAI.csv"
fileName_canopyLAI = "canopy_LAI.csv"
fileName_understoryLAI = "understory_LAI.csv"
folder_Input_Files = "input_files/"
filePath = "/Users/Gauss/Dropbox/Share_with_MAC/PhD_Studies/Model_work/Experiment_1_Alice_spring_seasonal" \
      "/materials/"

def makeInput_Python():
    return True


def getWholeString(allData, countFlag):
    string = ""

    if countFlag:
        string = str(len(allData)) + " "

    for data in allData:
        string += str(data) + " "

    return string


def fill_LAI_DATA(fileName, LAIName):

    with open(filePath + fileName) as file:
        fileContent = csv.reader(file)
        next(fileContent)
        LAI_DATA["DATE"] = []

        for rowContent in fileContent:
            date = str(rowContent[0])
            date = date.replace("-", "_")

            LAI = float(rowContent[1])

            LAI_DATA.get("DATE").append(date)
            LAI_DATA.get(LAIName).append(LAI)

    file.close()

def writeNormalFile(file):
    file.write(str(input_parameters.get("number_photon")) + " /* number_photon */\n")
    file.write(str(input_parameters.get("atm_type")) + " /* atm_type */\n")
    file.write(str(input_parameters.get("solar_angle")) + " ")
    file.write(str(input_parameters.get("solar_elevation")) + " /* solar_angle, solar_elevation */\n")

    tempData = getWholeString(input_parameters.get("zenith_angle"), True)
    file.write(tempData + " /* zenith_angle */\n")

    tempData = getWholeString(input_parameters.get("azimuth_angle"), True)
    file.write(tempData + " /* azimuth_angle */\n")

    file.write(str(input_parameters.get("integration_mode")) + " /* integration_mode */\n")
    file.write(str(input_parameters.get("single_wavelength")) + " /* single_wavelength */\n")
    file.write(str(input_parameters.get("atmosphere_type")) + " /* atmosphere_type */\n")
    file.write(str(input_parameters.get("aerosol_type")) + " /* aerosol_type */\n")
    file.write(str(input_parameters.get("AOT")) + " /* AOT */\n")
    file.write(str(input_parameters.get("cloud_type")) + " /* cloud_type */\n")
    file.write(str(input_parameters.get("surface_type")) + " /* surface_type */\n")
    file.write(str(input_parameters.get("calculation_mode")) + " /* calculation_mode */\n")

    tempData = getWholeString(input_parameters.get("BRF_zenith_angles"), True)
    file.write(tempData + " /* BRF_zenith_angles */\n")

    tempData = getWholeString(input_parameters.get("BRF_azimuth_angles"), True)
    file.write(tempData + " /* BRF_azimuth_angles */\n")

    file.write(str(input_parameters.get("tree_species")) + " /* tree_species */\n")

    tempData = getWholeString(input_parameters.get("optical_parameters"), False)
    file.write(tempData + " /* optical_parameters */\n")

    tempData = getWholeString(input_parameters.get("leaf_area_density"), False)
    file.write(tempData + " /* leaf_area_density */\n")

    file.write(str(input_parameters.get("forest_floor_LAI")) + " /* forest_floor_LAI */\n")

    tempData = getWholeString(input_parameters.get("branch_area_density"), False)
    file.write(tempData + " /* branch_area_density */\n")

    tempData = getWholeString(input_parameters.get("sbar"), False)
    file.write(tempData + " /* sbar */\n")

def writeNadirFile(file):
    file.write(str(input_parameters.get("number_photon")) + " /* number_photon */\n")
    file.write(str(input_parameters.get("atm_type")) + " /* atm_type */\n")
    file.write(str(input_parameters.get("solar_angle")) + " ")
    file.write(str(input_parameters.get("solar_elevation")) + " /* solar_angle, solar_elevation */\n")

    tempData = getWholeString(input_parameters.get("zenith_angle"), True)
    file.write(tempData + " /* zenith_angle */\n")

    tempData = getWholeString(input_parameters.get("azimuth_angle"), True)
    file.write(tempData + " /* azimuth_angle */\n")

    file.write(str(input_parameters.get("integration_mode")) + " /* integration_mode */\n")
    file.write(str(input_parameters.get("single_wavelength")) + " /* single_wavelength */\n")
    file.write(str(input_parameters.get("atmosphere_type")) + " /* atmosphere_type */\n")
    file.write(str(input_parameters.get("aerosol_type")) + " /* aerosol_type */\n")
    file.write(str(input_parameters.get("AOT")) + " /* AOT */\n")
    file.write(str(input_parameters.get("cloud_type")) + " /* cloud_type */\n")
    file.write(str(input_parameters.get("surface_type")) + " /* surface_type */\n")
    file.write(str(input_parameters.get("calculation_mode")) + " /* calculation_mode */\n")

    # tempData = getWholeString(input_parameters.get("BRF_zenith_angles"), True)
    # file.write(tempData + " /* BRF_zenith_angles */\n")
    #
    # tempData = getWholeString(input_parameters.get("BRF_azimuth_angles"), True)
    # file.write(tempData + " /* BRF_azimuth_angles */\n")

    file.write(str(input_parameters.get("tree_species")) + " /* tree_species */\n")

    tempData = getWholeString(input_parameters.get("optical_parameters"), False)
    file.write(tempData + " /* optical_parameters */\n")

    tempData = getWholeString(input_parameters.get("leaf_area_density"), False)
    file.write(tempData + " /* leaf_area_density */\n")

    file.write(str(input_parameters.get("forest_floor_LAI")) + " /* forest_floor_LAI */\n")

    tempData = getWholeString(input_parameters.get("branch_area_density"), False)
    file.write(tempData + " /* branch_area_density */\n")

    tempData = getWholeString(input_parameters.get("sbar"), False)
    file.write(tempData + " /* sbar */\n")

def writeShortWavelengthFile(file):
    file.write(str(input_parameters.get("number_photon")) + " /* number_photon */\n")
    file.write(str(input_parameters.get("atm_type")) + " /* atm_type */\n")
    file.write(str(input_parameters.get("solar_angle")) + " ")
    file.write(str(input_parameters.get("solar_elevation")) + " /* solar_angle, solar_elevation */\n")

    tempData = getWholeString(input_parameters.get("zenith_angle"), True)
    file.write(tempData + " /* zenith_angle */\n")

    tempData = getWholeString(input_parameters.get("azimuth_angle"), True)
    file.write(tempData + " /* azimuth_angle */\n")

    file.write(str(input_parameters.get("integration_mode")) + " /* integration_mode */\n")
    # file.write(str(input_parameters.get("single_wavelength")) + " /* single_wavelength */\n")
    file.write(str(input_parameters.get("atmosphere_type")) + " /* atmosphere_type */\n")
    file.write(str(input_parameters.get("aerosol_type")) + " /* aerosol_type */\n")
    file.write(str(input_parameters.get("AOT")) + " /* AOT */\n")
    file.write(str(input_parameters.get("cloud_type")) + " /* cloud_type */\n")
    file.write(str(input_parameters.get("surface_type")) + " /* surface_type */\n")
    file.write(str(input_parameters.get("calculation_mode")) + " /* calculation_mode */\n")

    # tempData = getWholeString(input_parameters.get("BRF_zenith_angles"), True)
    # file.write(tempData + " /* BRF_zenith_angles */\n")
    #
    # tempData = getWholeString(input_parameters.get("BRF_azimuth_angles"), True)
    # file.write(tempData + " /* BRF_azimuth_angles */\n")

    file.write(str(input_parameters.get("tree_species")) + " /* tree_species */\n")

    tempData = getWholeString(input_parameters.get("optical_parameters"), False)
    file.write(tempData + " /* optical_parameters */\n")

    tempData = getWholeString(input_parameters.get("optical_parameters_NIR"), False)
    file.write(tempData + " /* optical_parameters_NIR */\n")

    tempData = getWholeString(input_parameters.get("leaf_area_density"), False)
    file.write(tempData + " /* leaf_area_density */\n")

    file.write(str(input_parameters.get("forest_floor_LAI")) + " /* forest_floor_LAI */\n")

    tempData = getWholeString(input_parameters.get("branch_area_density"), False)
    file.write(tempData + " /* branch_area_density */\n")

    tempData = getWholeString(input_parameters.get("sbar"), False)
    file.write(tempData + " /* sbar */\n")

def makeInput_Fortran():

    # file LAI data into dataframe
    fill_LAI_DATA(fileName_totalLAI, "LAI_total")
    fill_LAI_DATA(fileName_canopyLAI, "LAI_canopy")
    fill_LAI_DATA(fileName_understoryLAI, "LAI_understory")

    # for band in BANDS_SETTING:
    #     i = 0
    #     for iDate in LAI_DATA.get("DATE"):
    #         inputFileName = "init_" + band + "_" + iDate + ".txt"
    #         print(inputFileName)
    #         input_parameters["single_wavelength"] = BANDS_SETTING.get(band)
    #         input_parameters["leaf_area_density"] = []
    #         input_parameters["leaf_area_density"].append(float(LAI_DATA.get("LAI_canopy")[i]))
    #         input_parameters["leaf_area_density"].append(float(LAI_DATA.get("LAI_understory")[i]))
    #
    #         i += 1
    #
    #         with open(filePath + folder_Input_Files + inputFileName, 'w') as file:
    #             writeShortWavelengthFile(file)
    #
    #         file.close()

    i = 0
    for iDate in LAI_DATA.get("DATE"):
        inputFileName = "init_" + iDate + ".txt"
        print(inputFileName)
        month = int(iDate.split("_")[1])
        print(month)
        if month in SUMMER:
            input_parameters["atmosphere_type"] = 2
        else:
            input_parameters["atmosphere_type"] = 3
        #input_parameters["single_wavelength"] = BANDS_SETTING.get(band)
        input_parameters["leaf_area_density"] = []
        input_parameters["leaf_area_density"].append(float(LAI_DATA.get("LAI_canopy")[i]))
        input_parameters["leaf_area_density"].append(float(LAI_DATA.get("LAI_understory")[i]))

        i += 1

        with open(filePath + folder_Input_Files + inputFileName, 'w') as file:
            writeShortWavelengthFile(file)

        file.close()

makeInput_Fortran()