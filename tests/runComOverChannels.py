import json
import os
import string
import sys
from numpy import empty
import scipy
import matlab.engine
import pandas as pd
import numpy as np
# If the user defined an environment variable called IEEE8023COM, then the path to COM should be taken from there
comDir = os.environ.get('IEEE8032COM')
if comDir is None:
    # If the user did not provide an environment variable, hardcode it here
    comDir = os.getcwd().replace("\\","/") #"c:/users/omer/com/"
# If the user defined an environment variable called IEEE8023CHANNELS then the path to channel data should be taken from there
channelsDir = os.environ.get('IEEE8032CHANNELS')
if channelsDir is None:
    # If the user did not provide an environment variable, hardcode it here
    channelsDir = os.getcwd().replace("\\","/") + "/channels" #c:/users/omer/channels/"

def exprimental():
    workDir = "C:/Users/Omer/channels/weaver_3dj_elec_02_230831/C2M_models/OSFP_X_PCB_3in_25C/ingress/"
    dfs = pd.read_excel("C:/Users/Omer/channels/weaver_3dj_elec_02_230831/C2M_models/OSFP_X_PCB_3in_25C/ingress/configuration.xlsx")
    print(dfs)
    systemDictionary = dfs.to_dict('records')
    if np.count_nonzero(dfs['channel'] == 'THRU') != 1:
        raise ValueError("The configuration supported at the moment is exactly one thru file per excel sheet")
    numberOfFext = np.count_nonzero(dfs['channel'] == 'FEXT')
    numberOfNext = np.count_nonzero(dfs['channel'] == 'FEXT')
    nextFiles = []
    fextFiles = []
    thruFile = empty
    # Get files and sort into fext, next and thru lists
    for entry in systemDictionary:
        if entry['channel'] == 'FEXT':
            fextFiles.append(str(workDir + entry['file']))
        elif entry['channel'] == 'NEXT':
            nextFiles.append(str(workDir + entry['file']))
        elif entry['channel'] == 'THRU':
            thruFile = str(workDir + entry['file'])
        else:
            pass
    filesAsList = [thruFile] + fextFiles + nextFiles
    filesAsString = thruFile
    for file in fextFiles:
        filesAsString = filesAsString + "\," + file
    for file in nextFiles:
        filesAsString = filesAsString + "','" + file
    #filesAsString = filesAsString + "'"
    print(filesAsString)
    print(numberOfFext)
    print(numberOfNext)
    #thruFile = "C:/Users/Omer/channels/weaver_3dj_elec_02_230831/C2M_models/OSFP_X_PCB_3in_25C/ingress/" + 
    resultsDir = "c:/users/omer/com/results/testScript/"
    # Check if directory for saving results exists, and if not create it. We're doing this now, because there is no point in running COM if we can't save the results.
    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)
    # start a MATLAB engine
    engine1 = matlab.engine.start_matlab()
    # Add COM directory to MATLAB path
    engine1.addpath(comDir)
    # Call COM with some parameters
    #'C:\Users\Omer\Downloads\mellitz_3df_02_2211\config_sheets_100G\config_com_ieee8023_93a=3ck_SA _TP0V_08_17_2022.xlsx',0,0,['C:\Users\Omer\channels\akinwale_3df_01_2209\85ohms\C2M_PCB_85ohms_10dB_202208016_v2\C2M_PCB_85ohms_10dB_202208016_v2_thru1.s4p']
    #'C:/Users/Omer/Downloads/mellitz_3df_02_2211/config_sheets_100G/config_com_ieee8023_93a=3ck_SA _TP0V_08_17_2022.xlsx',0,0
    #thruFilePath = 'C:/Users/Omer/channels/akinwale_3df_01_2209/85ohms/C2M_PCB_85ohms_10dB_202208016_v2/C2M_PCB_85ohms_10dB_202208016_v2_thru1.s4p'
    from pathlib import Path
    filesAsPath = Path(filesAsString)
    com = engine1
    #results = engine1.com_ieee8023_93a("C:\Users\Omer\COM\config_sheets_100G\config_com_ieee8023_93a=3ck_SA _TP0V_08_17_2022.xlsx",2,1, "C:\Users\Omer\coms\tests\THRU.s4p","C:/Users/Omer/com/tests/FEXT1.s4p","C:/Users/Omer/com/tests/FEXT2.s4p","C:/Users/Omer/com/tests/NEXT1.s4p")
    #This doesn't work, possibly because of \ as opposed to / results = engine1.com_ieee8023_93a_390('C:\Users\Omer\COM\Config_spreadsheets_200G_exploratory\config_com_ieee8023_93a=df_200G_PAM4_fr55_C2M_TP1a_11_2022.xlsx', 2,1,'C:\Users\Omer\COM\tests\KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_THRU.s4p', 'C:\Users\Omer\COM\tests\KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT1.s4p', 'C:\Users\Omer\COM\tests\KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT2.s4p','C:\Users\Omer\COM\tests\KR_1mCabledBP_TP0TP5_31p4dB_PCBHost_9p8dB_NEXT1.s4p' )
    # This one works - results = engine1.com_ieee8023_93a_390('C:/Users/Omer/COM/Config_spreadsheets_200G_exploratory/config_com_ieee8023_93a=df_200G_PAM4_fr55_C2M_TP1a_11_2022.xlsx', 2,1,'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_THRU.s4p', 'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT1.s4p', 'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT2.s4p','C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_31p4dB_PCBHost_9p8dB_NEXT1.s4p' )
    stringTest = "engine1.com_ieee8023_93a_390('C:/Users/Omer/COM/Config_spreadsheets_200G_exploratory/config_com_ieee8023_93a=df_200G_PAM4_fr55_C2M_TP1a_11_2022.xlsx', 2,1,'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_THRU.s4p', 'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT1.s4p', 'C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_19p3dB_PCBHost_3p8dB_FEXT2.s4p','C:/Users/Omer/COM/tests/KR_1mCabledBP_TP0TP5_31p4dB_PCBHost_9p8dB_NEXT1.s4p' )"
    eval(stringTest)
    #results = engine1.com_ieee8023_93a_390('C:/Users/Omer/COM/Config_spreadsheets_200G_exploratory/config_com_ieee8023_93a=df_200G_PAM4_fr55_C2M_TP1a_11_2022.xlsx', 2,1,stringTest )
    #scipy.io.savemat(fullResultsName, results)

    # Don't forget to turn off the engine
    engine1.quit()
    return
def runTests(testPath):
    jsonFiles = []
    for root, dirs, files in os.walk(testPath):
        if not dirs:
            for file in files:
                jsonFiles.append(root + "/" + file)
    engine1 = matlab.engine.start_matlab()
    engine1.addpath("c:/users/omer/com/")
    for jsonTestFile in jsonFiles:
        with open(jsonTestFile) as fid:
            data = json.load(fid)
            nextString = "','".join(data['NEXT'])
            fextString = "','".join(data['FEXT'])
            thruString = data['THRU'][0]
            excelSheetString = data['excelSheet']
            testString = "engine1.com_ieee8023_93a_390("+"'" + excelSheetString + "'" + "," + str(len(data['FEXT'])) + "," + str(len(data['NEXT'])) + "," + "'" + thruString +"' , '" + nextString + "' , '"+ fextString + "')"
            print(testString)
            eval(testString)
    engine1.quit()
    return
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type = str, help = "Path to json test file.")
    
    args = parser.parse_args()
    testPath = args.path.replace("\\","/")
    runTests(testPath)
    