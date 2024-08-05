"""
Author: Omer Sella, with help from Adee Ran and Rich Mellitz
Simple script to create json files that will be used to define COM tests
"""
import os
import json
import fnmatch
from pathlib import Path
comDir = os.environ.get('IEEE8032COM')
if comDir is None:
    # If the user did not provide an environment variable, use the parent directory
    comDir = os.path.dirname(os.getcwd()).replace("\\","/")
# If the user defined an environment variable called IEEE8023CHANNELS then the path to channel data should be taken from there
channelsDir = os.environ.get('IEEE8032CHANNELS')
if channelsDir is None:
    # If the user did not provide an environment variable, hardcode it here
    channelsDir = comDir + "/channels/"

xlsSheets100 = [
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_KR_08_17_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_CR_CA_08_17_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_162_ERL_HOST_10_26_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_120G_ERL_MODULE_10_26_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_120G_ERL_HOST_10_26_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_120g_C2M_tp1a_08_17_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA_120F_C2C_08_17_2022",
    comDir + "/config_sheets_100G/" + "config_com_ieee8023_93a=3ck_SA _TP0V_08_17_2022"
  ]
xlsSheets200 = [
    comDir + "/Config_spreadsheets_200G_exploratory/" + "config_com_ieee8023_93a=df_200G_PAM4_fr55_C2M_TP1a_11_2022",
    comDir + "/Config_spreadsheets_200G_exploratory/" + "config_com_ieee8023_93a=df_200G_PAM4_RCos_C2C_11_2022",
    comDir + "/Config_spreadsheets_200G_exploratory/" + "config_com_ieee8023_93a=df_200G_PAM4_RCos_C2M_TP1a_11_2022",
    comDir + "/Config_spreadsheets_200G_exploratory/" + "config_com_ieee8023_93a=df_200G_PAM4_RCos_CAKR_11_2022",
    comDir + "/Config_spreadsheets_200G_exploratory/" + "config_com_ieee8023_93a=df_200G_PAM4_RCos_Txpre_C2M_TP1a_11_2022"
  ]

def getThruFextNext(pathToDir):
    """
    No safety. 
    Thru list could be empty, or multiple files. 
    A file could be through AND FEXT AND NEXT if it's name suggests it.
    The files are assumed to be in this path
    """
    thru = []
    FEXT = []
    NEXT = []
    for file in os.listdir(pathToDir):
        if fnmatch.fnmatch(file, '*thru*'):
            thru.append(file)
        if fnmatch.fnmatch(file, '*through*'):
            thru.append(file)
        if fnmatch.fnmatch(file, '*fext*'):
            FEXT.append(file)
        if fnmatch.fnmatch(file, '*next*'):
            NEXT.append(file)
    return thru, NEXT, FEXT

def makeJsonForDir(pathToDir, pathToLeaveJson, pathToExcelSheet):
    THRU, NEXT, FEXT = getThruFextNext(pathToDir)
    data = {}
    data["excelSheet"] = pathToExcelSheet
    data["THRU"] = THRU
    data["NEXT"] = NEXT
    data["FEXT"] = FEXT
    jsonFileWithPath = pathToLeaveJson + '/systemConfiguration.json'
    print(jsonFileWithPath)
    with open(jsonFileWithPath, 'w', encoding='utf-8') as fid:
        json.dump(data, fid, ensure_ascii=False, indent=4)
    return

def crawl(pathToBegin, testFilesPath, xlsSheet = xlsSheets200[0]):
    # For personal use only:
    print("Inside crawler, path to begin is :" + pathToBegin)
    createJsonTree = True 
    createResultsTree = False 
    depth = 5
    endPaths = []
    jsonFilePaths = []
    resultPaths = []
    for root, dirs, files in os.walk(pathToBegin):
        newRoot = root.replace("\\","/")
        print(newRoot)
        if not dirs:
            suffixPath = newRoot.split("com")[1]
            newPath = testFilesPath + "/" + suffixPath #"c:/users/omer/com/tests/jsonFiles/channels/" + "/".join(path[depth:])
            resultPath = "c:/users/omer/com/tests/results/channels/" + "/" + suffixPath # .join(path[depth:])        
            if createJsonTree:
                Path(newPath).mkdir(parents=True, exist_ok = True)
            if createResultsTree:
                Path(resultPath).mkdir(parents=True, exist_ok = True)
            jsonFilePaths.append(newPath)
            resultPaths.append(resultPath)
            endPaths.append(newRoot)
            print("---new root")
            print(newRoot)
            print("---new path")
            print(newPath)
            
            makeJsonForDir(newRoot, newPath, xlsSheets200[0])
    return endPaths, jsonFilePaths, resultPaths
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--crawlPath', type = str, help = "Path to begin crawling.")
    parser.add_argument('--xlsSheet', type=str, help = "Excel sheet name.")
    parser.add_argument('--jsonPath', type=str, help = "Path to leave json files.")
    args = parser.parse_args()
    _,_,_ = crawl(args.crawlPath, args.jsonPath, args.xlsSheet)