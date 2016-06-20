#!/usr/bin/env python

import os
import subprocess
import argparse
import yaml

def subprocessCall(cmd):
    print(cmd)
    subprocess.call(cmd)

def GetFullTrainNumber(SearchPath, TrainName, TrainNumber):
    output = subprocess.check_output(["alien_find", "-d", "-l", "1", SearchPath, TrainName+"/"+str(TrainNumber)+"_*"], universal_newlines=True)
    #print("alien_find", "-d", "-l", "1", SearchPath, TrainName+"/"+str(TrainNumber)+"_*")
    #print(output)
    i = output.rfind(str(TrainNumber)+"_")
    j = len(str(TrainNumber)) + 14 + i
    FullTrainNumber = output[i:j]
    return FullTrainNumber

def PtHardBinMerging(LocalPath, Datasets, TrainName, TrainNumbers, MinPtHardBin, MaxPtHardBin, Year, AliPhysicsVersion, TestMode, GridTestMode):
    FileList = []
    FinalFileList = []

    for Dataset,TrainNumber in zip(sorted(Datasets.iterkeys()),TrainNumbers):
        
        AlienPath="/alice/sim/"+str(Year)+"/"+Dataset
        FirtRun = Datasets[Dataset][0]
        SearchPath = "{0}/{1}".format(AlienPath, FirtRun)
        FullTrainNumber = GetFullTrainNumber(SearchPath, TrainName, TrainNumber)
            
        dest = "/alice/cern.ch/user/s/saiola/{0}/{1}".format(TrainName, FullTrainNumber)
        
        subprocessCall(["alien_mkdir", "-p", dest])
        
        subprocessCall(["alien_rm", "{0}/grid_merge_train_pt_hard.sh".format(dest)])
        subprocessCall(["alien_rm", "{0}/MergeFiles.C".format(dest)])
        
        subprocessCall(["alien_cp", "./grid_merge_train_pt_hard.sh", "alien://{0}".format(dest)])
        subprocessCall(["alien_cp", "./MergeFiles.C", "alien://{0}".format(dest)])
        
        executableFile = "{0}/grid_merge_train_pt_hard.sh".format(dest)
        macroFile = "{0}/MergeFiles.C".format(dest)
        
        for PtHardBin in range(MinPtHardBin, MaxPtHardBin+1):
            
            localDest = "{0}/{1}/{2}".format(LocalPath, Dataset, PtHardBin)
            if not os.path.isdir(localDest):
                print "Creating directory "+localDest
                os.makedirs(localDest)
                
                
            alienFileList = "{0}/{1}/fileList.txt".format(dest, PtHardBin)
            alienOutputFile = "{0}/{1}".format(dest, PtHardBin)
            
            fileExists = True
            try:
                subprocess.check_call(["alien_ls", "{0}/AnalysisResults.root".format(alienOutputFile)], universal_newlines=True)
            except subprocess.CalledProcessError:
                fileExists = False   
            
            if fileExists:
                print("Output file alien://{0}/AnalysisResults.root already exists. Skipping.".format(alienOutputFile))
                continue
            
            fileList = "{0}/fileList.txt".format(localDest)
            
            f = open(fileList, 'w')
            
            jdlInputData = ""
            nruns = 0
            for Run in Datasets[Dataset]:
                AlienFile = "{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root ".format(AlienPath, Run, PtHardBin, TrainName, FullTrainNumber)
                
                fileExists = True
                
                try:
                    subprocess.check_call(["alien_ls", AlienFile], universal_newlines=True)
                except subprocess.CalledProcessError:
                    fileExists = False

                AlienFile = "alien://{0}".format(AlienFile)

                if not fileExists:
                    print("File '{0}' does not exist exist. Skipping.".format(AlienFile))
                    continue
            
                print("Adding alien location '{0}'".format(AlienFile))
                f.write(AlienFile)
            
                if jdlInputData:
                    jdlInputData += ",\n"
            
                jdlInputData +="\"LF:{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root,nodownload\"".format(AlienPath, Run, PtHardBin, TrainName, FullTrainNumber)
                
                nruns += 1
                    
            f.close()
            
            if nruns == 0:
                print("No input files found for pt,hard bin {0} of train {1}".format(PtHardBin, TrainNumber))
                continue

            subprocessCall(["alien_mkdir", "{0}/{1}".format(dest, PtHardBin)])
            
            subprocessCall(["alien_rm", "{0}/{1}/grid_merge_train_pt_hard.jdl".format(dest, PtHardBin)])
            subprocessCall(["alien_rm", "{0}/{1}/fileList.txt".format(dest, PtHardBin)])
            
            subprocessCall(["alien_cp", fileList, "alien://{0}/{1}".format(dest, PtHardBin)])
            
            jdlContent = "# This is the startup script \n\
Executable = \"{executable}\"; \n\
# Time after which the job is killed (500 min.) \n\
TTL = \"3600\"; \n\
OutputDir = \"{dest}/{pt_hard}\"; \n\
Output = {{ \n\
\"AnalysisResults*.root\" \n\
}}; \n\
Arguments = \"\"; \n\
Packages = {{ \n\
\"VO_ALICE@AliPhysics::{aliphysics}\", \n\
\"VO_ALICE@APISCONFIG::V1.1x\", \n\
\"VO_ALICE@jemalloc::v3.6.0\" \n\
}}; \n\
# JDL variables \n\
JDLVariables = \n\
{{ \n\
\"Packages\", \n\
\"OutputDir\" \n\
}}; \n\
# List of input files to be uploaded to workers \n\
InputFile = {{ \n\
\"LF:/{dest}/MergeFiles.C\", \n\
\"LF:/{dest}/{pt_hard}/fileList.txt\" \n\
}}; \n\
InputData = {{{inputData}}}; \n\
".format(executable=executableFile, dest=dest, inputData=jdlInputData, aliphysics=AliPhysicsVersion, pt_hard=PtHardBin)
            
            localJdlFile = "{0}/grid_merge_train_pt_hard.jdl".format(localDest)
            
            f = open(localJdlFile, 'w')
            f.write(jdlContent)
            f.close()
            
            subprocessCall(["alien_cp", localJdlFile, "alien://{0}/{1}".format(dest, PtHardBin)])
            
            jdlFile = "alien://{0}/{1}/grid_merge_train_pt_hard.jdl".format(dest, PtHardBin)
            
            if TestMode:
                subprocessCall(["chmod", "+x", "./grid_merge_train_pt_hard.sh"])
                subprocessCall(["./grid_merge_train_pt_hard.sh", alienFileList])
                return
            else:
                subprocessCall(["alien_submit", jdlFile])
                if GridTestMode:
                    return

    print "Done."

    subprocessCall(["ls", LocalPath])

def StartMerging(TrainNumbers, config, AliPhysicsVersion, TestMode, GridTestMode):
    try:
        rootPath=subprocess.check_output(["which", "root"]).rstrip()
        alirootPath=subprocess.check_output(["which", "aliroot"]).rstrip()
        alienPath=subprocess.check_output(["which", "alien-token-info"]).rstrip()
    except subprocess.CalledProcessError:
        print "Environment is not configured correctly!"
        exit()

    print "Root: "+rootPath
    print "AliRoot: "+alirootPath
    print "Alien: "+alienPath

    if "JETRESULTS" in os.environ:
        JetResults=os.environ["JETRESULTS"]
    else:
        JetResults="."

    try:
        print "Token info disabled"
        #tokenInfo=subprocess.check_output(["alien-token-info"])
    except subprocess.CalledProcessError:
        print "Alien token not available. Creating a token for you..."
        try:
            #tokenInit=subprocess.check_output(["alien-token-init", "saiola"], shell=True)
            print "Token init disabled"
        except subprocess.CalledProcessError:
            print "Error: could not create the token!"
            exit()

    Datasets = config["datasets"]

    if (len(Datasets)!=len(TrainNumbers)):
        print "The number of datasets {0} must be the same as the number of trains {1}.".format(len(Datasets),len(TrainNumbers))
        print "Datasets are"
        print Datasets
        print "Trains are"
        print TrainNumbers
        exit()

    LocalPath="{0}/{1}".format(JetResults, config["train"])

    for TrainNumber in TrainNumbers:
        LocalPath+="_"+str(TrainNumber)

    print "Train: "+config["train"]
    print "Local path: "+LocalPath
    print "Train numbers are: "
    print TrainNumbers
    
    if not os.path.isdir(LocalPath):
        print "Creating directory "+LocalPath
        os.makedirs(LocalPath)
        
    if config["pt_hard_bins"]:
        PtHardBinMerging(LocalPath, Datasets, config["train"], TrainNumbers, config["min_pt_hard_bin"], config["max_pt_hard_bin"], 
                         config["year"], AliPhysicsVersion, TestMode, GridTestMode)
    else:
        print("Error! This is only for pt hard binned productions! Fix YAML file.")

if __name__ == '__main__':
    # FinalMergeLocal.py executed as script
    
    parser = argparse.ArgumentParser(description='Local final merging for LEGO train results.')
    parser.add_argument('trainNumber', metavar='trainNumber',
                        help='Train numbers to be downloaded and merged')
    parser.add_argument('--aliphysics', metavar='vXXX',
                        help='AliPhysics version')
    parser.add_argument('--config', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--test', action='store_const',
                        default=False, const=True,
                        help='Test mode')
    parser.add_argument('--test-grid', action='store_const',
                        default=False, const=True,
                        help='Test grid mode')
    args = parser.parse_args()

    TrainNumberList = args.trainNumber.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range)-1])+1))
    
    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    StartMerging(TrainNumbers, config, args.aliphysics, args.test, args.test_grid)
