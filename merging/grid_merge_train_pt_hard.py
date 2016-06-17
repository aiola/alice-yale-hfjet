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

def PtHardBinMerging(LocalPath, Datasets, TrainName, TrainNumbers, MinPtHardBin, MaxPtHardBin, Year, AliPhysicsVersion):
    FileList = []
    FinalFileList = []

    for Dataset,TrainNumber in zip(sorted(Datasets.iterkeys()),TrainNumbers):
        
        AlienPath="/alice/sim/"+str(Year)+"/"+Dataset
        FirtRun = Datasets[Dataset][0]
        SearchPath = "{0}/{1}".format(AlienPath, FirtRun)
        FullTrainNumber = GetFullTrainNumber(SearchPath, TrainName, TrainNumber)
        AlienPath = "alien://" + AlienPath
            
        dest = "/alice/cern.ch/user/s/saiola/{0}/{1}".format(TrainName, FullTrainNumber)
        
        subprocessCall(["alien_mkdir", "-p", dest])
        subprocessCall(["alien_cp", "./grid_merge_train_pt_hard.sh", "alien://{0}".format(dest)])
        subprocessCall(["alien_cp", "./grid_merge_train_pt_hard.jdl", "alien://{0}".format(dest)])
        
        jdlFile = "alien:///{0}/grid_merge_train_pt_hard.jdl".format(dest)
        executableFile = "/{0}/grid_merge_train_pt_hard.sh".format(dest)
        
        for PtHardBin in range(MinPtHardBin, MaxPtHardBin+1):
        
            dest = "/alice/cern.ch/user/s/saiola/{0}/{1}/{2}".format(TrainName, FullTrainNumber, PtHardBin)
            
            localDest = "{0}/{1}/{2}".format(LocalPath, Dataset, PtHardBin)
            if not os.path.isdir(localDest):
                print "Creating directory "+localDest
                os.makedirs(localDest)
            
            fileList = "{0}/fileList.txt".format(localDest)
            
            f = open(fileList, 'w')
            
            for Run in Datasets[Dataset]:
                AlienFile="{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root ".format(AlienPath, Run, PtHardBin, TrainName, FullTrainNumber)
                print "Adding alien location '"+AlienFile
                f.write(AlienFile)
                
            f.close()
            
            alienFileList = "{0}/fileList.txt".format(dest)
            alienOutputFile = "{0}/AnalysisResults.root".format(dest)

            subprocessCall(["alien_mkdir", dest])
            subprocessCall(["alien_cp", fileList, "alien://{0}".format(dest)])
            subprocessCall(["alien_submit", jdlFile, alienFileList, alienOutputFile, executableFile, AliPhysicsVersion])

    print "Done."

    subprocessCall(["ls", LocalPath])

def StartMerging(TrainNumbers, config, AliPhysicsVersion):
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
        PtHardBinMerging(LocalPath, Datasets, config["train"], TrainNumbers, config["min_pt_hard_bin"], config["max_pt_hard_bin"], config["year"], AliPhysicsVersion)
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
    args = parser.parse_args()

    TrainNumberList = args.trainNumber.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range)-1])+1))
    
    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    StartMerging(TrainNumbers, config, args.aliphysics)
