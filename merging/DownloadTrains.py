#!/usr/bin/env python

import os
import subprocess
import sys
import shutil
import argparse
import yaml

import MergeFiles
import ScaleResults

def GetFullTrainNumber(SearchPath, TrainNumber):
    output = subprocess.check_output(["alien_find", "-d", "-l", "1", SearchPath, "{0}_20".format(TrainNumber)], universal_newlines=True)
    #print(output)
    i = output.rfind("{0}_".format(TrainNumber))
    j = len(str(TrainNumber)) + 14 + i
    FullTrainNumber = output[i:j]
    return FullTrainNumber

def GetMergeLists(SearchPath):
    output = subprocess.check_output(["alien_find", "-d", SearchPath, "merge"], universal_newlines=True)
    print output
    mergeLists = output.splitlines()
    resList = []
    for mergeList in mergeLists:
        start = mergeList.rfind("merge")
        stop = mergeList.rfind("/")
        if stop <= start:
            continue
        resList.append(mergeList[start:stop])
        
    return set(resList)

def StartDownload(LocalPath, Datasets, TrainNumbers, TrainName, Overwrite):
    FileName = "AnalysisResults.root"

    for Dataset,TrainNumber in zip(Datasets,TrainNumbers):
        SearchPath="/alice/cern.ch/user/a/alitrain/PWGJE/"+TrainName
        FullTrainNumber = GetFullTrainNumber(SearchPath, TrainNumber)
        mergeLists = GetMergeLists("{0}/{1}".format(SearchPath,FullTrainNumber))
        for mergeList in mergeLists:
            AlienPath = "alien://{0}/{1}/{2}/{3}".format(SearchPath, FullTrainNumber, mergeList, FileName)
            DestPath = "{0}/{1}/{2}".format(LocalPath, Dataset, mergeList)
            os.makedirs(DestPath)
            DestPath += "/{0}".format(FileName)
            print "Copying from alien location '{0}' to local location '{1}'".format(AlienPath, DestPath)
            subprocess.call(["alien_cp", AlienPath, DestPath])

def main(TrainNumbers, config, Overwrite=0):
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

    Datasets=config["datasets"]

    if (len(Datasets)!=len(TrainNumbers)):
        print "The number of datasets {0} must be the same as the number of trains {1}.".format(len(Datasets),len(TrainNumbers))
        print "Datasets are"
        print Datasets
        print "Trains are"
        print TrainNumbers
        exit()

    LocalPath="{0}/{1}".format(JetResults, config["train"])

    for TrainNumber in TrainNumbers:
        LocalPath = "{0}_{1}".format(LocalPath, TrainNumber)

    print "Train: "+config["train"]
    print "Local path: "+LocalPath
    print "Overwrite mode: {0}".format(Overwrite)
    print "Train numbers are: "
    print TrainNumbers
    
    if not os.path.isdir(LocalPath):
        print "Creating directory "+LocalPath
        os.makedirs(LocalPath)
 
    StartDownload(LocalPath, Datasets, TrainNumbers, config["train"], Overwrite)

if __name__ == '__main__':
    # FinalMergeLocal.py executed as script
    
    parser = argparse.ArgumentParser(description='Local final merging for MC production in pT hard bins.')
    parser.add_argument('trainNumber', metavar='trainNumber',
                        help='Train numbers to be downloaded and merged. Use ":" to define a range, and "," for a list')
    parser.add_argument('--yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--overwrite', metavar='overwrite',
                        default=0,
                        help='Overwrite level [0-4]. 0 = no overwrite')
    args = parser.parse_args()

    TrainNumberList = args.trainNumber.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range)-1])+1))
        
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(TrainNumbers, config, args.overwrite)
