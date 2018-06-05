#!/usr/bin/env python

import os
import subprocess
import sys
import shutil
import argparse
import yaml

import MergeFiles
import ScaleResults

def GetFullTrainNumber(SearchPath, TrainName, TrainNumber):
    TrainNumbers = GetFullTrainNumbers(SearchPath, TrainName, TrainNumber)
    if len(TrainNumbers) != 1:
        print("Error: was expecting 1 train number and I found {}".format(len(TrainNumbers)))
        print(TrainNumbers)
        exit(1)
    return TrainNumbers[0]

def GetFullTrainNumbers(SearchPath, TrainName, TrainNumber):
    cmd = ["alien_find", "-d", SearchPath, "{0}/{1}_201*/lego_train.jdl".format(TrainName, TrainNumber)]
    print(cmd)
    output = subprocess.check_output(cmd, universal_newlines=True)
    print(output)
    lines = output.splitlines()
    FullTrainNumbers = []
    for line in lines:
        if not line: continue
        i = line.rfind("{0}_201".format(TrainNumber))
        if i < 0: continue
        j = line.find("/", i + 1)
        FullTrainNumber = line[i:j]
        FullTrainNumbers.append(FullTrainNumber)
    return FullTrainNumbers

def GetMergeLists(SearchPath):
    cmd = ["alien_find", "-d", SearchPath, "merge"]
    print(cmd)
    output = subprocess.check_output(cmd, universal_newlines=True)
    print output
    mergeLists = output.splitlines()
    resList = []
    for mergeList in mergeLists:
        mergeList = mergeList.strip()
        if not mergeList.endswith("/"):
            continue

        start = mergeList.rfind("merge")
        stop = mergeList.rfind("/")
        if stop <= start:
            continue
        resList.append(mergeList[start:stop])

    return set(resList)

def StartDownloadMetadataset(LocalPath, Datasets, Children, TrainNumber, TrainName, Overwrite, FileName, DryRun):
    SearchPath = "/alice/cern.ch/user/a/alitrain/PWGJE/"
    FullTrainNumbers = GetFullTrainNumbers(SearchPath, TrainName, TrainNumber)
    for FullTrainNumber in FullTrainNumbers:
        print("The full train number is '{}'". format(FullTrainNumber))
        ichild = FullTrainNumber.rfind("_", 0, FullTrainNumber.rfind("_") - 1)
        if ichild < 6:
            child = ""
        else:
            child = FullTrainNumber[ichild+1:]
        print("The child name is '{}'".format(child))
        Dataset = Datasets[Children.index(child)]
        print("The datatset name is '{}'".format(Dataset))
        mergeLists = GetMergeLists("{0}/{1}/{2}/".format(SearchPath, TrainName, FullTrainNumber))
        for mergeList in mergeLists:
            AlienPath = "alien://{0}/{1}/{2}/{3}/{4}".format(SearchPath, TrainName, FullTrainNumber, mergeList, FileName)
            DestPath = "{}/{}_{}/{}".format(LocalPath, Dataset, FullTrainNumber, mergeList)
            if not DryRun and not os.path.isdir(DestPath): os.makedirs(DestPath)
            DestPath += "/{0}".format(FileName)
            print "Copying from alien location '{0}' to local location '{1}'".format(AlienPath, DestPath)
            if not DryRun: subprocess.call(["alien_cp", AlienPath, DestPath])

def StartDownload(LocalPath, Datasets, TrainNumbers, TrainName, Overwrite, FileName, DryRun):
    SearchPath = "/alice/cern.ch/user/a/alitrain/PWGJE/"
    for Dataset, TrainNumber in zip(Datasets, TrainNumbers):
        FullTrainNumber = GetFullTrainNumber(SearchPath, TrainName, TrainNumber)
        print("The full train number is '{}'". format(FullTrainNumber))
        print("The datatset name is '{}'".format(Dataset))
        mergeLists = GetMergeLists("{0}/{1}/{2}/".format(SearchPath, TrainName, FullTrainNumber))
        for mergeList in mergeLists:
            AlienPath = "alien://{0}/{1}/{2}/{3}/{4}".format(SearchPath, TrainName, FullTrainNumber, mergeList, FileName)
            DestPath = "{}/{}_{}/{}".format(LocalPath, Dataset, FullTrainNumber, mergeList)
            if not DryRun and not os.path.isdir(DestPath): os.makedirs(DestPath)
            DestPath += "/{0}".format(FileName)
            print "Copying from alien location '{0}' to local location '{1}'".format(AlienPath, DestPath)
            if not DryRun: subprocess.call(["alien_cp", AlienPath, DestPath])

def main(TrainNumbers, Config, FileName, Overwrite, DryRun):
    try:
        rootPath = subprocess.check_output(["which", "root"]).rstrip()
        alirootPath = subprocess.check_output(["which", "aliroot"]).rstrip()
        alienPath = subprocess.check_output(["which", "alien-token-info"]).rstrip()
    except subprocess.CalledProcessError:
        print "Environment is not configured correctly!"
        exit()

    print "Root: " + rootPath
    print "AliRoot: " + alirootPath
    print "Alien: " + alienPath

    if "JETRESULTS" in os.environ:
        JetResults = os.environ["JETRESULTS"]
    else:
        JetResults = "."

    try:
        print "Token info disabled"
        # tokenInfo=subprocess.check_output(["alien-token-info"])
    except subprocess.CalledProcessError:
        print "Alien token not available. Creating a token for you..."
        try:
            # tokenInit=subprocess.check_output(["alien-token-init", "saiola"], shell=True)
            print "Token init disabled"
        except subprocess.CalledProcessError:
            print "Error: could not create the token!"
            exit()

    Datasets = Config["datasets"]

    Metadataset = False
    if "meta_dataset" in Config and Config["meta_dataset"]:
        Metadataset = True

    if Metadataset and len(TrainNumbers) != 1:
        print("For a metadataset expected one and only one train number!")
        exit()

    if not Metadataset and (len(Datasets) != len(TrainNumbers)):
        print "The number of datasets {0} must be the same as the number of trains {1}.".format(len(Datasets), len(TrainNumbers))
        print "Datasets are"
        print Datasets
        print "Trains are"
        print TrainNumbers
        exit()

    LocalPath = "{0}/{1}".format(JetResults, Config["train"])

    for TrainNumber in TrainNumbers:
        LocalPath = "{0}_{1}".format(LocalPath, TrainNumber)

    print "Train: " + Config["train"]
    print "Local path: " + LocalPath
    print "Overwrite mode: {0}".format(Overwrite)
    print "Train numbers are: "
    print TrainNumbers

    if not DryRun and not os.path.isdir(LocalPath):
        print "Creating directory " + LocalPath
        os.makedirs(LocalPath)

    if Metadataset:
        StartDownloadMetadataset(LocalPath, Datasets, Config["children"], TrainNumbers[0], Config["train"], Overwrite, FileName, DryRun)
    else:
        StartDownload(LocalPath, Datasets, TrainNumbers, Config["train"], Overwrite, FileName, DryRun)

if __name__ == '__main__':
    # FinalMergeLocal.py executed as script

    parser = argparse.ArgumentParser(description='Script to download results from LEGO trains.')
    parser.add_argument('trainNumber', metavar='trainNumber',
                        help='Train numbers to be downloaded and merged. Use ":" to define a range, and "," for a list')
    parser.add_argument('--yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--overwrite', metavar='overwrite',
                        default=0,
                        help='Overwrite level [0-4]. 0 = no overwrite')
    parser.add_argument('--file-name', metavar='AnalysisResults.root',
                        default='AnalysisResults.root',
                        help='Output file name')
    parser.add_argument('--dry-run', action='store_const',
                        default=False, const=True,
                        help='Dry run')
    args = parser.parse_args()

    trainNumberList = args.trainNumber.split(",")
    trainNumbers = []
    for trainNumberRange in trainNumberList:
        myrange = trainNumberRange.split(":")
        trainNumbers.extend(range(int(myrange[0]), int(myrange[len(myrange) - 1]) + 1))

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(trainNumbers, config, args.file_name, args.overwrite, args.dry_run)
