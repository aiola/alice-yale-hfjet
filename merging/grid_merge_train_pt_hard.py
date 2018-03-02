#!/usr/bin/env python

import os
import subprocess
import argparse
import yaml


def AlienDelete(fileName):
    if fileName.find("alien://") == -1:
        fname = fileName
    else:
        fname = fileName[8:]

    subprocessCall(["alien_rm", fname])


def AlienFileExists(fileName):
    if fileName.find("alien://") == -1:
        fname = fileName
    else:
        fname = fileName[8:]

    fileExists = True
    try:
        subprocessCheckCall(["alien_ls", fname])
    except subprocess.CalledProcessError:
        fileExists = False

    return fileExists


def AlienCopy(source, destination, attempts=3, overwrite=False):
    i = 0
    fileExists = False

    if AlienFileExists(destination):
        if overwrite:
            AlienDelete(destination)
        else:
            return True

    if destination.find("alien://") == -1:
        dest = "alien://{0}".format(destination)
    else:
        dest = destination

    while True:
        subprocessCall(["alien_cp", source, dest])
        i += 1
        fileExists = AlienFileExists(destination)
        if fileExists:
            break
        if i >= attempts:
            print("After {0} attempts I could not copy {1} to {2}".format(i, source, dest))
            break

    return fileExists


def subprocessCall(cmd):
    print(cmd)
    return  subprocess.call(cmd)


def subprocessCheckCall(cmd):
    print(cmd)
    return subprocess.check_call(cmd)


def subprocessCheckOutput(cmd):
    print(cmd)
    return subprocess.check_output(cmd, universal_newlines=True)


def GetFullTrainNumber(SearchPath, TrainName, TrainNumber):
    output = subprocessCheckOutput(["alien_find", "-d", "-l", "2", SearchPath, TrainName + "/" + str(TrainNumber) + "_*"])
    i = output.rfind(str(TrainNumber) + "_")
    j = len(str(TrainNumber)) + 14 + i
    FullTrainNumber = output[i:j]
    return FullTrainNumber


def PtHardBinMerging(LocalPath, Datasets, TrainName, TrainNumbers, MinPtHardBin, MaxPtHardBin, Year, InvertedScheme, AliPhysicsVersion, TestMode, GridTestMode, GridUpdate):
    for Dataset, TrainNumber in zip(sorted(Datasets), TrainNumbers):

        AlienPath = "/alice/sim/" + str(Year) + "/" + Dataset
        FullTrainNumber = GetFullTrainNumber(AlienPath, TrainName, TrainNumber)

        dest = "/alice/cern.ch/user/s/saiola/{0}/{1}".format(TrainName, FullTrainNumber)

        subprocessCall(["alien_mkdir", "-p", dest])

        executableFile = "{0}/grid_merge_train_pt_hard.sh".format(dest)
        macroFile = "{0}/MergeFiles.C".format(dest)
        validationScript = "{0}/grid_merge_train_pt_hard_validation.sh".format(dest)

        AlienCopy("./grid_merge_train_pt_hard.sh", "alien://{0}".format(executableFile), 3, GridUpdate)
        AlienCopy("./MergeFiles.C", "alien://{0}".format(macroFile), 3, GridUpdate)
        AlienCopy("./grid_merge_train_pt_hard_validation.sh", "alien://{0}".format(validationScript), 3, GridUpdate)

        for PtHardBin in range(MinPtHardBin, MaxPtHardBin + 1):

            localDest = "{0}/{1}/{2}".format(LocalPath, Dataset, PtHardBin)
            if not os.path.isdir(localDest):
                print "Creating directory " + localDest
                os.makedirs(localDest)

            alienOutputFile = "{0}/{1}/output".format(dest, PtHardBin)

            if AlienFileExists("{0}/AnalysisResults.root".format(alienOutputFile)):
                print("Output file alien://{0}/AnalysisResults.root already exists. Skipping.".format(alienOutputFile))
                continue

            AlienDelete(alienOutputFile)

            if InvertedScheme:
                SearchPath = "/{0}/*/PWGJE/{1}/{2}/AnalysisResults.root".format(PtHardBin, TrainName, FullTrainNumber)
            else:
                SearchPath = "/*/{0}/PWGJE/{1}/{2}/AnalysisResults.root".format(PtHardBin, TrainName, FullTrainNumber)

            xmlColl = subprocessCheckOutput(["alien_find", "-x", "merge_files.xml", AlienPath, SearchPath])

            subprocessCall(["alien_mkdir", "{0}/{1}".format(dest, PtHardBin)])

            jdlContent = "# This is the startup script \n\
Executable = \"{executable}\"; \n\
# Time after which the job is killed (500 min.) \n\
TTL = \"14400\"; \n\
OutputDir = \"{dest}/{pt_hard}/output\"; \n\
Output = {{ \n\
\"AnalysisResults*.root\", \n\
\"std*\" \n\
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
\"LF:{dest}/MergeFiles.C\", \n\
\"LF:{dest}/{pt_hard}/merge_files.xml\" \n\
}}; \n\
InputDataList = \"merge_files.xml\"; \n\
InputDataListFormat = \"xml-single\"; \n\
InputDataCollection = {{ \n\
\"LF:{dest}/{pt_hard}/merge_files.xml,nodownload\" \n\
}}; \n\
Validationcommand = \"{validationScript}\"; \n\
".format(executable=executableFile, dest=dest, aliphysics=AliPhysicsVersion, pt_hard=PtHardBin, validationScript=validationScript)

            localJdlFile = "{0}/grid_merge_train_pt_hard.jdl".format(localDest)

            f = open(localJdlFile, 'w')
            f.write(jdlContent)
            f.close()

            localXmlFile = "{0}/merge_files.xml".format(localDest)

            f = open(localXmlFile, 'w')
            f.write(xmlColl)
            f.close()

            jdlFile = "{0}/{1}/grid_merge_train_pt_hard.jdl".format(dest, PtHardBin)

            AlienCopy(localJdlFile, "alien://{0}".format(jdlFile), 3, GridUpdate)
            AlienCopy(localXmlFile, "alien://{0}/{1}/merge_files.xml".format(dest, PtHardBin), 3, GridUpdate)

            if TestMode:
                subprocessCall(["chmod", "+x", "./grid_merge_train_pt_hard.sh"])
                return
            else:
                subprocessCall(["alien_submit", "alien://{0}".format(jdlFile)])
                if GridTestMode:
                    return

    print "Done."

    subprocessCall(["ls", LocalPath])


def StartMerging(TrainNumbers, config, AliPhysicsVersion, TestMode, GridTestMode, GridUpdate):
    try:
        rootPath = subprocess.check_output(["which", "root"]).rstrip()
        alirootPath = subprocess.check_output(["which", "aliroot"]).rstrip()
        alienPath = subprocess.check_output(["which", "alien-token-info"]).rstrip()
    except subprocess.CalledProcessError:
        print "Environment is not configured correctly!"
        exit()

    if not AliPhysicsVersion:
        print("Please provide the AliPhysics version!")
        exit(1)

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

    Datasets = config["datasets"]

    if (len(Datasets) != len(TrainNumbers)):
        print "The number of datasets {0} must be the same as the number of trains {1}.".format(len(Datasets), len(TrainNumbers))
        print "Datasets are"
        print Datasets
        print "Trains are"
        print TrainNumbers
        exit()

    LocalPath = "{0}/{1}".format(JetResults, config["train"])

    for TrainNumber in TrainNumbers:
        LocalPath += "_" + str(TrainNumber)

    print "Train: " + config["train"]
    print "Local path: " + LocalPath
    print "Train numbers are: "
    print TrainNumbers

    if not os.path.isdir(LocalPath):
        print "Creating directory " + LocalPath
        os.makedirs(LocalPath)

    if "inverted_scheme" in config: InvertedScheme = config["inverted_scheme"]
    else: InvertedScheme = False

    if config["pt_hard_bins"]:
        PtHardBinMerging(LocalPath, Datasets, config["train"], TrainNumbers, config["min_pt_hard_bin"], config["max_pt_hard_bin"],
                         config["year"], InvertedScheme, AliPhysicsVersion, TestMode, GridTestMode, GridUpdate)
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
    parser.add_argument("--update", action='store_const',
                        default=False, const=True,
                        help='Update all scripts and macros on the grid.')
    args = parser.parse_args()

    TrainNumberList = args.trainNumber.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range) - 1]) + 1))

    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    StartMerging(TrainNumbers, config, args.aliphysics, args.test, args.test_grid, args.update)
