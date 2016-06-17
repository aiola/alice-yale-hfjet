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
    output = subprocess.check_output(["alien_find", "-d", "-l", "1", SearchPath, TrainName+"/"+str(TrainNumber)+"_*"], universal_newlines=True)
    #print("alien_find", "-d", "-l", "1", SearchPath, TrainName+"/"+str(TrainNumber)+"_*")
    #print(output)
    i = output.rfind(str(TrainNumber)+"_")
    j = len(str(TrainNumber)) + 14 + i
    FullTrainNumber = output[i:j]
    return FullTrainNumber

def RegularMerging(LocalPath, Datasets, TrainName, TrainNumbers, Overwrite, Year, MC, Pass, jetTypes, triggers, skipList, acceptList, NoDownload):
    FileList = []
    FinalFileList = []
    
    FileName = "AnalysisResults_TriggerQA.root"

    for Dataset,TrainNumber in zip(sorted(Datasets.keys()),TrainNumbers):    
        FileList[:] = []
        
        if not NoDownload:
            if MC:
                FirtRun = Datasets[Dataset][0]
                AlienPath="/alice/sim/"+Year+"/"+Dataset
            else:
                FirtRun = "000" + RunLists[Dataset][0]
                AlienPath="/alice/data/"+Year+"/"+Dataset
                
            SearchPath = "{0}/{1}/{2}".format(AlienPath, FirtRun, Pass)
            FullTrainNumber = GetFullTrainNumber(SearchPath, TrainName, TrainNumber)
            AlienPath = "alien://" + AlienPath
            
        for Run in RunLists[Dataset]:
            if not MC:
                Run = "000" + Run
            dest="{0}/{1}/{2}".format(LocalPath, Dataset, Run)
            if not os.path.isdir(dest):
                print "Creating directory "+dest
                os.makedirs(dest)

            dest+="/AnalysisResults.root"
            
            if not NoDownload:
                if os.path.exists(dest) and Overwrite > 4:
                    print "Deleting file: "+dest
                    os.remove(dest)
    
                if not os.path.exists(dest):
                    AlienFile="{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root".format(AlienPath, Run, Pass, TrainName, FullTrainNumber)
                    print "Copying from alien location '"+AlienFile+"' to local location '"+dest+"'"
                    subprocess.call(["alien_cp", AlienFile, dest])

            if os.path.exists(dest):
                FileList.append(dest)

        dest="{0}/{1}/{2}".format(LocalPath, Dataset, FileName)
        if os.path.exists(dest) and Overwrite > 3:
            print "Deleting file "+dest
            os.remove(dest)

        if not os.path.exists(dest):
            print "Merging for dataset: {0}. Total number of files is {1}".format(Dataset, len(FileList))
            MergeFiles.MergeFiles(dest, FileList, skipList, acceptList)

        FinalFileList.append(dest)

    dest="{0}/{1}".format(LocalPath, FileName)
    if os.path.exists(dest) and Overwrite > 1:
        print "Deleting file "+dest
        os.remove(dest)

    if not os.path.exists(dest):
        MergeFiles.MergeFiles(dest, FinalFileList, skipList, acceptList, 2)

    print "Done."

    subprocess.call(["ls", LocalPath])

def PtHardBinMerging(LocalPath, Datasets, TrainName, TrainNumbers, Overwrite, MinPtHardBin, MaxPtHardBin, Year):
    FileList = []
    FinalFileList = []

    for PtHardBin in range(MinPtHardBin, MaxPtHardBin+1):
        FileList[:] = []
        for Dataset,TrainNumber in zip(sorted(Datasets.iterkeys()),TrainNumbers):
            AlienPath="/alice/sim/"+str(Year)+"/"+Dataset
            FirtRun = Datasets[Dataset][0]
            SearchPath = "{0}/{1}".format(AlienPath, FirtRun)
            FullTrainNumber = GetFullTrainNumber(SearchPath, TrainName, TrainNumber)
            AlienPath = "alien://" + AlienPath
            
            for Run in Datasets[Dataset]:
                dest="{0}/{1}/{2}/{3}".format(LocalPath, Dataset, PtHardBin, Run)
                if not os.path.isdir(dest):
                    print "Creating directory "+dest
                    os.makedirs(dest)

                dest+="/AnalysisResults.root"
                if os.path.exists(dest) and Overwrite > 4:
                    print "Deleting file: "+dest
                    os.remove(dest)

                if not os.path.exists(dest):
                    AlienFile="{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root".format(AlienPath, Run, PtHardBin, TrainName, FullTrainNumber)
                    print "Copying from alien location '"+AlienFile+"' to local location '"+dest+"'"
                    subprocess.call(["alien_cp", AlienFile, dest])

                if os.path.exists(dest):
                    FileList.append(dest)

            dest="{0}/{1}/{2}/AnalysisResults.root".format(LocalPath, Dataset, PtHardBin)
            if os.path.exists(dest) and Overwrite>3:
                print "Deleting file "+dest
                os.remove(dest)

            if not os.path.exists(dest):
                print "Merging for pT hard bin: {0}. Total number of files is {1}".format(PtHardBin, len(FileList))
                MergeFiles.MergeFiles(dest, FileList)

    print "Done."

    subprocess.call(["ls", LocalPath])

def StartMerging(TrainNumbers, config, Overwrite=0, NoDownload=False):
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
    print "Overwrite mode: {0}".format(Overwrite)
    print "Train numbers are: "
    print TrainNumbers
    
    if not os.path.isdir(LocalPath):
        print "Creating directory "+LocalPath
        os.makedirs(LocalPath)
        
    if config["pt_hard_bins"]:
        PtHardBinMerging(LocalPath, Datasets, config["train"], TrainNumbers, Overwrite, config["min_pt_hard_bin"], config["max_pt_hard_bin"], config["year"])
    else:
        RegularMerging(LocalPath, Datasets, config["train"], TrainNumbers, Overwrite, config["year"], config["MC"], config["pass"], 
                       config["jet_types"], config["triggers"], config["skip_list"], config["accept_list"], NoDownload)

if __name__ == '__main__':
    # FinalMergeLocal.py executed as script
    
    parser = argparse.ArgumentParser(description='Local final merging for LEGO train results.')
    parser.add_argument('trainNumber', metavar='trainNumber',
                        help='Train numbers to be downloaded and merged')
    parser.add_argument('--overwrite', metavar='overwrite',
                        default=0,
                        help='Overwrite level [0-4]. 0 = no overwrite')
    parser.add_argument('--config', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--no-download', action='store_const',
                        default=False, const=True,
                        help='Do not download')
    args = parser.parse_args()

    TrainNumberList = args.trainNumber.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range)-1])+1))
    
    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    StartMerging(TrainNumbers, config, args.overwrite, args.no_download)
