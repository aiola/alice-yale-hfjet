#!/usr/bin/env python

import os
import subprocess
import sys
import shutil

import MergeFiles
import ScaleResults

def CreateRunLists():
    
    RunLists=dict(LHC15i2b=["114786", "114798", "114918", "114920", "114924", "114930", "114931", "115186", "115193", "115310", "115318", "115322", "115328", "115335", "115345", "115393", "115399", "115401", "115414", "115521", "116079", "116081", "116102", "116288", "116402", "116403", "116562", "116571", "116574", "116643", "116645", "117048", "117050", "117052", "117053", "117059", "117060", "117063", "117092", "117099", "117109", "117112", "117116", "117220", "117222"], LHC15i2c=["118506", "118507", "118512", "118518", "118556", "118558", "118560", "118561", "119159", "119161", "119163", "119841", "119842", "119844", "119845", "119846", "119849", "119853", "119856", "119859", "119862", "120067", "120069", "120072", "120073", "120076", "120079", "120244", "120503", "120504", "120505", "120616", "120617", "120671", "120741", "120750", "120758", "120820", "120821", "120822", "120823", "120824", "120825", "120829", "121039", "121040"], LHC15i2d=["118506", "118507", "118512", "118518", "118556", "118558", "118560", "118561", "119159", "119161", "119163", "119841", "119842", "119844", "119845", "119846", "119849", "119853", "119856", "119859", "119862", "120067", "120069", "120072", "120073", "120076", "120079", "120244", "120503", "120504", "120505", "120616", "120617", "120671", "120741", "120750", "120758", "120820", "120821", "120822", "120823", "120824", "120825", "120829", "121039", "121040"], LHC15i2e=["128366", "128452", "128486", "128494", "128495", "128498", "128503", "128504", "128505", "128506", "128582", "128590", "128592", "128594", "128596", "128605", "128609", "128611", "128615", "128621", "128677", "128678", "128777", "128778", "128819", "128820", "128823", "128824", "128833", "128834", "128835", "128836", "128843", "128850", "128853", "128855", "128913", "129042", "129512", "129513", "129514", "129515", "129516", "129519", "129520", "129521", "129523", "129524", "129525", "129527", "129528", "129536", "129540", "129586", "129587", "129599", "129639", "129641", "129647", "129650", "129651", "129652", "129653", "129659", "129666", "129723", "129725", "129726", "129729", "129734", "129735", "129736", "129738", "129742", "129744", "129959", "129960", "129961", "129962", "129966", "129983", "130149", "130151", "130157", "130158", "130168", "130172", "130178", "130342", "130343", "130354", "130356", "130358", "130360", "130375", "130479", "130480", "130481", "130517", "130519", "130520", "130524", "130526", "130601", "130608", "130609", "130620", "130621", "130623", "130628", "130696", "130704", "130793", "130795", "130798", "130799", "130802", "130803", "130804", "130834", "130840", "130842", "130844", "130847", "130848", "130850"])

    #RunLists=dict(LHC15i2b=["114786", "114798"], LHC15i2c=["118506", "118507"], LHC15i2d=["118506", "118507"], LHC15i2e=["128366", "128452"])
    
    return RunLists

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

if len(sys.argv)>1:
    TrainNumbers=sys.argv[1].split(":")
else:
    print sys.argv[0]+" [trainNumber1:trainNumber2:...] [overwrite=0] [year=2015] [dataset=LHC15i2b]  [trainName=Jets_EMC_pp_MC]"
    exit()

StrippedTrainNumbers=[]
    
for i in range(len(TrainNumbers)):
    StrippedTrainNumbers[i:]=[TrainNumbers[i][:3]]

if len(sys.argv) > 2:
    Overwrite = int(sys.argv[2])
else:
    Overwrite = 0

if len(sys.argv)>3:
    Year = sys.argv[3]
else:
    Year = "2015"

if len(sys.argv)>4:
    UserDataset=sys.argv[4]
else:
    UserDataset="LHC15i2b"

if len(sys.argv)>5:
    TrainName=sys.argv[5]
else:
    TrainName="Jets_EMC_pp_MC"

RunLists=CreateRunLists()

Datasets=[]

if UserDataset=="LHC15i2x":
    print "Will work on LHC15i2{b,c,d,e}."
    Datasets[0:]=["LHC15i2b"]
    Datasets[1:]=["LHC15i2c"]
    Datasets[2:]=["LHC15i2d"]
    Datasets[3:]=["LHC15i2e"]
else:
    if UserDataset in RunLists:
        print "Dataset "+UserDataset+" was chosen."
        Datasets[0:]=[UserDataset]
        print RunLists[UserDataset]
    else:
        print "Dataset "+UserDataset+" is not defined."
        exit()

MinPtHardBin=0
MaxPtHardBin=8

if (len(Datasets)!=len(TrainNumbers)):
    print "The number of datasets {0} must be the same as the number of trains {1}.".format(len(Datasets),len(TrainNumbers))
    print "Datasets are"
    print Datasets
    print "Trains are"
    print TrainNumbers
    exit()

LocalPath="{0}/{1}".format(JetResults, TrainName)

for StrippedTrainNumber in StrippedTrainNumbers:
    LocalPath+="_"+StrippedTrainNumber

print "Train: "+TrainName
print "Local path: "+LocalPath
print "Overwrite mode: {0}".format(Overwrite)
print "Train numbers are: "
print TrainNumbers

if not os.path.isdir(LocalPath):
    print "Creating directory "+LocalPath
    os.makedirs(LocalPath)

if os.path.exists("./"+UserDataset+".xsec.root"):
    shutil.copy("./"+UserDataset+".xsec.root", LocalPath+"/"+UserDataset+".xsec.root")

FileList = []
FinalFileList = []
    
for PtHardBin in range(MinPtHardBin, MaxPtHardBin+1):
    FileList[:] = []
    for Dataset,TrainNumber in zip(Datasets,TrainNumbers):
        AlienPath="alien:///alice/sim/"+Year+"/"+Dataset
        for Run in RunLists[Dataset]:
            dest="{0}/{1}/{2}".format(LocalPath, PtHardBin, Run)
            if not os.path.isdir(dest):
                print "Creating directory "+dest
                os.makedirs(dest)
                
            dest+="/AnalysisResults.root"
            if os.path.exists(dest) and Overwrite > 4:
                print "Deleting file: "+dest
                os.remove(dest)

            if not os.path.exists(dest):
                AlienFile="{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root".format(AlienPath, Run, PtHardBin, TrainName, TrainNumber)
                print "Copying from alien location '"+AlienFile+"' to local location '"+dest+"'"
                #subprocess.call(["alien_cp", AlienFile, dest])

            if os.path.exists(dest):
                FileList.append(dest)

    dest="{0}/{1}/AnalysisResults.root".format(LocalPath, PtHardBin)
    if os.path.exists(dest) and Overwrite>3:
        print "Deleting file "+dest
        os.remove(dest)
 
    if not os.path.exists(dest):
      print "Merging for pT hard bin: {0}. Total number of files is {1}".format(PtHardBin, len(FileList))
      MergeFiles.MergeFiles(dest, FileList)

    dest="{0}/{1}/ScaledResults.root".format(LocalPath, PtHardBin)
    if os.path.exists(dest) and Overwrite > 2:
        print "Removing file "+dest
        os.remove(dest)
 
    if not os.path.exists(dest):
      print "Scaling for pT hard bin: {0}".format(PtHardBin)
      ScaleResults.ScaleResults(LocalPath, PtHardBin, Dataset)
  
    if os.path.exists(dest) and PtHardBin > 0:
      FinalFileList.append(dest)

dest="{0}/AnalysisResults.root".format(LocalPath)
if os.path.exists(dest) and Overwrite > 1:
    print "Deleting file "+dest
    os.remove(dest)

if not os.path.exists(dest):
    MergeFiles.MergeFiles(dest, FinalFileList, 2)
     
print "Done."

subprocess.call(["ls", LocalPath])
