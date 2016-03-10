#!/usr/bin/env python

import os
import subprocess
import sys
import shutil
import argparse

import MergeFiles
import ScaleResults

def CreateRunLists():
    
    RunLists=dict(LHC15i2b=["114786", "114798", "114918", "114920", "114924", "114930", "114931", "115186", "115193", "115310", "115318", "115322", 
                            "115328", "115335", "115345", "115393", "115399", "115401", "115414", "115521", "116079", "116081", "116102", "116288", 
                            "116402", "116403", "116562", "116571", "116574", "116643", "116645", "117048", "117050", "117052", "117053", "117059", 
                            "117060", "117063", "117092", "117099", "117109", "117112", "117116", "117220", "117222"], 
                  LHC15i2c=["118506", "118507", "118512", "118518", "118556", "118558", "118560", "118561", "119159", "119161", "119163", "119841", 
                            "119842", "119844", "119845", "119846", "119849", "119853", "119856", "119859", "119862", "120067", "120069", "120072", 
                            "120073", "120076", "120079", "120244", "120503", "120504", "120505", "120616", "120617", "120671", "120741", "120750", 
                            "120758", "120820", "120821", "120822", "120823", "120824", "120825", "120829", "121039", "121040"], 
                  LHC15i2d=["118506", "118507", "118512", "118518", "118556", "118558", "118560", "118561", "119159", "119161", "119163", "119841", 
                            "119842", "119844", "119845", "119846", "119849", "119853", "119856", "119859", "119862", "120067", "120069", "120072", 
                            "120073", "120076", "120079", "120244", "120503", "120504", "120505", "120616", "120617", "120671", "120741", "120750", 
                            "120758", "120820", "120821", "120822", "120823", "120824", "120825", "120829", "121039", "121040"], 
                  LHC15i2e=["128366", "128452", "128486", "128494", "128495", "128498", "128503", "128504", "128505", "128506", "128582", "128590", 
                            "128592", "128594", "128596", "128605", "128609", "128611", "128615", "128621", "128677", "128678", "128777", "128778", 
                            "128819", "128820", "128823", "128824", "128833", "128834", "128835", "128836", "128843", "128850", "128853", "128855", 
                            "128913", "129042", "129512", "129513", "129514", "129515", "129516", "129519", "129520", "129521", "129523", "129524", 
                            "129525", "129527", "129528", "129536", "129540", "129586", "129587", "129599", "129639", "129641", "129647", "129650", 
                            "129651", "129652", "129653", "129659", "129666", "129723", "129725", "129726", "129729", "129734", "129735", "129736", 
                            "129738", "129742", "129744", "129959", "129960", "129961", "129962", "129966", "129983", "130149", "130151", "130157", 
                            "130158", "130168", "130172", "130178", "130342", "130343", "130354", "130356", "130358", "130360", "130375", "130479", 
                            "130480", "130481", "130517", "130519", "130520", "130524", "130526", "130601", "130608", "130609", "130620", "130621", 
                            "130623", "130628", "130696", "130704", "130793", "130795", "130798", "130799", "130802", "130803", "130804", "130834", 
                            "130840", "130842", "130844", "130847", "130848", "130850"], 
                  LHC16a1c=["179796", "179802", "179803", "179858", "179859", "179916", "179917", "179918", "179919", "179920", "180000", "180042", 
                            "180044", "180127", "180129", "180130", "180131", "180132", "180133", "180230", "180500", "180501", "180510", "180515", 
                            "180517", "180561", "180562", "180564", "180567", "180569", "180716", "180717", "180719", "180720", "182017", "182018", 
                            "182022", "182023", "182106", "182110", "182111", "182207", "182289", "182295", "182297", "182299", "182300", "182302", 
                            "182322", "182323", "182324", "182325", "182624", "182635", "182684", "182687", "182691", "182692"],
                  
                  LHC12a=["177182", "177180", "177173", "177167", "177160", "177157", "177148", "177011", "176929", "176927", "176926", "176924", 
                          "176859", "176854", "176849", "176753", "176752", "176749", "176730", "176715", "176707", "176704", "176701", "176661"],
                  LHC12b=["178030", "178026", "178025", "178024", "178018", "177942", "177869", "177864", "177861", "177860", "177858", "177810", 
                          "177805", "177804", "177799", "177798", "177682", "177681", "177680", "177679", "177671", "177624", "177620", "177612", 
                          "177597", "177592", "177580"],
                  LHC12c=["182692", "182691", "182684", "180501", "180500", "180230", "180133", "180132", "180130", "180129", "180127", "180044", 
                          "180042", "180000", "179920", "179919", "179918", "179917", "179916", "179859", "179858", "179803", "179802", "179796", 
                          "179639", "179638", "179621", "179618"],
                  LHC12d=["186320", "186319", "186208", "186205", "186167", "186165", "186164", "186163", "186083", "185375", "185371", "185363", 
                          "185362", "185361", "185360", "185359", "185356", "185351", "185350", "185349", "185303", "185302", "185300", "185299", 
                          "185296", "185293", "185292", "185291", "185289", "185284", "185282", "185221", "185217", "185208", "185206", "185203", 
                          "185198", "185196", "185189", "185164", "185160", "185157", "185134", "185132", "185127", "185126", "185116", "185031", 
                          "185029", "184990", "184988", "184987", "184968", "184967", "184964", "184938", "184933", "184786", "184784", "184687", 
                          "184682", "184678", "184673", "184371", "184215", "184209", "184208", "184188", "184138", "184137", "184135", "184132", 
                          "184131", "184127", "183916", "183913"],
                  LHC12e=["186389", "186388", "186387", "186386", "186365"],
                  LHC12f=["188101", "187796", "187791", "187785", "187749", "187739", "187698", "187656", "187627", "187624", "187623", "187510", 
                          "187489", "187488", "187487", "187341", "187340", "187339", "187203", "187202", "187152", "187151", "187150", "187149", 
                          "187146", "187143", "186992", "186990", "186969", "186966", "186939", "186938", "186937", "186851", "186811", "186692", 
                          "186690", "186689", "186688", "186668"],
                  LHC12g=["188503", "188500", "188499", "188447", "188446", "188444", "188443", "188442"],
                  LHC12h=["192732", "192731", "192729", "192535", "192534", "192510", "192505", "192499", "192492", "192468", "192461", "192453", 
                          "192417", "192415", "192349", "192348", "192347", "192246", "192205", "192202", "192201", "192200", "192199", "192197", 
                          "192194", "192177", "192174", "192172", "192141", "192140", "192136", "192128", "192075", "192073", "192072", "192004", 
                          "191248", "191247", "191245", "191244", "191242", "191234", "191232", "191231", "191229", "191227", "190307", "190305", 
                          "190303", "190216", "190215", "190214", "190212", "190210", "190150"],
                  LHC12i=["193151", "193094", "193051", "193049", "193047", "193014", "193011", "193010", "193008", "193007", "193005", "193004", 
                          "192824", "192822", "192820", "192779", "192778", "192775", "192772"]
                  )

    
    
    return RunLists

def GetFullTrainNumber(SearchPath, TrainNumber):
    output = subprocess.check_output(["alien_find", "-d", "-l", "1", SearchPath, TrainNumber+"_*"], universal_newlines=True)
    #print(output)
    i = output.rfind(TrainNumber+"_")
    j = len(str(TrainNumber)) + 14 + i
    FullTrainNumber = output[i:j]
    return FullTrainNumber

def RegularMerging(UserDataset, LocalPath, Datasets, RunLists, TrainNumbers, TrainName, Overwrite, Year, MC, Pass):
    
    jetTypes = ["Charged", "Full"]
    triggers = ["AnyINT", "EMC", "EMCEGA", "EMCEJE"]
    skipList = []
    
    for trigger in triggers:
        skipList.append("AliEmcalTriggerQATaskPP_{0}_histos".format(trigger))
            
    for jetType in jetTypes:
        for trigger in triggers:
            skipList.append("AliAnalysisTaskDmesonJets_{0}_{1}_histos".format(jetType,trigger))
    
    FileList = []
    FinalFileList = []

    for Dataset,TrainNumber in zip(Datasets,TrainNumbers):    
        FileList[:] = []
        if MC:
            FirtRun = RunLists[Dataset][0]
            AlienPath="/alice/sim/"+Year+"/"+Dataset
        else:
            FirtRun = "000" + RunLists[Dataset][0]
            AlienPath="/alice/data/"+Year+"/"+Dataset
            
        SearchPath = "{0}/{1}/{2}/PWGJE/{3}/".format(AlienPath, FirtRun, Pass, TrainName)
        FullTrainNumber = GetFullTrainNumber(SearchPath, TrainNumber)
        
        AlienPath = "alien://" + AlienPath
            
        for Run in RunLists[Dataset]:
            if not MC:
                Run = "000" + Run
            dest="{0}/{1}/{2}".format(LocalPath, Dataset, Run)
            if not os.path.isdir(dest):
                print "Creating directory "+dest
                os.makedirs(dest)

            dest+="/AnalysisResults.root"
            if os.path.exists(dest) and Overwrite > 4:
                print "Deleting file: "+dest
                os.remove(dest)

            if not os.path.exists(dest):
                AlienFile="{0}/{1}/{2}/PWGJE/{3}/{4}/AnalysisResults.root".format(AlienPath, Run, Pass, TrainName, FullTrainNumber)
                print "Copying from alien location '"+AlienFile+"' to local location '"+dest+"'"
                subprocess.call(["alien_cp", AlienFile, dest])

            if os.path.exists(dest):
                FileList.append(dest)

        dest="{0}/{1}/AnalysisResults.root".format(LocalPath, Dataset)
        if os.path.exists(dest) and Overwrite > 3:
            print "Deleting file "+dest
            os.remove(dest)

        if not os.path.exists(dest):
            print "Merging for dataset: {0}. Total number of files is {1}".format(Dataset, len(FileList))
            MergeFiles.MergeFiles(dest, FileList, skipList)

    FinalFileList.append(dest)

    dest="{0}/AnalysisResults.root".format(LocalPath)
    if os.path.exists(dest) and Overwrite > 1:
        print "Deleting file "+dest
        os.remove(dest)

    if not os.path.exists(dest):
        MergeFiles.MergeFiles(dest, FinalFileList, skipList, 2)

    print "Done."

    subprocess.call(["ls", LocalPath])

def PtHardBinMerging(UserDataset, LocalPath, Datasets, TrainNumbers, Overwrite, Year):
    MinPtHardBin=0
    MaxPtHardBin=8
    
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
                    subprocess.call(["alien_cp", AlienFile, dest])

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

def StartMerging(TrainNumbers, Overwrite=0, Year="2015", UserDataset="LHC15i2b", TrainName="Jets_EMC_pp_MC", PtHardMerging=False, MC=False):
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

    StrippedTrainNumbers=[]
            
    for i in range(len(TrainNumbers)):
        StrippedTrainNumbers[i:]=[TrainNumbers[i][:3]]

    RunLists=CreateRunLists()

    Datasets=[]

    if UserDataset=="LHC15i2x":
        print "Will work on LHC15i2{b,c,d,e}."
        Datasets[0:]=["LHC15i2b"]
        Datasets[1:]=["LHC15i2c"]
        Datasets[2:]=["LHC15i2d"]
        Datasets[3:]=["LHC15i2e"]
    elif UserDataset=="LHC12x":
        print "Will work on LHC12{a,b,c,d,e,f,g,h,i}."
        Datasets[0:]=["LHC12a"]
        Datasets[1:]=["LHC12b"]
        Datasets[2:]=["LHC12c"]
        Datasets[3:]=["LHC12d"]
        Datasets[4:]=["LHC12e"]
        Datasets[5:]=["LHC12f"]
        Datasets[6:]=["LHC12g"]
        Datasets[7:]=["LHC12h"]
        Datasets[8:]=["LHC12i"]
    else:
        if UserDataset in RunLists:
            print "Dataset "+UserDataset+" was chosen."
            Datasets[0:]=[UserDataset]
            print RunLists[UserDataset]
        else:
            print "Dataset "+UserDataset+" is not defined."
            exit()

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
        
    if PtHardMerging:
        PtHardBinMerging(UserDataset, LocalPath, Datasets, TrainNumbers, Overwrite, Year)
    else:
        RegularMerging(UserDataset, LocalPath, Datasets, RunLists, TrainNumbers, TrainName, Overwrite, Year, MC, "pass2")

if __name__ == '__main__':
    # FinalMergeLocal.py executed as script
    
    parser = argparse.ArgumentParser(description='Local final merging for MC production in pT hard bins.')
    parser.add_argument('trainNumber', metavar='trainNumber',
                        help='Train numbers to be downloaded and merged')
    parser.add_argument('--overwrite', metavar='overwrite',
                        default=0,
                        help='Overwrite level [0-4]. 0 = no overwrite')
    parser.add_argument('--year', metavar='year',
                        default='2015',
                        help='MC production year')
    parser.add_argument('--dataset', metavar='dataset',
                        default='LHC15i2b',
                        help='MC production name')
    parser.add_argument('--trainName', metavar='trainName',
                        default='Jets_EMC_pp_MC',
                        help='Train name')
    parser.add_argument('--pthard', action='store_const',
                        default=False, const=True,
                        help='Pt hard bin merging mode')
    parser.add_argument('--MC', action='store_const',
                        default=False, const=True,
                        help='MC mode')
    args = parser.parse_args()

    TrainNumbers=args.trainNumber.split(":")

    StartMerging(TrainNumbers, args.overwrite, args.year, args.dataset, args.trainName, args.pthard, args.MC)
