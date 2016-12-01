#!/usr/bin/env python
#python script to compare data analysis results

import argparse
import yaml
import IPython
import ROOT
import subprocess

globalList = []

def main(configs):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    baseline = None
    baselineName = None

    files = dict()
    for config in configs:
        fname = "{input_path}/{train}/{analysis}.root".format(input_path=config["input_path"], train=config["train"], analysis=config["name"])
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Skipping file {0}, since it was not open successfully".format(fname))
            continue
        if not baseline:
            baseline = file
            baselineName = config["name"]
        else:
            files[config["name"]] = file

    if len(files) < 1:
        print(files)
        print("Less than 2 analysis to compare, quitting...")
        exit(0)

    mylist = CompareObjects(baseline, files)
    outputFileName = "_".join([config["name"] for config in configs])
    outputFileName = "{0}/Comparison_{1}.root".format(configs[0]["input_path"], outputFileName)
    outputFile = ROOT.TFile(outputFileName, "recreate")
    outputFile.cd()
    rlist = GenerateRootList(mylist)
    for robj in rlist:
        robj.Write(robj.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()
    print("Comparison results stored in {0}".format(outputFileName))

def GenerateRootList(mylist):
    if len(mylist) < 1:
        return None
    rlist = ROOT.TList()
    for name,obj in mylist.iteritems():
        robj = None
        if isinstance(obj, dict):
            robj = GenerateRootList(obj)
            if robj:
                robj.SetName(name)            
        elif isinstance(obj, ROOT.TObject):
            robj = obj
        else:
            print("Unexpected object type")
            print(obj)
        if robj:
            rlist.Add(robj)
    return rlist

def CompareObjects(baseline, inputObjects):
    print("Comparing object {0}".format(baseline.GetName()))
    mylist = dict()
    if len(inputObjects) < 1:
        print("Nothing to compare with, returning")
        return
    if isinstance(baseline, ROOT.TDirectory):
        root_keys = baseline.GetListOfKeys()
        for root_key in root_keys:    
            baselineObj = baseline.Get(root_key.GetName())
            objects = dict()
            for name,file in inputObjects.iteritems():
                obj = file.Get(root_key.GetName())
                if not obj:
                    print("{0} not found for analysis {1}".format())
                    continue
                objects[name] = obj
            mylist[root_key.GetName()] = CompareObjects(baselineObj, objects)

    elif isinstance(baseline, ROOT.TList): 
        for baselineObj in baseline:
            objects = dict()
            for name,rlist in inputObjects.iteritems():
                obj = rlist.FindObject(baselineObj.GetName())
                if not obj:
                    print("{0} not found for analysis {1}".format(baselineObj.GetName(), name))
                    continue
                objects[name] = obj
            mylist[baselineObj.GetName()] = CompareObjects(baselineObj, objects)

    elif isinstance(baseline, ROOT.TH1):
        for name, hist in inputObjects.iteritems():
            ratio = hist.Clone("{0}_ratio".format(name))
            ratio.Divide(baseline)
            mylist[ratio.GetName()] = ratio

    return mylist

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='D meson jet analysis for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml', nargs='+',
                        help='YAML configuration file')

    args = parser.parse_args()

    configs = []
    for yamlFileName in args.yaml:
        f = open(yamlFileName, 'r')
        config = yaml.load(f)
        f.close()
        configs.append(config)

    main(configs)

    IPython.embed()
