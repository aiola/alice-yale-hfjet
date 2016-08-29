#!/usr/bin/env python
#python script to do some checks on the SB analysis 

import argparse
import yaml
import IPython
import ROOT

globalList = []

def main(config):
    file = OpenFile(config)
    signalOnly = LoadHistograms("D0_kSignalOnly_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    histos = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    subtracted = SubtractSignal(histos, signalOnly)
    CompareBackground(histos, subtracted)

def OpenFile(config):
    ROOT.TH1.AddDirectory(0)
    file = ROOT.TFile.Open("{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"]))
    return file

def LoadHistograms(lname, file):
    rlist = file.Get(lname)
    rSBlist = rlist.FindObject("SideBandAnalysis")
    result = dict()
    for h in rSBlist:
        result[h.GetName()] = h
    return result

def SubtractSignal(histos, signalOnly):
    result = dict()
    for hname, h in signalOnly.iteritems():
        if "SideBandWindow" in hname:
            continue
        hname2 = hname.replace("D0_kSignalOnly", "D0")
        hnewname = hname2 + "_Subtracted"
        hnew = histos[hname2].Clone(hnewname)
        hnew.SetTitle(hnewname)
        hnew.Add(h, -1)
        result[hnewname] = hnew
    return result

def CompareBackground(histos, subtracted):
    for hsub in subtracted.itervalues():
        c = ROOT.TCanvas(hsub.GetName(), hsub.GetTitle())
        c.cd()
        globalList.append(c)
        hsubcopy = hsub.DrawCopy()
        hsubcopy.SetMarkerStyle(ROOT.kFullCircle)
        hsubcopy.SetMarkerSize(0.9)
        hsubcopy.SetMarkerColor(ROOT.kBlue+2)
        hsubcopy.SetLineColor(ROOT.kBlue+2)
        globalList.append(hsubcopy)
        hname = hsub.GetName().replace("_Subtracted", "")
        hname = hname.replace("SignalWindow", "SideBandWindow")
        hcopy = histos[hname].DrawCopy("same")
        hcopy.SetMarkerStyle(ROOT.kFullSquare)
        hcopy.SetMarkerSize(0.9)
        hcopy.SetMarkerColor(ROOT.kRed+2)
        hcopy.SetLineColor(ROOT.kRed+2)
        if hcopy.GetMaximum() > hsubcopy.GetMaximum():
            hsubcopy.SetMaximum(hcopy.GetMaximum()*1.2)
        if hcopy.GetMinimum() < hsubcopy.GetMinimum():
            if hcopy.GetMinimum() > 0:
                hsubcopy.SetMinimum(0)
            else:
                hsubcopy.SetMinimum(hcopy.GetMinimum()*1.2)
        globalList.append(hcopy)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Side Band analysis.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)
    
    IPython.embed()