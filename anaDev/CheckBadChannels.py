#!/usr/local/bin/python
#python script to test Jet D meson analysis

import argparse
import ROOT
import IPython

globalList = []

def GetFEEBadChannels(run, geom):
    contBC = ROOT.AliOADBContainer("")
    fbad = ROOT.TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root","read")
    if not fbad or fbad.IsZombie():
        print("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root was not found")
        exit(1)
    else:
        fbad.Close()
    
    contBC.InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root","AliEMCALBadChannels")
  
    arrayBC = contBC.GetObject(run);
    if not arrayBC:
        print("No external hot channel set for run number: {0}".format(runBC))
        exit(2)
    
    BClist = []
    nSupMod = geom.GetNumberOfSuperModules()
    for i in range(0,nSupMod):
        h = arrayBC.FindObject("EMCALBadChannelMap_Mod{0}".format(i));

        if not h: 
            print("Can not get EMCALBadChannelMap_Mod{0}".format(i))
        else:
            BClist.append(h)   
            
    return BClist

def GetListOfFEEBadChannels(FEEbc, geom):
    sm = 0
    bclist = []
    for h in FEEbc:
        for x in range(1, h.GetNbinsX()):
            for y in range(1, h.GetNbinsY()):
                status = h.GetBinContent(x,y)
                if status != 0:
                    absId = geom.GetAbsCellIdFromCellIndexes(sm, y-1, x-1)
                    bclist.append(absId)
        sm += 1
    return bclist

def GetListOfTRUBadChannels(TRUbadChannelsFile):
    bclist = []
    f = open(TRUbadChannelsFile, 'r')
    for line in f:
        bclist.append(int(line))
    return bclist

def DrawFEEBadChannels(FEEbcList, geom):
    h = ROOT.TH2C("FEEBadChannels", "FEEBadChannels;eta;phi", 96, 0, 96, 240, 0, 240)
    for absId in FEEbcList:
        nSupMod = ROOT.Long()
        nMod = ROOT.Long()
        mEta = ROOT.Long()
        mPhi = ROOT.Long()
        smeta = ROOT.Long()
        smphi = ROOT.Long()
        geom.GetCellIndex(absId, nSupMod, nMod, mPhi, mEta)
        geom.GetCellPhiEtaIndexInSModule(nSupMod, nMod, mPhi, mEta, smphi, smeta)
        eta = smeta + (nSupMod % 2) * 48
        phi = smphi + (nSupMod / 2) * 24
        if nSupMod == 13 or nSupMod == 15 or nSupMod == 17:
            eta += 16
        h.SetBinContent(eta, phi, 1)
    
    canvas = ROOT.TCanvas("FEEBadChannels","FEEBadChannels", 600, 1000)
    canvas.cd()
    h.Draw("colz")
    globalList.append(h)
    globalList.append(canvas)
    return canvas
    
def DrawTRUBadChannels(TRUbcList, geom):
    h = ROOT.TH2C("TRUBadChannels", "TRUBadChannels;eta;phi", 48, 0, 48, 120, 0, 120)
    for absId in TRUbcList:
        eta = ROOT.Long()
        phi = ROOT.Long()
        geom.GetPositionInEMCALFromAbsFastORIndex(absId, eta, phi)
        h.SetBinContent(eta, phi, 1)
    
    canvas = ROOT.TCanvas("TRUBadChannels","TRUBadChannels", 600, 1000)
    canvas.cd()
    h.Draw("colz")
    globalList.append(h)
    globalList.append(canvas)
    return canvas
        
def CheckBadChannels(TRUbadChannelsFile, run):
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    print "Getting geometry for run {0}".format(run)
    geom = ROOT.AliEMCALGeometry.GetInstanceFromRunNumber(run)
    FEEbc = GetFEEBadChannels(run, geom)
    FEEbcList = GetListOfFEEBadChannels(FEEbc, geom)
    #TRUbcList = GetListOfTRUBadChannels(TRUbadChannelsFile)
    canvas = DrawFEEBadChannels(FEEbcList, geom)
    canvas.SaveAs("FEEBadChannels_{0}.pdf".format(run))
    #canvas = DrawTRUBadChannels(TRUbcList, geom)
    #canvas.SaveAs("TRUBadChannels_{0}.pdf".format(run))
    
if __name__ == '__main__':
    # CheckBadChannels.py executed as script
    
    parser = argparse.ArgumentParser(description='Check bad channels.')
    parser.add_argument('-r', '--run',
                        type=int,
                        help='Run (e.g. 237673)')
    parser.add_argument('--TRU',
                        help='TRU bad channel map file')
    args = parser.parse_args()
    
    CheckBadChannels(args.TRU, args.run)
    
    IPython.embed()
