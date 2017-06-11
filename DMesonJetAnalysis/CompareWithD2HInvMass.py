#!/usr/local/bin/python

import ROOT
import math
import IPython
import argparse
import DMesonJetUtils

globalList = []

DPTBins = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 36, 9999];

def GetD2HInvMass(file, dpt_bins):
    # listname = "PWG3_D2H_D0InvMass_MC_wtopo_p4/coutputmassD0Mass_MC_wtopo_p40100"
    listname = "PWG3_D2H_D0InvMass/coutputmassD0Mass0100"
    hlist = DMesonJetUtils.GetObject(file, listname)
    if not hlist:
        print("Could not get list {}".format(listname))
        exit(1)
    hSgn = []
    hRfl = []
    for i in range(0, len(dpt_bins) - 1):
        h = hlist.FindObject("histSgn_{}".format(i))
        h.Rebin(6)
        hSgn.append(h)
        h = hlist.FindObject("histRfl_{}".format(i))
        h.Rebin(6)
        hRfl.append(h)
    return hSgn, hRfl

def ProjectMyTree(tree, dpt_bins, hname):
    histos = [ROOT.TH1F("{}_{}".format(hname, i), "{}_{};M (GeV/#it{{c}}^{{2}});arb. units".format(hname, i), 100, 1.6248, 2.2248) for i in range(0, len(dpt_bins) - 1)]
    for i, obj in enumerate(tree):
        if i % 1000 == 0:
            print("Event {}".format(i))
        invMass = obj.DmesonJet.fInvMass
        pt = obj.DmesonJet.fPt
        for iPtBin, ptLim in enumerate(dpt_bins[1:]):
            if ptLim > pt: break
        # print("pt = {}, bin = {}".format(pt, iPtBin))
        h = histos[iPtBin]
        h.Fill(invMass)
    for h in histos: h.Sumw2()
    return histos

def GetMyInvMass(file, dpt_bins):
    treename = "AliAnalysisTaskDmesonJets_D0_kSignalOnly"
    tree = DMesonJetUtils.GetObject(file, treename)
    hSgn = ProjectMyTree(tree, dpt_bins, "hSgn")

    treename = "AliAnalysisTaskDmesonJets_D0_WrongPID"
    tree = DMesonJetUtils.GetObject(file, treename)
    hRfl = ProjectMyTree(tree, dpt_bins, "hRfl")

    return hSgn, hRfl

def GetWeight(file):
    tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
    hTrials = tlist.FindObject("fHistTrialsVsPtHardNoSel")
    hXsec = tlist.FindObject("fHistXsectionVsPtHardNoSel")
    trials = hTrials.Integral()
    xsec = hXsec.GetMean(2)
    w = xsec / trials
    return w

def single_file_analysis(filename, dpt_bins, filename_d2h=None):
    print("Opening file {0}".format(filename))
    file = ROOT.TFile(filename)

    if filename_d2h:
        print("Opening file {0}".format(filename_d2h))
        file_d2h = ROOT.TFile(filename_d2h)
    else:
        file_d2h = file

    hSgn_D2H_histos, hRfl_D2H_histos = GetD2HInvMass(file_d2h, dpt_bins)
    hSgn_histos, hRfl_histos = GetMyInvMass(file, dpt_bins)

    globalList.extend(hSgn_D2H_histos)
    globalList.extend(hRfl_D2H_histos)
    globalList.extend(hSgn_histos)
    globalList.extend(hRfl_histos)

    for i, (minPt, maxPt, hSgn_D2H, hSgn, hRfl_D2H, hRfl) in enumerate(zip(dpt_bins[:-1], dpt_bins[1:], hSgn_D2H_histos, hSgn_histos, hRfl_D2H_histos, hRfl_histos)):
        cname = "cSgn_{}_{}".format(int(minPt * 10), int(maxPt * 10))
        cSgn = ROOT.TCanvas(cname, cname)
        globalList.append(cSgn)
        cSgn.cd()
        hSgn_D2H.Draw("hist")
        hSgn_D2H.SetMarkerSize(0.9)
        hSgn_D2H.SetMarkerColor(ROOT.kRed)
        hSgn_D2H.SetMarkerStyle(ROOT.kOpenCircle)
        hSgn_D2H.SetLineColor(ROOT.kRed)
        hSgn.Draw("same hist")
        hSgn.SetMarkerSize(0.9)
        hSgn.SetMarkerColor(ROOT.kBlue)
        hSgn.SetMarkerStyle(ROOT.kOpenSquare)
        hSgn.SetLineColor(ROOT.kBlue)

        cname = "cRfl_{}_{}".format(int(minPt * 10), int(maxPt * 10))
        cRfl = ROOT.TCanvas(cname, cname)
        globalList.append(cRfl)
        cRfl.cd()
        hRfl_D2H.Draw("hist")
        hRfl_D2H.SetMarkerSize(0.9)
        hRfl_D2H.SetMarkerColor(ROOT.kRed)
        hRfl_D2H.SetMarkerStyle(ROOT.kOpenCircle)
        hRfl_D2H.SetLineColor(ROOT.kRed)
        hRfl.Draw("same hist")
        hRfl.SetMarkerSize(0.9)
        hRfl.SetMarkerColor(ROOT.kBlue)
        hRfl.SetMarkerStyle(ROOT.kOpenSquare)
        hRfl.SetLineColor(ROOT.kBlue)

        sgn_D2H = hSgn_D2H.Integral()
        sgn = hSgn.Integral()
        rfl_D2H = hRfl_D2H.Integral()
        rfl = hRfl.Integral()

        if sgn_D2H > 0:
            ros_D2H = rfl_D2H / sgn_D2H
        else:
            ros_D2H = -1
        if sgn > 0:
            ros = rfl / sgn
        else:
            ros = -1

        print("Pt = [{},{}] R/S = {} - D2H R/S = {}".format(minPt, maxPt, ros, ros_D2H))

def pthard_analysis(path, periods, pt_hard_bins, dpt_bins):
    hSgn_D2H_histos = [ROOT.TH1F("hSgn_D2H_{}".format(i), "hSgn_D2H_{};M (GeV/#it{{c}}^{{2}});arb. units".format(i), 100, 1.6248, 2.2248) for i in range(0, len(dpt_bins) - 1)]
    hRfl_D2H_histos = [ROOT.TH1F("hRfl_D2H_{}".format(i), "hRfl_D2H_{};M (GeV/#it{{c}}^{{2}});arb. units".format(i), 100, 1.6248, 2.2248) for i in range(0, len(dpt_bins) - 1)]
    hSgn_histos = [ROOT.TH1F("hSgn_{}".format(i), "hSgn_{};M (GeV/#it{{c}}^{{2}});arb. units".format(i), 100, 1.6248, 2.2248) for i in range(0, len(dpt_bins) - 1)]
    hRfl_histos = [ROOT.TH1F("hRfl_{}".format(i), "hRfl_{};M (GeV/#it{{c}}^{{2}});arb. units".format(i), 100, 1.6248, 2.2248) for i in range(0, len(dpt_bins) - 1)]

    for pt_hard in pt_hard_bins:
        for period in periods:
            fileName = "{}/{}/{}/AnalysisResults.root".format(path, period, pt_hard)
            print("Opening file {0}".format(fileName))
            file = ROOT.TFile(fileName)

            hSgn_D2H_pthard_histos, hRfl_D2H_pthard_histos = GetD2HInvMass(file, dpt_bins)
            hSgn_pthard_histos, hRfl_pthard_histos = GetMyInvMass(file, dpt_bins)

            w = GetWeight(file)

            def addH(totHistos, partHistos):
                for tot, part in zip(totHistos, partHistos): tot.Add(part)

            addH(hSgn_D2H_histos, hSgn_D2H_pthard_histos)
            addH(hRfl_D2H_histos, hRfl_D2H_pthard_histos)
            addH(hSgn_histos, hSgn_pthard_histos)
            addH(hRfl_histos, hRfl_pthard_histos)

    globalList.extend(hSgn_D2H_histos)
    globalList.extend(hRfl_D2H_histos)
    globalList.extend(hSgn_histos)
    globalList.extend(hRfl_histos)

    for i, (minPt, maxPt, hSgn_D2H, hSgn, hRfl_D2H, hRfl) in enumerate(zip(dpt_bins[:-1], dpt_bins[1:], hSgn_D2H_histos, hSgn_histos, hRfl_D2H_histos, hRfl_histos)):
        cname = "cSgn_{}_{}".format(int(minPt * 10), int(maxPt * 10))
        cSgn = ROOT.TCanvas(cname, cname)
        globalList.append(cSgn)
        cSgn.cd()
        hSgn_D2H.Draw("hist")
        hSgn_D2H.SetMarkerSize(0.9)
        hSgn_D2H.SetMarkerColor(ROOT.kRed)
        hSgn_D2H.SetMarkerStyle(ROOT.kOpenCircle)
        hSgn_D2H.SetLineColor(ROOT.kRed)
        hSgn.Draw("same hist")
        hSgn.SetMarkerSize(0.9)
        hSgn.SetMarkerColor(ROOT.kBlue)
        hSgn.SetMarkerStyle(ROOT.kOpenSquare)
        hSgn.SetLineColor(ROOT.kBlue)

        cname = "cRfl_{}_{}".format(int(minPt * 10), int(maxPt * 10))
        cRfl = ROOT.TCanvas(cname, cname)
        globalList.append(cRfl)
        cRfl.cd()
        hRfl_D2H.Draw("hist")
        hRfl_D2H.SetMarkerSize(0.9)
        hRfl_D2H.SetMarkerColor(ROOT.kRed)
        hRfl_D2H.SetMarkerStyle(ROOT.kOpenCircle)
        hRfl_D2H.SetLineColor(ROOT.kRed)
        hRfl.Draw("same hist")
        hRfl.SetMarkerSize(0.9)
        hRfl.SetMarkerColor(ROOT.kBlue)
        hRfl.SetMarkerStyle(ROOT.kOpenSquare)
        hRfl.SetLineColor(ROOT.kBlue)

        sgn_D2H = hSgn_D2H.Integral()
        sgn = hSgn.Integral()
        rfl_D2H = hRfl_D2H.Integral()
        rfl = hRfl.Integral()

        ros_D2H = rfl_D2H / sgn_D2H
        ros = rfl / sgn

        print("Pt = [{},{}] R/S = {} - D2H R/S = {}".format(minPt, maxPt, ros, ros_D2H))

def ParseTrainNumbers(train):
    TrainNumberList = train.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range) - 1]) + 1))
    TrainNumbersLabel = "_".join([str(i) for i in TrainNumbers])
    return TrainNumbersLabel

def main(train, test, pthard, train_name, train_d2h, train_name_d2h):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if test:
        fileName = "../anaDev/AnalysisResults.root"
        single_file_analysis(fileName, DPTBins)
    else:
        TrainNumbersLabel = ParseTrainNumbers(train)

        if pthard:
            trainno = int(train)
            path = "/Volumes/DATA/ALICE/JetResults/{}_{}".format(train_name, TrainNumbersLabel)
            periods = ["LHC15i2b", "LHC15i2c", "LHC15i2d", "LHC15i2e"]
            pt_hard_bins = range(1, 9)
            pthard_analysis(path, periods, pt_hard_bins, DPTBins)
        else:
            fileName = "/Volumes/DATA/ALICE/JetResults/{}_{}/AnalysisResults.root".format(train_name, TrainNumbersLabel)
            fileName_d2h = None
            if train_d2h or train_name_d2h:
                if train_d2h:
                    TrainNumbersLabel_d2h = ParseTrainNumbers(train_d2h)
                else:
                    TrainNumbersLabel_d2h = TrainNumbersLabel
                if not train_name_d2h:
                    train_name_d2h = train_name
                fileName_d2h = "/Volumes/DATA/ALICE/JetResults/{}_{}/AnalysisResults.root".format(train_name_d2h, TrainNumbersLabel_d2h)
            if fileName_d2h == fileName: fileName_d2h = None

            single_file_analysis(fileName, DPTBins, fileName_d2h)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compare invariant mass with D2H.')
    parser.add_argument('--train',
                        help='Train number(s) (use : for a range and , for a list)')
    parser.add_argument('--test', action='store_const',
                        default=False, const=True,
                        help='Ignore train number to use ../anaDev/AnalysisResults.root')
    parser.add_argument('--pthard', action='store_const',
                        default=False, const=True,
                        help='Pt hard bin production')
    parser.add_argument('--train-name',
                        default="Jets_EMC_pp_MC",
                        help='Train name')
    parser.add_argument('--train-d2h',
                        default=None,
                        help='Train number(s) (use : for a range and , for a list)')
    parser.add_argument('--train-name-d2h',
                        default=None,
                        help='Train name')
    args = parser.parse_args()

    main(args.train, args.test, args.pthard, args.train_name, args.train_d2h, args.train_name_d2h)

    IPython.embed()
