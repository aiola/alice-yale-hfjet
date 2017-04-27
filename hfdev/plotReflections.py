#!/usr/bin/env python
# python script to plot reflection studies on the decay simulator results

import argparse
import ROOT
import IPython

globalList = []

def main(fname, ptMin, ptMax):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file '{}'".format(fname))
        exit(1)

    tree = file.Get("DecaySimulation")
    if not tree:
        print("Could not get tree!")
        exit(1)

    ReflDiffVsDPt = ROOT.TH2D("ReflDiffVsDPt", "ReflDiffVsDPt;#it{p}_{T,D} (GeV/#it{c});m_{refl} - m (GeV/#it{c}^{2})", ptMax - ptMin, ptMin, ptMax, 1500, -1.3, 3.2)
    globalList.append(ReflDiffVsDPt)

    ReflDiffVsDPtHighPt = ROOT.TH2D("ReflDiffVsDPtHighPt", "ReflDiffVsDPtHighPt;#it{p}_{T,D} (GeV/#it{c});#Delta_{approx}^{2} (GeV/#it{c}^{2})", ptMax - ptMin, ptMin, ptMax, 1500, -1.3, 3.2)
    globalList.append(ReflDiffVsDPtHighPt)

    for dmeson in tree:
        mass1 = dmeson.Daughter1.M()
        mass2 = dmeson.Daughter2.M()
        reflDaugh1 = ROOT.TLorentzVector(dmeson.Daughter1)
        reflDaugh1.SetPtEtaPhiM(dmeson.Daughter1.Pt(), dmeson.Daughter1.Eta(), dmeson.Daughter1.Phi(), mass2)
        reflDaugh2 = ROOT.TLorentzVector(dmeson.Daughter2)
        reflDaugh2.SetPtEtaPhiM(dmeson.Daughter2.Pt(), dmeson.Daughter2.Eta(), dmeson.Daughter2.Phi(), mass1)
        refl = reflDaugh1 + reflDaugh2
        massRefl = refl.M()
        mass = dmeson.Mother.M()
        pt = dmeson.Mother.Pt()
        ReflDiffVsDPt.Fill(pt, massRefl - mass)
        ReflDiffVsDPtHighPt.Fill(pt, (mass1 * mass1 - mass2 * mass2) * (dmeson.Daughter2.P() / dmeson.Daughter1.P() - dmeson.Daughter1.P() / dmeson.Daughter2.P()))

    D0mass = 1.86484

    c = ROOT.TCanvas("ReflDiffVsDPt", "ReflDiffVsDPt")
    globalList.append(c)
    c.cd()
    ReflDiffVsDPt.Draw("colz")
    line1 = ROOT.TLine(ptMin, 1.715 - D0mass, ptMax, 1.715 - D0mass)
    globalList.append(line1)
    line1.SetLineColor(ROOT.kRed)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(ptMin, 2.015 - D0mass, ptMax, 2.015 - D0mass)
    globalList.append(line2)
    line2.SetLineColor(ROOT.kRed)
    line2.SetLineWidth(2)
    line2.Draw()

    c = ROOT.TCanvas("ReflDiffVsDPtHighPt", "ReflDiffVsDPtHighPt")
    globalList.append(c)
    c.cd()
    ReflDiffVsDPtHighPt.Draw("colz")
    line1 = ROOT.TLine(ptMin, 1.715 - D0mass, ptMax, 1.715 - D0mass)
    globalList.append(line1)
    line1.SetLineColor(ROOT.kRed)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(ptMin, 2.015 - D0mass, ptMax, 2.015 - D0mass)
    globalList.append(line2)
    line2.SetLineColor(ROOT.kRed)
    line2.SetLineWidth(2)
    line2.Draw()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reflection studies on the decay simulator results.')
    parser.add_argument('filename', metavar='simulation.root',
                        help='ROOT file containing the results of the decay simulator')
    parser.add_argument('--pt-min', metavar='PTMIN',
                        default=0, type=int,
                        help='Minimum pt of the simulation')
    parser.add_argument('--pt-max', metavar='PTMAX',
                        default=100, type=int,
                        help='Maximum pt of the simulation')
    args = parser.parse_args()

    main(args.filename, args.pt_min, args.pt_max)

    IPython.embed()
