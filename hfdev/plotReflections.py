#!/usr/bin/env python
# python script to plot reflection studies on the decay simulator results

import argparse
import ROOT
import IPython
import math

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

    histograms = []

    ReflDiffVsDPt = ROOT.TH2D("ReflDiffVsDPt", "ReflDiffVsDPt;#it{p}_{T,D} (GeV/#it{c});m_{refl} - m (GeV/#it{c}^{2})", ptMax - ptMin, ptMin, ptMax, 1500, -1.3, 3.2)
    globalList.append(ReflDiffVsDPt)
    histograms.append(ReflDiffVsDPt)

    ReflDiffVsDPtHighPt = ROOT.TH2D("ReflDiffVsDPtHighPt", "ReflDiffVsDPtHighPt;#it{p}_{T,D} (GeV/#it{c});#Delta_{approx}^{2} (GeV/#it{c}^{2})", ptMax - ptMin, ptMin, ptMax, 1500, -1.3, 3.2)
    globalList.append(ReflDiffVsDPtHighPt)
    histograms.append(ReflDiffVsDPtHighPt)

    ApproxInvMass = ROOT.TH2D("ApproxInvMass", "ApproxInvMass;#it{p}_{T,D} (GeV/#it{c});#it{m}_{K}^{2} + #it{m}_{#pi}^{2} + 2#it{p}_{K}#it{p}_{#pi}(1 - cos(#theta)) (GeV/#it{c}^{2})^{2}", ptMax - ptMin, ptMin, ptMax, 400, 1.0, 4.0)
    globalList.append(ApproxInvMass)
    histograms.append(ApproxInvMass)

    DaughtersPt = ROOT.TH2D("DaughtersPt", "DaughtersPt;#it{p}_{K} (GeV/#it{c});#it{p}_{#pi} (GeV/#it{c})", int(ptMax * 1.5), 0, ptMax * 1.5, int(ptMax * 1.5), 0, ptMax * 1.5)
    globalList.append(DaughtersPt)
    histograms.append(DaughtersPt)

    ApproxErrorInvMass = ROOT.TH1D("ApproxErrorInvMass", "ApproxErrorInvMass;2(#sqrt{(#it{p}_{K}^{2} + #it{m}_{K}^{2})(#it{p}_{#pi}^{2} + #it{m}_{#pi}^{2})} - #it{p}_{K}#it{p}_{#pi}) (GeV/#it{c}^{2})^{2}", 1600, 0, 4.0)
    globalList.append(ApproxErrorInvMass)
    histograms.append(ApproxErrorInvMass)

    ApproxErrorInvMassVsPtRatio = ROOT.TH2D("ApproxErrorInvMassVsPtRatio", "ApproxErrorInvMassVsPtRatio;#it{p}_{K}#it{m_{#pi}} / #it{p}_{#pi}#it{m_{K}};2(#sqrt{(#it{p}_{K}^{2} + #it{m}_{K}^{2})(#it{p}_{#pi}^{2} + #it{m}_{#pi}^{2})} - #it{p}_{K}#it{p}_{#pi}) (GeV/#it{c}^{2})^{2}", 2000, 0, 100, 2000, 0, 4.0)
    globalList.append(ApproxErrorInvMassVsPtRatio)
    histograms.append(ApproxErrorInvMassVsPtRatio)

    ApproxErrorInvMassVsDP = ROOT.TH2D("ApproxErrorInvMassVsDP", "ApproxErrorInvMassVsDP;#it{p}_{D} (GeV/#it{c});2(#sqrt{(#it{p}_{K}^{2} + #it{m}_{K}^{2})(#it{p}_{#pi}^{2} + #it{m}_{#pi}^{2})} - #it{p}_{K}#it{p}_{#pi}) (GeV/#it{c}^{2})^{2}", ptMax - ptMin, ptMin, ptMax, 2000, 0, 4.0)
    globalList.append(ApproxErrorInvMassVsDP)
    histograms.append(ApproxErrorInvMassVsDP)

    D0mass = 1.86484
    piMass = 0.139570
    Kmass = 0.493677

    minE = 100

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
        theta = dmeson.Daughter1.Vect().Angle(dmeson.Daughter2.Vect())
        approxMass = mass1 ** 2 + mass2 ** 2 + 2 * dmeson.Daughter1.P() * dmeson.Daughter2.P() * (1 - math.cos(theta))
        ApproxInvMass.Fill(pt, approxMass)
        approxE = D0mass ** 2 - approxMass
        if approxE < minE:
            minE = approxE
        ApproxErrorInvMass.Fill(approxE)
        DaughtersPt.Fill(dmeson.Daughter2.P(), dmeson.Daughter1.P())
        ApproxErrorInvMassVsPtRatio.Fill(dmeson.Daughter2.P() * dmeson.Daughter1.M() / dmeson.Daughter2.M() / dmeson.Daughter1.P(), approxE)
        ApproxErrorInvMassVsDP.Fill(dmeson.Mother.P(), approxE)

    print("Minimum approx error {}".format(minE))

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
    c.SaveAs("{}.pdf".format(c.GetName()))

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
    c.SaveAs("{}.pdf".format(c.GetName()))

    c = ROOT.TCanvas("ApproxInvMass", "ApproxInvMass")
    globalList.append(c)
    c.cd()
    ApproxInvMass.Draw("colz")
    line1 = ROOT.TLine(ptMin, D0mass ** 2, ptMax, D0mass ** 2)
    globalList.append(line1)
    line1.SetLineColor(ROOT.kRed)
    line1.SetLineWidth(2)
    line1.Draw()
    c.SaveAs("{}.pdf".format(c.GetName()))

    c = ROOT.TCanvas("DaughtersPt", "DaughtersPt")
    globalList.append(c)
    c.cd()
    DaughtersPt.Draw("colz")
    c.SetLogz()
    line1 = ROOT.TLine(0, 0, ptMax * 1.5, ptMax * 1.5 / piMass * Kmass)
    globalList.append(line1)
    line1.SetLineColor(ROOT.kRed)
    line1.SetLineWidth(2)
    line1.Draw()
    c.SaveAs("{}.pdf".format(c.GetName()))

    c = ROOT.TCanvas("ApproxErrorInvMass", "ApproxErrorInvMass")
    globalList.append(c)
    c.cd()
    ApproxErrorInvMass.Draw()
    c.SaveAs("{}.pdf".format(c.GetName()))

    c = ROOT.TCanvas("ApproxErrorInvMassVsPtRatio", "ApproxErrorInvMassVsPtRatio")
    globalList.append(c)
    c.cd()
    ApproxErrorInvMassVsPtRatio.Draw("colz")
    ApproxErrorInvMassVsPtRatio.GetXaxis().SetRangeUser(0, 10)
    ApproxErrorInvMassVsPtRatio.GetYaxis().SetRangeUser(0, 1)
    line1 = ROOT.TLine(0, 2 * piMass * Kmass, 10, 2 * piMass * Kmass)
    globalList.append(line1)
    line1.SetLineColor(ROOT.kRed)
    line1.SetLineWidth(2)
    line1.Draw()
    c.SaveAs("{}.pdf".format(c.GetName()))

    c = ROOT.TCanvas("ApproxErrorInvMassVsDP", "ApproxErrorInvMassVsDP")
    globalList.append(c)
    c.cd()
    ApproxErrorInvMassVsDP.Draw("colz")
    c.SaveAs("{}.pdf".format(c.GetName()))

    fname_out = fname.replace(".root", "_hists.root")
    file_out = ROOT.TFile(fname_out, "recreate")
    file_out.cd()
    for h in histograms:
        h.Write()
    file_out.Close()

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
