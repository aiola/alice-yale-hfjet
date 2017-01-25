#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import argparse
import math

globalList = []

def PlotReflections(config, var, reflFit):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    fname = "{input_path}/{train}/{ana}/reflTemp/{ana}_{var}_fitted_{fit}.root".format(input_path=config["input_path"], train=config["train"], ana=config["name"], var=var, fit=reflFit)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    signHist = dict()
    reflHist = dict()
    reflFitHist = dict()

    for k in file.GetListOfKeys():
        h = file.Get(k.GetName())
        if "histSgn" in h.GetName():
            i = int(h.GetName()[-1])
            signHist[i] = h
        if "histRflFitted" in h.GetName():
            i = int(h.GetName()[-1])
            reflFitHist[i] = h
        if "histoRfl" in h.GetName():
            i = int(h.GetName()[-1])
            reflHist[i] = h

    globalList.extend(signHist)
    globalList.extend(reflHist)

    n = len(reflHist)
    ncols = int(math.ceil((math.sqrt(n))))
    nrows = int(math.floor((n / ncols)))

    cname = "ReflectionTemplates_{var}_{fit}".format(var=var, fit=reflFit)
    canvas = ROOT.TCanvas(cname, cname, 1200, 800)
    globalList.append(canvas)
    canvas.Divide(ncols, nrows)

    for i, (hSig, hRef, hRefFit) in enumerate(zip(signHist.itervalues(), reflHist.itervalues(), reflFitHist.itervalues())):
        pad = canvas.cd(i + 1)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.22)
        pad.SetRightMargin(0.02)
        pad.SetTopMargin(0.13)
        pad.SetBottomMargin(0.15)
        h = hRef.DrawCopy("axis")
        h.GetYaxis().SetRangeUser(0, 0.13)
        h.GetYaxis().SetTitle("arb. units")
        h.GetXaxis().SetTitleFont(43)
        h.GetXaxis().SetTitleOffset(2.6)
        h.GetXaxis().SetTitleSize(19)
        h.GetXaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelOffset(0.009)
        h.GetXaxis().SetLabelSize(18)
        h.GetYaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleOffset(4.5)
        h.GetYaxis().SetTitleSize(19)
        h.GetYaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelOffset(0.009)
        h.GetYaxis().SetLabelSize(23)
        globalList.append(h)

        hRef_copy = hRef.DrawCopy("p0 x0 same")
        globalList.append(hRef_copy)
        hRef_copy.SetLineColor(ROOT.kBlue + 2)
        hRef_copy.SetMarkerColor(ROOT.kBlue + 2)
        hRef_copy.SetMarkerStyle(ROOT.kFullCircle)
        hRef_copy.SetMarkerSize(1.0)

#         hSig_copy = hSig.DrawCopy("p0 x0 same")
#         globalList.append(hSig_copy)
#         hSig_copy.SetLineColor(ROOT.kGray)
#         hSig_copy.SetMarkerColor(ROOT.kGray)
#         hSig_copy.SetMarkerStyle(ROOT.kOpenCircle)
#         hSig_copy.SetMarkerSize(1.2)

        hRefFit_copy = hRefFit.DrawCopy("p0 x0 same")
        globalList.append(hRefFit_copy)
        hRefFit_copy.SetLineColor(ROOT.kGray)
        hRefFit_copy.SetMarkerColor(ROOT.kGray)
        hRefFit_copy.SetMarkerStyle(ROOT.kOpenCircle)
        hRefFit_copy.SetMarkerSize(1.2)

        refIntErr = ROOT.Double(0)
        sigIntErr = ROOT.Double(0)
        refInt = hRef.IntegralAndError(1, hRef.GetNbinsX(), refIntErr)
        sigInt = hSig.IntegralAndError(1, hSig.GetNbinsX(), sigIntErr)

        ros = refInt / sigInt
        rosErr = math.sqrt((refIntErr / refInt) ** 2 + (sigIntErr / sigInt) ** 2) * ros
        # rosAfterFit = hRefFit.Integral() / hSig.Integral()
        rosPave = ROOT.TPaveText(0.61, 0.68, 0.96, 0.78, "NB NDC")
        rosPave.SetBorderSize(0)
        rosPave.SetFillStyle(0)
        rosPave.SetTextFont(43)
        rosPave.SetTextSize(18)
        rosPave.SetTextAlign(31)
        rosPave.AddText("R/S = {0:.3f} #pm {1:.3f}".format(ros, rosErr))
        # rosPave.AddText("R/S (after fit) = {0:.3f}".format(rosAfterFit))
        rosPave.Draw()
        globalList.append(rosPave)

        binTitle = hRef.GetTitle()
        binTitle = binTitle[binTitle.index(":") + 2:]
        htitle = ROOT.TPaveText(0.15, 0.88, 0.95, 0.97, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(20)
        htitle.SetTextAlign(21)
        htitle.AddText(binTitle)
        htitle.Draw()
        globalList.append(htitle)

    canvas.cd(1)
    leg = ROOT.TLegend(0.25, 0.65, 0.58, 0.81, "", "NB NDC")
    globalList.append(leg)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(18)
    leg.SetTextAlign(13)
    leg.AddEntry(hRef_copy, "Reflections", "pe")
    leg.AddEntry(hRef_copy.GetListOfFunctions()[0], "Fit", "l")
    leg.AddEntry(hRefFit_copy, "Histogram from fit", "pe")
    leg.Draw()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot reflections.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    PlotReflections(config, "DPt", "DoubleGaus")
    PlotReflections(config, "JetPt", "DoubleGaus")

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{input_path}/{train}/{ana}/reflTemp/{name}.pdf".format(input_path=config["input_path"], train=config["train"], ana=config["name"], name=obj.GetName()))

    IPython.embed()
