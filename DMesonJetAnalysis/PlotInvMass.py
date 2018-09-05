#!/usr/bin/env python

import math
import argparse
import subprocess
import yaml
import IPython
import ROOT
import DMesonJetUtils

globalList = []

def DrawFitResults(mass_fitter, w):
    if mass_fitter is None or not mass_fitter.FitSuccessfull():
        return

    mass_fitter.Draw("same")

    if w:
        chi2Text = mass_fitter.GetChisquareWString().Data()
        signifText = mass_fitter.GetSignificanceWString().Data()
    else:
        chi2Text = mass_fitter.GetChisquareString().Data()
        signifText = mass_fitter.GetSignificanceString().Data()

    paveSig = ROOT.TPaveText(0.165, 0.732, 0.490, 0.92, "NB NDC")
    globalList.append(paveSig)
    paveSig.SetBorderSize(0)
    paveSig.SetFillStyle(0)
    paveSig.SetTextFont(43)
    paveSig.SetTextSize(14)
    paveSig.SetTextAlign(13)
    paveSig.AddText("{0}, {1}".format(mass_fitter.GetSignalString().Data(),
                                        mass_fitter.GetBackgroundString().Data()))
    paveSig.AddText("{0}, {1}".format(mass_fitter.GetSignalOverBackgroundString().Data(),
                                        signifText))
    paveSig.AddText(chi2Text)
    paveSig.Draw()

    paveFit = ROOT.TPaveText(0.47, 0.54, 0.955, 0.80, "NB NDC")
    globalList.append(paveFit)
    paveFit.SetBorderSize(0)
    paveFit.SetFillStyle(0)
    paveFit.SetTextFont(43)
    paveFit.SetTextSize(14)
    paveFit.SetTextAlign(33)

    paveFit.AddText(mass_fitter.GetSignalMeanString().Data())
    paveFit.AddText(mass_fitter.GetSignalWidthString().Data())
    paveFit.AddText(mass_fitter.GetBkgPar1String().Data())
    paveFit.AddText(mass_fitter.GetTotalEntriesString().Data())
    paveFit.Draw()

    if mass_fitter.GetReflOverSign() > 0:
        paveRos = ROOT.TPaveText(0.165, 0.51, 0.490, 0.77, "NB NDC")
        globalList.append(paveRos)
        paveRos.SetBorderSize(0)
        paveRos.SetFillStyle(0)
        paveRos.SetTextFont(43)
        paveRos.SetTextSize(14)
        paveRos.SetTextAlign(23)

        paveRos.AddText(mass_fitter.GetReflOverSignString().Data())
        paveRos.Draw()

def PlotInvMassDistr(pad, inv_mass_histo, inv_mass_fitter, title_offset, w):
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.05)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.13)
    h = inv_mass_histo.DrawCopy()
    globalList.append(h)
    if w:
        h.SetMaximum(h.GetMaximum() * 2.2)
    else:
        h.SetMaximum(h.GetMaximum() * 2.0)
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(19)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetTitleOffset(title_offset)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(19)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleOffset(title_offset+0.3)
    htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(18)
    title = h.GetTitle()
    title = title[title.index(":")+2:]
    htitle.AddText(title)
    htitle.Draw()
    globalList.append(htitle)
    DrawFitResults(inv_mass_fitter, w)

def InvMassPlot(config, var, kincuts, index_list):
    if "efficiency" in config["name"]:
        w = True
    else:
        w = False

    file_name = "{}/{}/{}.root".format(config["input_path"], config["train"], config["name"])
    root_file = ROOT.TFile(file_name)

    if not root_file or root_file.IsZombie():
        print("Could not open file '{}'".format(file_name))
        exit(1)
    
    jet_def = "Charged_R040_pt_scheme"
    list_name = "AnyINT_D0_D0toKpiCuts/{jet_def}/D0_D0toKpiCuts_{jet_def}_{var}Bins_{kincuts}".format(var=var, kincuts=kincuts, jet_def=jet_def)
    root_list = DMesonJetUtils.GetObject(root_file, list_name)

    if not root_list:
        print("Could not find list '{}' in file '{}'".format(list_name, file_name))
        exit(1)

    n = len(index_list)
    ncols = int(math.ceil((math.sqrt(n))))
    nrows = int(math.floor((n / ncols)))
    if ncols*nrows < n:
        ncols += 1
    
    if nrows > 1:
        title_offset = 2
    else:
        title_offset = 1

    cname = "InvMass_{}_{}_Thesis".format(var, kincuts)
    canvas = ROOT.TCanvas(cname, cname, ncols*390, nrows*390)
    globalList.append(canvas)
    canvas.Divide(ncols, nrows)

    for ibin, bin in enumerate(index_list):
        inv_mass_histo_name = "InvMass_AnyINT_D0_D0toKpiCuts_{var}_{bin}".format(var=var, bin=bin)
        inv_mass_histo = root_list.FindObject(inv_mass_histo_name)
        if not inv_mass_histo:
            print("Could not find histo '{}'".format(inv_mass_histo_name))
            exit(1)

        inv_mass_fitter_name = "InvMass_AnyINT_D0_D0toKpiCuts_{var}_{bin}_fitter_DoubleGaus".format(var=var, bin=bin)
        inv_mass_fitter = root_list.FindObject(inv_mass_fitter_name)
        if not inv_mass_fitter:
            print("Could not find fitter '{}'".format(inv_mass_fitter_name))
            exit(1)
        pad = canvas.cd(ibin+1)
        PlotInvMassDistr(pad, inv_mass_histo, inv_mass_fitter, title_offset, w)

def main(config, var, kincuts, index_list):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    InvMassPlot(config, var, kincuts, index_list)

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            file_name = "{}/{}/{}/{}".format(config["input_path"], config["train"], config["name"], obj.GetName())
            obj.SaveAs("{}.pdf".format(file_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot reflections.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--plot', nargs='+', type=str)
    parser.add_argument('--var', type=str, default="DPt")
    parser.add_argument('--kincuts', type=str, default="JetPt_5_30")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.var, args.kincuts, args.plot)

    IPython.embed()
