#!/usr/bin/env python

import argparse
import subprocess
import yaml
import IPython
import ROOT
import DMesonJetUtils

globalList = []

def InvMassFitResultsPlot(config, var, kincuts):
    file_name = "{}/{}/{}.root".format(config["input_path"], config["train"], config["name"])
    root_file = ROOT.TFile(file_name)

    if not root_file or root_file.IsZombie():
        print("Could not open file '{}'".format(file_name))
        exit(1)
    
    jet_def = "Charged_R040_pt_scheme"
    list_name = "AnyINT_D0_D0toKpiCuts/{jet_def}/D0_D0toKpiCuts_{jet_def}_{var}Spectrum_{kincuts}_InvMassFit_DoubleGaus".format(var=var, kincuts=kincuts, jet_def=jet_def)
    root_list = DMesonJetUtils.GetObject(root_file, list_name)

    if not root_list:
        print("Could not find list '{}' in file '{}'".format(list_name, file_name))
        exit(1)

    signal_name = "D0_D0toKpiCuts_{jet_def}_{var}Spectrum_{kincuts}_InvMassFit_DoubleGaus".format(var=var, kincuts=kincuts, jet_def=jet_def)
    bkg_name = "{}_Bkg".format(signal_name)
    significance_name = "{}_Significance".format(signal_name)
    mass_name = "{}_Mass".format(signal_name)
    mass_width_name = "{}_MassWidth".format(signal_name)

    signal_histo = root_list.FindObject(signal_name)
    if not signal_histo:
        print("Could not find histo '{}'".format(signal_name))
        exit(1)
    bkg_histo = root_list.FindObject(bkg_name)
    if not bkg_histo:
        print("Could not find histo '{}'".format(bkg_name))
        exit(1)
    significance_histo = root_list.FindObject(significance_name)
    if not significance_histo:
        print("Could not find histo '{}'".format(significance_name))
        exit(1)
    mass_histo = root_list.FindObject(mass_name)
    if not mass_histo:
        print("Could not find histo '{}'".format(mass_name))
        exit(1)
    mass_width_histo = root_list.FindObject(mass_width_name)
    if not mass_width_histo:
        print("Could not find histo '{}'".format(mass_width_name))
        exit(1)

    PlotSignal(signal_histo, bkg_histo, significance_histo, var, kincuts)
    PlotMassPars(mass_histo, mass_width_histo, var, kincuts)

def PlotMassPars(mass_histo, mass_width_histo, var, kincuts):
    cname = "MassPars_{}_{}_Thesis".format(var, kincuts)
    canvas = ROOT.TCanvas(cname, cname, 500, 500)
    globalList.append(canvas)

    canvas = ROOT.TCanvas(cname, cname, 500, 500)
    globalList.append(canvas)
    canvas.Divide(1,2)

    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.50, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.15)
    padMain.SetRightMargin(0.05)
    padMain.SetTopMargin(0.05)
    padMain.SetTicks(1, 1)

    padBottom = canvas.cd(2)
    padBottom.SetPad(0, 0., 1, 0.50)
    padBottom.SetTopMargin(0)
    padBottom.SetBottomMargin(0.30)
    padBottom.SetLeftMargin(0.15)
    padBottom.SetRightMargin(0.05)
    padBottom.SetTicks(1, 1)

    padMain.cd()
    mass_histo_copy = mass_histo.DrawCopy()
    globalList.append(mass_histo_copy)
    mass_histo_copy.GetYaxis().SetRangeUser(1.852, 1.876)
    mass_histo_copy.GetYaxis().SetTitleFont(43)
    mass_histo_copy.GetYaxis().SetTitleSize(18)
    mass_histo_copy.GetYaxis().SetTitleOffset(1.7)
    mass_histo_copy.GetYaxis().SetLabelFont(43)
    mass_histo_copy.GetYaxis().SetLabelSize(16)
    mass_histo_copy.GetYaxis().SetLabelOffset(0.007)
    mass_histo_copy.SetMarkerColor(ROOT.kGreen + 2)
    mass_histo_copy.SetMarkerStyle(ROOT.kOpenCircle)
    mass_histo_copy.SetMarkerSize(1.0)
    mass_histo_copy.SetLineColor(ROOT.kGreen + 2)

    d0mass = 1.86484
    d0mass_line = ROOT.TLine(mass_histo_copy.GetXaxis().GetXmin(), d0mass, mass_histo_copy.GetXaxis().GetXmax(), d0mass)
    globalList.append(d0mass_line)
    d0mass_line.SetLineWidth(2)
    d0mass_line.SetLineColor(ROOT.kBlack)
    d0mass_line.SetLineStyle(2)
    d0mass_line.Draw()

    padBottom.cd()
    mass_width_histo_copy = mass_width_histo.DrawCopy()
    globalList.append(mass_width_histo_copy)
    mass_width_histo_copy.GetYaxis().SetRangeUser(0, 0.029)
    mass_width_histo_copy.GetYaxis().SetNdivisions(509)
    mass_width_histo_copy.GetYaxis().SetTitleFont(43)
    mass_width_histo_copy.GetYaxis().SetTitleSize(18)
    mass_width_histo_copy.GetYaxis().SetTitleOffset(1.7)
    mass_width_histo_copy.GetYaxis().SetLabelFont(43)
    mass_width_histo_copy.GetYaxis().SetLabelSize(16)
    mass_width_histo_copy.GetYaxis().SetLabelOffset(0.007)
    mass_width_histo_copy.GetXaxis().SetTitleFont(43)
    mass_width_histo_copy.GetXaxis().SetTitleSize(18)
    mass_width_histo_copy.GetXaxis().SetTitleOffset(3.0)
    mass_width_histo_copy.GetXaxis().SetLabelFont(43)
    mass_width_histo_copy.GetXaxis().SetLabelSize(16)
    mass_width_histo_copy.GetXaxis().SetLabelOffset(0.007)
    mass_width_histo_copy.SetMarkerColor(ROOT.kRed + 2)
    mass_width_histo_copy.SetMarkerStyle(ROOT.kOpenSquare)
    mass_width_histo_copy.SetMarkerSize(1.0)
    mass_width_histo_copy.SetLineColor(ROOT.kRed + 2)

def PlotSignal(signal_histo, _, significance_histo, var, kincuts):
    cname = "Signal_{}_{}_Thesis".format(var, kincuts)
    canvas = ROOT.TCanvas(cname, cname, 500, 500)
    globalList.append(canvas)
    canvas.Divide(1,2)

    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.30, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.14)
    padMain.SetRightMargin(0.05)
    padMain.SetTopMargin(0.05)
    padMain.SetTicks(1, 1)

    padBottom = canvas.cd(2)
    padBottom.SetPad(0, 0., 1, 0.30)
    padBottom.SetTopMargin(0)
    padBottom.SetBottomMargin(0.30)
    padBottom.SetLeftMargin(0.14)
    padBottom.SetRightMargin(0.05)
    padBottom.SetTicks(1, 1)

    padMain.cd()
    signal_histo_copy = signal_histo.DrawCopy()
    globalList.append(signal_histo_copy)
    signal_histo_copy.SetMinimum(0.001)
    signal_histo_copy.GetYaxis().SetTitleFont(43)
    signal_histo_copy.GetYaxis().SetTitleSize(18)
    signal_histo_copy.GetYaxis().SetTitleOffset(1.5)
    signal_histo_copy.GetYaxis().SetLabelFont(43)
    signal_histo_copy.GetYaxis().SetLabelSize(16)
    signal_histo_copy.GetYaxis().SetLabelOffset(0.007)
    signal_histo_copy.SetMarkerColor(ROOT.kBlue + 2)
    signal_histo_copy.SetMarkerStyle(ROOT.kFullCircle)
    signal_histo_copy.SetMarkerSize(1.0)
    signal_histo_copy.SetLineColor(ROOT.kBlue + 2)

    # bkg_histo_copy = bkg_histo.DrawCopy("same")
    # globalList.append(bkg_histo_copy)
    # bkg_histo_copy.SetMarkerColor(ROOT.kRed + 2)
    # bkg_histo_copy.SetMarkerStyle(ROOT.kOpenCircle)
    # bkg_histo_copy.SetMarkerSize(1.0)
    # bkg_histo_copy.SetLineColor(ROOT.kRed + 2)
    # bkg_histo_copy.SetLineWidth(2)

    padBottom.cd()
    significance_histo_copy = significance_histo.DrawCopy()
    globalList.append(significance_histo_copy)
    significance_histo_copy.GetYaxis().SetRangeUser(0, 19.9)
    significance_histo_copy.GetYaxis().SetNdivisions(509)
    significance_histo_copy.GetYaxis().SetTitleFont(43)
    significance_histo_copy.GetYaxis().SetTitleSize(18)
    significance_histo_copy.GetYaxis().SetTitleOffset(1.5)
    significance_histo_copy.GetYaxis().SetLabelFont(43)
    significance_histo_copy.GetYaxis().SetLabelSize(16)
    significance_histo_copy.GetYaxis().SetLabelOffset(0.007)
    significance_histo_copy.GetXaxis().SetTitleFont(43)
    significance_histo_copy.GetXaxis().SetTitleSize(18)
    significance_histo_copy.GetXaxis().SetTitleOffset(3.0)
    significance_histo_copy.GetXaxis().SetLabelFont(43)
    significance_histo_copy.GetXaxis().SetLabelSize(16)
    significance_histo_copy.GetXaxis().SetLabelOffset(0.007)
    significance_histo_copy.SetFillColor(ROOT.kCyan - 10)
    significance_histo_copy.SetLineColor(ROOT.kBlue + 2)
    significance_histo_copy.SetLineWidth(2)

def main(config, var, kincuts):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    InvMassFitResultsPlot(config, var, kincuts)

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            file_name = "{}/{}/{}/{}".format(config["input_path"], config["train"], config["name"], obj.GetName())
            obj.SaveAs("{}.pdf".format(file_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot reflections.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--var', type=str, default="DPt")
    parser.add_argument('--kincuts', type=str, default="JetPt_5_30")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.var, args.kincuts)

    IPython.embed()
