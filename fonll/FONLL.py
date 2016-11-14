#!/usr/bin/env python
#python script to plot FONLL calculations

import argparse
import IPython
import ROOT
import array
import math

globalList = []

class FONLL:
    def __init__(self, fonll_file):
        self.name = fonll_file
        self.n = 0
        self.values = dict()
        self.OpenFONLL(fonll_file)

    def OpenFONLL(self, fonll_file):
        with open(fonll_file, "r") as myfile:
            lines = myfile.read().splitlines()

        vars = None
        prevLine = None
        for line in lines:
            line = line.strip()
            #print(line)
            if line[0] == "#":
                prevLine = line
                continue
            if not vars:
                prevLine = "  {0}".format(prevLine[1:])
                vars = [prevLine[i:11+i].strip().strip("#") for i in range(0, len(prevLine), 11)]
                print("Variables taken from {0}".format(vars))
                for var in vars:
                    self.values[var] = []
            valuesStr = line.split(" ")
            values = [float(v) for v in valuesStr]
            self.n += 1
            for name,val in zip(vars, values):
                self.values[name].append(val)
            prevLine = line
    
    def GenerateGraph(self, x, y, yerrup, yerrdown):
        g = ROOT.TGraphAsymmErrors(self.n, array.array('d', self.values[x]), array.array('d', self.values[y]), ROOT.nullptr, ROOT.nullptr, 
                                    array.array('d', [v-down for down,v in zip(self.values[yerrdown], self.values[y])]),
                                    array.array('d', [up-v for up,v in zip(self.values[yerrup], self.values[y])]))
        return g

class MCGEN:
    def __init__(self, name, title, file_name, spectrum_name):
        self.file_name = file_name
        self.spectrum_name = spectrum_name
        self.name = name
        self.title = title
        self.xsec = 1e9/2
        self.LoadSpectrum()

    def LoadSpectrum(self):
        file = ROOT.TFile.Open(self.file_name)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(self.file_name))
            return
        spectrum = file.Get(self.spectrum_name)
        if not spectrum:
            print("Could not get histogram {0} from file {1}".format(self.spectrum_name, self.file_name))
            return
        self.spectrum = spectrum.Clone(self.name)
        self.spectrum.SetTitle(self.title)
        self.spectrum.Scale(self.xsec)

def MakeUniform(g,h):
    xval = []
    xerr = []
    yval = []
    yerr = []
    ibin = 1
    while ibin <= h.GetNbinsX():
        xval.append(h.GetXaxis().GetBinCenter(ibin))
        xerr.append((h.GetXaxis().GetBinUpEdge(ibin)-h.GetXaxis().GetBinLowEdge(ibin))/2)
        yval.append(h.GetBinContent(ibin))
        yerr.append(h.GetBinError(ibin))
        ibin += 1

    h_new = ROOT.TGraphErrors(len(xval), array.array('d', xval), array.array('d', yval), 
                              array.array('d', xerr), array.array('d', yerr))

    gxval = []
    gyval = []
    gyerrup = []
    gyerrdown = []
    i = 0
    while i < g.GetN():
        if g.GetX()[i] < 9:
            gxval.append(g.GetX()[i])
            gyval.append(g.GetY()[i])
            gyerrup.append(g.GetEYhigh()[i])
            gyerrdown.append(g.GetEYlow()[i])
            i = i+1
        elif g.GetX()[i] < 16:
            gxval.append((g.GetX()[i] + g.GetX()[i+1]) / 2)
            gyval.append((g.GetY()[i] + g.GetY()[i+1]) / 2)
            gyerrup.append(math.sqrt(g.GetEYhigh()[i]**2+g.GetEYhigh()[i+1]**2)/2)
            gyerrdown.append(math.sqrt(g.GetEYlow()[i]**2+g.GetEYlow()[i+1]**2)/2)
            i = i+2
        else:
            gxval.append((g.GetX()[i] + g.GetX()[i+1] + g.GetX()[i+2] + g.GetX()[i+3]) / 4)
            gyval.append((g.GetY()[i] + g.GetY()[i+1] + g.GetY()[i+2] + g.GetY()[i+3]) / 4)
            gyerrup.append(math.sqrt(g.GetEYhigh()[i]**2+g.GetEYhigh()[i+1]**2+g.GetEYhigh()[i+2]**2+g.GetEYhigh()[i+3]**2)/4)
            gyerrdown.append(math.sqrt(g.GetEYlow()[i]**2+g.GetEYlow()[i+1]**2+g.GetEYlow()[i+2]**2+g.GetEYlow()[i+3]**2)/4)
            i = i+4

    g_new = ROOT.TGraphAsymmErrors(len(gxval), array.array('d', gxval), array.array('d', gyval), 
                              ROOT.nullptr, ROOT.nullptr, array.array('d', gyerrdown), array.array('d', gyerrup))

    return g_new,h_new

def MakeRatio(g, h):
    i = 0
    j = 0
    ratio = []
    ratioerrup = []
    ratioerrdown = []
    while i < g.GetN() and j < h.GetN():
        if g.GetX()[i] < h.GetX()[j]:
            i += 1
            continue
        if g.GetX()[i] > h.GetX()[j]:
            j += 1
            continue
        ratio.append(h.GetX()[j] / g.GetX()[j])

def main(fonll_file, spectrum, gen, proc, ts, compare):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    FONLL_ntuple = FONLL(fonll_file)
    g = FONLL_ntuple.GenerateGraph("pt", "central", "max", "min")
    globalList.append(g)

    c = ROOT.TCanvas(fonll_file,fonll_file)
    globalList.append(c)
    c.SetLogy()
    c.cd()

    if not compare:
        print("FONLL")
        g.SetMarkerStyle(ROOT.kFullCircle)
        g.SetMarkerSize(0.9)
        g.SetMarkerColor(ROOT.kBlue+2)
        g.SetLineColor(ROOT.kBlue+2)
        g.SetFillColor(ROOT.kCyan+2)
        g.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        g.GetYaxis().SetTitle("d#sigma / d#it{p}_{T} [pb (GeV/#it{c})^{-1}]")
        g.Draw("A3")
    elif compare == "fastsim":
        print("Compare with fastsim")
        file_name = " /Volumes/DATA/ALICE/JetResults/FastSim_{gen}_{proc}_{ts}/stage_1/output/FastSimAnalysis_{gen}_{proc}_{ts}.root".format(gen=gen, proc=proc, ts=ts)
        MCGEN_data = MCGEN("POWHEG_7TeV", "POWHEG_7TeV", file_name, spectrum)
        h = MCGEN_data.spectrum
        globalList.append(h)

        (g_new, h_new) = MakeUniform(g, h)
        globalList.append(g_new)
        globalList.append(h_new)

        g_new.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        g_new.GetYaxis().SetTitle("d#sigma / d#it{p}_{T} [pb (GeV/#it{c})^{-1}]")
        g_new.SetMarkerStyle(ROOT.kFullCircle)
        g_new.SetMarkerSize(0.9)
        g_new.SetMarkerColor(ROOT.kBlue+2)
        g_new.SetLineColor(ROOT.kBlue+2)
        g_new.SetFillColor(ROOT.kCyan+1)
        h_new.SetMarkerStyle(ROOT.kOpenCircle)
        h_new.SetMarkerSize(0.9)
        h_new.SetMarkerColor(ROOT.kRed+2)
        h_new.SetLineColor(ROOT.kRed+2)
        g_new.Draw("A3")
        h_new.Draw("P")

        MakeRatio(g_new, h_new)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Comparison between FONLL and POWHEG.')
    parser.add_argument('fonll_file', metavar='fonll.dat',
                        help='Data file with FONLL calculations')
    parser.add_argument("--compare", metavar='COMP',
                        default=None)
    parser.add_argument('--spectrum', metavar='spectrum',
                        default="D0_MCTruth_D_Pt_Spectrum_Normalized")
    parser.add_argument('--gen', metavar='GEN',
                        default="powheg")
    parser.add_argument('--proc', metavar='PROC',
                        default="charm")
    parser.add_argument('--ts', metavar='TS',
                        default="local")
    args = parser.parse_args()

    main(args.fonll_file, args.spectrum, args.gen, args.proc, args.ts, args.compare)

    IPython.embed()
