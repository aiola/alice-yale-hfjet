#!/usr/bin/env python
#python script to plot FONLL calculations

import argparse
import IPython
import ROOT
import array
import math

globalList = []

class FONLL:
    def __init__(self, fonll_file, gname="D0Kpiprediction"):
        self.name = fonll_file
        self.n = 0
        if ".root" in fonll_file:
            self.LoadD2HGraph(fonll_file, gname)
        else:
            self.values = dict()
            self.OpenFONLL(fonll_file)
            self.GenerateGraph("pt", "central", "max", "min")

    def LoadD2HGraph(self, filename, gname):
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(filename))
            exit(1)
        g = file.Get(gname)
        if not g:
            print("Could not get graph {0}".format(gname))
            exit(1)
        self.n = g.GetN()
        self.fGraph = g.Clone()

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
        self.fGraph = ROOT.TGraphAsymmErrors(self.n, array.array('d', self.values[x]), array.array('d', self.values[y]), ROOT.nullptr, ROOT.nullptr, 
                                             array.array('d', [v-down for down,v in zip(self.values[yerrdown], self.values[y])]),
                                             array.array('d', [up-v for up,v in zip(self.values[yerrup], self.values[y])]))

class MCGEN:
    def __init__(self, name, title, file_name, spectrum_name):
        self.file_name = file_name
        self.spectrum_name = spectrum_name
        self.name = name
        self.title = title
        self.xsec = 1e9 # mb -> pb
        self.xsec /= 2 # normalize to particle (antiparticle)
        #self.xsec *= 3.93e-2 # decay D0->Kpi
        self.LoadSpectrum()

    def LoadSpectrum(self):
        file = ROOT.TFile.Open(self.file_name)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(self.file_name))
            return
        mesonlistname = "D0_MCTruth"
        mesonlist = file.Get(mesonlistname)
        if not mesonlist:
            print("Could not get list {0} from file {1}".format(mesonlistname, self.file_name))
            return
        spectrum = mesonlist.FindObject(self.spectrum_name)
        if not spectrum:
            print("Could not get histogram {0} from list {1} in file {2}".format(self.spectrum_name, mesonlistname, self.file_name))
            return
        self.spectrum = spectrum.Clone(self.name)
        self.spectrum.SetTitle(self.title)
        self.spectrum.Scale(self.xsec)

def GenerateGraphFromHist(h):
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
    return h_new

def MakeUniform(g,h_new, scaling=1):
    gxval = []
    gyval = []
    gyerrup = []
    gyerrdown = []
    i = 0
    j = 0
    myXval = 0
    myYval = 0
    myerrup2 = 0
    myerrdown2 = 0
    nsum = 0
    while i < g.GetN() and j < h_new.GetN():
        if nsum == 0 or myXval/nsum <= h_new.GetX()[j]:
            print("Adding {0}".format(g.GetX()[i]))
            myXval += g.GetX()[i]
            myYval += g.GetY()[i]
            myerrup2 += g.GetEYhigh()[i]**2
            myerrdown2 += g.GetEYlow()[i]**2
            nsum += 1
        if nsum > 0 and myXval/nsum - h_new.GetX()[j] > 1e-6:
            print("Something wrong...")
        if nsum > 0 and math.fabs(myXval/nsum - h_new.GetX()[j]) < 1e-6:
            j += 1
            print("Bin X {0}".format(myXval/nsum))
            gxval.append(round(myXval/nsum,2))
            gyval.append(myYval/nsum*scaling)
            gyerrup.append(math.sqrt(myerrup2)/nsum*scaling)
            gyerrdown.append(math.sqrt(myerrdown2)/nsum*scaling)
            myXval = 0
            myYval = 0
            myerrup2 = 0
            myerrdown2 = 0
            nsum = 0
        i += 1

    print(gxval)
    g_new = ROOT.TGraphAsymmErrors(len(gxval), array.array('d', gxval), array.array('d', gyval), 
                                   ROOT.nullptr, ROOT.nullptr, array.array('d', gyerrdown), array.array('d', gyerrup))

    return g_new

def MakeRatio(g, h):
    i = 0
    j = 0
    xval = []
    xerr = []
    ratio = []
    ratioerr = []
    unc_up = []
    unc_down = []
    while i < g.GetN() and j < h.GetN():
        if g.GetX()[i] < h.GetX()[j]:
            i += 1
            continue
        if g.GetX()[i] > h.GetX()[j]:
            j += 1
            continue
        xval.append(h.GetX()[j])
        xerr.append((h.GetErrorX(j)))
        ratio.append(h.GetY()[j] / g.GetY()[j])
        ratioerr.append(h.GetErrorY(j) / g.GetY()[j])
        unc_up.append(g.GetEYhigh()[j] / g.GetY()[j])
        unc_down.append(g.GetEYlow()[j] / g.GetY()[j])
        i += 1
        j += 1
    ratio_h = ROOT.TGraphErrors(len(xval), array.array('d', xval), array.array('d', ratio), 
                                array.array('d', xerr), array.array('d', ratioerr))
    unc_g = ROOT.TGraphAsymmErrors(len(xval), array.array('d', xval), array.array('d', [1]*len(xval)), 
                                   array.array('d', xerr), array.array('d', xerr), 
                                   array.array('d', unc_down), array.array('d', unc_up))
    
    return unc_g, ratio_h
    

def main(fonll_file, spectrum, gen, proc, ts, compare, fonll_file_2):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    FONLL_ntuple = FONLL(fonll_file)
    g = FONLL_ntuple.fGraph
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
        g.SetFillColor(ROOT.kCyan+1)
        g.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        g.GetYaxis().SetTitle("d#sigma / d#it{p}_{T} [pb (GeV/#it{c})^{-1}]")
        g.Draw("A3")
    elif compare == "fastsim":
        print("Compare with fastsim")
        file_name = " /Volumes/DATA/ALICE/JetResults/FastSim_{gen}_{proc}_{ts}/stage_1/output/FastSimAnalysis_{gen}_{proc}_{ts}.root".format(gen=gen, proc=proc, ts=ts)
        MCGEN_data = MCGEN("POWHEG_7TeV", "POWHEG_7TeV", file_name, spectrum)
        h = MCGEN_data.spectrum
        globalList.append(h)

        h_new = GenerateGraphFromHist(h)
        g_new = MakeUniform(g, h_new)
        globalList.append(g_new)
        globalList.append(h_new)

        g_new.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
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

        leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.85)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(g_new, "FONLL", "f")
        leg.AddEntry(h_new, "POWHEG+PYTHIA", "pe")
        leg.Draw()
        globalList.append(leg)

        (unc_g, ratio_h) = MakeRatio(g_new, h_new)
        print("Ratio ok")
        globalList.append(unc_g)
        globalList.append(ratio_h)
        cRatio = ROOT.TCanvas("{0}_ratio".format(fonll_file),"{0}_ratio".format(fonll_file))
        globalList.append(cRatio)
        cRatio.cd()

        unc_g.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
        unc_g.GetYaxis().SetTitle("POWHEG+PYTHIA / FONLL")
        unc_g.SetMarkerStyle(ROOT.kFullCircle)
        unc_g.SetMarkerSize(0.9)
        unc_g.SetMarkerColor(ROOT.kBlue+2)
        unc_g.SetLineColor(ROOT.kBlue+2)
        unc_g.SetFillColor(ROOT.kCyan+1)
        ratio_h.SetMarkerStyle(ROOT.kOpenCircle)
        ratio_h.SetMarkerSize(0.9)
        ratio_h.SetMarkerColor(ROOT.kRed+2)
        ratio_h.SetLineColor(ROOT.kRed+2)
        unc_g.Draw("A3")
        ratio_h.Draw("P")

    elif compare == "fonll":
        print("Compare with FONLL")
        FONLL_ntuple2 = FONLL(fonll_file_2)
        g2 = FONLL_ntuple2.fGraph
        globalList.append(g2)

        g_new = MakeUniform(g, g2, 1./3.93e-2)
        globalList.append(g_new)

        g_new.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
        g_new.GetYaxis().SetTitle("d#sigma / d#it{p}_{T} [pb (GeV/#it{c})^{-1}]")
        g_new.SetMarkerStyle(ROOT.kFullCircle)
        g_new.SetMarkerSize(0.9)
        g_new.SetMarkerColor(ROOT.kBlue+2)
        g_new.SetLineColor(ROOT.kBlue+2)
        g_new.SetFillColor(ROOT.kCyan+1)
        g2.SetMarkerStyle(ROOT.kOpenCircle)
        g2.SetMarkerSize(0.9)
        g2.SetMarkerColor(ROOT.kRed+2)
        g2.SetLineColor(ROOT.kRed+2)
        g_new.Draw("AP")
        g2.Draw("P")

        leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.85)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(g_new, "FONLL 1", "pe")
        leg.AddEntry(g2, "FONLL 2", "pe")
        leg.Draw()
        globalList.append(leg)

        (unc_g, ratio_g2) = MakeRatio(g_new, g2)
        print("Ratio ok")
        globalList.append(unc_g)
        globalList.append(ratio_g2)
        cRatio = ROOT.TCanvas("{0}_ratio".format(fonll_file),"{0}_ratio".format(fonll_file))
        globalList.append(cRatio)
        cRatio.cd()

        unc_g.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
        unc_g.GetYaxis().SetTitle("FONLL 2 / FONLL 1")
        unc_g.SetMarkerStyle(ROOT.kFullCircle)
        unc_g.SetMarkerSize(0.9)
        unc_g.SetMarkerColor(ROOT.kBlue+2)
        unc_g.SetLineColor(ROOT.kBlue+2)
        unc_g.SetFillColor(ROOT.kCyan+1)
        ratio_g2.SetMarkerStyle(ROOT.kOpenCircle)
        ratio_g2.SetMarkerSize(0.9)
        ratio_g2.SetMarkerColor(ROOT.kRed+2)
        ratio_g2.SetLineColor(ROOT.kRed+2)
        #unc_g.Draw("APX")
        ratio_g2.Draw("APX")
        

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
    parser.add_argument('--fonll2', metavar='fonll.dat',
                        default=None)
    args = parser.parse_args()

    main(args.fonll_file, args.spectrum, args.gen, args.proc, args.ts, args.compare, args.fonll2)

    IPython.embed()
