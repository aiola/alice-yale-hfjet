#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import DMesonJetCompare

globalList = []

class SpectrumDef:
    def __init__(self, _name, _title, _styles):
        self.name = _name
        self.title = _title
        self.styles = _styles

class SpectraSet:
    def __init__(self, _name, _title):
        self.name = _name
        self.title = _title
        self.spectra = []

    def add(self, s):
        self.spectra.append(s)

class QuarkSetting:
    def __init__(self, _name):
        self.name = _name
        self.histos = dict()

def main(jet_type, jet_radius, ppdata, pPbdata):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    pPb_comparison(jet_type, jet_radius, ppdata, pPbdata)

class DataSpectrumDef:
    def __init__(self, spectrumName, unfoldingMethod, unfoldingReg, unfoldingPrior):
        self.fSpectrumName = spectrumName
        self.fUnfoldingMethod = unfoldingMethod
        self.fUnfoldingReg = unfoldingReg
        self.fUnfoldingPrior = unfoldingPrior
        self.fCrossSection = 62.3  # mb CINT1
        self.fBranchingRatio = 0.0393  # D0->Kpi
        self.fAntiParticleNorm = 2.0

    def LoadSpectrum(self, file):
        spectrumList = file.Get(self.fSpectrumName)
        if not spectrumList:
            print("Could not get list {0} from file {1}".format(self.fSpectrumName, file.GetName()))
            file.ls()
            exit(1)
        else:
            print("List {0} loaded from file {1}".format(self.fSpectrumName, file.GetName()))
        methodList = spectrumList.FindObject(self.fUnfoldingMethod)
        if not methodList:
            print("Could not get list {0} from list {1}".format(self.fUnfoldingMethod, self.fSpectrumName))
            spectrumList.Print()
            exit(1)
        else:
            print("List {0} loaded from list {1}".format(self.fUnfoldingMethod, self.fSpectrumName))
        hname = "_".join([self.fSpectrumName, "UnfoldedSpectrum", self.fUnfoldingMethod, self.fUnfoldingReg, self.fUnfoldingPrior])
        h = methodList.FindObject(hname)
        if not h:
            print("Could not get histogram {0} from list {1}".format(hname, self.fUnfoldingMethod))
            methodList.Print()
            exit(1)
        else:
            print("Histogram {0} loaded from list {1}".format(hname, self.fUnfoldingMethod))
        self.fHistogram = h.Clone("{0}_copy".format(h.GetName()))
        self.fEvents = spectrumList.FindObject("Events")
        if not h:
            print("Could not get histogram {0} from list {1}".format("Events", spectrumList.GetName()))
            exit(1)
        else:
            print("Histogram {0} loaded from list {1}".format("Events", spectrumList.GetName()))
        self.fNumberOfEvents = self.fEvents.GetBinContent(1)
        self.fNormalizedHistogram = self.fHistogram.Clone("{0}_Normalized".format(hname))
        self.fNormalizedHistogram.Scale(self.fCrossSection / (self.fNumberOfEvents * self.fBranchingRatio * self.fAntiParticleNorm), "width")

def GetpPbSpectrum(pPbdata):
    file = ROOT.TFile(pPbdata)
    h = file.Get("hJetPtSB_unf")
    h.GetXaxis().SetTitle("#it{p}_{T, ch jet} (GeV/#it{c})")
    h.Scale(0.5)
    return h

def pPb_comparison(jet_type, jet_radius, ppdata, pPbdata):
    rootPath = "/Volumes/DATA/ALICE/JetResults"

    dataFileName = "{0}/{1}/{1}.root".format(rootPath, ppdata)
    dataFile = ROOT.TFile(dataFileName)
    if not dataFile or dataFile.IsZombie():
        print("Could not open file {0}".format(dataFileName))
        exit(1)
    else:
        print("File {0} successfully open".format(dataFileName))

    spectra = []
    spectrum = DataSpectrumDef("SideBand_DPt_30", "Bayes", "Reg4", "PriorResponseTruth")
    spectrum.LoadSpectrum(dataFile)
    spectrum.fNormalizedHistogram.SetTitle("pp, #sqrt{#it{s}} = 7 TeV")
    globalList.append(spectrum.fNormalizedHistogram)
    spectrum.fNormalizedHistogram.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    spectrum.fNormalizedHistogram.GetXaxis().SetTitle("#it{p}_{T, ch jet} (GeV/#it{c})")

    cname = "_".join(["pPbComparison", jet_type, jet_radius] + [ppdata])
    ratioAxis = "#frac{d^{2}#sigma_{p-Pb}}{d#it{p}_{T}d#eta} / #frac{d^{2}#sigma_{pp}}{d#it{p}_{T}d#eta}"
    pPbSpectrum = GetpPbSpectrum(pPbdata)
    pPbSpectrum.SetTitle("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV")
    globalList.append(pPbSpectrum)
    histos = []
    histos.append(pPbSpectrum)

    globalList.extend(histos)
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fYaxisRatio = ratioAxis
    r = comp.CompareSpectra(spectrum.fNormalizedHistogram, histos)
    for obj in r:
        globalList.append(obj)
        if isinstance(obj, ROOT.TH1) and "Ratio" in obj.GetName():
            obj.GetYaxis().SetRangeUser(0, 4)
    for obj in r:
        if isinstance(obj, ROOT.TCanvas):
            if "Ratio" in obj.GetName():
                obj.SetLeftMargin(0.12)
                obj.SetGridy()
            obj.SaveAs("{0}/{1}.pdf".format(rootPath, obj.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down.')
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--pp', metavar='LHC10_Train823_LHC15i2_Train961_efficiency',
                        default=None)
    parser.add_argument('--pPb', metavar='CorrectedJetSpectra_ppBinningSB_pPb_Barbara_20161213.root',
                        default=None)
    args = parser.parse_args()

    main(args.jet_type, args.jet_radius, args.pp, args.pPb)

    IPython.embed()
