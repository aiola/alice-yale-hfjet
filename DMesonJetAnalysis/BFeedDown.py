#!/usr/bin/env python
#python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils

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
    def __init__(self,_name):
        self.name = _name
        self.histos = dict()

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    generators = ["powheg"]
    for gen in generators:
        analysis_for_generator(gen)

def analysis_for_generator(gen):
    rootPath = "/Volumes/DATA/ALICE/JetResults"
    charmQuark = QuarkSetting("charm")
    beautyQuark = QuarkSetting("beauty")
    quarks = dict()
    quarks["charm"] = charmQuark
    quarks["beauty"] = beautyQuark

    # change this as needed
    charmQuark.ts = "1478104073"
    beautyQuark.ts = "1478104389"    

    for quark in quarks.itervalues():
        quark.path = "{0}/FastSim_{1}_{2}_{3}/stage_2/output/001".format(rootPath, gen, quark.name, quark.ts)
        quark.filename = "{0}/FastSimAnalysis_{1}_{2}_{3}.root".format(quark.path, gen, quark.name, quark.ts)

    ptD = SpectraSet("BFeedDownVsPtD_{0}_{1}_{2}".format(gen, charmQuark.ts, beautyQuark.ts), "B feed-down vs #it{p}_{T,D}")
    ptD.add(SpectrumDef("D0_MCTruth_D_Pt_Spectrum", "#it{p}_{T,ch jet} > 0", dict(colors=[ROOT.kBlue+2,ROOT.kGreen+2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_D_Pt_Spectrum_JetPt_2", "#it{p}_{T,ch jet} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed+2,ROOT.kOrange+2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_D_Pt_Spectrum_JetPt_5", "#it{p}_{T,ch jet} > 5 GeV/#it{c}", dict(colors=[ROOT.kAzure+2,ROOT.kCyan+2], markers=[ROOT.kOpenDiamond, ROOT.kOpenDiamond], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_D_Pt_Spectrum_JetPt_8", "#it{p}_{T,ch jet} > 8 GeV/#it{c}", dict(colors=[ROOT.kMagenta+2,ROOT.kPink+2], markers=[ROOT.kOpenStar, ROOT.kOpenStar], lines=[None, None])))

    ptJet = SpectraSet("BFeedDownVsPtJet_{0}_{1}_{2}".format(gen, charmQuark.ts, beautyQuark.ts), "B feed-down vs #it{p}_{T,ch jet}")
    ptJet.add(SpectrumDef("D0_MCTruth_D_Tagged_Jet_Pt_Spectrum", "#it{p}_{T,D} > 0", dict(colors=[ROOT.kBlue+2,ROOT.kGreen+2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptJet.add(SpectrumDef("D0_MCTruth_D_Tagged_Jet_Pt_Spectrum_PtD_2", "#it{p}_{T,D} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed+2,ROOT.kOrange+2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    spectraSets = [ptD, ptJet]

    for spectraSet in spectraSets:
        cname = "BFeedDown_{0}_{1}".format(gen, spectraSet.name)
        opt = ""
        optRatio = ""
        c = None
        cRatio = None
        leg = None
        legRatio = None
        for spectrumDef in spectraSet.spectra:
            histos = []
            for quark in quarks.itervalues():
                h = GetSpectrum(quark.filename, spectrumDef.name)
                h.SetName("{0}_{1}".format(quark.name, h.GetName()))
                h.SetTitle("{0} #rightarrow D^{{0}}, {1}".format(quark.name[0], spectrumDef.title))
                histos.append(h)
                quark.histos[spectrumDef] = h
            r = DMesonJetUtils.CompareSpectra(histos[0], histos[1:], cname, opt, optRatio, "ratio", "logy", "lineary", c, cRatio, leg, legRatio, spectrumDef.styles)
            opt = "same"
            optRatio = "same"
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TLegend):
                    if "Ratio" in obj.GetName():
                        legRatio = obj
                    else:
                        leg = obj
                elif isinstance(obj, ROOT.TCanvas):
                    if "Ratio" in obj.GetName():
                        cRatio = obj
                    else:
                        c = obj
        c.SaveAs("{0}/{1}.pdf".format(rootPath, c.GetName()))
        cRatio.SaveAs("{0}/{1}.pdf".format(rootPath,cRatio.GetName()))

def GetSpectrum(filename, spectrum):
    file = ROOT.TFile(filename)
    rlist = file.Get(spectrum)
    h = rlist.FindObject(spectrum).Clone("{0}_copy".format(spectrum))
    h.Scale(1, "width")
    h.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    file.Close()
    return h

if __name__ == '__main__':
    main()
    
    IPython.embed()
