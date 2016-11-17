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

def main(charm_ts, beauty_ts, jet_type, jet_radius):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    generators = ["powheg"]
    for gen in generators:
        analysis_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius)

def analysis_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius):
    rootPath = "/Volumes/DATA/ALICE/JetResults"
    charmQuark = QuarkSetting("charm")
    beautyQuark = QuarkSetting("beauty")
    quarks = dict()
    quarks["charm"] = charmQuark
    quarks["beauty"] = beautyQuark

    charmQuark.ts = charm_ts
    beautyQuark.ts = beauty_ts 

    for quark in quarks.itervalues():
        quark.path = "{0}/FastSim_{1}_{2}_{3}/stage_1/output".format(rootPath, gen, quark.name, quark.ts)
        quark.filename = "{0}/FastSimAnalysis_{1}_{2}_{3}.root".format(quark.path, gen, quark.name, quark.ts)

    ptD = SpectraSet("BFeedDownVsPtD_{0}_{1}_{2}".format(gen, charmQuark.ts, beautyQuark.ts), "B feed-down vs #it{p}_{T,D}")
    ptD.add(SpectrumDef("D0_MCTruth_{0}_{1}_D_Pt_Spectrum_JetPt_0_Normalized".format(jet_type, jet_radius), "#it{p}_{T,ch jet} > 0", dict(colors=[ROOT.kBlue+2,ROOT.kGreen+2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_{0}_{1}_D_Pt_Spectrum_JetPt_2_Normalized".format(jet_type, jet_radius), "#it{p}_{T,ch jet} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed+2,ROOT.kOrange+2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_{0}_{1}_D_Pt_Spectrum_JetPt_5_Normalized".format(jet_type, jet_radius), "#it{p}_{T,ch jet} > 5 GeV/#it{c}", dict(colors=[ROOT.kAzure+2,ROOT.kCyan+2], markers=[ROOT.kOpenDiamond, ROOT.kOpenDiamond], lines=[None, None])))
    ptD.add(SpectrumDef("D0_MCTruth_{0}_{1}_D_Pt_Spectrum_JetPt_8_Normalized".format(jet_type, jet_radius), "#it{p}_{T,ch jet} > 8 GeV/#it{c}", dict(colors=[ROOT.kMagenta+2,ROOT.kPink+2], markers=[ROOT.kOpenStar, ROOT.kOpenStar], lines=[None, None])))

    ptJet = SpectraSet("BFeedDownVsPtJet_{0}_{1}_{2}".format(gen, charmQuark.ts, beautyQuark.ts), "B feed-down vs #it{p}_{T,ch jet}")
    ptJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Pt_Spectrum_PtD_0_Normalized".format(jet_type, jet_radius), "#it{p}_{T,D} > 0", dict(colors=[ROOT.kBlue+2,ROOT.kGreen+2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Pt_Spectrum_PtD_2_Normalized".format(jet_type, jet_radius), "#it{p}_{T,D} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed+2,ROOT.kOrange+2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    
    zJet = SpectraSet("BFeedDownVsZ_{0}_{1}_{2}".format(gen, charmQuark.ts, beautyQuark.ts), "B feed-down vs #it{z}_{||,D}^{ch jet}")
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_0_5_Normalized".format(jet_type, jet_radius), "0 < #it{p}_{T,D} < 5 GeV/#it{c}", dict(colors=[ROOT.kBlue+2,ROOT.kGreen+2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_5_10_Normalized".format(jet_type, jet_radius), "5 < #it{p}_{T,D} < 10 GeV/#it{c}", dict(colors=[ROOT.kRed+2,ROOT.kOrange+2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_10_15_Normalized".format(jet_type, jet_radius), "10 < #it{p}_{T,D} < 15 GeV/#it{c}", dict(colors=[ROOT.kAzure+2,ROOT.kCyan+2], markers=[ROOT.kOpenDiamond, ROOT.kOpenDiamond], lines=[None, None])))
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_15_20_Normalized".format(jet_type, jet_radius), "15 < #it{p}_{T,D} < 20 GeV/#it{c}", dict(colors=[ROOT.kMagenta+2,ROOT.kPink+2], markers=[ROOT.kOpenStar, ROOT.kOpenStar], lines=[None, None])))
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_20_25_Normalized".format(jet_type, jet_radius), "20 < #it{p}_{T,D} < 25 GeV/#it{c}", dict(colors=[ROOT.kTeal+2,ROOT.kSpring+2], markers=[ROOT.kStar, ROOT.kStar], lines=[None, None])))
    zJet.add(SpectrumDef("D0_MCTruth_{0}_{1}_Jet_Z_Spectrum_PtJet_25_Normalized".format(jet_type, jet_radius), "#it{p}_{T,D} > 25 GeV/#it{c}", dict(colors=[ROOT.kViolet+2,ROOT.kYellow+2], markers=[ROOT.kOpenCross, ROOT.kOpenCross], lines=[None, None])))
    
    spectraSets = [ptD, ptJet, zJet]

    ratioAxis = "({0} #rightarrow D^{{0}}) / ({1} #rightarrow D^{{0}})".format(quarks.values()[1].name[0], quarks.values()[0].name[0])

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
                h = GetSpectrum(quark.filename, jet_type, jet_radius, spectrumDef.name)
                h.SetName("{0}_{1}".format(quark.name, h.GetName()))
                h.SetTitle("{0} #rightarrow D^{{0}}, {1}".format(quark.name[0], spectrumDef.title))
                histos.append(h)
                globalList.append(h)
                quark.histos[spectrumDef] = h
            r = DMesonJetUtils.CompareSpectra(histos[0], histos[1:], cname, opt, optRatio, ratioAxis, "logy", "lineary", c, cRatio, leg, legRatio, spectrumDef.styles)
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

def GetSpectrum(filename, jet_type, jet_radius, spectrum):
    file = ROOT.TFile(filename)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(filename))
        return
    mesonlistname = "D0_MCTruth"
    mesonlist = file.Get(mesonlistname)
    if not mesonlist:
        print("Could not get list {0} from file {1}".format(mesonlistname, filename))
        return
    jetlistname = "{0}_{1}".format(jet_type, jet_radius)
    jetlist = mesonlist.FindObject(jetlistname)
    if not jetlist:
        print("Could not get list {0} from list {1} in file {2}".format(jetlistname, mesonlistname, filename))
        return
    h = jetlist.FindObject(spectrum)
    if not h:
        print("Could not find object {0} in list {1}/{2} in file {3}".format(spectrum, mesonlistname, jetlistname, filename))
        exit(1)
    h_copy = h.Clone("{0}_copy".format(spectrum))
    file.Close()
    return h_copy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down.')
    parser.add_argument('--charm', metavar='CHARM',
                        default="local")
    parser.add_argument('--beauty', metavar='CHARM',
                        default="local")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    args = parser.parse_args()

    main(args.charm, args.beauty, args.jet_type, args.jet_radius)

    IPython.embed()
