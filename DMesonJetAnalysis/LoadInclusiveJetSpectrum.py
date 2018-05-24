# Load inclusive jet spectrum

import numpy
import ROOT
import DMesonJetUtils

def GetCrossSection():
    sigmapp = 73.2  # published x-section in mbarn
    triggerEff = 62.2 / 73.2  # eff of MB OR trigger to be applied to data
    bin0CorrData = 0.91  # ratio event after vertexNcontributors cut to evts after phys selection

    fname = "../obusch/outData_spec_Bayes_combPtH_rebinned.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat_old = file.Get("hSpecComb")
    if not hStat_old:
        print("Cannot get inclusive cross section with statistical uncertainty!")
        exit(1)
    file.Close()

    hStat_old.Scale(sigmapp * triggerEff)
    hStat_old.Scale(bin0CorrData)

    jetptbins = [5, 6, 8, 10, 14, 20, 30]
    hStat = DMesonJetUtils.Rebin1D_fromBins(hStat_old, "{}_rebinned".format(hStat_old.GetName()), len(jetptbins) - 1, numpy.array(jetptbins, dtype=numpy.float64))

    fname = "../obusch/outData_spec_Bayes_combPtH.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat_orig = file.Get("hSpecComb")
    if not hStat_orig:
        print("Cannot get inclusive cross section with statistical uncertainty (not rebinned)!")
        exit(1)
    file.Close()

    hStat_orig.Scale(sigmapp * triggerEff)
    hStat_orig.Scale(bin0CorrData)

    fname = "../obusch/outSys_tot_spec.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hSyst_orig = file.Get("grSysErrTotNoUE")
    if not hSyst_orig:
        print("Cannot get inclusive cross section with systematic uncertainty!")
        exit(1)
    file.Close()

    hSyst = ROOT.TGraphErrors(hStat.GetNbinsX())

    for i in range(1, hStat.GetNbinsX() + 1):
        x = hStat.GetXaxis().GetBinCenter(i)
        y = hStat.GetBinContent(i)
        xerr = hStat.GetXaxis().GetBinWidth(i) / 2
        print("hStat Bin {}, x in [{}, {}], y = {}".format(i, x - xerr, x + xerr, y))
        yerr_rel = 0.
        y_check = 0.
        for iorig in range(0, hSyst_orig.GetN()):
            if hSyst_orig.GetX()[iorig] < hStat.GetXaxis().GetBinLowEdge(i): continue
            if hSyst_orig.GetX()[iorig] > hStat.GetXaxis().GetBinUpEdge(i): break
            yerr_rel += hSyst_orig.GetErrorY(iorig) * hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            y_check += hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            print("hSyst_orig Bin {}, x = {}, err_y_rel = {}".format(iorig, hSyst_orig.GetX()[iorig], hSyst_orig.GetErrorY(iorig)))
            print("hStat_orig Bin {}, x = {}, y = {}, err_y_rel = {}".format(iorig, hStat_orig.GetXaxis().GetBinCenter(iorig + 1), hStat_orig.GetBinContent(iorig + 1), hStat_orig.GetBinError(iorig + 1)))
        yerr_rel /= y_check
        y_check /= hStat.GetXaxis().GetBinWidth(i)
        yerr = yerr_rel * y
        hSyst.SetPoint(i - 1, x, y)
        hSyst.SetPointError(i - 1, xerr, yerr)
        print("Bin {} Tot Error rel {}".format(i, yerr_rel))
        print("Bin {}, tot y {}, check {}\n".format(i, y, y_check))

    return hStat, hSyst