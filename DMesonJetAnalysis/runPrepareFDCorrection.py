#!/usr/bin/env python
#python script to prepare B feed-down correction file

import argparse
import IPython
import ROOT
import array
import math

def main(ts, bResponse, cResponse):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    suffix = "powheg_beauty_{0}".format(ts)
    name = "FastSim_{0}".format(suffix)
    input_path = "/Volumes/DATA/ALICE/JetResults"
    input_file_name = "{0}/FastSim_{1}/stage_1/output/FastSimAnalysis_{1}.root".format(input_path, suffix)
    output_file_name = "{0}/FastSim_{1}/stage_1/output/BFeedDown_{1}.root".format(input_path, suffix)

    FDhistogram_jetpt_orig, FDhistogram_dpt_orig = LoadFDHistogram(input_file_name)
    histograms = []
    histograms += PrepareFDhist_dpt(FDhistogram_dpt_orig, input_path, bResponse, cResponse)
    histograms += PrepareFDhist_jetpt(FDhistogram_jetpt_orig, input_path, bResponse, cResponse)
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
        print("Could not open output file {0}".format(output_file_name))
        exit(1)
    for h in histograms:
        h.Write()
    outputFile.Close()
    print("Done: {0}.".format(output_file_name))

def PrepareFDhist_dpt(FDhistogram_old, input_path, bResponse, cResponse):
    print("Preparing D pt FD histograms")
    result = []

    dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 16, 30]

    print("Rebinning the FD original histogram")
    FDhistogram_orig = FDhistogram_old.Rebin(len(dptbins)-1, FDhistogram_old.GetName(), array.array('d', dptbins))
    result.append(FDhistogram_orig)

    print("Opening the b->D0 response matrix file")
    bResponseFile = OpenResponseFile(input_path, bResponse, False)
    print("Loading the response matrix")
    temp, bEfficiency_dpt = LoadResponse(bResponseFile, "JetPtDPtSpectrum", "_NoJet", "b", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
    result.append(bEfficiency_dpt)

    print("Applying the b->D0 reconstruction efficiency")
    FDhistogram_dpt = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
    FDhistogram_dpt.Multiply(bEfficiency_dpt)
    result.append(FDhistogram_dpt)

    print("Opening the c->D0 response matrix file")
    cResponseFile = OpenResponseFile(input_path, cResponse, False)
    print("Loading the response matrix")
    temp, cEfficiency_dpt = LoadResponse(cResponseFile, "JetPtDPtSpectrum", "_NoJet", "c", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
    result.append(cEfficiency_dpt)

    print("Applying the correction for the c->D0 reconstruction efficiency")
    FDhistogram_efficiency_dpt = FDhistogram_dpt.Clone("{0}_efficiency".format(FDhistogram_dpt.GetName()))
    FDhistogram_efficiency_dpt.Divide(cEfficiency_dpt)
    result.append(FDhistogram_efficiency_dpt)

    return result

def PrepareFDhist_jetpt(FDhistogram_old, input_path, bResponse, cResponse):
    print("Preparing jet pt FD histograms")
    result = []

    jetptbins = [5, 6, 8, 10, 14, 20, 30]

    print("Rebinning the FD original histogram")
    FDhistogram_orig = Rebin2D(FDhistogram_old, FDhistogram_old.GetName(), len(jetptbins)-1, array.array('d', jetptbins), FDhistogram_old.GetNbinsY(), FDhistogram_old.GetYaxis().GetXbins().GetArray())
    result.append(FDhistogram_orig)
    
    print("Opening the b->D0 response matrix file")
    bResponseFile = OpenResponseFile(input_path, bResponse, False)
    print("Loading the response matrix")
    temp, bEfficiency_dpt = LoadResponse(bResponseFile, "JetPtDPtSpectrum", "_JetPt_500_3000", "b", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
    result.append(bEfficiency_dpt)

    print("Applying the b->D0 reconstruction efficiency")
    FDhistogram = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
    ApplyEfficiency(FDhistogram, bEfficiency_dpt, False)
    result.append(FDhistogram)

    print("Opening the c->D0 response matrix file")
    cResponseFile = OpenResponseFile(input_path, cResponse, False)
    print("Loading the response matrix")
    temp, cEfficiency_dpt = LoadResponse(cResponseFile, "JetPtDPtSpectrum", "_JetPt_500_3000", "c", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
    result.append(cEfficiency_dpt)

    print("Applying the correction for the c->D0 reconstruction efficiency")
    FDhistogram_efficiency = FDhistogram.Clone("{0}_efficiency".format(FDhistogram.GetName()))
    ApplyEfficiency(FDhistogram_efficiency, cEfficiency_dpt, True)
    result.append(FDhistogram_efficiency)
    
    print("Opening the b->D0 response matrix file (w/ efficiency)")
    bResponseFile_efficiency = OpenResponseFile(input_path, bResponse, True)

    for ptd in range(2,5):
        print("Loading the response matrix for jet pt resolution correction (w/o efficiency) for pt > {0} GeV/c".format(ptd))
        bDetectorResponse_jetpt, temp = LoadResponse(bResponseFile, "JetPtSpectrum_DPt_{0}".format(ptd*10), "", "b", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
        result.append(bDetectorResponse_jetpt)

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency)".format(ptd))
        FDhistogram_jetpt = FDhistogram.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram.GetName(),ptd*10))
        result.append(FDhistogram_jetpt)

        print("Applying the jet pt response matrix")
        FDhistogram_jetpt_detector = ApplyResponse(FDhistogram_jetpt, bDetectorResponse_jetpt)
        result.append(FDhistogram_jetpt_detector)
        
        print("Loading the response matrix for jet pt resolution correction (w/ efficiency) for pt > {0} GeV/c".format(ptd))
        bDetectorResponse_efficiency_jetpt, temp = LoadResponse(bResponseFile_efficiency, "JetPtSpectrum_DPt_{0}".format(ptd*10), "", "b_cEff", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
        result.append(bDetectorResponse_efficiency_jetpt)

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/ efficiency)".format(ptd))
        FDhistogram_efficiency_jetpt = FDhistogram_efficiency.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram_efficiency.GetName(),ptd*10), FDhistogram_efficiency.GetYaxis().FindBin(ptd), FDhistogram_efficiency.GetNbinsY()+1)
        result.append(FDhistogram_efficiency_jetpt)

        print("Applying the jet pt response matrix")
        FDhistogram_efficiency_jetpt_detector = ApplyResponse(FDhistogram_efficiency_jetpt, bDetectorResponse_efficiency_jetpt)
        result.append(FDhistogram_efficiency_jetpt_detector)

    return result

def ApplyEfficiency(hist, efficiency, reverse):
    print("Applying efficiency {0} to histogram {1}".format(efficiency.GetName(), hist.GetName()))
    for ybin in range(0, hist.GetNbinsY()+2):
        if reverse and efficiency.GetBinContent(ybin)>0:
            print("Dividing by {0}".format(efficiency.GetBinContent(ybin)))
            eff = 1./efficiency.GetBinContent(ybin)
        else:
            print("Multiplying by {0}".format(efficiency.GetBinContent(ybin)))
            eff = efficiency.GetBinContent(ybin)
        for xbin in range(0, hist.GetNbinsX()+2):
            hist.SetBinContent(xbin, ybin, hist.GetBinContent(xbin, ybin)*eff)
            hist.SetBinError(xbin, ybin, hist.GetBinError(xbin, ybin)*eff)

def ApplyResponse(histogram, response):
    result = ROOT.TH1D("{0}_detector".format(histogram.GetName()), "{0} detector".format(histogram.GetTitle()), response.GetNbinsX(), response.GetXaxis().GetXbins().GetArray())
    result.GetXaxis().SetTitle("#it{p}_{T, ch jet}^{det} (GeV/#it{c})")
    result.GetYaxis().SetTitle(histogram.GetYaxis().GetTitle())
    result.Sumw2()
    for xbin in range(0, response.GetNbinsX()+2):
        content = 0
        error = 0
        sumw = 0
        for ybin in range(0, response.GetNbinsY()+2):
            w = response.GetBinContent(xbin, ybin)
            sumw += w
            content += w*histogram.GetBinContent(ybin)
            error += w*histogram.GetBinError(ybin)
        if sumw <= 0:
            continue
        content /= sumw
        error /= sumw
        result.SetBinContent(xbin, content)
        result.SetBinError(xbin, error)
    return result

def OpenResponseFile(input_path, responseTrain, efficiency):
    if efficiency:
        file_name = "{0}/Jets_EMC_pp_MC_{1}_{2}_{3}_{4}/LHC15i2_Train{1}_efficiency.root".format(input_path, responseTrain, responseTrain+1, responseTrain+2, responseTrain+3)
    else:
        file_name = "{0}/Jets_EMC_pp_MC_{1}_{2}_{3}_{4}/LHC15i2_Train{1}.root".format(input_path, responseTrain, responseTrain+1, responseTrain+2, responseTrain+3)
    file = ROOT.TFile(file_name)
    print("Response file: {0}".format(file_name))
    if not file or file.IsZombie():
        print("Could not open input file {0}".format(file_name))
        exit(1)
    return file

def LoadResponse(responseFile, spectrumName, suffix_in, suffix_out, detAxisNbins, detAxisBin):
    rlistName = "D0_Jet_AKTChargedR040_pt_scheme_{0}".format(spectrumName)
    rlist = responseFile.Get(rlistName)
    if not rlist:
        print("Could not get list {0}".format(rlistName))
        exit(1)
    truthName = "D0_Jet_AKTChargedR040_pt_scheme_{0}_Truth{1}".format(spectrumName, suffix_in)
    truth = rlist.FindObject(truthName)
    if not truth:
        print("Could not get histogram {0}".format(truthName))
        exit(1)
    truth.SetName("{0}_{1}".format(truthName, suffix_out))
    truth_coarse = truth.Rebin(detAxisNbins, "{0}_coarse_{1}".format(truthName, suffix_out), detAxisBin)
    respName = "D0_Jet_AKTChargedR040_pt_scheme_{0}_DetectorResponse{1}".format(spectrumName, suffix_in)
    resp = rlist.FindObject(respName)
    if not resp:
        print("Could not get histogram {0}".format(respName))
        exit(1)
    resp.SetName("{0}_{1}".format(respName, suffix_out))
    resp_coarse = Rebin2D(resp, "{0}_coarse_{1}".format(respName, suffix_out), detAxisNbins, detAxisBin, detAxisNbins, detAxisBin)
    eff_coarse = resp_coarse.ProjectionY(truth_coarse.GetName().replace("Truth", "Efficiency"))
    eff_coarse.Divide(truth_coarse)
    return resp_coarse, eff_coarse

def Rebin2D(hist, name, nbinsX, binsX, nbinsY, binsY):
    result = ROOT.TH2D(name, name, nbinsX, binsX, nbinsY, binsY)
    result.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    result.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    result.GetZaxis().SetTitle(hist.GetZaxis().GetTitle())
    result.Sumw2()
    for xbin in range(0, hist.GetNbinsX()+2):
        xbinRes = result.GetXaxis().FindBin(hist.GetXaxis().GetBinCenter(xbin))
        for ybin in range(0, hist.GetNbinsY()+2):
            ybinRes = result.GetYaxis().FindBin(hist.GetYaxis().GetBinCenter(ybin))
            result.SetBinContent(xbinRes, ybinRes, result.GetBinContent(xbinRes,ybinRes)+hist.GetBinContent(xbin,ybin))
            result.SetBinError(xbinRes, ybinRes, math.sqrt(result.GetBinError(xbinRes,ybinRes)**2+hist.GetBinError(xbin,ybin)**2))
    return result

def LoadFDHistogram(file_name):
    file = ROOT.TFile(file_name)
    if not file or file.IsZombie():
        print("Could not open input file {0}".format(file_name))
        exit(1)
    else:
        print("File {0} open".format(file_name))
    dlist = file.Get("D0_MCTruth")
    jlist = dlist.FindObject("Charged_R040")
    slist = jlist.FindObject("D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum")
    jetpt_hist = slist.FindObject("D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum")
    if not jetpt_hist:
        print("Could not find FD histogram (jet pt)!")
        slist.ls()
        exit(1)
    slist = dlist.FindObject("D0_MCTruth_D_Pt_Spectrum")
    dpt_hist = slist.FindObject("D0_MCTruth_D_Pt_Spectrum")
    if not dpt_hist:
        print("Could not find FD histogram (d pt)!")
        slist.ls()
        exit(1)

    return jetpt_hist, dpt_hist

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('--charm-response', metavar='charm',
                        default=None, type=int)
    parser.add_argument('--beauty-response', metavar='beauty',
                        default=None, type=int)
    parser.add_argument('--ts', metavar='TS',
                        default=None)
    args = parser.parse_args()
    main(args.ts, args.beauty_response, args.charm_response)

    IPython.embed()
