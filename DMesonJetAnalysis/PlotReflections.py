#!/usr/bin/env python

import argparse
import math
import yaml
import IPython
import ROOT

globalList = []

InvMassRange = [1.72, 2.014]

def PlotReflections(config, index_list, var, reflFit):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    fname = "{input_path}/{train}/{ana}/reflTemp/{ana}_{var}_fitted_{fit}.root".format(input_path=config["input_path"], train=config["train"], ana=config["name"], var=var, fit=reflFit)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        return

    signHist = dict()
    reflHist = dict()
    reflFitHist = dict()

    for k in file.GetListOfKeys():
        h = file.Get(k.GetName())
        if "histSgn" in h.GetName():
            i = int(h.GetName()[-1])
            signHist[i] = h
            print(h.GetName())
        if "histRflFitted" in h.GetName():
            i = int(h.GetName()[-1])
            reflFitHist[i] = h
            print(h.GetName())
        if "histoRfl" in h.GetName():
            i = int(h.GetName()[-1])
            reflHist[i] = h
            print(h.GetName())

    globalList.extend(signHist)
    globalList.extend(reflHist)

    if len(index_list) == 0:
        n = len(reflHist)
    else:
        n = len(index_list)
    ncols = int(math.ceil((math.sqrt(n))))
    nrows = int(math.floor((n / ncols)))

    cname = "ReflectionTemplates_{var}_{fit}".format(var=var, fit=reflFit)
    canvas = ROOT.TCanvas(cname, cname, ncols*350, nrows*350)
    globalList.append(canvas)
    canvas.Divide(ncols, nrows)

    icanvas = 1
    for i, (hSig, hRef, hRefFit) in enumerate(zip(signHist.itervalues(), reflHist.itervalues(), reflFitHist.itervalues())):
        if len(index_list) > 0 and not i in index_list:
            continue
        pad = canvas.cd(icanvas)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.22)
        pad.SetRightMargin(0.02)
        pad.SetTopMargin(0.13)
        pad.SetBottomMargin(0.15)
        h = hRef.DrawCopy("axis")
        h.GetYaxis().SetRangeUser(0, 0.16)
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

        refIntErr = ROOT.Double(0)
        sigIntErr = ROOT.Double(0)

        refInt = hRef.IntegralAndError(hRef.GetXaxis().FindBin(InvMassRange[0]), hRef.GetXaxis().FindBin(InvMassRange[1]), refIntErr)
        sigInt = hSig.IntegralAndError(hSig.GetXaxis().FindBin(InvMassRange[0]), hSig.GetXaxis().FindBin(InvMassRange[1]), sigIntErr)
        ros = refInt / sigInt
        rosErr = math.sqrt((refIntErr / refInt) ** 2 + (sigIntErr / sigInt) ** 2) * ros

        print("Bin {0}: ros = {1:.3f} #pm {2:.3f}".format(i, ros, rosErr))

        refInt = hRef.IntegralAndError(1, hRef.GetNbinsX(), refIntErr)
        sigInt = hSig.IntegralAndError(1, hSig.GetNbinsX(), sigIntErr)
        ros = refInt / sigInt
        rosErr = math.sqrt((refIntErr / refInt) ** 2 + (sigIntErr / sigInt) ** 2) * ros

        print("Bin {0} (full range): ros = {1:.3f} #pm {2:.3f}".format(i, ros, rosErr))

        rosPave = ROOT.TPaveText(0.25, 0.77, 0.57, 0.87, "NB NDC")
        rosPave.SetBorderSize(0)
        rosPave.SetFillStyle(0)
        rosPave.SetTextFont(43)
        rosPave.SetTextSize(18)
        rosPave.SetTextAlign(13)
        rosPave.AddText("#it{{R}}_{{f}} = {0:.2f} #pm {1:.2f}".format(ros, rosErr))
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

        icanvas += 1

    canvas.cd(1)
    leg = ROOT.TLegend(0.25, 0.61, 0.58, 0.75, "", "NB NDC")
    globalList.append(leg)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(18)
    leg.SetTextAlign(13)
    leg.AddEntry(hRef_copy, "Reflections", "pe")
    leg.AddEntry(hRef_copy.GetListOfFunctions()[0], "Fit", "l")
    leg.Draw()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot reflections.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--plot', nargs='*', type=int)
    parser.add_argument('--var', type=str, default="DPt")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    if args.var == "DPt":
        name = "D0toKpiCuts_DPt_Charged_R040_DPtBins_JetPt_5_30"
    elif args.var == "JetPt":
        name = "D0toKpiCuts_JetPt_Charged_R040_JetPtBins_DPt_30"
    else:
        print("Error: var '{}' not known!".format(args.var))
        exit(1)
    
    if args.plot:
        index_list = args.plot
    else:
        index_list = []

    PlotReflections(yconfig, index_list, name, "DoubleGaus")

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{input_path}/{train}/{ana}/reflTemp/{name}.pdf".format(input_path=yconfig["input_path"], train=yconfig["train"], ana=yconfig["name"], name=obj.GetName()))

    IPython.embed()
