#!/usr/bin/env python
# python script to test Jet D meson analysis

import argparse
import ROOT

def MakeRDHFCuts(dmesons, period, recopass, cent):
    cutlist = []

    for dmeson in dmesons:
        print("Processing {}".format(dmeson))
        if period == "LHC10" and recopass == "pass2" and cent == -1:
            if dmeson == "D0":
                cuts = ROOT.AliRDHFCutsD0toKpi()
            elif dmeson == "DStar":
                cuts = ROOT.AliRDHFCutsDStartoKpipi()
            else:
                print("D meson '{0}' not valid.".format(dmeson))
                return
            cuts.SetStandardCutsPP2010()
        elif period == "LHC10" and recopass == "pass4" and cent == -1:
            if dmeson == "D0":
                ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include")
                ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include")
                ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_pp2010pass4_arXiv_1702_00766_CTerrevoli.C+g")
                cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_arXiv_1702_00766_CTerrevoli(False)
            elif dmeson == "DStar":
                fnameInput = "DStartoKpipiCuts_pp2010pass4_arXiv_1702_00766_AGrelli.root"
                fileInput = ROOT.TFile(fnameInput)
                if not fileInput or fileInput.IsZombie():
                    print("Could not find file {}".format(fnameInput))
                    exit(1)
                cuts = fileInput.Get("DStartoKpipiCuts")
                fileInput.Close()
            else:
                print("D meson '{0}' not valid.".format(dmeson))
                return
        elif period == "LHC15o" and recopass == "pass1" and (cent == 0 or cent == 1):
            if dmeson == "D0":
                ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include")
                ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include")
                ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_PbPb2015_Cent0_10_Raa_WIP_XPeng.C+g")
                cuts = ROOT.MakeD0toKpiCuts_PbPb2015_Cent0_10_Raa_WIP_XPeng(False)
            elif dmeson == "DStar":
                fnameInput = "DStartoKpiCuts_PbPb2015_Cent0_10_Raa_WIP_SJaelani.root"
                fileInput = ROOT.TFile(fnameInput)
                if not fileInput or fileInput.IsZombie():
                    print("Could not find file {}".format(fnameInput))
                    exit(1)
                cuts = fileInput.Get("DStartoKpipiCuts")
                fileInput.Close()
            else:
                print("D meson '{0}' not valid.".format(dmeson))
                return
        elif period == "LHC15o" and recopass == "pass1" and cent == 2:
            if dmeson == "D0":
                ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include")
                ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include")
                ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_PbPb2015_Cent30_50_Raa_QM17_XPeng.C+g")
                cuts = ROOT.MakeD0toKpiCuts_PbPb2015_Cent30_50_Raa_QM17_XPeng(False)
            elif dmeson == "DStar":
                fnameInput = "DStartoKpiCuts_PbPb2015_Cent0_10_Raa_WIP_SJaelani.root"
                fileInput = ROOT.TFile(fnameInput)
                if not fileInput or fileInput.IsZombie():
                    print("Could not find file {}".format(fnameInput))
                    exit(1)
                cuts = fileInput.Get("DStartoKpipiCuts")
                fileInput.Close()
            else:
                print("D meson '{0}' not valid.".format(dmeson))
                return
        else:
            print("Period {0}, pass {1}, cent {2} not valid".format(period, recopass, cent))
            return

        cuts.SetUsePhysicsSelection(False)
        cuts.SetTriggerClass("", "")
        cuts.SetMinCentrality(0)
        cuts.SetMaxCentrality(100)
        if dmeson == "D0":
            cuts.SetName("D0toKpiCuts")
        elif dmeson == "Dstar":
            cuts.SetName("DStartoKpipiCuts")
        cutlist.append(cuts)

    fname = "RDHFCuts_{}_{}".format(period, recopass)
    if cent >= 0: fname += "_Cent{}".format(cent)
    fname += ".root"
    print("Opening file {} for writing".format(fname))
    file = ROOT.TFile(fname, "recreate")
    if not file or file.IsZombie():
        print("Could not create file '{0}'".format(fname))
        return
    file.cd()
    for cuts in cutlist:
        cuts.Write()
    file.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates RDHF cuts.')
    parser.add_argument('dmesons',
                        default="D0", nargs='*',
                        help='D meson (D0, Dstar, Dplus, Ds')
    parser.add_argument('--cent',
                        default=-1, type=int,
                        help='Centrality class')
    parser.add_argument('--period',
                        help='Period (e.g. LHC10)')
    parser.add_argument('--recopass',
                        default="pass1",
                        help='Reco pass (pass1)')
    args = parser.parse_args()

    MakeRDHFCuts(args.dmesons, args.period, args.recopass, args.cent)
