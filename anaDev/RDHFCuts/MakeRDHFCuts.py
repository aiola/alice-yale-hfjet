#!/usr/bin/env python
# python script to test Jet D meson analysis

import argparse
import ROOT

def MakeRDHFCuts(fname, dmesons, period, recopass):
    cutlist = []

    for dmeson in dmesons:
        print("Processing {}".format(dmeson))
        if period == "LHC10" and recopass == "pass2":
            if dmeson == "D0":
                cuts = ROOT.AliRDHFCutsD0toKpi()
            elif dmeson == "Dstar":
                cuts = ROOT.AliRDHFCutsDStartoKpipi()
            else:
                print("D meson '{0}' not valid.".format(dmeson))
                return
            cuts.SetStandardCutsPP2010()
            cuts.GetPidHF().SetOldPid(False)
        elif period == "LHC10" and recopass == "pass4":
            if dmeson == "D0":
                ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include")
                ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include")
                ROOT.gROOT.LoadMacro("Make2010pp_pass4_cuts.C+g")
                cuts = ROOT.Make2010pp_pass4_cuts(False)
            elif dmeson == "Dstar":
                fnameInput = "DStartoKpipiCuts_pp2010pass4_AGrelli_20170314.root"
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
            print("Period {0}, pass {1} not valid".format(period, recopass))
            return

        cuts.SetUsePhysicsSelection(False)
        cuts.SetTriggerClass("", "")
        if dmeson == "D0":
            cuts.SetName("D0toKpiCuts")
        elif dmeson == "Dstar":
            cuts.SetName("DStartoKpipiCuts")
        cutlist.append(cuts)

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
    parser.add_argument('-o',
                        default="cuts.root",
                        help='Output file name')
    parser.add_argument('--period',
                        default="LHC10",
                        help='Period (LHC10)')
    parser.add_argument('--recopass',
                        default="pass2",
                        help='Reco pass (pass2)')
    args = parser.parse_args()

    MakeRDHFCuts(args.o, args.dmesons, args.period, args.recopass)
