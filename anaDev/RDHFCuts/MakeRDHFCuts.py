#!/usr/bin/env python
# python script to test Jet D meson analysis

import argparse
import ROOT
import array

def MakeRDHFCuts(dmesons, period, recopass, cent):
    ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include")
    ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include")

    ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_pp2010pass4_arXiv_1702_00766_CTerrevoli.C+g")
    ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_pp2010pass4_Variations.C+g")
    ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_PbPb2015_Cent0_10_Raa_WIP_XPeng.C+g")
    ROOT.gROOT.LoadMacro("MakeD0toKpiCuts_PbPb2015_Cent30_50_Raa_QM17_XPeng.C+g")

    cutlist = []
    for dmeson in dmesons:
        if dmeson == "D0" and period == "LHC10" and recopass == "pass4":
            cut_strenghts = ["standard", "loosest_nopid", "loosest_pid", "LoosePointingLoosed0d0", "LoosePointing", "Loosed0d0",
                             "TopoOnlyNSigma1", "TopoOnlyNSigma2", "TopoOnlyNSigma3", "TopoOnlyNSigma4",
                             "D0JetOptimLowJetPtv1", "D0JetOptimHighJetPtv1", 
                             "D0JetOptimLowJetPtv2", "D0JetOptimHighJetPtv2",
                             "D0JetOptimLowJetPtv3", "D0JetOptimHighJetPtv3"]
        else:
            cut_strenghts = ["standard"]
        for cut_strenght in cut_strenghts:
            print("Processing {} {}".format(dmeson, cut_strenght))
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
                    cuts = MakeD0toKpiCuts_pp2010pass4(cut_strenght)
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
            elif period == "LHC16k" and recopass == "pass1" and cent == -1:
                if dmeson == "D0":
                    fnameInput = "D0toKpiCuts_pp2016_SQM17_SCostanza.root"
                    fileInput = ROOT.TFile(fnameInput)
                    if not fileInput or fileInput.IsZombie():
                        print("Could not find file {}".format(fnameInput))
                        exit(1)
                    cuts = fileInput.Get("D0toKpiCuts")
                    fileInput.Close()
                else:
                    print("D meson '{0}' not valid.".format(dmeson))
                    return
            elif period == "LHC12a" and recopass == "pass2" and cent == -1:
                if dmeson == "D0":
                    fnameInput = "D0toKpiCuts_pp2012_SQM17_SCostanza.root"
                    fileInput = ROOT.TFile(fnameInput)
                    if not fileInput or fileInput.IsZombie():
                        print("Could not find file {}".format(fnameInput))
                        exit(1)
                    cuts = fileInput.Get("D0toKpiCuts")
                    fileInput.Close()
                else:
                    print("D meson '{0}' not valid.".format(dmeson))
                    return
            elif period == "LHC15o" and recopass == "pass1" and (cent == 0 or cent == 1):
                if dmeson == "D0":
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
                    cuts = ROOT.MakeD0toKpiCuts_PbPb2015_Cent30_50_Raa_QM17_XPeng(False)
                elif dmeson == "DStar":
                    fnameInput = "DStartoKpiCuts_PbPb2015_Cent30_50_Raa_QM17_SJaelani.root"
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

            if dmeson == "D0":
                cuts_name = "D0toKpiCuts"
            elif dmeson == "DStar":
                cuts_name = "DStartoKpipiCuts"
            if not cut_strenght == "standard":
                cuts_name += "_{}".format(cut_strenght)

            cuts.SetUsePhysicsSelection(False)
            cuts.SetTriggerClass("", "")
            cuts.SetMinCentrality(0)
            cuts.SetMaxCentrality(100)
            cuts.SetName(cuts_name)
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

def MakeD0toKpiCuts_pp2010pass4(cut_strenght):
    esdTrackCuts = ROOT.AliESDtrackCuts("AliESDtrackCuts", "default")
    esdTrackCuts.SetRequireSigmaToVertex(False)  # not present in the AOD filtering
    esdTrackCuts.SetRequireTPCRefit(True)
    esdTrackCuts.SetMinNClustersTPC(70)  # 50 in the AOD filtering
    esdTrackCuts.SetRequireITSRefit(True)
    esdTrackCuts.SetClusterRequirementITS(ROOT.AliESDtrackCuts.kSPD, ROOT.AliESDtrackCuts.kAny)
    esdTrackCuts.SetMinDCAToVertexXY(0.)
    esdTrackCuts.SetPtRange(0.3, 1.e10)
    esdTrackCuts.SetEtaRange(-0.8, +0.8)
    esdTrackCuts.SetMinRatioCrossedRowsOverFindableClustersTPC(0.8)  # not present in the AOD filtering

    pidObj = ROOT.AliAODPidHF()
    mode = 1
    plims = array.array('d', [0.6, 0.8])  # TPC limits in momentum [GeV/c]
    compat = True  # effective only for this mode
    asym = True
    sigmas = array.array('d', [2., 1., 0., 3., 0.])  # to be checked and to be modified with new implementation of setters by Rossella
    pidObj.SetAsym(asym)  # if you want to use the asymmetric bands in TPC
    pidObj.SetMatch(mode)
    pidObj.SetPLimit(plims, len(plims))
    pidObj.SetSigma(sigmas)
    pidObj.SetCompat(compat)
    pidObj.SetTPC(True)
    pidObj.SetTOF(True)
    pidObj.SetPCompatTOF(1.5)
    pidObj.SetSigmaForTPCCompat(3.)
    pidObj.SetSigmaForTOFCompat(3.)
    pidObj.SetOldPid(False)

    if cut_strenght == "standard":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_arXiv_1702_00766_CTerrevoli(False)
    elif cut_strenght == "loosest_nopid":
        # Same as AOD filtering in AliPhysics/PWGHF/vertexingHF/ConfigVertexingHF.C
        cuts = ROOT.AliRDHFCutsD0toKpi()
        cuts.SetName("D0toKpiCuts")
        cuts.SetTitle("Cuts for D0 analysis (loosest)")
        cuts.SetSelectCandTrackSPDFirst(True, 3.)
        # PILE UP REJECTION
        cuts.SetOptPileup(ROOT.AliRDHFCuts.kRejectPileupEvent)
        # EVENT CUTS
        cuts.SetMinVtxContr(1)
        # MAX Z-VERTEX CUT
        cuts.SetMaxVtxZ(10.)

        cutsArrayD0toKpi = array.array('f', [0.3, 999999., 1.1, 0., 0., 999999., 999999., 999999., 0., -1, 0.])
        cuts.SetCuts(11, cutsArrayD0toKpi)
        cuts.AddTrackCuts(esdTrackCuts)

        # No PID cuts in the AOD filtering
        cuts.SetUsePID(False)

        # Do not recalculate the vertex
        cuts.SetRemoveDaughtersFromPrim(True)  # activate for pp
        cuts.SetUseDefaultPID(False)
        cuts.SetLowPt(False)
        cuts.PrintAll()
    elif cut_strenght == "loosest_pid":
        # Same as AOD filtering in AliPhysics/PWGHF/vertexingHF/ConfigVertexingHF.C
        cuts = ROOT.AliRDHFCutsD0toKpi()
        cuts.SetName("D0toKpiCuts")
        cuts.SetTitle("Cuts for D0 analysis (loosest)")
        cuts.SetSelectCandTrackSPDFirst(True, 3.)
        # PILE UP REJECTION
        cuts.SetOptPileup(ROOT.AliRDHFCuts.kRejectPileupEvent)
        # EVENT CUTS
        cuts.SetMinVtxContr(1)
        # MAX Z-VERTEX CUT
        cuts.SetMaxVtxZ(10.)

        cutsArrayD0toKpi = array.array('f', [0.3, 999999., 1.1, 0., 0., 999999., 999999., 999999., 0., -1, 0.])
        cuts.SetCuts(11, cutsArrayD0toKpi)
        cuts.AddTrackCuts(esdTrackCuts)

        # No PID cuts in the AOD filtering, following comes from standard pass4 cuts (C. Terrevoli)
        cuts.SetUsePID(True)
        cuts.SetPidHF(pidObj)

        # Do not recalculate the vertex
        cuts.SetRemoveDaughtersFromPrim(True)  # activate for pp
        cuts.SetUseDefaultPID(False)
        cuts.SetLowPt(False)
        cuts.PrintAll()
    elif cut_strenght == "LoosePointingLoosed0d0":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kLoosePointingLoosed0d0, True)
    elif cut_strenght == "Loosed0d0":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kLoosed0d0, True)
    elif cut_strenght == "LoosePointing":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kLoosePointing, True)
    elif cut_strenght == "TopoOnlyNSigma1":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kTopoOnlyNSigma1, False)
    elif cut_strenght == "TopoOnlyNSigma2":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kTopoOnlyNSigma2, False)
    elif cut_strenght == "TopoOnlyNSigma3":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kTopoOnlyNSigma3, False)
    elif cut_strenght == "TopoOnlyNSigma4":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kTopoOnlyNSigma4, False)
    elif cut_strenght == "D0JetOptimLowJetPtv1":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimLowJetPtv1, True)
    elif cut_strenght == "D0JetOptimHighJetPtv1":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimHighJetPtv1, True)
    elif cut_strenght == "D0JetOptimLowJetPtv2":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimLowJetPtv2, True)
    elif cut_strenght == "D0JetOptimHighJetPtv2":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimHighJetPtv2, True)
    elif cut_strenght == "D0JetOptimLowJetPtv3":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimLowJetPtv3, True)
    elif cut_strenght == "D0JetOptimHighJetPtv3":
        cuts = ROOT.MakeD0toKpiCuts_pp2010pass4_Variations(ROOT.kD0JetOptimHighJetPtv3, True)
    else:
        print("Cut '{}' not known!".format(cut_strenght))
        cuts = None
    return cuts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates RDHF cuts.')
    parser.add_argument('dmesons',
                        default="D0", nargs='*',
                        help='D meson (D0, DStar, DPlus, Ds)')
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
