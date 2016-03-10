import ROOT

def LoadDJetCorrClasses():
    ROOT.gSystem.Load("DJetCorr.so")
    #gROOT.LoadMacro("HistoStyler.cxx+g")
    #gROOT.LoadMacro("MassFitter.cxx+g")
    #gROOT.LoadMacro("DJetCorrAnalysisParams.cxx+g")
    #gROOT.LoadMacro("DJetCorrBase.cxx+g")
    #gROOT.LoadMacro("DJetCorrAnalysis.cxx+g")
    #gROOT.LoadMacro("DJetCorrResponse.cxx+g")
    #gROOT.LoadMacro("DJetCorrUnfold.cxx+g")