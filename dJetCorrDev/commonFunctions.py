import ROOT
import subprocess

def LoadDJetCorrClasses():
    subprocess.call("make")
    ROOT.gSystem.Load("DJetCorr.so")
    