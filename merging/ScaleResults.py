#!/usr/local/bin/python

import ROOT

from ROOT import gROOT

def ScaleResults(localPath, ptHard, dataset, type="", exclude=""):
    gROOT.LoadMacro("ScaleResults.C+g")
    ROOT.ScaleResults(localPath, ptHard, dataset, type, exclude)
