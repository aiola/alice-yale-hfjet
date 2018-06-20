#!/usr/bin/env python

import ROOT
import IPython
import LoadInclusiveJetSpectrum

globalList = []

def main():
    ROOT.TH1.AddDirectory(False)
    LoadInclusiveJetSpectrum.GetCrossSection([5, 15, 30])

if __name__ == '__main__':
    main()

    IPython.embed()
