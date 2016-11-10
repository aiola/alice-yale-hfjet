#!/usr/bin/env python
#python script to plot FONLL calculations

import argparse
import IPython
import ROOT
import array

globalList = []

class FONLL:
    def __init__(self, fonll_file):
        self.name = fonll_file
        self.n = 0
        self.values = dict()
        self.OpenFONLL(fonll_file)

    def OpenFONLL(self, fonll_file):
        with open(fonll_file, "r") as myfile:
            lines = myfile.read().splitlines()

        vars = None
        prevLine = None
        for line in lines:
            line = line.strip()
            print(line)
            if line[0] == "#":
                prevLine = line
                continue
            if not vars:
                prevLine = "  {0}".format(prevLine[1:])
                vars = [prevLine[i:11+i].strip().strip("#") for i in range(0, len(prevLine), 11)]
                print("Variables taken from {0}".format(vars))
                for var in vars:
                    self.values[var] = []
            valuesStr = line.split(" ")
            values = [float(v) for v in valuesStr]
            self.n += 1
            for name,val in zip(vars, values):
                self.values[name].append(val)
            prevLine = line
    
    def GenerateGraph(self, x, y, yerrup, yerrdown):
        g = ROOT.TGraphAsymmErrors(self.n, array.array('d', self.values[x]), array.array('d', self.values[y]), ROOT.nullptr, ROOT.nullptr, 
                                    array.array('d', [v-down for down,v in zip(self.values[yerrdown], self.values[y])]),
                                    array.array('d', [up-v for up,v in zip(self.values[yerrup], self.values[y])]))
        return g

def main(fonll_file, data, spectrum):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    ntuple = FONLL(fonll_file)
    g = ntuple.GenerateGraph("pt", "central", "max", "min")
    globalList.append(g)

    c = ROOT.TCanvas(fonll_file,fonll_file)
    globalList.append(c)
    c.cd()
    g.Draw()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Comparison between FONLL and POWHEG.')
    parser.add_argument('fonll_file', metavar='fonll.dat',
                        help='Data file with FONLL calculations')
    parser.add_argument('--data', metavar='data.root',
                        default=None)
    parser.add_argument('--spectrum', metavar='spectrum',
                        default="D0_MCTruth_D_Pt_Spectrum_Normalized")
    args = parser.parse_args()

    main(args.fonll_file, args.data, args.spectrum)

    IPython.embed()
