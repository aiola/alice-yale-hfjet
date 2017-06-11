#!/usr/local/bin/python
# python script to run a simple decay simulator

import argparse
import ROOT
import numpy

globalList = []

def main(ptMin, ptMax, n, fname):
  ROOT.gROOT.LoadMacro("DecaySimulator.cxx+g")

  simulator = ROOT.DecaySimulator()
  if ptMin < ptMax:
      ptBins = range(ptMin, ptMax + 1)
      if ptBins[0] == 0: ptBins[0] = 0.5
      print("Simulating the following pt bins: {}".format(ptBins))
      simulator.SetPtBins(len(ptBins), numpy.array(ptBins, dtype=numpy.float64))
  simulator.StartSimulation(n)
  simulator.SaveTree(fname)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simple decay simulator.')
    parser.add_argument('--pt-min', metavar='PTMIN',
                        default=0, type=int,
                        help='Minimum pt of the simulation')
    parser.add_argument('--pt-max', metavar='PTMAX',
                        default=0, type=int,
                        help='Maximum pt of the simulation')
    parser.add_argument('-n', metavar='N',
                        default=10000,
                        help='Number of decays')
    parser.add_argument('-o', metavar='simulation.root',
                        default='simulation.root',
                        help='ROOT file containing the results of the decay simulator')
    args = parser.parse_args()

    main(args.pt_min, args.pt_max, args.n, args.o)
