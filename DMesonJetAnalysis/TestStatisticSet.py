#!/usr/bin/env python
#python script to test the StatisticSet class

import ROOT
import IPython
import DMesonJetUtils

def TestStatisticSet():
    stat = DMesonJetUtils.StatisticSet()
    stat.Fill(0.010,1.5)
    stat.Fill(0.06,1.5)
    stat.Fill(0.08,1.5)
    stat.PrintSummary()

if __name__ == '__main__':

    TestStatisticSet()
    
    IPython.embed()
