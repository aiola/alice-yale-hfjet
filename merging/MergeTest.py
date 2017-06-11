#!/usr/local/bin/python

import argparse
import ROOT
import MergeFiles
import subprocess
import os

def MergeTest(fileList, destPath):
    ROOT.gSystem.Load("libCGAL")
    
    skipList = []
    f = open(fileList, 'r')
    i = 0
    mergeList = []
    for line in f:
        line = line.rstrip()
        line = line.lstrip()
        print line
        dest = line
        if line.startswith("alien://"):
            dest = "{0}/AnalysisResults_{1}.root".format(destPath, i)
            if not os.path.exists(dest):
                subprocess.call(["alien_cp", line, dest])
        mergeList.append(dest)
        i += 1
    output = "{0}/AnalysisResults.root".format(destPath, i)
    MergeFiles.MergeFiles(output, mergeList, skipList, 100)
            
if __name__ == '__main__':
    # FinalMergeLocal.py executed as script
    
    parser = argparse.ArgumentParser(description='Local final merging for MC production in pT hard bins.')
    parser.add_argument('fileList', metavar='fileList',
                        help='File list to be merged. If names start with alien:// it will first download them for alien')
    parser.add_argument('destPath', metavar='destPath',
                        help='Destination path')
    
    args = parser.parse_args()

    MergeTest(args.fileList, args.destPath)
