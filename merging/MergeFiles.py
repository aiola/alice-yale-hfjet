#!/usr/bin/env python

from ROOT import TFileMerger

def MergeFiles(output, fileList, skipList=[], acceptList=[], n=20):
    merger = TFileMerger(False)
    merger.OutputFile(output);
    merger.SetMaxOpenedFiles(n);

    print "Total number of files is {0}".format(len(fileList))
    
    for fileName in fileList:
        print "Adding file {0}".format(fileName)
        merger.AddFile(fileName)
        
    mode = TFileMerger.kAllIncremental

    if len(skipList) > 0:
        mode = mode | TFileMerger.kSkipListed
        if (len(acceptList) > 0):
            print("Accept list is being ignored!!!")
        for skipObject in skipList:
            merger.AddObjectNames(skipObject)
            
    elif len(acceptList) > 0:
        mode = mode | TFileMerger.kAcceptListed
        for acceptObject in acceptList:
            merger.AddObjectNames(acceptObject)
  
    merger.PrintFiles("");
    r = merger.PartialMerge(mode);

    if not r:
        print "Merge error!"

    return r
