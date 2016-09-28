#!/usr/bin/env python
#python program to perform a D meson jet unfolding

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
import collections

globalList = []

class DMesonJetUnfoldingEngine:
    def __init__(self):
        return
    
class DMesonJetUnfolding:
    def __init__(self, name):
        self.fName = name
        self.fUnfoldingEngine = []
        self.fCanvases = []
        
    def StartUnfolding(self, config):
        return
