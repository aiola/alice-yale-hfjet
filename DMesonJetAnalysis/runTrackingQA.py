#!/usr/bin/env python
# python script to project TTree containing inclusive jet spectra

import argparse
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import numpy
import yaml
from enum import Enum

globalList = []


def ResetTHnFilters(hn):
    for iaxis in range(0, hn.GetNdimensions()):
        hn.GetAxis(iaxis).SetRange(0, -1)


class Variable:

    def __init__(self, name, title, units, logy, tracks_axis, particles_axis, matched_tracks_axis, matched_particles_axis):
        self.fName = name
        self.fTitle = title
        self.fUnits = units
        self.fTracksAxis = tracks_axis
        self.fParticlesAxis = particles_axis
        self.fMatchedTracksAxis = matched_tracks_axis
        self.fMatchedParticlesAxis = matched_particles_axis
        self.fLogy = logy


class Histogram:

    def __init__(self, variable, matched):
        self.fVariable = variable
        self.fMatched = matched
        self.fHistograms = dict()
        self.fNormHistograms = dict()
        self.fTrackTypes = { 1 : "Global tracks", 2 : "Constrained tracks w/ ITS refit", 3 : "Constrained tracks w/o ITS refit"}

    def MergeHistograms(self):
        if len(self.fHistograms) < 2: return
        htot = None
        for h in self.fHistograms.itervalues():
            if htot:
                htot.Add(h)
            else:
                htot = h.Clone(h.GetName()[:-2])
        self.fHistograms["all"] = htot

    def NormalizeHistogram(self, events):
        for tt, h in self.fHistograms.iteritems():
            hnorm = h.Clone("{}_Normalized".format(h.GetName()))
            hnorm.Scale(1.0 / events, "width")
            if self.fVariable.fUnits:
                yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}} {units}^{{-1}}".format(var=self.fVariable.fTitle, units=self.fVariable.fUnits)
            else:
                yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fVariable.fTitle)
            self.fNormHistograms[tt] = hnorm

    def GetFullHistogram(self, normalized=True):
        if normalized:
            return self.fNormHistograms["all"]
        else:
            return self.fHistograms["all"]

    def GetPartialHistograms(self, normalized=True):
        if normalized:
            l = self.fNormHistograms
        else:
            l = self.fHistograms
        return [obj for k, obj in l.iteritems() if not k == "all"]


class TrackHistogram(Histogram):

    def Add(self, hn, track_type_axis):
        print("Adding histogram '{}' (variable '{}')...".format(hn.GetName(), self.fVariable.fName))
        for track_type, track_type_title in self.fTrackTypes.iteritems():
            hn.GetAxis(track_type_axis).SetRange(track_type, track_type)
            if self.fMatched:
                h = hn.Projection(self.fVariable.fMatchedTracksAxis, "A")
                hname = "MatchedTracks{}_{}".format(self.fVariable.fName, track_type)
            else:
                h = hn.Projection(self.fVariable.fTracksAxis, "A")
                hname = "Tracks{}_{}".format(self.fVariable.fName, track_type)
            hn.GetAxis(track_type_axis).SetRange(0, -1)
            xaxis_title = self.fVariable.fTitle
            if self.fVariable.fUnits: xaxis_title += " " + self.fVariable.fUnits
            h.SetName(hname)
            h.SetTitle(track_type_title)
            h.GetXaxis().SetTitle(xaxis_title)
            h.GetYaxis().SetTitle("counts")
            h.Sumw2()
            if track_type in self.fHistograms:
                self.fHistograms[track_type].Add(h)
            else:
                self.fHistograms[track_type] = h


class ParticleHistogram(Histogram):

    def Add(self, hn, track_type_axis):
        if self.fMatched:
            track_types = self.fTrackTypes
        else:
            track_types = { "all" : "Hybrid tracks" }
        for track_type, track_type_title in track_types.iteritems():
            if self.fMatched:
                hn.GetAxis(track_type_axis).SetRange(track_type, track_type)
                h = hn.Projection(self.fVariable.fMatchedParticlesAxis, "A")
                hn.GetAxis(track_type_axis).SetRange(0, -1)
                hname = "MatchedParticles{}_{}".format(self.fVariable.fName, track_type)
            else:
                hn.GetAxis(track_type_axis).SetRange(2, 2)  # findable
                h = hn.Projection(self.fVariable.fParticlesAxis, "A")
                hn.GetAxis(track_type_axis).SetRange(0, -1)
                hname = "Particles{}".format(self.fVariable.fName)
            xaxis_title = self.fVariable.fTitle
            if self.fVariable.fUnits: xaxis_title += " " + self.fVariable.fUnits
            h.SetName(hname)
            h.SetTitle(track_type_title)
            h.GetXaxis().SetTitle(xaxis_title)
            h.GetYaxis().SetTitle("counts")
            h.Sumw2()
            if track_type in self.fHistograms:
                self.fHistograms[track_type].Add(h)
            else:
                self.fHistograms[track_type] = h


class TrackingQA:

    def __init__(self, input_path, train):
        self.fInputPath = input_path
        self.fTrain = train
        self.fTaskName = "AliEmcalTrackingQATask_histos"
        self.fHnTracksName = "fTracks"
        self.fHnParticlesName = "fParticlesPhysPrim"
        self.fHnMatchedName = "fParticlesMatched"
        self.fVariables = []
        self.fVariables.append(Variable("Pt", "#it{p}_{T}", "(GeV/#it{c})", True, 0, 0, 3, 0))
        self.fVariables.append(Variable("Eta", "#eta", "", True, 1, 1, 4, 1))
        self.fVariables.append(Variable("Phi", "#phi", "", True, 2, 2, 5, 2))
        self.fTrackTypeAxis = 4
        self.fTrackMomResAxis = 5
        self.fParticleFindableAxis = 4
        self.MatchedTrackTypeAxis = 6
        self.MatchedTrackMomRelDiffAxis = 7
        self.fTrackHistograms = dict()
        self.fParticleHistograms = dict()
        self.fMatchedTracksHistograms = dict()
        self.fMatchedParticleHistograms = dict()
        self.fHnTracks = []
        self.fHnParticles = []
        self.fHnMatched = []
        self.fNEvents = []
        self.fTotEvents = 0

    def ExtractPeriod(self, fname):
        train_index = fname.find(self.fTrain)
        index1 = fname.find("/", train_index)
        index2 = fname.find("/", index1 + 1)
        period = fname[index1 + 1: index2]
        return period

    def LoadHns(self):
        file_name = "AnalysisResults.root"
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)
        print("Looking for file {0} in path {1}".format(file_name, path))
        file_names = DMesonJetUtils.find_file(path, file_name)

        for fname in file_names:
            print("Opening file '{}'".format(fname))
            period = self.ExtractPeriod(fname)

            file = ROOT.TFile(fname, "read")

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnTracksName))
            h.SetName("{}_{}".format(h.GetName(), period))
            self.fHnTracks.append(h)

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnParticlesName))
            h.SetName("{}_{}".format(h.GetName(), period))
            self.fHnParticles.append(h)

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnMatchedName))
            h.SetName("{}_{}".format(h.GetName(), period))
            self.fHnMatched.append(h)

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, "fHistEventCount"))
            events = h.GetBinContent(1)
            self.fNEvents.append(events)
            self.fTotEvents += events

    def ProjectHistograms(self):
        for var in self.fVariables:
            track_histogram = TrackHistogram(var, False)
            particle_histogram = ParticleHistogram(var, False)
            matched_track_histogram = TrackHistogram(var, True)
            matched_particle_histogram = ParticleHistogram(var, True)
            for hn in self.fHnTracks:
                track_histogram.Add(hn, self.fTrackTypeAxis)
            for hn in self.fHnParticles:
                particle_histogram.Add(hn, self.fParticleFindableAxis)
            for hn in self.fHnMatched:
                matched_track_histogram.Add(hn, self.MatchedTrackTypeAxis)
                matched_particle_histogram.Add(hn, self.MatchedTrackTypeAxis)
            track_histogram.MergeHistograms()
            particle_histogram.MergeHistograms()
            matched_track_histogram.MergeHistograms()
            matched_particle_histogram.MergeHistograms()
            track_histogram.NormalizeHistogram(self.fTotEvents)
            particle_histogram.NormalizeHistogram(self.fTotEvents)
            matched_track_histogram.NormalizeHistogram(self.fTotEvents)
            matched_particle_histogram.NormalizeHistogram(self.fTotEvents)
            self.fTrackHistograms[var.fName] = track_histogram
            self.fParticleHistograms[var.fName] = particle_histogram
            self.fMatchedTracksHistograms[var.fName] = matched_track_histogram
            self.fMatchedParticleHistograms[var.fName] = matched_particle_histogram

    def Plot(self):
        self.PlotSpectra("Tracks", self.fTrackHistograms)
        self.PlotSpectra("MatchedTracks", self.fMatchedTracksHistograms)
        self.PlotSpectra("MatchedParticles", self.fMatchedParticleHistograms)

    def PlotSpectra(self, name, histograms):
        for var in self.fVariables:
            comp_name = "{}_{}".format(name, var.fName)
            comp = DMesonJetCompare.DMesonJetCompare(comp_name)
            comp.fDoRatioPlot = "logy"
            if var.fLogy:
                comp.fDoSpectraPlot = "logy"
            else:
                comp.fDoSpectraPlot = "lineary"
            results = comp.CompareSpectra(histograms[var.fName].GetFullHistogram(), histograms[var.fName].GetPartialHistograms())
            for r in results: globalList.append(r)


def main(input_path, train):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    qa = TrackingQA(input_path, train)
    globalList.append(qa)
    qa.LoadHns()
    qa.ProjectHistograms()
    qa.Plot()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking QA.')
    parser.add_argument('train', metavar='train')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Volumes/DATA/ALICE/JetResults")
    args = parser.parse_args()

    main(args.input_path, args.train)

    IPython.embed()
