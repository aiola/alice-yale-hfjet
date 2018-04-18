#!/usr/bin/env python
# python script analyze the PID cuts

import ROOT
print("ROOT {} imported".format(ROOT.gROOT.GetVersionInt()))
ROOT.PyConfig.IgnoreCommandLineOptions = True

import argparse
import yaml
import IPython
import os

import DMesonJetUtils

globalList = []


def LoadHistograms(input_path, train, file_name, task_name, trigger, d_meson, d_meson_suffix):
    histogram_names = ["fHistNAcceptedDmesonsPt0", "fHistNAcceptedDmesonsPt1", "fHistNAcceptedDmesonsPt2", "fHistNAcceptedDmesonsPt3", "fHistNAcceptedDmesonsPt4",
                       "fHistNTotAcceptedDmesonsPt0", "fHistNTotAcceptedDmesonsPt1", "fHistNTotAcceptedDmesonsPt2", "fHistNTotAcceptedDmesonsPt3", "fHistNTotAcceptedDmesonsPt4"]
    if trigger:
        histogram_path = "{task_name}_{trigger}_histos/histos{task_name}_{trigger}/{d_meson}_{d_meson_suffix}".format(task_name=task_name, trigger=trigger, d_meson=d_meson, d_meson_suffix=d_meson_suffix)
    else:
        histogram_path = "{task_name}_histos/histos{task_name}/{d_meson}_{d_meson_suffix}".format(task_name=task_name, d_meson=d_meson, d_meson_suffix=d_meson_suffix)
    histos = dict()
    for hname in histogram_names:
        histos[hname] = []

    path = "{0}/{1}".format(input_path, train)

    print("Looking for file {0} in path {1}".format(file_name, path))
    file_names = DMesonJetUtils.find_file(path, file_name)

    for fname in file_names:
        file = ROOT.TFile(fname, "read")
        for hname in histogram_names:
            h = DMesonJetUtils.GetObject(file, "{}/{}".format(histogram_path, hname))
            histos[hname].append(h)
    return histos


def MergeHistograms(histogram_list):
    hmerged = histogram_list[0].Clone("{}_Merged".format(histogram_list[0].GetName()))
    for h in histogram_list[1:]:
        hmerged.Add(h)
    return hmerged


def ProcessPID(histos, histos_wrong_pid):
    canvas = ROOT.TCanvas("PIDSummary", "PIDSummary")
    globalList.append(canvas)
    colors = [ROOT.kBlack, ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kOrange + 2]
    legend = ROOT.TLegend(0.15, 0.53, 0.47, 0.80)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    globalList.append(legend)
    table_title = "Minimum \ptd\ (\GeVc) & \Dzero\ candidates & \Dzerobar\ candidates & No-PID candidates"
    if histos_wrong_pid: table_title += "& Wrong-PID candidates"
    table_title += " \\\\ \hline"
    print(table_title)
    for ptcut in range(0, 5):
        histogram_list = histos["fHistNTotAcceptedDmesonsPt{}".format(ptcut)]
        hmerged = MergeHistograms(histogram_list)
        globalList.append(hmerged)
        D0_count = hmerged.GetBinContent(hmerged.GetXaxis().FindBin("D"))
        D0bar_count = hmerged.GetBinContent(hmerged.GetXaxis().FindBin("Anti-D"))
        either_count = hmerged.GetBinContent(hmerged.GetXaxis().FindBin("Both"))
        tot_count = D0_count + D0bar_count + either_count
        D0_frac = D0_count / tot_count
        D0bar_frac = D0bar_count / tot_count
        either_frac = either_count / tot_count
        table_line1 = "${}$ & ${:.0f}$ & ${:.0f}$ & ${:.0f}$".format(ptcut, D0_count, D0bar_count, either_count)
        table_line2 = "fraction: & ${}$ & ${}$ & ${}$".format(sci_notation(D0_frac, 2), sci_notation(D0bar_frac, 2), sci_notation(either_frac, 2))
        hmerged.SetLineColor(colors[ptcut])
        hmerged.SetLineWidth(2)
        hmerged.SetLineStyle(ptcut)
        if ptcut == 0:
            hmerged.Draw()
            hmerged.GetYaxis().SetRangeUser(0, 500000)
        else:
            hmerged.Draw("same")
        legend.AddEntry(hmerged, "#it{{p}}_{{T,D}} > {} GeV/#it{{c}}".format(ptcut))

        if histos_wrong_pid:
            histogram_list_wrong_pid = histos_wrong_pid["fHistNTotAcceptedDmesonsPt{}".format(ptcut)]
            hmerged_wrong_pid = MergeHistograms(histogram_list_wrong_pid)
            wrong_pid_count = hmerged_wrong_pid.GetBinContent(1) + hmerged_wrong_pid.GetBinContent(2)
            wrong_pid_frac = wrong_pid_count / tot_count
            table_line1 += " & ${:.0f}$".format(wrong_pid_count)
            table_line2 += " & ${}$".format(sci_notation(wrong_pid_frac, 2))

        table_line1 += " \\\\"
        table_line2 += " \\\\ \hline"
        print(table_line1)
        print(table_line2)

    legend.Draw()
    return [canvas]


def sci_notation(number, sig_fig):
  ret_string = "{0:.{1:d}e}".format(number, sig_fig)
  a, b = ret_string.split("e")
  b = int(b)  # removed leading "+" and strips leading zeros too.
  return a + " \\times 10^{" + str(b) + "}"


def ProcessDMesonCount(histos):
    canvas = ROOT.TCanvas("DMesonCountSummary", "DMesonCountSummary")
    canvas.SetLogy()
    globalList.append(canvas)
    colors = [ROOT.kBlack, ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kOrange + 2]
    legend = ROOT.TLegend(0.59, 0.57, 0.91, 0.84)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    globalList.append(legend)
    max_cand = 4
    table_title = "Minimum \ptd\ (\GeVc)"
    for n in range(0, max_cand + 1):
        table_title = "{} & {} candidates".format(table_title, n)
    table_title = "{} \\\\ \hline".format(table_title)
    print(table_title)
    for ptcut in range(0, 5):
        histogram_list = histos["fHistNAcceptedDmesonsPt{}".format(ptcut)]
        hmerged = MergeHistograms(histogram_list)
        globalList.append(hmerged)
        n_cand = dict()
        tot_events = 0
        frac_n_cand = dict()
        for n in range(0, max_cand + 1):
            n_cand[n] = hmerged.GetBinContent(n + 1)
            tot_events += hmerged.GetBinContent(n + 1)
        for n in range(0, max_cand + 1):
            frac_n_cand[n] = n_cand[n] / tot_events

        table_line1 = "${}$".format(ptcut)
        table_line2 = "fraction: "
        for n in range(0, max_cand + 1):
            table_line1 += " & ${:.0f}$".format(n_cand[n])
            table_line2 += " & ${}$".format(sci_notation(frac_n_cand[n], 2))

        table_line1 += " \\\\"
        table_line2 += " \\\\ \hline"

        print(table_line1)
        print(table_line2)

        hmerged.SetLineColor(colors[ptcut])
        hmerged.SetLineWidth(2)
        hmerged.SetLineStyle(ptcut)
        if ptcut == 0:
            hmerged.Draw()
            hmerged.GetXaxis().SetRangeUser(-0.5, 5)
            hmerged.GetYaxis().SetRangeUser(1e-1, 1e9)
        else:
            hmerged.Draw("same")
        legend.AddEntry(hmerged, "#it{{p}}_{{T,D}} > {} GeV/#it{{c}}".format(ptcut))

    legend.Draw()
    return [canvas]


def main(config, format):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    histos = LoadHistograms(config["input_path"], config["train"], config["file_name"], config["task_name"], config["analysis"][0]["trigger"][0], "D0", "D0toKpiCuts")

    if config["monte_carlo"]:
        histos_wrong_pid = LoadHistograms(config["input_path"], config["train"], config["file_name"], config["task_name"], config["analysis"][0]["trigger"][0], "D0_OnlyWrongPIDAccepted", "D0toKpiCuts")
    else:
        histos_wrong_pid = None

    if not histos or len(histos) == 0:
        print("No histograms were loaded!")
        exit(1)

    canvases = []
    canvases.extend(ProcessPID(histos, histos_wrong_pid))
    canvases.extend(ProcessDMesonCount(histos))

    output_path = "{}/{}/CountDMesons".format(config["input_path"], config["train"])
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for c in canvases:
        c.SaveAs("{}/{}.{}".format(output_path, c.GetName(), format))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analysis of the topological cuts.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.format)

    IPython.embed()
