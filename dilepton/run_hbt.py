import os
import sys
import shutil
import math
import numpy as np
import yaml
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory, TH1F, TH2F, THnSparseD, TMath, TF1
import sys
sys.path.append("../common/");
from hbt_utils import get_cf_1d, get_cf_3d

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

config = "";
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
print(config)
system = config["common"]["system"];
energy = config["common"]["energy"]; #float
suffix = args.suffix;

period = "";
is_mc = False;
if "mc" in args.type:
    is_mc = True;
    period = config["period_mc"];
else:
    is_mc = False;
    period = config["period_data"];
pass_number = config["pass_number"];
#__________________________________________________________________
def run_hbt_1d(filename, tasks, arr_qinv, arr_qabs, arr_kt):
    print(sys._getframe().f_code.co_name);
    print("filename = ", filename);
    print("arr_qinv = ", arr_qinv);
    print("arr_qabs = ", arr_qabs);
    print("arr_kt = ", arr_kt);
    __delta = 1e-3;

    rootfile = TFile.Open(filename, "READ");
    outname = "";
    if "pp" in system:
        outname = "photon_hbt_1d_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    else:
        outname = "photon_hbt_1d_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)

    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    nkt = len(arr_kt);
    for itask in range(0, len(tasks)):
        taskname = tasks[itask]["name"];
        cent_min = tasks[itask]["cent_min"];
        cent_max = tasks[itask]["cent_min"];
        print("analyzing ... ", taskname, cent_min, cent_max); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        #rootdir_event.ls();
        #rootdir_pair .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();

        hs_1d_same  = rootdir_pair.Get("same/hs_1d") .Clone("hs_1d_same")
        hs_1d_mix   = rootdir_pair.Get("mix/hs_1d")  .Clone("hs_1d_mix")
        hs_1d_same .Sumw2();
        hs_1d_mix  .Sumw2();
        outdir.WriteTObject(hs_1d_same);
        outdir.WriteTObject(hs_1d_mix);
        #hDeltaEtaDeltaPhi_same = rootdir_pair.Get("same/hDeltaEtaDeltaPhi") .Clone("hDeltaEtaDeltaPhi_same")
        #hDeltaEtaDeltaPhi_mix  = rootdir_pair.Get("mix/hDeltaEtaDeltaPhi") .Clone("hDeltaEtaDeltaPhi_mix")
        #outdir.WriteTObject(hDeltaEtaDeltaPhi_same);
        #outdir.WriteTObject(hDeltaEtaDeltaPhi_mix);

        h1lambda_qinv = TH1F("h1lambda_qinv", "correlation strength", nkt-1, arr_kt);
        h1lambda_qinv.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1lambda_qinv.SetYTitle("#it{#lambda}_{inv}");
        h1lambda_qabs = TH1F("h1lambda_qabs", "correlation strength", nkt-1, arr_kt);
        h1lambda_qabs.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1lambda_qabs.SetYTitle("#it{#lambda}_{LCMS}");

        h1R_qinv = TH1F("h1R_qinv", "radius of homogeneous region", nkt-1, arr_kt);
        h1R_qinv.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R_qinv.SetYTitle("#it{R}_{inv} (fm)");
        h1R_qabs = TH1F("h1R_qabs", "radius of homogeneous region", nkt-1, arr_kt);
        h1R_qabs.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R_qabs.SetYTitle("#it{R}_{LCMS} (fm)");

        for ikt in range(0, nkt-1):
            kt1 = arr_kt[ikt];
            kt2 = arr_kt[ikt+1];
            bin1 = hs_1d_same.GetAxis(1).FindBin(kt1 + __delta);
            bin2 = hs_1d_same.GetAxis(1).FindBin(kt2 - __delta);
            hs_1d_same.GetAxis(1).SetRange(bin1, bin2);
            hs_1d_mix .GetAxis(1).SetRange(bin1, bin2);

            h1qinv_same = hs_1d_same.Projection(0);
            h1qinv_same.SetName("h1qinv_same_kt{0:d}".format(ikt));
            h1qinv_same.SetTitle("h1qinv_same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1qinv_mix  = hs_1d_mix.Projection(0);
            h1qinv_mix.SetName("h1qinv_mix_kt{0:d}".format(ikt));
            h1qinv_mix.SetTitle("h1qinv_mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1qabs_same = hs_1d_same.Projection(0);
            h1qabs_same.SetName("h1qabs_same_kt{0:d}".format(ikt));
            h1qabs_same.SetTitle("h1qabs_same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1qabs_mix  = hs_1d_mix.Projection(0);
            h1qabs_mix.SetName("h1qabs_mix_kt{0:d}".format(ikt));
            h1qabs_mix.SetTitle("h1qabs_mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));

            #h1qinv_same.RebinX(2);
            #h1qinv_mix .RebinX(2);
            #h1qabs_same.RebinX(2);
            #h1qabs_mix .RebinX(2);

            if h1qinv_same.GetEntries() < 0.5 or h1qabs_same.GetEntries() < 0.5 or h1qinv_mix.GetEntries() < 0.5 or h1qabs_mix.GetEntries() < 0.5:
                continue;

            [h1qinv_same, h1qinv_mix, h1cf_qinv, f1cf_qinv] = get_cf_1d(h1qinv_same, h1qinv_mix, 0.0, 0.03);
            h1qinv_same.SetName("h1qinv_same_kt{0:d}".format(ikt));
            h1qinv_same.SetTitle("h1qinv_same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1qinv_mix.SetName("h1qinv_mix_kt{0:d}".format(ikt));
            h1qinv_mix.SetTitle("h1qinv_mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1cf_qinv.SetName("h1cf_qinv_kt{0:d}".format(ikt));
            h1cf_qinv.SetTitle("h1cf_qinv: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            f1cf_qinv.SetName("f1cf_qinv_kt{0}".format(ikt));
            outdir.WriteTObject(h1qinv_same);
            outdir.WriteTObject(h1qinv_mix);
            outdir.WriteTObject(h1cf_qinv);
            outdir.WriteTObject(f1cf_qinv);

            [h1qabs_same, h1qabs_mix, h1cf_qabs, f1cf_qabs] = get_cf_1d(h1qabs_same, h1qabs_mix, 0.0, 0.06);
            h1qabs_same.SetName("h1qabs_same_kt{0:d}".format(ikt));
            h1qabs_same.SetTitle("h1qabs_same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1qabs_mix.SetName("h1qabs_mix_kt{0:d}".format(ikt));
            h1qabs_mix.SetTitle("h1qabs_mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1cf_qabs.SetName("h1cf_qabs_kt{0:d}".format(ikt));
            h1cf_qabs.SetTitle("h1cf_qabs: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            f1cf_qabs.SetName("f1cf_qabs_kt{0}".format(ikt));
            outdir.WriteTObject(h1qabs_same);
            outdir.WriteTObject(h1qabs_mix);
            outdir.WriteTObject(h1cf_qabs);
            outdir.WriteTObject(f1cf_qabs);

            h1lambda_qinv.SetBinContent(ikt+1, f1cf_qinv.GetParameter(0));
            h1lambda_qinv.SetBinError(ikt+1, f1cf_qinv.GetParError(0));
            h1lambda_qabs.SetBinContent(ikt+1, f1cf_qabs.GetParameter(0));
            h1lambda_qabs.SetBinError(ikt+1, f1cf_qabs.GetParError(0));

            h1R_qinv.SetBinContent(ikt+1, f1cf_qinv.GetParameter(1));
            h1R_qinv.SetBinError(ikt+1, f1cf_qinv.GetParError(1));
            h1R_qabs.SetBinContent(ikt+1, f1cf_qabs.GetParameter(1));
            h1R_qabs.SetBinError(ikt+1, f1cf_qabs.GetParError(1));

        outdir.WriteTObject(h1lambda_qinv);
        outdir.WriteTObject(h1R_qinv);
        outdir.WriteTObject(h1lambda_qabs);
        outdir.WriteTObject(h1R_qabs);
    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
def run_hbt_3d(filename, tasks, arr_qout, arr_qside, arr_qlong, arr_kt):
    print(sys._getframe().f_code.co_name);
    print("filename = ", filename);
    print("arr_qout = ", arr_qout);
    print("arr_qside = ", arr_qside);
    print("arr_qlong = ", arr_qlong);
    print("arr_kt = ", arr_kt);
    __delta = 1e-3;

    rootfile = TFile.Open(filename, "READ");
    outname = "";
    if "pp" in system:
        outname = "photon_hbt_3d_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    else:
        outname = "photon_hbt_3d_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)

    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    nkt = len(arr_kt);
    for itask in range(0, len(tasks)):
        taskname = tasks[itask]["name"];
        cent_min = tasks[itask]["cent_min"];
        cent_max = tasks[itask]["cent_min"];
        print("analyzing ... ", taskname, cent_min, cent_max); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        #rootdir_event.ls();
        #rootdir_pair .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();

        hs_3d_same  = rootdir_pair.Get("same/hs_3d").Clone("hs_3d_same")
        hs_3d_mix   = rootdir_pair.Get("mix/hs_3d") .Clone("hs_3d_mix")
        hs_3d_same .Sumw2();
        hs_3d_mix  .Sumw2();
        outdir.WriteTObject(hs_3d_same);
        outdir.WriteTObject(hs_3d_mix);
        #hDeltaEtaDeltaPhi_same = rootdir_pair.Get("same/hDeltaEtaDeltaPhi") .Clone("hDeltaEtaDeltaPhi_same")
        #hDeltaEtaDeltaPhi_mix  = rootdir_pair.Get("mix/hDeltaEtaDeltaPhi") .Clone("hDeltaEtaDeltaPhi_mix")
        #outdir.WriteTObject(hDeltaEtaDeltaPhi_same);
        #outdir.WriteTObject(hDeltaEtaDeltaPhi_mix);

        h1lambda = TH1F("h1lambda", "correlation strength", nkt-1, arr_kt);
        h1lambda.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1lambda.SetYTitle("#it{#lambda}");

        h1R_out = TH1F("h1R_out", "Rout", nkt-1, arr_kt);
        h1R_out.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R_out.SetYTitle("#it{R}_{out} (fm)");
        h1R_side = TH1F("h1R_side", "Rside", nkt-1, arr_kt);
        h1R_side.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R_side.SetYTitle("#it{R}_{side} (fm)");
        h1R_long = TH1F("h1R_long", "Rlong", nkt-1, arr_kt);
        h1R_long.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R_long.SetYTitle("#it{R}_{long} (fm)");

        for ikt in range(0, nkt-1):
            kt1 = arr_kt[ikt];
            kt2 = arr_kt[ikt+1];
            bin1 = hs_3d_same.GetAxis(3).FindBin(kt1 + __delta);
            bin2 = hs_3d_same.GetAxis(3).FindBin(kt2 - __delta);
            hs_3d_same.GetAxis(3).SetRange(bin1, bin2);
            hs_3d_mix .GetAxis(3).SetRange(bin1, bin2);

            h3same = hs_3d_same.Projection(0, 1, 2);
            h3same.SetName("h3same_kt{0:d}".format(ikt));
            h3same.SetTitle("h3same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h3mix  = hs_3d_mix.Projection(0, 1, 2);
            h3mix.SetName("h3mix_kt{0:d}".format(ikt));
            h3mix.SetTitle("h3mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));

            if h3same.GetEntries() < 0.5 or h3mix.GetEntries() < 0.5:
                continue;
            [h3same, h3mix, h3cf, f3cf] = get_cf_3d(h3same, h3mix, -0.06, 0.06);
            h3same.SetName("h3same_kt{0:d}".format(ikt));
            h3same.SetTitle("h3same: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h3mix.SetName("h3mix_kt{0:d}".format(ikt));
            h3mix.SetTitle("h3mix: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h3cf.SetName("h3cf_kt{0:d}".format(ikt));
            h3cf.SetTitle("h3cf: {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            f3cf.SetName("f3cf_kt{0:d}".format(ikt));
            outdir.WriteTObject(h3same);
            outdir.WriteTObject(h3mix);
            outdir.WriteTObject(h3cf);
            outdir.WriteTObject(f3cf);

            h1lambda.SetBinContent(ikt+1, f3cf.GetParameter(0));
            h1lambda.SetBinError(ikt+1, f3cf.GetParError(0));
            h1R_out.SetBinContent(ikt+1, f3cf.GetParameter(1));
            h1R_out.SetBinError(ikt+1, f3cf.GetParError(1));
            h1R_side.SetBinContent(ikt+1, f3cf.GetParameter(2));
            h1R_side.SetBinError(ikt+1, f3cf.GetParError(2));
            h1R_long.SetBinContent(ikt+1, f3cf.GetParameter(3));
            h1R_long.SetBinError(ikt+1, f3cf.GetParError(3));

        outdir.WriteTObject(h1lambda);
        outdir.WriteTObject(h1R_out);
        outdir.WriteTObject(h1R_side);
        outdir.WriteTObject(h1R_long);
    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
if __name__ == "__main__":
    filename = args.input;

    arr_kt = np.array(config["common"]["kt_bin"], dtype=float);
    arr_qinv = np.array(config["common"]["qinv_bin"], dtype=float);
    arr_qabs = np.array(config["common"]["qabs_bin"], dtype=float);
    arr_qout = np.array(config["common"]["qout_bin"], dtype=float);
    arr_qside = np.array(config["common"]["qside_bin"], dtype=float);
    arr_qlong = np.array(config["common"]["qlong_bin"], dtype=float);

    tasks = config["data"]["tasks"];
    #print(tasks);
    if config["common"]["do_1d"]:
        run_hbt_1d(filename, tasks, arr_qinv, arr_qabs, arr_kt);
    if config["common"]["do_3d"]:
        run_hbt_3d(filename, tasks, arr_qout, arr_qside, arr_qlong, arr_kt);

