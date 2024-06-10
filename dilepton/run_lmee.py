import os
import sys
import shutil
import numpy as np
import yaml
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory
from dilepton_analyzer import DileptonAnalyzer

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();


config = "";
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
#print(config);
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
def run_mee_ptee_dcaee(filename, tasknames, arr_mee, arr_ptee, arr_dcaee):
    print(sys._getframe().f_code.co_name);
    print("filename = ", filename);
    print("arr_mee = ", arr_mee);
    print("arr_ptee = ", arr_ptee);
    print("arr_dcaee = ", arr_dcaee);

    rootfile = TFile.Open(filename, "READ");
    outname = "";
    if "pp" in system:
        outname = "mee_ptee_dcaee_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    else:
        outname = "mee_ptee_dcaee_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");
    #for itask in range(0, len(tasknames)):
    for itask in range(0, 1):
        taskname = tasknames[itask];
        print("analyzing ... ", taskname); 

        rootdir_task = rootfile.Get(taskname);
        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        ana = DileptonAnalyzer(rootdir_event, rootdir_pair);
        ana.set_arr_3d(arr_mee, arr_ptee, arr_dcaee);
        outdir = ana.analyze3d(0, 1, 2, "#it{m}_{ee}", "#it{p}_{T,ee}", "DCA_{ee}^{3D}", "(GeV/#it{c}^{2})", "(GeV/#it{c})", "(#it{#sigma})");
        outdir.SetName(taskname);
        outfile.WriteTObject(outdir);
        outdir.Clear();

#    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
if __name__ == "__main__":
    filename = args.input;
 

    arr_mee = np.array(config["common"]["mee_bin"],dtype=float);
    arr_ptee = np.array(config["common"]["ptee_bin"],dtype=float);
    arr_dcaee = np.array(config["common"]["dcaee_bin"],dtype=float);

    if is_mc :
        tasknames = config["mc"]["tasknames"]
        #print(tasknames);
    else:
        tasknames = config["data"]["tasknames"]
        #print(tasknames);
        run_mee_ptee_dcaee(filename, tasknames, arr_mee, arr_ptee, arr_dcaee);


