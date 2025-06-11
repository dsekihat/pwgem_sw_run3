import os
import sys
import shutil
import math
import numpy as np
import yaml
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory, TH1F, TH2F, THnSparseD
from dilepton_utils import get_yield_1d, get_flow
import sys
sys.path.append("../common/");

#__________________________________________________________
def extract_dca_template(filename, taskname, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();

    rootdire_pair_rec = rootdire.Get("Pair");
    rootdire_pair_rec.ls();

    arr_idx = np.array([0, 1, 10], dtype=np.int32); #mee, pTee, DCAee

    parnames_org = [
        "sm/Photon/hs",
        "sm/PromptPi0/hs",
        "sm/NonPromptPi0/hs",
        "sm/Eta/hs",
        "sm/EtaPrime/hs",
        "sm/Rho/hs",
        "sm/Omega/hs",
        "sm/Phi/hs",
        "sm/PromptJPsi/hs",
        "sm/NonPromptJPsi/hs",
        "sm/PromptPsi2S/hs",
        "sm/NonPromptPsi2S/hs",
        "ccbar/c2l_c2l/hadron_hadron/hs",
        "bbbar/b2l_b2l/hadron_hadron/hs",
        "bbbar/b2c2l_b2c2l/hadron_hadron/hs",
        "bbbar/b2c2l_b2l_sameb/hadron_hadron/hs",
        "bbbar/b2c2l_b2l_diffb/hadron_hadron/hs",
    ];
    parnames = [
        "photon",
        "promptpi0",
        "nonpromptpi0",
        "eta",
        "etaprime",
        "rho",
        "omega",
        "phi",
        "promptjpsi",
        "nonpromptjpsi",
        "promptpsi2s",
        "nonpromptpsi2s",
        "c2l_c2l",
        "b2l_b2l",
        "b2c2l_b2c2l",
        "b2c2l_b2l_sameb",
        "b2c2l_b2l_diffb",
    ];

    h1z = rootdire.Get("Event/after/hZvtx");
    nev = h1z.GetEntries();
    print("nev = {0:e}".format(nev));
    outname = "output_template{0}.root".format(suffix)
    outfile = TFile(outname, "RECREATE");
    for ip in range(0, len(parnames_org)):
        hs_org = rootdire_pair_rec.Get(parnames_org[ip]);
        hs_3d = hs_org.ProjectionND(len(arr_idx), arr_idx, "")
        hs_3d.SetName("hs_template_{0}".format(parnames[ip]));
        hs_3d.Scale(1/nev);
        outfile.WriteTObject(hs_3d);

    rootfile.Close();
    outfile.Close();

#__________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_286035.root";
    #taskname = "dielectron-mc_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_3050_LHC24g8";
    #extract_dca_template(filename, taskname, 30, 50, suffix);


    #filename = "AnalysisResults_HL_289057.root";
    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dca3d";
    #suffix = "_pp_13.6TeV_dca3d_LHC23k4g_HL_289057";
    #extract_dca_template(filename, taskname, suffix);

    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dcaxy";
    #suffix = "_pp_13.6TeV_dcaxy_LHC23k4g_HL_289057";
    #extract_dca_template(filename, taskname, suffix);

    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dcaz";
    #suffix = "_pp_13.6TeV_dcaz_LHC23k4g_HL_289057";
    #extract_dca_template(filename, taskname, suffix);

    #filename = "AnalysisResults_HL_290262.root";
    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dca3d";
    #suffix = "_pp_13.6TeV_dca3d_LHC23k4g_HL_290262";
    #extract_dca_template(filename, taskname, suffix);

    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dcaxy";
    #suffix = "_pp_13.6TeV_dcaxy_LHC23k4g_HL_290262";
    #extract_dca_template(filename, taskname, suffix);

    #taskname = "dielectron-mc_TPChadrejorTOFreq_minpt200_dcaz";
    #suffix = "_pp_13.6TeV_dcaz_LHC23k4g_HL_290262";
    #extract_dca_template(filename, taskname, suffix);


    #filename = "AnalysisResults_HL_389477.root";
    #taskname = "dimuon-mc_global_tight_minpt300";
    #suffix = "_pp_13.6TeV_dcaxy_LHC24f3c_fix_HL_389477";
    #extract_dca_template(filename, taskname, suffix);

    filename = "AnalysisResults_HL_394914.root"; # LHC23k4g
    taskname = "dimuon-mc";
    suffix = "_pp_13.6TeV_dcaxy_LHC23k4g_HL_394914";
    extract_dca_template(filename, taskname, suffix);


#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
