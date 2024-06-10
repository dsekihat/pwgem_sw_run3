import numpy as np
import pandas as pd

import ROOT
from ROOT import TH1D, TGraph, TGraphErrors, TGraphAsymmErrors, TFile, TDirectory
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from signal_extractor import get_R_factor, get_corrected_bkg, get_significance

class DileptonAnalyzer:
    def __init__(self, rootdir_event, rootdir_pair):
        self.rootdir_event = rootdir_event;
        self.rootdir_pair  = rootdir_pair;
        #self.hs_uls_same = hs_uls_same.Clone("{0}_clone".format(hs_uls_same.GetName()));
        #self.hs_lspp_same = hs_lspp_same.Clone("{0}_clone".format(hs_lspp_same.GetName()));
        #self.hs_lsmm_same = hs_lsmm_same.Clone("{0}_clone".format(hs_lsmm_same.GetName()));
        #self.hs_uls_mix = hs_uls_mix.Clone("{0}_clone".format(hs_uls_mix.GetName()));
        #self.hs_lspp_mix = hs_lspp_mix.Clone("{0}_clone".format(hs_lspp_mix.GetName()));
        #self.hs_lsmm_mix = hs_lsmm_mix.Clone("{0}_clone".format(hs_lsmm_mix.GetName()));
        self.arr0 = np.array([0, 1], dtype=float);
        self.arr1 = np.array([0, 1], dtype=float);
        self.arr2 = np.array([0, 1], dtype=float);
        self.arr3 = np.array([0, 1], dtype=float);

    def set_arr_3d(self, arr0, arr1, arr2):
        self.arr0 = arr0;
        self.arr1 = arr1;
        self.arr2 = arr2;

    def set_arr_2d(self, arr0, arr1):
        self.arr0 = arr0;
        self.arr1 = arr1;

    def analyze3d(self, axis_id0=0, axis_id1=1, axis_id2=2, var0="#it{m}_{ee}", var1="#it{p}_{T,ee}", var2="DCA_{ee}^{3D}", unit0="(GeV/#it{c}^{2})", unit1="(GeV/#it{c})", unit2="(#it{#sigma})"):
        print("test 3d");
        outdir = TDirectory("tmp", "tmp");
        h1 = TH1D("h1", "h1", 10,0,10);
        outdir.Add(h1);
        return outdir;

    def analyze2d(self, axis_id0=0, axis_id1=1, title0="", title1=""):
        print("test 2d");
