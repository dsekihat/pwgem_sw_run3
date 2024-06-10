import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd

import ROOT
from ROOT import TH1D, TGraph, TGraphErrors, TGraphAsymmErrors, TFile, TDirectory
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from signal_extractor import get_R_factor, get_corrected_bkg, get_corrected_bkg_simple, get_significance
from histo_manager import rebin_histogram, get_bkg_subtracted

class DileptonAnalyzer:
    def __init__(self, h1_uls_same, h1_lspp_same, h1_lsmm_same, h1_uls_mix, h1_lspp_mix, h1_lsmm_mix, arr):
        self.h1_uls_same  = h1_uls_same.Clone("{0}_clone".format(h1_uls_same.GetName()));
        self.h1_lspp_same = h1_lspp_same.Clone("{0}_clone".format(h1_lspp_same.GetName()));
        self.h1_lsmm_same = h1_lsmm_same.Clone("{0}_clone".format(h1_lsmm_same.GetName()));

        if h1_uls_mix is not None:
            self.h1_uls_mix = h1_uls_mix.Clone("{0}_clone".format(h1_uls_mix.GetName()));
            self.h1_lspp_mix = h1_lspp_mix.Clone("{0}_clone".format(h1_lspp_mix.GetName()));
            self.h1_lsmm_mix = h1_lsmm_mix.Clone("{0}_clone".format(h1_lsmm_mix.GetName()));
        else:
            self.h1_uls_mix  = None;
            self.h1_lspp_mix = None;
            self.h1_lsmm_mix = None;
        self.arr = arr;

    def run(self):
        h1_uls_same_rebin  = rebin_histogram(self.h1_uls_same , self.arr, False, False);
        h1_lspp_same_rebin = rebin_histogram(self.h1_lspp_same, self.arr, False, False);
        h1_lsmm_same_rebin = rebin_histogram(self.h1_lsmm_same, self.arr, False, False);
        if self.h1_uls_mix is not None:
            h1_uls_mix_rebin  = rebin_histogram(self.h1_uls_mix , self.arr, False, False);
            h1_lspp_mix_rebin = rebin_histogram(self.h1_lspp_mix, self.arr, False, False);
            h1_lsmm_mix_rebin = rebin_histogram(self.h1_lsmm_mix, self.arr, False, False);
            h1r = get_R_factor(h1_uls_mix_rebin, None, h1_lspp_mix_rebin, h1_lsmm_mix_rebin);
            h1bkg = get_corrected_bkg(h1r, h1_lspp_same_rebin, h1_lsmm_same_rebin);
            h1sig = get_bkg_subtracted(h1_uls_same_rebin, h1bkg);
            return [h1_uls_same_rebin, h1_lspp_same_rebin, h1_lsmm_same_rebin, h1_uls_mix_rebin, h1_lspp_mix_rebin, h1_lsmm_mix_rebin, h1r, h1bkg, h1sig];
        else:
            h1r = None;
            h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1_lspp_same_rebin, h1_lsmm_same_rebin);
            h1sig = get_bkg_subtracted(h1_uls_same_rebin, h1bkg);
            return [h1_uls_same_rebin, h1_lspp_same_rebin, h1_lsmm_same_rebin, None, None, None, h1r, h1bkg, h1sig];
