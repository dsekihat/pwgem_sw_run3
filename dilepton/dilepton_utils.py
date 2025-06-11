import math
import ctypes
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd

from signal_extractor import get_R_factor, get_corrected_bkg, get_corrected_bkg_simple, get_significance, get_Rn, get_R_factor_1number
from histo_manager import rebin_histogram, get_bkg_subtracted, rebin_profile

#__________________________________________________________
def get_yield_1d(h1_uls_same, h1_lspp_same, h1_lsmm_same, h1_uls_mix, h1_lspp_mix, h1_lsmm_mix, arr):
    isR1 = False;
    if h1_uls_mix is None or h1_lspp_mix is None or h1_lsmm_mix is None:
        isR1 = True; 

    h1_uls_same_rebin  = rebin_histogram(h1_uls_same , arr, False, False);
    h1_lspp_same_rebin = rebin_histogram(h1_lspp_same, arr, False, False);
    h1_lsmm_same_rebin = rebin_histogram(h1_lsmm_same, arr, False, False);

    if isR1:
        h1r = None;
        h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1_lspp_same_rebin, h1_lsmm_same_rebin);
        h1sig = get_bkg_subtracted(h1_uls_same_rebin, h1bkg);
        return [h1_uls_same_rebin, h1_lspp_same_rebin, h1_lsmm_same_rebin, None, None, None, h1r, h1bkg, h1sig];
    else:
        h1_uls_mix_rebin  = rebin_histogram(h1_uls_mix , arr, False, False);
        h1_lspp_mix_rebin = rebin_histogram(h1_lspp_mix, arr, False, False);
        h1_lsmm_mix_rebin = rebin_histogram(h1_lsmm_mix, arr, False, False);
        h1r = get_R_factor(h1_uls_mix_rebin, None, h1_lspp_mix_rebin, h1_lsmm_mix_rebin);
        h1bkg = get_corrected_bkg(h1r, h1_lspp_same_rebin, h1_lsmm_same_rebin);
        h1sig = get_bkg_subtracted(h1_uls_same_rebin, h1bkg);
        return [h1_uls_same_rebin, h1_lspp_same_rebin, h1_lsmm_same_rebin, h1_uls_mix_rebin, h1_lspp_mix_rebin, h1_lsmm_mix_rebin, h1r, h1bkg, h1sig];
#__________________________________________________________
def get_flow_old(h1_uls_same_org, h1_lspp_same_org, h1_lsmm_same_org, h1_uls_mix_org, h1_lspp_mix_org, h1_lsmm_mix_org,
        h1prf_uls_uq_same_org,
        h1prf_uls_uq_cosdphi_leg1_mix_org , 
        h1prf_uls_2uq1uq2cosdphi12_mix_org,
        arr, Rn):

    # reference : https://arxiv.org/pdf/2005.14518
    [h1_uls_same, h1_lspp_same, h1_lsmm_same, h1_uls_mix, h1_lspp_mix, h1_lsmm_mix, h1R, h1bkg, h1sig] = get_yield_1d(h1_uls_same_org, h1_lspp_same_org, h1_lsmm_same_org, h1_uls_mix_org, h1_lspp_mix_org, h1_lsmm_mix_org, arr); # for yield

    h1prf_uls_uq_same = rebin_profile(h1prf_uls_uq_same_org, arr);
    h1vn_uls_same = h1prf_uls_uq_same.ProjectionX("h1vn_uls_same");
    #h1vn_uls_same.Scale(1./Rn);
    h1vn_uls_same.SetYTitle("v_{n}^{S+B}");

    h1prf_uls_uq_cosdphi_leg1_mix  = rebin_profile(h1prf_uls_uq_cosdphi_leg1_mix_org , arr);
    h1prf_uls_2uq1uq2cosdphi12_mix = rebin_profile(h1prf_uls_2uq1uq2cosdphi12_mix_org, arr);

    h1uq_cosdphi_uls_leg1_mix = h1prf_uls_uq_cosdphi_leg1_mix .ProjectionX("h1uq_cosdphi_uls_leg1_mix"     );
    h1_2uq1uq2cosdphi12_mix   = h1prf_uls_2uq1uq2cosdphi12_mix.ProjectionX("h1_2uq1uq2cosdphi12_mix"     );

    h1vn_uls_leg1_mix = h1uq_cosdphi_uls_leg1_mix.Clone("h1vn_uls_leg1_mix");
    #h1vn_uls_leg1_mix.Scale(1./Rn);
    h1vn_bkg = h1vn_uls_leg1_mix.Clone("h1vn_bkg");

    h1_2v1v2cosdphi12_mix = h1_2uq1uq2cosdphi12_mix.Clone("h1_2v1v2cosdphi12_mix"); # modification of bkg shape due to single lepton flow
    #h1_2v1v2cosdphi12_mix.Scale(1./Rn/Rn);

    h1_1_2v1v2cosdphi12_mix = h1_2v1v2cosdphi12_mix.Clone("h1_1_2v1v2cosdphi12_mix");
    for i in range(0, h1_1_2v1v2cosdphi12_mix.GetNbinsX()):
        y_mod = h1_2v1v2cosdphi12_mix.GetBinContent(i+1);
        y_mod_err = h1_2v1v2cosdphi12_mix.GetBinError(i+1);
        h1_1_2v1v2cosdphi12_mix.SetBinContent(i+1, y_mod + 1);
        h1_1_2v1v2cosdphi12_mix.SetBinError(i+1, y_mod_err);

    #x_min = arr[0];
    #x_max = arr[-1];
    #bin1 = h1bkg.FindBin(arr[0] + 1e-3);
    #bin2 = h1bkg.FindBin(arr[-1] - 1e-3);
    #nbkg_uls_err = ctypes.c_double(0);
    #nbkg_uls = h1bkg.IntegralAndError(bin1, bin2, nbkg_uls_err, ""); # h1bkg is NOT differential yield

    #nbkg_uls_mix_err = ctypes.c_double(0);
    #nbkg_uls_mix = h1_uls_mix.IntegralAndError(bin1, bin2, nbkg_uls_mix_err, ""); # h1_uls_mix is NOT differential yield
    #h1_uls_mix.Scale(nbkg_uls/nbkg_uls_mix); # normalized to nbkg_uls in the same event
    #h1bkg_flow = h1bkg.Clone("h1bkg_flow");
    #h1bkg_flow.Divide(h1_uls_mix); #h1bkg_flow represents modification of bkg shape because of single lepton flow.

    h1vn_bkg.Divide(h1_1_2v1v2cosdphi12_mix);
    #h1vn_bkg.Divide(h1bkg_flow);
    h1vn_bkg.SetYTitle("v_{n}^{B}");
    ##h1vn_bkg.Scale(Rn);

    h1vn_bkg_corr = h1vn_bkg.Clone("h1vn_bkg_corr");
    h1vn_bkg_corr.SetYTitle("#frac{B}{S+B} v_{n}^{B}");

    h1vn_sig = h1vn_uls_same.Clone("h1vn_sig"); # only for binning.
    h1vn_sig.Reset();
    h1vn_sig.SetYTitle("v_{n}^{S}");
    
    h1frac_sig = h1sig.Clone("h1frac_sig");
    h1frac_sig.Reset();
    h1frac_sig.Divide(h1sig, h1_uls_same, 1., 1., "B");
    for i in range(0, h1vn_sig.GetNbinsX()):
        frac_sig = h1frac_sig.GetBinContent(i+1);
        if frac_sig < 1e-9:
            continue;

        #frac_sig_err = h1frac_sig.GetBinError(i+1);
        frac_sig_err = math.sqrt(math.pow(h1bkg.GetBinContent(i+1) * h1_uls_same.GetBinError(i+1),2)/math.pow(h1_uls_same.GetBinContent(i+1),4) + math.pow(h1bkg.GetBinError(i+1)/h1_uls_same.GetBinContent(i+1),2) );
        frac_bkg = 1 - frac_sig;
        frac_bkg_err = frac_sig_err;
        frac_sig_err = 0;
        frac_bkg_err = 0;

        vn_sb = h1vn_uls_same.GetBinContent(i+1);
        vn_sb_err = h1vn_uls_same.GetBinError(i+1);
        vn_b = h1vn_bkg.GetBinContent(i+1);
        vn_b_err = h1vn_bkg.GetBinError(i+1);
        h1vn_bkg_corr.SetBinContent(i+1, vn_b * frac_bkg);
        h1vn_bkg_corr.SetBinError(i+1, vn_b_err);

        vn_s = (vn_sb - frac_bkg * vn_b) / frac_sig;
        vn_s_err = math.sqrt(math.pow(vn_sb_err/frac_sig, 2) + math.pow(frac_bkg*vn_b_err/frac_sig, 2));
        #vn_s_err = math.sqrt(math.pow(vn_sb_err/frac_sig, 2) + math.pow(frac_bkg*vn_b_err/frac_sig, 2) + math.pow(vn_sb * frac_sig_err,2)/math.pow(frac_sig,4) + math.pow(vn_b * frac_bkg_err/frac_sig,2));

        h1vn_sig.SetBinContent(i+1, vn_s);
        h1vn_sig.SetBinError(i+1, vn_s_err);
    return [h1vn_uls_same, h1uq_cosdphi_uls_leg1_mix, h1_1_2v1v2cosdphi12_mix, h1frac_sig, h1vn_bkg, h1vn_bkg_corr, h1vn_sig];
#__________________________________________________________
def get_flow(h1uq_uls_same, h1uq_lspp_same, h1uq_lsmm_same, frac_sig, frac_sig_err, R, R_err, Rn):
    h1uq_uls_same.ResetStats();
    h1uq_lspp_same.ResetStats();
    h1uq_lsmm_same.ResetStats();

    vn_tot = h1uq_uls_same.GetMean();
    vn_tot_err = h1uq_uls_same.GetMeanError();

    h1uq_bkg = get_corrected_bkg_simple(1, 0, h1uq_lspp_same, h1uq_lsmm_same);
    h1uq_bkg.ResetStats();
    vn_bkg = h1uq_bkg.GetMean();
    vn_bkg_err = h1uq_bkg.GetMeanError();

    #print("vn_tot: ", h1uq_uls_same.Integral(), vn_tot, vn_tot_err);
    #print("vn_bkg: ", h1uq_bkg.Integral(), vn_bkg, vn_bkg_err);

    frac_bkg = 1 - frac_sig;
    frac_bkg_err = frac_sig_err;
    if frac_sig > 1e-6:
        vn_sig = (vn_tot - frac_bkg * vn_bkg) / frac_sig;
        vn_sig_err = math.sqrt(math.pow(vn_tot_err/frac_sig, 2) + math.pow(frac_bkg*vn_bkg_err/frac_sig, 2));
        #vn_sig_err = math.sqrt(math.pow(vn_tot_err/frac_sig, 2) + math.pow(frac_bkg*vn_bkg_err/frac_sig, 2) + math.pow(vn_tot * frac_sig_err,2)/math.pow(frac_sig,4) + math.pow(vn_bkg * frac_bkg_err/frac_sig,2));
    else:
        vn_sig     = 0;
        vn_sig_err = 0;

    return [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err, h1uq_bkg];
    #return [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err, h1uq_bkg, h1uq_sig];
#__________________________________________________________
def extract_flow(hs_uls_same, hs_lspp_same, hs_lsmm_same, hs_uls_mix, hs_lspp_mix, hs_lsmm_mix, Rn, m_min, m_max, pt_min, pt_max, dca_min, dca_max):
    __delta = 1e-3;
    #set range for mass axis
    bin_m0 = hs_uls_same.GetAxis(0).FindBin(m_min + __delta);    
    bin_m1 = hs_uls_same.GetAxis(0).FindBin(m_max - __delta);    
    hs_uls_same.GetAxis(0).SetRange(bin_m0, bin_m1);
    hs_lspp_same.GetAxis(0).SetRange(bin_m0, bin_m1);
    hs_lsmm_same.GetAxis(0).SetRange(bin_m0, bin_m1);
    hs_uls_mix.GetAxis(0).SetRange(bin_m0, bin_m1);
    hs_lspp_mix.GetAxis(0).SetRange(bin_m0, bin_m1);
    hs_lsmm_mix.GetAxis(0).SetRange(bin_m0, bin_m1);

    #set range for pair pt axis
    bin_pt0 = hs_uls_same.GetAxis(1).FindBin(pt_min + __delta);    
    bin_pt1 = hs_uls_same.GetAxis(1).FindBin(pt_max - __delta);    
    hs_uls_same.GetAxis(1).SetRange(bin_pt0, bin_pt1);
    hs_lspp_same.GetAxis(1).SetRange(bin_pt0, bin_pt1);
    hs_lsmm_same.GetAxis(1).SetRange(bin_pt0, bin_pt1);
    hs_uls_mix.GetAxis(1).SetRange(bin_pt0, bin_pt1);
    hs_lspp_mix.GetAxis(1).SetRange(bin_pt0, bin_pt1);
    hs_lsmm_mix.GetAxis(1).SetRange(bin_pt0, bin_pt1);

    #set range for pair dca axis
    bin_dca0 = hs_uls_same.GetAxis(2).FindBin(dca_min + __delta);    
    bin_dca1 = hs_uls_same.GetAxis(2).FindBin(dca_max - __delta);    
    hs_uls_same.GetAxis(2).SetRange(bin_dca0, bin_dca1);
    hs_lspp_same.GetAxis(2).SetRange(bin_dca0, bin_dca1);
    hs_lsmm_same.GetAxis(2).SetRange(bin_dca0, bin_dca1);
    hs_uls_mix.GetAxis(2).SetRange(bin_dca0, bin_dca1);
    hs_lspp_mix.GetAxis(2).SetRange(bin_dca0, bin_dca1);
    hs_lsmm_mix.GetAxis(2).SetRange(bin_dca0, bin_dca1);

    #hs_same has 4D, hs_mix has 3D

    h1m_uls_mix_tmp = hs_uls_mix.Projection(0);
    h1m_uls_mix_tmp .SetName("h1m_uls_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}_tmp".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    h1m_lspp_mix_tmp = hs_lspp_mix.Projection(0);
    h1m_lspp_mix_tmp.SetName("h1m_lspp_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}_tmp".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    h1m_lsmm_mix_tmp = hs_lsmm_mix.Projection(0);
    h1m_lsmm_mix_tmp.SetName("h1m_lsmm_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}_tmp".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    #print("h1m_uls_mix_tmp.GetNbinsX() = ", h1m_uls_mix_tmp.GetNbinsX());
    n_uls_mix_err  = ctypes.c_double(0);
    n_lspp_mix_err = ctypes.c_double(0);
    n_lsmm_mix_err = ctypes.c_double(0);
    n_uls_mix  = h1m_uls_mix_tmp .IntegralAndError(1, h1m_uls_mix_tmp.GetNbinsX() , n_uls_mix_err, "");
    n_lspp_mix = h1m_lspp_mix_tmp.IntegralAndError(1, h1m_lspp_mix_tmp.GetNbinsX(), n_uls_mix_err, "");
    n_lsmm_mix = h1m_lsmm_mix_tmp.IntegralAndError(1, h1m_lsmm_mix_tmp.GetNbinsX(), n_uls_mix_err, "");

    [R, R_err] = get_R_factor_1number(n_uls_mix, n_lspp_mix, n_lsmm_mix, n_uls_mix_err.value, n_lspp_mix_err.value, n_lsmm_mix_err.value);
    #print("R factor = ", R, R_err);

    #n_uls_same  = 0;
    #n_lspp_same = 0;
    #n_lsmm_same = 0;
    #for im in range(bin_m0, bin_m1 + 1):
    #    for ipt in range(bin_pt0, bin_pt1 + 1):
    #        for idca in range(bin_dca0, bin_dca1 + 1):
    #            for iuq in range(1, hs_uls_same.GetAxis(3).GetNbins() + 1):
    #                n_uls_same      += hs_uls_same .GetBinContent(np.array([im, ipt, idca, iuq], dtype=np.int32));
    #                n_lspp_same     += hs_lspp_same.GetBinContent(np.array([im, ipt, idca, iuq], dtype=np.int32));
    #                n_lsmm_same     += hs_lspp_same.GetBinContent(np.array([im, ipt, idca, iuq], dtype=np.int32));
    #n_uls_same_err  = math.sqrt(n_uls_same );
    #n_lspp_same_err = math.sqrt(n_lspp_same);
    #n_lsmm_same_err = math.sqrt(n_lsmm_same);
    #n_bkg = 2.0 * R * math.sqrt(n_lspp_same * n_lsmm_same);
    #n_sig = n_uls_same - n_bkg;
    #print("yield = ", n_uls_same, n_bkg, n_sig);

    h1uq_uls_same = hs_uls_same.Projection(3);
    h1uq_uls_same .SetName("h1us_uls_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    h1uq_lspp_same = hs_lspp_same.Projection(3);
    h1uq_lspp_same.SetName("h1us_lspp_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    h1uq_lsmm_same = hs_lsmm_same.Projection(3);
    h1uq_lsmm_same.SetName("h1us_lsmm_same_m{0:d}_{1:d}_pt{2:d}_{3:d}_dca{4:d}_{5:d}".format(bin_m0, bin_m1, bin_pt0, bin_pt1, bin_dca0, bin_dca1));
    h1uq_uls_same .RebinX(2);
    h1uq_lspp_same.RebinX(2);
    h1uq_lsmm_same.RebinX(2);
    h1uq_uls_same.ResetStats();
    h1uq_lspp_same.ResetStats();
    h1uq_lsmm_same.ResetStats();

    n_uls_same_err  = ctypes.c_double(0);
    n_lspp_same_err = ctypes.c_double(0);
    n_lsmm_same_err = ctypes.c_double(0);
    n_uls_same  = h1uq_uls_same .IntegralAndError(1, h1uq_uls_same .GetNbinsX(), n_uls_same_err , "");
    n_lspp_same = h1uq_lspp_same.IntegralAndError(1, h1uq_lspp_same.GetNbinsX(), n_lspp_same_err, "");
    n_lsmm_same = h1uq_lsmm_same.IntegralAndError(1, h1uq_lsmm_same.GetNbinsX(), n_lsmm_same_err, "");

    n_bkg = 2.0 * R * math.sqrt(n_lspp_same * n_lsmm_same);
    n_sig = n_uls_same - n_bkg;

    frac_sig = n_sig/n_uls_same;
    frac_bkg = 1.0 - frac_sig;
    print("yield = ", n_uls_same, n_bkg, n_sig, frac_sig, frac_bkg);

    vn_tot = h1uq_uls_same.GetMean() / Rn;
    vn_tot_err = h1uq_uls_same.GetMeanError() / Rn;

    h1uq_bkg = get_corrected_bkg_simple(R, R_err, h1uq_lspp_same, h1uq_lsmm_same);
    h1uq_bkg.ResetStats();
    vn_bkg = h1uq_bkg.GetMean() / Rn;
    vn_bkg_err = h1uq_bkg.GetMeanError() / Rn;

    if frac_sig > 1e-6:
        vn_sig = (vn_tot - frac_bkg * vn_bkg)/frac_sig;
        vn_sig_err = math.sqrt(math.pow(vn_tot_err/frac_sig, 2) + math.pow(frac_bkg*vn_bkg_err/frac_sig, 2));
    else:
        vn_sig = 0;
        vn_sig_err = 0;
    return [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err]
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________

