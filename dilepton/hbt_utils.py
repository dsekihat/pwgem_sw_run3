import math
import ctypes
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
from ROOT import TF1, TMath, TF3

#__________________________________________________________
def get_cf_1d(h1same_org, h1mix_org, fit_min, fit_max):
    h1same = h1same_org.Clone("h1same_org");
    h1mix = h1mix_org.Clone("h1mix_org");
    #npair_same = h1same.Integral(1, h1same.GetNbinsX(), "");
    #npair_mix  = h1mix .Integral(1, h1mix.GetNbinsX(), "");
    npair_same = h1same.GetEntries();
    npair_mix  = h1mix .GetEntries();
    print(npair_same, npair_mix);
    h1same.Scale(1/npair_same);
    h1mix .Scale(1/npair_mix);

    h1cf = h1same.Clone("h1cf");
    h1cf.Reset();
    h1cf.Divide(h1same, h1mix, 1., 1., "G");
    h1cf.SetYTitle("C_{2}");

    __sf = TMath.Hbar() * TMath.C() / TMath.Qe() /1e+9 / 1e-15; #GeV x fm
    f1cf = TF1("f1cf", "1 + [0] * exp(-[1]*[1]*x*x/[2]/[2])", 0, 0.1);
    f1cf.SetNpx(1000);
    f1cf.SetParameter(0, 0.1);
    f1cf.SetParameter(1, 10.0);
    f1cf.FixParameter(2, __sf);
    f1cf.SetParLimits(0, 0, 10);
    f1cf.SetParLimits(1, 1e-2, 1e+2);
    #f1cf.SetParNames("#lambda", "R (fm)", "#hbar");
    f1cf.SetParNames("#lambda", "R", "#hbar");
    h1cf.Fit(f1cf, "SME" ,"" ,fit_min, fit_max);
    return [h1same, h1mix, h1cf, f1cf];
#__________________________________________________________
#__________________________________________________________
def get_cf_3d(h3same_org, h3mix_org, fit_min, fit_max):
    h3same = h3same_org.Clone("h3same_org");
    h3mix = h3mix_org.Clone("h3mix_org");
    #npair_same = h3same.Integral(1, h3same.GetNbinsX(), 1, h3same.GetNbinsY(), 1, h3same.GetNbinsZ(), "");
    #npair_mix  = h3mix .Integral(1, h3mix.GetNbinsX() , 1, h3mix.GetNbinsY() , 1, h3mix.GetNbinsZ() , "");
    npair_same = h3same.GetEntries();
    npair_mix  = h3mix .GetEntries();
    print(npair_same, npair_mix);
    h3same.Scale(1/npair_same);
    h3mix .Scale(1/npair_mix);

    h3cf = h3same.Clone("h3cf");
    h3cf.Reset();
    h3cf.Divide(h3same, h3mix, 1., 1., "G");

    __sf = TMath.Hbar() * TMath.C() / TMath.Qe() /1e+9 / 1e-15; #GeV x fm
    f3cf = TF3("f3cf", "1 + [0] * exp(-[1]*[1]*x*x/[4]/[4] -[2]*[2]*y*y/[4]/[4] -[3]*[3]*z*z/[4]/[4])", fit_min, fit_max, fit_min, fit_max, fit_min, fit_max);
    f3cf.SetNpx(int(1000 * (fit_max - fit_min)));
    f3cf.SetNpy(int(1000 * (fit_max - fit_min)));
    f3cf.SetNpz(int(1000 * (fit_max - fit_min)));
    f3cf.SetParameter(0, 0.1);
    f3cf.SetParameter(1, 10.0);
    f3cf.SetParameter(2, 10.0);
    f3cf.SetParameter(3, 10.0);
    f3cf.FixParameter(4, __sf);
    f3cf.SetParLimits(0, 0, 10);
    f3cf.SetParLimits(1, 1e-2, 1e+2);
    f3cf.SetParLimits(2, 1e-2, 1e+2);
    f3cf.SetParLimits(3, 1e-2, 1e+2);
    #f3cf.SetParNames("#lambda", "R_{out} (fm)", "R_{side} (fm)", "R_{long} (fm)", "#hbar"); #space, underscore crash this code.
    f3cf.SetParNames("#lambda", "Rout", "Rside", "Rlong", "#hbar");
    h3cf.Fit(f3cf, "SMER" ,"");
    return [h3same, h3mix, h3cf, f3cf];
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
