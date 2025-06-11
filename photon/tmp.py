import numpy as np
import ctypes
import ROOT
from ROOT import TFile
from nm_fitter import NMFitter

filename = "AnalysisResults_HL_233027.root";
rootfile = TFile.Open(filename, "READ");
rootfile.ls();
roottask = rootfile.Get("pi0eta-to-gammagamma-pcmpcm");
roottask.ls();
dir_pair = roottask.Get("Pair");
hs_same = dir_pair.Get("same").Get("hs");
hs_mix  = dir_pair.Get("mix").Get("hs");

pt1 = 0.4;
pt2 = 0.5;

bin1 = hs_same.GetAxis(1).FindBin(pt1 + 1e-3);
bin2 = hs_same.GetAxis(1).FindBin(pt2 - 1e-3);

hs_same.GetAxis(1).SetRange(bin1, bin2);
hs_mix .GetAxis(1).SetRange(bin1, bin2);

h1same = hs_same.Projection(0);
h1mix = hs_mix.Projection(0);
h1same.RebinX(2);
h1mix .RebinX(2);

npair_same = h1same.GetEntries();
npair_mix  = h1mix.GetEntries();
h1mix_scaled = h1mix.Clone("h1mix_scaled_pt{0}".format(0));
h1mix_scaled.Scale(npair_same/npair_mix);


nmf = NMFitter(h1same, h1mix_scaled, "cb", "pol1"); 
nmf.set_parameters(0.13, 0.005, 12, 0.6, True, True);
fit_result = nmf.fit("SME", "", 0.04, 0.24);
h1sig = fit_result[1];
h1bkg = fit_result[2];
h1ratio = fit_result[3];
f1sig = fit_result[4];
f1bkg = fit_result[5];
f1total = fit_result[6];


mean = f1sig.GetParameter(1);
sigma = f1sig.GetParameter(2);

bin1 = h1sig.FindBin(mean - 3 * sigma);
bin2 = h1sig.FindBin(mean + 3 * sigma);
print(bin1, bin2);
ry_err = ctypes.c_double(0);
ry = h1sig.IntegralAndError(bin1, bin2, ry_err, "");
print(ry, ry_err.value);


h1bkg.Draw("");
h1same.Draw("same");
h1sig.Draw("same");
h1same.SetDirectory(0);
ROOT.SetOwnership(h1same, False);



