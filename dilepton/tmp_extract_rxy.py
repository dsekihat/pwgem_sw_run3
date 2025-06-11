import datetime
import sys
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D, gStyle
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);


filename = "AnalysisResults_HL_251833.root";
taskname = "pcm-qc";
rootfile = TFile.Open(filename, "READ");
roottask = rootfile.Get(taskname);
roottask.ls();

roottask_ev = roottask.Get("Event/after");
roottask_ev.ls();
h1z = roottask_ev.Get("hZvtx");
nev = h1z.GetEntries();
print(nev/1e+9);
roottask_v0 = roottask.Get("V0");
h2 = roottask_v0.Get("hGammaRxy");
h2.SetContour(1000);
h2.Scale(1/nev);

outfile = TFile("20241112_rxy_pp_13.6TeV.root", "RECREATE");
outfile.WriteTObject(h2);
outfile.Close();

rootfile.Close();
