import numpy as np
import pandas as pd
import ROOT
from ROOT import TProfile, TH2D
from ROOT import gPad

#h1 = TProfile("h1", "h1prf", 10, 0, 10);
#h1.Sumw2();
#ROOT.SetOwnership(h1, False);
#h1.Draw();
##h1.Fill(1, 1);
##h1.Fill(1, 2);
#h1.SetBinContent(1, 1.5);
#h1.SetBinError(1, 0.2);
#h1.SetBinEntries(1, 1);
#
##h1.Fill(2, 2);
##h1.Fill(2, 3);
#h1.SetBinContent(2, 2.5);
#h1.SetBinError(2, 0.3);
#h1.SetBinEntries(2, 1);
#
#print(h1.GetEntries(), h1.GetEffectiveEntries(), h1.GetBinEffectiveEntries(1));
#print(h1.GetBinContent(1), h1.GetBinError(1));

h2 = TH2D("h2", "h2;x;y", 10, 0, 10, 10, 0, 10);
h2.Sumw2();
ROOT.SetOwnership(h2, False);

h2.SetBinContent(1, 1, 1.5);
h2.SetBinError  (1, 1, 0.5);
h2.SetBinContent(2, 1, 0.5);
h2.SetBinError  (2, 1, 0.1);

h2.SetBinContent(1, 2, 2.0);
h2.SetBinError  (1, 2, 1.0);
h2.SetBinContent(2, 2, 2.5);
h2.SetBinError  (2, 2, 1.0);
#h2.Draw("colz");

h1 = h2.ProfileX("h1", 1, 2, "");
h1.Sumw2();
ROOT.SetOwnership(h1, False);
h1.Draw("E0");
print(h1.GetBinContent(1), h1.GetBinError(1));

gPad.Modified();
gPad.Update();
