import ROOT
from ROOT import TFile, TList

filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_253788.root";
taskname = "dielectron_v2_3050_TOFreq_minpt400";

rootfile = TFile.Open(filename, "READ");
roottask = rootfile.Get(taskname);
roottask.ls();


h1R2 = roottask.Get("h1R2");
h1R2.SetName("h1_R2_FT0M_BPos_BNeg");
outfile = TFile("20240905_sp_resolution.root", "RECREATE");
outlist = TList();
outlist.SetName("ccdb_object"); #Don't touch.
outlist.SetOwner(True);
outlist.Add(h1R2);
outfile.WriteTObject(outlist);
outfile.Close();
outlist.Clear();
