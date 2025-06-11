import numpy as np
import ROOT
from ROOT import TFile, TProfile3D

def rebin3d_profile(h3prf, arr_x, arr_y, arr_z):
    print("test");
    



if __name__ == "__main__":
    arr_mee = np.array([0, 0.14, 0.5, 1.1, 2.0, 2.7, 3.5], dtype=np.float32);
    arr_ptee = np.array([0, 8], dtype=np.float32);
    arr_dcaee = np.array([0, 10], dtype=np.float32);

    #filename = "AnalysisResults_test_dielectron_v2.root";
    #filename = "AnalysisResults_test_dimuon_v2.root";
    filename = "AnalysisResults_HL_239513.root";
    rootfile = TFile.Open(filename, "READ");
    #h3prf = rootfile.Get("dielectron-qc_3050_v2/Pair/same/uls/hPrfUQ");
    #h3prf = rootfile.Get("dielectron-qc_3050_v2/Pair/mix/uls/hPrfUQ_leg1");
    #h3prf = rootfile.Get("dielectron-qc_3050_v2/Pair/mix/uls/hPrfCosDPhi_leg1");

    h3prf = rootfile.Get("dimuon-qc_standalone_3050_v2/Pair/same/uls/hPrfUQ");
    #h3prf = rootfile.Get("dimuon-qc_standalone_3050_v2/Pair/mix/uls/hPrfUQ_leg1");
    #h3prf = rootfile.Get("dimuon-qc_standalone_3050_v2/Pair/mix/uls/hPrfCosDPhi_leg1");
    h3prf.SetDirectory(0);
    ROOT.SetOwnership(h3prf, False);
    #h3prf.Draw("colz");

    dca_min = 0.0;
    dca_max = 99.0;
    bin1 = h3prf.GetZaxis().FindBin(dca_min + 1e-3);
    bin2 = h3prf.GetZaxis().FindBin(dca_max - 1e-3);
    h3prf.GetZaxis().SetRange(bin1, bin2);

    h2prf = h3prf.Project3DProfile("yx");
    #h2prf.Draw("colz");

    h1prf = h2prf.ProfileX("h1prf", 2, 4, "");
    h1prf.Draw();
    #h3 = rebin3d_profile(h3prf, arr_mee, arr_ptee, arr_dcaee);
