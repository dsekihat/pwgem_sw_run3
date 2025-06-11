import datetime
import sys
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross

#__________________________________________________________
def extract_lf_cocktail(filename, taskname, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get("{0}/rec".format(taskname));
    #rootdire.ls();

    outname = "lmee_lf_cocktail" + suffix + ".root"
    outfile = TFile(outname, "RECREATE");
    parnames = ["sum", "pi0", "eta", "etaP", "rho", "omega", "phi"];
    print("input file = ", filename);
    print("output file = ", outname);

    for parname in parnames:
        print("extracting...", parname);
        if parname == "sum":
            hs = rootdire.Get("MeeVsPteeVsCos2DPhiRP");
        else:
            hs = rootdire.Get("{0}/MeeVsPteeVsCos2DPhiRP".format(parname));
        hs.SetName("hs_mee_ptee_v2ee_{0}".format(parname));
        hs.GetAxis(0).SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
        hs.GetAxis(1).SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
        hs.GetAxis(2).SetTitle("#it{v}_{2,ee}");
        hs.Scale(1/1e+6); #only temporary
        outfile.WriteTObject(hs);

    outfile.Close();
    rootfile.Close();
#__________________________________________________________
if __name__ == "__main__":
    filename = "AnalysisResults_HL_293890_3040.root";
    taskname = "em-lmee-lf-cocktail_minpt200";
    suffix = "_PbPb_5.36TeV_3040_minpt200_maxeta08";
    extract_lf_cocktail(filename, taskname, suffix);
    taskname = "em-lmee-lf-cocktail_minpt400";
    suffix = "_PbPb_5.36TeV_3040_minpt400_maxeta08";
    extract_lf_cocktail(filename, taskname, suffix);

    filename = "AnalysisResults_HL_293915_4050.root";
    taskname = "em-lmee-lf-cocktail_minpt200";
    suffix = "_PbPb_5.36TeV_4050_minpt200_maxeta08";
    extract_lf_cocktail(filename, taskname, suffix);
    taskname = "em-lmee-lf-cocktail_minpt400";
    suffix = "_PbPb_5.36TeV_4050_minpt400_maxeta08";
    extract_lf_cocktail(filename, taskname, suffix);

#__________________________________________________________
#__________________________________________________________
