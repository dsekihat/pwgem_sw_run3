import ROOT
from ROOT import TH1D, TH2D, TH3D, TGraph, TGraphErrors, TGraphAsymmErrors
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross

def make_common_style(g1, marker, size, color, width=2, fill=0, linestyle=1):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetLineStyle(linestyle);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);
