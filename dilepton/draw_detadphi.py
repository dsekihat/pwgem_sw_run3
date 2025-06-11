import datetime
import sys
import math
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D, TColor
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
gStyle.SetPalette(55);
#gStyle.SetOptStat(0);
#gStyle.SetOptTitle(0);

##gStyle.SetPalette(70); # kBlackBody
##TColor.InvertPalette();

#__________________________________________________________
def draw_ft0c_occupancy(filename, tasknames, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();

    list_h2 = [];
    for taskname in tasknames:
        rootdire_ev = rootfile.Get(taskname + "/Event/after");
        #rootdire_ev.ls();
        list_h2.append(rootdire_ev.Get("hCorrOccupancy"));

    h2 = list_h2[0].Clone("h2");
    h2.Sumw2();
    h2.SetContour(1000);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    for i in range(1, len(list_h2)):
        h2.Add(list_h2[i], 1);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.12, 0.1, 0.05);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0, 0, 200e+3, 20e+3);
    frame1.GetXaxis().SetTitle("FT0C occupancy");
    frame1.GetYaxis().SetTitle("Track occupancy");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetXaxis().SetMaxDigits(4);
    frame1.GetYaxis().SetMaxDigits(4);
    ROOT.SetOwnership(frame1, False);
    h2.Draw("colz,same");

    c1.Modified();
    c1.Update();
    palette = h2.GetListOfFunctions().FindObject("palette");
    #print(palette.GetX1NDC());
    #print(palette.GetX2NDC());
    #print(palette.GetY1NDC());
    #print(palette.GetY2NDC());
    #palette.SetX1NDC(0.9);
    #palette.SetX2NDC(0.99);
    palette.SetY1NDC(0.16);
    #palette.SetY2NDC(0.95);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_occupancy_correlation{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_occupancy_correlation{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_occupancy_correlation{1}.png".format(date, suffix));

#__________________________________________________________
def draw_dpt(filename, taskname, deta_min, deta_max, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    hs_uls_same  = rootdire.Get("Pair/same/uls/hsDeltaP");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hsDeltaP");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hsDeltaP");
    hs_uls_mix   = rootdire.Get("Pair/mix/uls/hsDeltaP");
    hs_lspp_mix  = rootdire.Get("Pair/mix/lspp/hsDeltaP");
    hs_lsmm_mix  = rootdire.Get("Pair/mix/lsmm/hsDeltaP");

    hs_uls_same .Sumw2();
    hs_lspp_same.Sumw2();
    hs_lsmm_same.Sumw2();
    hs_uls_mix  .Sumw2();
    hs_lspp_mix .Sumw2();
    hs_lsmm_mix .Sumw2();

    bin0_pos = hs_uls_same.GetAxis(1).FindBin(deta_min + 1e-6);
    bin1_pos = hs_uls_same.GetAxis(1).FindBin(999 - 1e-6);
    #bin1_pos = hs_uls_same.GetAxis(1).FindBin(deta_max - 1e-6);

    hs_uls_same.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lspp_same.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lsmm_same.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_uls_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lspp_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lsmm_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);

    h1_uls_same_pos  = hs_uls_same.Projection(0).Clone("h1same_pos");
    h1_lspp_same_pos = hs_lspp_same.Projection(0).Clone("h1same_pos");
    h1_lsmm_same_pos = hs_lsmm_same.Projection(0).Clone("h1same_pos");
    h1_uls_mix_pos   = hs_uls_mix.Projection(0).Clone("h1mix_pos");
    h1_lspp_mix_pos  = hs_lspp_mix.Projection(0).Clone("h1mix_pos");
    h1_lsmm_mix_pos  = hs_lsmm_mix.Projection(0).Clone("h1mix_pos");

    bin0_neg = hs_uls_same.GetAxis(1).FindBin(- 999 + 1e-6);
    bin1_neg = hs_uls_same.GetAxis(1).FindBin(- deta_min - 1e-6);
    #bin1_neg = hs_uls_same.GetAxis(1).FindBin(deta_max - 1e-6);

    hs_uls_same.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_lspp_same.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_lsmm_same.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_uls_mix.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_lspp_mix.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_lsmm_mix.GetAxis(1).SetRange(bin0_neg, bin1_neg);

    h1_uls_same_neg  = hs_uls_same.Projection(0).Clone("h1same_neg");
    h1_lspp_same_neg = hs_lspp_same.Projection(0).Clone("h1same_neg");
    h1_lsmm_same_neg = hs_lsmm_same.Projection(0).Clone("h1same_neg");
    h1_uls_mix_neg   = hs_uls_mix.Projection(0).Clone("h1mix_neg");
    h1_lspp_mix_neg  = hs_lspp_mix.Projection(0).Clone("h1mix_neg");
    h1_lsmm_mix_neg  = hs_lsmm_mix.Projection(0).Clone("h1mix_neg");

    h1_uls_same  = h1_uls_same_pos .Clone("h1_uls_same");
    h1_lspp_same = h1_lspp_same_pos.Clone("h1_lspp_same");
    h1_lsmm_same = h1_lsmm_same_pos.Clone("h1_lsmm_same");
    h1_uls_mix   = h1_uls_mix_pos  .Clone("h1_uls_mix");
    h1_lspp_mix  = h1_lspp_mix_pos .Clone("h1_lspp_mix");
    h1_lsmm_mix  = h1_lsmm_mix_pos .Clone("h1_lsmm_mix");

    h1_uls_same  .Add(h1_uls_same_pos , 1);
    h1_lspp_same .Add(h1_lspp_same_pos, 1);
    h1_lsmm_same .Add(h1_lsmm_same_pos, 1);
    h1_uls_mix   .Add(h1_uls_mix_pos  , 1);
    h1_lspp_mix  .Add(h1_lspp_mix_pos , 1);
    h1_lsmm_mix  .Add(h1_lsmm_mix_pos , 1);

    h1_uls_same  .RebinX(5);
    h1_lspp_same .RebinX(5);
    h1_lsmm_same .RebinX(5);
    h1_uls_mix   .RebinX(5);
    h1_lspp_mix  .RebinX(5);
    h1_lsmm_mix  .RebinX(5);

    sf_uls = h1_uls_same.GetEntries()/h1_uls_mix.GetEntries();
    h1_uls_mix.Scale(sf_uls);
    sf_lspp = h1_lspp_same.GetEntries()/h1_lspp_mix.GetEntries();
    h1_lspp_mix.Scale(sf_lspp);
    sf_lsmm = h1_lsmm_same.GetEntries()/h1_lsmm_mix.GetEntries();
    h1_lsmm_mix.Scale(sf_lsmm);

    h1_uls_ratio = h1_uls_same.Clone("h1_uls_ratio");
    h1_uls_ratio.Divide(h1_uls_mix);
    h1_uls_ratio.SetContour(1000);
    h1_uls_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_ratio, False);
    make_common_style(h1_uls_ratio, 24, 1.2, kBlack, 2, 0);

    h1_lspp_ratio = h1_lspp_same.Clone("h1_lspp_ratio");
    h1_lspp_ratio.Divide(h1_lspp_mix);
    h1_lspp_ratio.SetContour(1000);
    h1_lspp_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lspp_ratio, False);
    make_common_style(h1_lspp_ratio, 22, 1.5, kRed+1, 2, 0);

    h1_lsmm_ratio = h1_lsmm_same.Clone("h1_lsmm_ratio");
    h1_lsmm_ratio.Divide(h1_lsmm_mix);
    h1_lsmm_ratio.SetContour(1000);
    h1_lsmm_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lsmm_ratio, False);
    make_common_style(h1_lsmm_ratio, 23, 1.5, kBlue+1, 2, 0);

    h1_uls_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_ratio, False);
    h1_lspp_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lspp_ratio, False);
    h1_lsmm_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lsmm_ratio, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);

    frame1 = c1.DrawFrame(0, 0.98, 1, 1.05);
    frame1.GetXaxis().SetTitle("p_{T} asymmetry");
    frame1.GetYaxis().SetTitle("same/mix");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h1_uls_ratio .Draw("hist,E0same");
    h1_lspp_ratio.Draw("hist,E0same");
    h1_lsmm_ratio.Draw("hist,E0same");

    leg = TLegend(0.6, 0.7, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    #leg.SetHeader("{0:3.2f} < #Delta#eta < {1:3.2f}".format(deta_min, deta_max));
    leg.SetHeader("|#Delta#eta| > {0:3.2f}".format(deta_min));
    leg.AddEntry(h1_uls_ratio, "ULS", "P");
    leg.AddEntry(h1_lspp_ratio, "LS++", "P");
    leg.AddEntry(h1_lsmm_ratio, "LS--", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.eps".format(date, deta_min, deta_max, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.pdf".format(date, deta_min, deta_max, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.png".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_dpt_deta{1:3.2f}{2}.eps".format(date, deta_min, suffix));
    c1.SaveAs("{0}_dpt_deta{1:3.2f}{2}.pdf".format(date, deta_min, suffix));
    c1.SaveAs("{0}_dpt_deta{1:3.2f}{2}.png".format(date, deta_min, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dphi_yield(filename, taskname, sign, deta_min, deta_max, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    hs_uls_same  = rootdire.Get("Pair/same/{0}/hsDeltaP".format(sign));
    hs_uls_mix   = rootdire.Get("Pair/mix/{0}/hsDeltaP".format(sign));

    hs_uls_same .Sumw2();
    hs_uls_mix  .Sumw2();

    bin0_pos = hs_uls_same.GetAxis(1).FindBin(deta_min + 1e-6);
    bin1_pos = hs_uls_same.GetAxis(1).FindBin(999 - 1e-6);
    #bin1_pos = hs_uls_same.GetAxis(1).FindBin(deta_max - 1e-6);

    hs_uls_same.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_uls_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);

    h1_uls_same_pos  = hs_uls_same.Projection(2).Clone("h1same_pos");
    h1_uls_mix_pos   = hs_uls_mix.Projection(2).Clone("h1mix_pos");

    bin0_neg = hs_uls_same.GetAxis(1).FindBin(- 999 + 1e-6);
    bin1_neg = hs_uls_same.GetAxis(1).FindBin(- deta_min - 1e-6);
    #bin1_neg = hs_uls_same.GetAxis(1).FindBin(deta_max - 1e-6);

    hs_uls_same.GetAxis(1).SetRange(bin0_neg, bin1_neg);
    hs_uls_mix.GetAxis(1).SetRange(bin0_neg, bin1_neg);

    h1_uls_same_neg  = hs_uls_same.Projection(2).Clone("h1same_neg");
    h1_uls_mix_neg   = hs_uls_mix.Projection(2).Clone("h1mix_neg");

    h1_uls_same  = h1_uls_same_pos .Clone("h1_uls_same");
    h1_uls_mix   = h1_uls_mix_pos  .Clone("h1_uls_mix");

    h1_uls_same  .Add(h1_uls_same_pos , 1);
    h1_uls_mix   .Add(h1_uls_mix_pos  , 1);

    h1_uls_same.RebinX(1);
    h1_uls_mix .RebinX(1);

    h1_uls_same.Scale(1/h1_uls_same.GetEntries());
    h1_uls_mix.Scale(1/h1_uls_mix.GetEntries());
    h1_uls_same.Scale(1, "width");
    h1_uls_mix .Scale(1, "width");

    print(h1_uls_same.GetMaximum());

    h1_uls_same.SetDirectory(0);
    h1_uls_mix.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_same, False);
    ROOT.SetOwnership(h1_uls_mix , False);
    make_common_style(h1_uls_same, 20, 1.2, kRed+1 , 2, 0);
    make_common_style(h1_uls_mix , 20, 1.2, kBlue+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.14, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
    #c1.SetLogy(1);

    frame1 = c1.DrawFrame(-math.pi, 0.15, +math.pi, 0.18);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("dN/d#Delta#varphi (a.u.)");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h1_uls_same.Draw("hist,E0same");
    h1_uls_mix.Draw("hist,E0same");

    leg = TLegend(0.6, 0.7, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);

    if "uls" in sign:
        leg.SetHeader("ULS, |#Delta#eta| > {0:3.2f}".format(deta_min));
    if "lspp" in sign:
        leg.SetHeader("LS++, |#Delta#eta| > {0:3.2f}".format(deta_min));
    if "lsmm" in sign:
        leg.SetHeader("LS--, |#Delta#eta| > {0:3.2f}".format(deta_min));
    leg.AddEntry(h1_uls_same, "same", "P");
    leg.AddEntry(h1_uls_mix, "mixed", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.eps".format(date, deta_min, deta_max, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.pdf".format(date, deta_min, deta_max, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.png".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_dphi_{1}_yield_deta{2:3.2f}{3}.eps".format(date, sign, deta_min, suffix));
    c1.SaveAs("{0}_dphi_{1}_yield_deta{2:3.2f}{3}.pdf".format(date, sign, deta_min, suffix));
    c1.SaveAs("{0}_dphi_{1}_yield_deta{2:3.2f}{3}.png".format(date, sign, deta_min, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_rfactor_vs_dphi(filename, taskname, deta_min, deta_max, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    hs_uls_mix   = rootdire.Get("Pair/mix/uls/hsDeltaP");
    hs_lspp_mix  = rootdire.Get("Pair/mix/lspp/hsDeltaP");
    hs_lsmm_mix  = rootdire.Get("Pair/mix/lsmm/hsDeltaP");

    hs_uls_mix  .Sumw2();
    hs_lspp_mix .Sumw2();
    hs_lsmm_mix .Sumw2();

    bin0_pos = hs_uls_mix.GetAxis(1).FindBin(deta_min + 1e-6);
    bin1_pos = hs_uls_mix.GetAxis(1).FindBin(deta_max - 1e-6);

    hs_uls_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lspp_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);
    hs_lsmm_mix.GetAxis(1).SetRange(bin0_pos, bin1_pos);

    h1_uls_mix = hs_uls_mix .Projection(2).Clone("h1_uls_mix");
    h1_lspp_mix= hs_lspp_mix.Projection(2).Clone("h1_lspp_mix");
    h1_lsmm_mix= hs_lsmm_mix.Projection(2).Clone("h1_lsmm_mix");
    h1_uls_mix   .RebinX(1);
    h1_lspp_mix  .RebinX(1);
    h1_lsmm_mix  .RebinX(1);

    h1r = h1_uls_mix.Clone("h1r");
    h1r.Reset();
    for i in range(0, h1_uls_mix.GetNbinsX()):
        n_uls  = h1_uls_mix .GetBinContent(i+1);
        n_lspp = h1_lspp_mix.GetBinContent(i+1);
        n_lsmm = h1_lsmm_mix.GetBinContent(i+1);
        r = n_uls / (2 * math.sqrt(n_lspp * n_lsmm));
        h1r.SetBinContent(i+1, r);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);

    frame1 = c1.DrawFrame(-math.pi, 0.95, +math.pi, 1.05);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("R factor");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h1r.Draw("E0same");
    make_common_style(h1r, 24, 1.2, kBlack, 2, 0);

    leg = TLegend(0.6, 0.9, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.SetHeader("{0:3.2f} < #Delta#eta < {1:3.2f}".format(deta_min, deta_max));
    leg.Draw("");
    ROOT.SetOwnership(leg,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_rfactor_vs_dphi_deta{1:3.2f}_{1:3.2f}{2}.eps".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_rfactor_vs_dphi_deta{1:3.2f}_{1:3.2f}{2}.pdf".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_rfactor_vs_dphi_deta{1:3.2f}_{1:3.2f}{2}.png".format(date, deta_min, deta_max, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dphi(filename, taskname, deta_min, deta_max, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    #hs_uls_same  = rootdire.Get("Pair/same/uls/hsDeltaP");
    #hs_lspp_same = rootdire.Get("Pair/same/lspp/hsDeltaP");
    #hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hsDeltaP");
    #hs_uls_mix   = rootdire.Get("Pair/mix/uls/hsDeltaP");
    #hs_lspp_mix  = rootdire.Get("Pair/mix/lspp/hsDeltaP");
    #hs_lsmm_mix  = rootdire.Get("Pair/mix/lsmm/hsDeltaP");

    hs_uls_same  = rootdire.Get("Pair/same/uls/hDeltaEtaDeltaPhi");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hDeltaEtaDeltaPhi");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hDeltaEtaDeltaPhi");
    hs_uls_mix   = rootdire.Get("Pair/mix/uls/hDeltaEtaDeltaPhi");
    hs_lspp_mix  = rootdire.Get("Pair/mix/lspp/hDeltaEtaDeltaPhi");
    hs_lsmm_mix  = rootdire.Get("Pair/mix/lsmm/hDeltaEtaDeltaPhi");

    hs_uls_same .Sumw2();
    hs_lspp_same.Sumw2();
    hs_lsmm_same.Sumw2();
    hs_uls_mix  .Sumw2();
    hs_lspp_mix .Sumw2();
    hs_lsmm_mix .Sumw2();

    bin0_pos = hs_uls_same.GetYaxis().FindBin(deta_min + 1e-6);
    bin1_pos = hs_uls_same.GetYaxis().FindBin(deta_max - 1e-6);

    hs_uls_same .GetYaxis().SetRange(bin0_pos, bin1_pos);
    hs_lspp_same.GetYaxis().SetRange(bin0_pos, bin1_pos);
    hs_lsmm_same.GetYaxis().SetRange(bin0_pos, bin1_pos);
    hs_uls_mix  .GetYaxis().SetRange(bin0_pos, bin1_pos);
    hs_lspp_mix .GetYaxis().SetRange(bin0_pos, bin1_pos);
    hs_lsmm_mix .GetYaxis().SetRange(bin0_pos, bin1_pos);

    h1_uls_same_pos  = hs_uls_same .ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1same_pos");
    h1_lspp_same_pos = hs_lspp_same.ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1same_pos");
    h1_lsmm_same_pos = hs_lsmm_same.ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1same_pos");
    h1_uls_mix_pos   = hs_uls_mix  .ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1mix_pos");
    h1_lspp_mix_pos  = hs_lspp_mix .ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1mix_pos");
    h1_lsmm_mix_pos  = hs_lsmm_mix .ProjectionX("h1tmp", bin0_pos, bin1_pos, "").Clone("h1mix_pos");

    bin0_neg = hs_uls_same.GetYaxis().FindBin(- deta_max + 1e-6);
    bin1_neg = hs_uls_same.GetYaxis().FindBin(- deta_min - 1e-6);

    hs_uls_same .GetYaxis().SetRange(bin0_neg, bin1_neg);
    hs_lspp_same.GetYaxis().SetRange(bin0_neg, bin1_neg);
    hs_lsmm_same.GetYaxis().SetRange(bin0_neg, bin1_neg);
    hs_uls_mix  .GetYaxis().SetRange(bin0_neg, bin1_neg);
    hs_lspp_mix .GetYaxis().SetRange(bin0_neg, bin1_neg);
    hs_lsmm_mix .GetYaxis().SetRange(bin0_neg, bin1_neg);

    h1_uls_same_neg  = hs_uls_same .ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1same_neg");
    h1_lspp_same_neg = hs_lspp_same.ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1same_neg");
    h1_lsmm_same_neg = hs_lsmm_same.ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1same_neg");
    h1_uls_mix_neg   = hs_uls_mix  .ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1mix_neg");
    h1_lspp_mix_neg  = hs_lspp_mix .ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1mix_neg");
    h1_lsmm_mix_neg  = hs_lsmm_mix .ProjectionY("h1tmp", bin0_neg, bin1_neg, "").Clone("h1mix_neg");

    h1_uls_same  = h1_uls_same_pos .Clone("h1_uls_same");
    h1_lspp_same = h1_lspp_same_pos.Clone("h1_lspp_same");
    h1_lsmm_same = h1_lsmm_same_pos.Clone("h1_lsmm_same");
    h1_uls_mix   = h1_uls_mix_pos  .Clone("h1_uls_mix");
    h1_lspp_mix  = h1_lspp_mix_pos .Clone("h1_lspp_mix");
    h1_lsmm_mix  = h1_lsmm_mix_pos .Clone("h1_lsmm_mix");

    h1_uls_same  .Add(h1_uls_same_pos , 1);
    h1_lspp_same .Add(h1_lspp_same_pos, 1);
    h1_lsmm_same .Add(h1_lsmm_same_pos, 1);
    h1_uls_mix   .Add(h1_uls_mix_pos  , 1);
    h1_lspp_mix  .Add(h1_lspp_mix_pos , 1);
    h1_lsmm_mix  .Add(h1_lsmm_mix_pos , 1);

    h1_uls_same  .RebinX(2);
    h1_lspp_same .RebinX(2);
    h1_lsmm_same .RebinX(2);
    h1_uls_mix   .RebinX(2);
    h1_lspp_mix  .RebinX(2);
    h1_lsmm_mix  .RebinX(2);

    sf_uls = h1_uls_same.GetEntries()/h1_uls_mix.GetEntries();
    h1_uls_mix.Scale(sf_uls);
    sf_lspp = h1_lspp_same.GetEntries()/h1_lspp_mix.GetEntries();
    h1_lspp_mix.Scale(sf_lspp);
    sf_lsmm = h1_lsmm_same.GetEntries()/h1_lsmm_mix.GetEntries();
    h1_lsmm_mix.Scale(sf_lsmm);

    h1_uls_ratio = h1_uls_same.Clone("h1_uls_ratio");
    h1_uls_ratio.Divide(h1_uls_mix);
    h1_uls_ratio.SetContour(1000);
    h1_uls_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_ratio, False);
    make_common_style(h1_uls_ratio, 24, 1.2, kBlack, 2, 0);

    h1_lspp_ratio = h1_lspp_same.Clone("h1_lspp_ratio");
    h1_lspp_ratio.Divide(h1_lspp_mix);
    h1_lspp_ratio.SetContour(1000);
    h1_lspp_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lspp_ratio, False);
    make_common_style(h1_lspp_ratio, 22, 1.5, kRed+1, 2, 0);

    h1_lsmm_ratio = h1_lsmm_same.Clone("h1_lsmm_ratio");
    h1_lsmm_ratio.Divide(h1_lsmm_mix);
    h1_lsmm_ratio.SetContour(1000);
    h1_lsmm_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lsmm_ratio, False);
    make_common_style(h1_lsmm_ratio, 23, 1.5, kBlue+1, 2, 0);

    h1_uls_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_ratio, False);
    h1_lspp_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lspp_ratio, False);
    h1_lsmm_ratio.SetDirectory(0);
    ROOT.SetOwnership(h1_lsmm_ratio, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetGridx(1);
    #c1.SetGridy(1);

    frame1 = c1.DrawFrame(-math.pi, 0.98, +math.pi, 1.05);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("same/mix");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h1_uls_ratio .Draw("hist,E0same");
    h1_lspp_ratio.Draw("hist,E0same");
    h1_lsmm_ratio.Draw("hist,E0same");

    leg = TLegend(0.7, 0.7, 0.9, 0.90);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.SetHeader("{0:2.1f} < |#Delta#eta| < {1:2.1f}".format(deta_min, deta_max));
    #leg.SetHeader("|#Delta#eta| > {0:3.2f}".format(deta_min));
    leg.AddEntry(h1_uls_ratio, "ULS", "P");
    leg.AddEntry(h1_lspp_ratio, "LS++", "P");
    leg.AddEntry(h1_lsmm_ratio, "LS--", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    cen1 = 0;
    cen2 = 0;
    if "0090" in suffix:
        cen1 = 0;
        cen2 = 90;
    elif "1090" in suffix:
        cen1 = 10;
        cen2 = 90;
    elif "0010" in suffix:
        cen1 = 0;
        cen2 = 10;

    txt = TPaveText(0.15, 0.9, 0.45, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.eps".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.pdf".format(date, deta_min, deta_max, suffix));
    c1.SaveAs("{0}_dphi_deta{1:3.2f}_{2:3.2f}{3}.png".format(date, deta_min, deta_max, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}{2}.eps".format(date, deta_min, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}{2}.pdf".format(date, deta_min, suffix));
    #c1.SaveAs("{0}_dphi_deta{1:3.2f}{2}.png".format(date, deta_min, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_100(filename, taskname, phiv1, phiv2, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    sign = "uls";
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hMvsPhiV").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hMvsPhiV").Clone("h2mix");

    bin1 = h2same.GetXaxis().FindBin(phiv1 + 1e-6);
    bin2 = h2same.GetXaxis().FindBin(phiv2 - 1e-6);
    h1same = h2same.ProjectionY("h1same", bin1, bin2, "");
    h1mix  = h2mix .ProjectionY("h1mix" , bin1, bin2, "");
    #h1same.RebinX(1);
    #h1mix .RebinX(1);
    print(bin1, bin2);

    h1same.SetDirectory(0);
    h1mix .SetDirectory(0);
    ROOT.SetOwnership(h1same, False);
    ROOT.SetOwnership(h1mix , False);

    sf = h1same.GetEntries()/h1mix.GetEntries();
    h1mix.Scale(sf);

    h1ratio = h1same.Clone("h1ratio");
    h1ratio.Divide(h1mix);
    h1ratio.SetContour(1000);
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);
    make_common_style(h1ratio, 20, 1.2, kBlack, 2, 0);

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.04, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetGridx(1);
    #c1.SetGridy(1);
    #c1.SetLogy(1);

    frame1 = c1.DrawFrame(0., 0., 0.02, 2e+3);
    #frame1 = c1.DrawFrame(0., 0., 0.02, 30.0);
    frame1.GetXaxis().SetTitle("m_{ee} (GeV/c^{2})");
    frame1.GetYaxis().SetTitle("same/mix");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetXaxis().SetNdivisions(505, True);
    ROOT.SetOwnership(frame1, False);
    h1ratio.Draw("hist,E0same");

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_low_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.eps".format(date, sign, phiv1, phiv2, suffix));
    c1.SaveAs("{0}_{1}_low_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.pdf".format(date, sign, phiv1, phiv2, suffix));
    c1.SaveAs("{0}_{1}_low_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.png".format(date, sign, phiv1, phiv2, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee(filename, taskname, phiv1, phiv2, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    sign = "uls";
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hMvsPhiV").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hMvsPhiV").Clone("h2mix");

    bin1 = h2same.GetXaxis().FindBin(phiv1 + 1e-6);
    bin2 = h2same.GetXaxis().FindBin(phiv2 - 1e-6);
    h1same = h2same.ProjectionY("h1same", bin1, bin2, "");
    h1mix  = h2mix .ProjectionY("h1mix" , bin1, bin2, "");
    #h1same.RebinX(1);
    #h1mix .RebinX(1);
    print(bin1, bin2);

    h1same.SetDirectory(0);
    h1mix .SetDirectory(0);
    ROOT.SetOwnership(h1same, False);
    ROOT.SetOwnership(h1mix , False);

    sf = h1same.GetEntries()/h1mix.GetEntries();
    h1mix.Scale(sf);

    h1ratio = h1same.Clone("h1ratio");
    h1ratio.Divide(h1mix);
    h1ratio.SetContour(1000);
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);
    make_common_style(h1ratio, 20, 1.2, kBlack, 2, 0);

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetGridx(1);
    #c1.SetGridy(1);
    #c1.SetLogy(1);

    frame1 = c1.DrawFrame(0., 0.8, 1.0, 2.4);
    #frame1 = c1.DrawFrame(0., 0.6, 0.5, 3.0);
    frame1.GetXaxis().SetTitle("m_{ee} (GeV/c^{2})");
    frame1.GetYaxis().SetTitle("same/mix");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h1ratio.Draw("hist,E0same");

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.eps".format(date, sign, phiv1, phiv2, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.pdf".format(date, sign, phiv1, phiv2, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv{2:3.2f}_{3:3.2f}rad{4}.png".format(date, sign, phiv1, phiv2, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_meephiv_prop(filename, taskname, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hMvsPhiV_prop").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hMvsPhiV_prop").Clone("h2mix");

    h2same.Rebin2D(1, 1);
    h2mix .Rebin2D(1, 1);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.5);
    #h2ratio.SetMinimum(0.7);
    h2ratio.SetMaximum(1.5);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.1);

    #gStyle.SetPalette(55);
    #gStyle.SetPalette(52); #kGreyScale
    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    #gStyle.SetPalette(107); # kVisibleSpectrum

    #gStyle.SetPalette(70); # kBlackBody
    #TColor.InvertPalette();
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.13, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0., 0, math.pi, 1.0);
    frame1.GetXaxis().SetTitle("#varphi_{V} (rad.)");
    frame1.GetYaxis().SetTitle("m_{ee} (GeV/c^{2})");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_mee_phiv_prop{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv_prop{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv_prop{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_opangle(filename, taskname, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hMvsOpAng").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hMvsOpAng").Clone("h2mix");

    h2same.Rebin2D(1, 1);
    h2mix .Rebin2D(1, 1);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.9);
    h2ratio.SetMaximum(1.1);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.4);

    #gStyle.SetPalette(55);
    #gStyle.SetPalette(52); #kGreyScale
    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    #gStyle.SetPalette(107); # kVisibleSpectrum

    #gStyle.SetPalette(70); # kBlackBody
    #TColor.InvertPalette();
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.15, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0., 0, math.pi, 1.0);
    frame1.GetXaxis().SetTitle("#omega_{ee} (rad.)");
    frame1.GetYaxis().SetTitle("m_{ee} (GeV/c^{2})");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_mee_opang{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_opang{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_opang{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_meephiv(filename, taskname, cen1, cen2, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hMvsPhiV").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hMvsPhiV").Clone("h2mix");

    h2same.Rebin2D(1, 1);
    h2mix .Rebin2D(1, 1);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.9);
    h2ratio.SetMaximum(1.1);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.4);

    #gStyle.SetPalette(55);
    #gStyle.SetPalette(52); #kGreyScale
    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    #gStyle.SetPalette(107); # kVisibleSpectrum

    #gStyle.SetPalette(70); # kBlackBody
    #TColor.InvertPalette();
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.15, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0., 0, math.pi, 1.0);
    frame1.GetXaxis().SetTitle("#varphi_{V} (rad.)");
    frame1.GetYaxis().SetTitle("m_{ee} (GeV/c^{2})");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    txt = TPaveText(0.15, 0.85, 0.45, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);

    energy = 13.6;
    if "13.6" in suffix:
        energy = 13.6;
    elif "5.36" in suffix:
        energy = 5.36;

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{{#it{{s}}}} = {0:3.2f} TeV".format(energy));
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = {2:3.2f} TeV".format(cen1, cen2, energy));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = {2:3.2f} TeV".format(cen1, cen2, energy));
    else:
        txt.AddText("pp at #sqrt{{#it{{s}}}} = {0:3.2f} TeV".format(energy));

    if "uls" in sign:
        txt.AddText("ULS");
    elif "lspp" in sign:
        txt.AddText("LS++");
    elif "lsmm" in sign:
        txt.AddText("LS#minus#minus");

    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #f1 = TF1("f1", "[0] + [1] * x", 0, math.pi);
    #f1.SetNpx(1000);
    #f1.SetParameters(-3.7, 1.40);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_mee_phiv{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_mee_phiv{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_dzrdphi_geom(filename, taskname, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hGeomDeltaZRDeltaPhi").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hGeomDeltaZRDeltaPhi").Clone("h2mix");

    h2same.Rebin2D(1, 1);
    h2mix .Rebin2D(1, 1);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.9);
    h2ratio.SetMaximum(2.0);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.5);

    #gStyle.SetPalette(70); # kBlackBody
    #TColor.InvertPalette();

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.16, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-100, -20, 100, +20);
    frame1.GetXaxis().SetTitle("r#Delta#varphi (cm)");
    frame1.GetYaxis().SetTitle("#Deltaz (cm)");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.2);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    min_dz = 4.0;
    min_rdphi = 10.0;
    f1up = TF1("f1up", "[1] * sqrt(1 - x*x/[0]/[0])", -min_rdphi, +min_rdphi);
    f1up.SetNpx(8000);
    f1up.SetParameters(min_rdphi, min_dz);
    f1up.Draw("same");
    ROOT.SetOwnership(f1up, False);
    f1down = TF1("f1down", "-[1] * sqrt(1 - x*x/[0]/[0])", -min_rdphi, +min_rdphi);
    f1down.SetNpx(8000);
    f1down.SetParameters(min_rdphi, min_dz);
    f1down.Draw("same");
    ROOT.SetOwnership(f1down, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_geom_dz_rdphi{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_geom_dz_rdphi{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_geom_dz_rdphi{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_detadphi_geom(filename, taskname, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hGeomDeltaEtaDeltaPhi").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hGeomDeltaEtaDeltaPhi").Clone("h2mix");

    #h2same.Rebin2D(2, 2);
    #h2mix .Rebin2D(2, 2);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.8);
    h2ratio.SetMaximum(1.2);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.5);

    #gStyle.SetPalette(70); # kBlackBody
    #TColor.InvertPalette();

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.16, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-math.pi/2, -0.5, +math.pi/2, +0.5);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("#Delta#eta");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    min_deta = 0.1;
    min_dphi = 0.4;

    f1up = TF1("f1up", "[1] * sqrt(1 - x*x/[0]/[0])", -0.4, +0.4);
    f1up.SetNpx(8000);
    f1up.SetParameters(min_dphi, min_deta);
    f1up.Draw("same");
    ROOT.SetOwnership(f1up, False);

    f1down = TF1("f1down", "-[1] * sqrt(1 - x*x/[0]/[0])", -0.4, +0.4);
    f1down.SetNpx(8000);
    f1down.SetParameters(min_dphi, min_deta);
    f1down.Draw("same");
    ROOT.SetOwnership(f1down, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_geom_deta_dphi{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_geom_deta_dphi{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_geom_deta_dphi{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_detadphi(filename, taskname, cen1, cen2, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();

    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    #rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    h2same = rootdire_pair_same.Get("hDeltaEtaDeltaPhi").Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hDeltaEtaDeltaPhi").Clone("h2mix");
    #h2same = rootdire_pair_same.Get("hsDeltaP").Projection(1,2).Clone("h2same");
    #h2mix  = rootdire_pair_mix .Get("hsDeltaP").Projection(1,2).Clone("h2mix");

    h2same.Rebin2D(1, 1);
    h2mix .Rebin2D(1, 1);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.9);
    h2ratio.SetMaximum(1.1);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.5);

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.16, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-math.pi, -2, +math.pi, +2);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("#Delta#eta");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    txt = TPaveText(0.15, 0.85, 0.45, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);

    energy = 13.6;
    if "13.6" in suffix:
        energy = 13.6;
    elif "5.36" in suffix:
        energy = 5.36;

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{{#it{{s}}}} = {0:3.2f} TeV".format(energy));
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = {2:3.2f} TeV".format(cen1, cen2, energy));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}%, p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = {2:3.2f} TeV".format(cen1, cen2, energy));
    else:
        txt.AddText("pp at #sqrt{{#it{{s}}}} = {0:3.2f} TeV".format(energy));

    if "uls" in sign:
        txt.AddText("ULS");
    elif "lspp" in sign:
        txt.AddText("LS++");
    elif "lsmm" in sign:
        txt.AddText("LS#minus#minus");

    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #min_dphi = 100;
    #min_deta = 0.02;
    #f1up = TF1("f1up", "[1] * sqrt(1 - x*x/[0]/[0])", -min_dphi, +min_dphi);
    #f1up.SetNpx(1000);
    #f1up.SetParameters(min_dphi, min_deta);
    #f1up.Draw("same");
    #ROOT.SetOwnership(f1up, False);
    #f1down = TF1("f1down", "-[1] * sqrt(1 - x*x/[0]/[0])", -min_dphi, +min_dphi);
    #f1down.SetNpx(1000);
    #f1down.SetParameters(min_dphi, min_deta);
    #f1down.Draw("same");
    #ROOT.SetOwnership(f1down, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_deta_dphi{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_deta_dphi{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_deta_dphi{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dptdphi(filename, taskname, sign, suffix):
    rootfile = TFile(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    rootdire_pair_same = rootdire.Get("Pair/same/{0}".format(sign));
    rootdire_pair_mix = rootdire.Get("Pair/mix/{0}".format(sign));
    #rootdire_pair_same.ls();
    #rootdire_pair_mix.ls();

    #h2same = rootdire_pair_same.Get("hDeltaEtaDeltaPhi").Clone("h2same");
    #h2mix  = rootdire_pair_mix .Get("hDeltaEtaDeltaPhi").Clone("h2mix");
    h2same = rootdire_pair_same.Get("hsDeltaP").Projection(0,2).Clone("h2same");
    h2mix  = rootdire_pair_mix .Get("hsDeltaP").Projection(0,2).Clone("h2mix");

    #h2same.Rebin2D(2, 2);
    #h2mix .Rebin2D(2, 2);
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix , False);

    sf = h2same.GetEntries()/h2mix.GetEntries();
    h2mix.Scale(sf);

    h2ratio = h2same.Clone("h2ratio");
    h2ratio.Divide(h2mix);
    h2ratio.SetContour(1000);
    h2ratio.SetDirectory(0);
    ROOT.SetOwnership(h2ratio, False);
    h2ratio.SetMinimum(0.8);
    h2ratio.SetMaximum(1.2);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.5);

    gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.16, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-math.pi/2, 0., math.pi/2, 1.0);
    frame1.GetXaxis().SetTitle("#Delta#varphi (rad.)");
    frame1.GetYaxis().SetTitle("#Deltap_{T} (GeV/c)");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    h2ratio.Draw("colz,same");

    #f1 = TF1("f1", "[1] * sqrt(1 - x*x/[0]/[0])", -0.25, +0.25);
    #f1.SetNpx(1000);
    #f1.SetParameters(0.25, 0.08);
    #f1.Draw("same");
    #ROOT.SetOwnership(f1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_dpt_dphi{2}.eps".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_dpt_dphi{2}.pdf".format(date, sign, suffix));
    c1.SaveAs("{0}_{1}_dpt_dphi{2}.png".format(date, sign, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_279546.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #filename = "AnalysisResults_HL_280731.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "";
    #draw_detadphi(filename, taskname, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_tpcshared04";
    #draw_detadphi(filename, taskname, suffix);

    #filename = "AnalysisResults_HL_280891.root";
    #taskname = "dielectron_TPChadrejorTOFreq_minpt200";
    #suffix = "_pp_13.6TeV";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    #filename = "AnalysisResults_HL_281955.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_wo_TPCshared03_HL_281955";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    #filename = "AnalysisResults_HL_281955.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400";
    #suffix = "_w_TPCshared03_HL_281955";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    #filename = "AnalysisResults_HL_281955.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_TPCshared03_minpt400";
    #suffix = "";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);


    #filename = "AnalysisResults_HL_284904.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_TPCshared03";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    #filename = "AnalysisResults_HL_285248.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_TPCshared03";
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    filename = "AnalysisResults_HL_286598.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_FT0Coccupancy0_10000";
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, -0.5, -0.1, suffix);
#    draw_dphi(filename, taskname, -0.5, -0.1, suffix);
#    draw_dphi(filename, taskname, -0.5, -0.1, suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#    draw_dptdphi(filename, taskname, "uls", suffix);
#    draw_dptdphi(filename, taskname, "lspp", suffix);
#    draw_dptdphi(filename, taskname, "lsmm", suffix);

    #taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400";
    #suffix = "_FT0Coccupancy10000_20000";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

#    taskname = "dielectron_occupancy20000_30000_TPChadrejorTOFreq_minpt400";
#    suffix = "_FT0Coccupancy20000_30000";
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#
#    taskname = "dielectron_occupancy30000_40000_TPChadrejorTOFreq_minpt400";
#    suffix = "_FT0Coccupancy30000_40000";
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#
#
    taskname = "dielectron_occupancy40000_50000_TPChadrejorTOFreq_minpt400";
    suffix = "_FT0Coccupancy40000_50000";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


#    filename = "AnalysisResults_HL_286598.root";
#    tasknames = [
#        "single-electron-qc_occupancy0_10000_TPChadrejorTOFreq_minpt400",
#        "single-electron-qc_occupancy10000_20000_TPChadrejorTOFreq_minpt400",
#        "single-electron-qc_occupancy20000_30000_TPChadrejorTOFreq_minpt400",
#        "single-electron-qc_occupancy30000_40000_TPChadrejorTOFreq_minpt400",
#        "single-electron-qc_occupancy40000_50000_TPChadrejorTOFreq_minpt400",
#        "single-electron-qc_occupancy50000_99999_TPChadrejorTOFreq_minpt400",
#    ];
#    suffix = "_PbPb_2023_apass4";
#    draw_ft0c_occupancy(filename, tasknames, suffix);



    #filename = "AnalysisResults_HL_287694.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_FT0Coccupancy0_10000_HL_287694";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dptdphi(filename, taskname, "uls", suffix);
    #draw_dptdphi(filename, taskname, "lspp", suffix);
    #draw_dptdphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);

    #taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400";
    #suffix = "_FT0Coccupancy10000_20000_HL_287694";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dptdphi(filename, taskname, "uls", suffix);
    #draw_dptdphi(filename, taskname, "lspp", suffix);
    #draw_dptdphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);

    filename = "AnalysisResults_HL_288465.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_FT0Coccupancy0_10000_TPCshared07_HL_288465";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dptdphi(filename, taskname, "uls", suffix);
    #draw_dptdphi(filename, taskname, "lspp", suffix);
    #draw_dptdphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);


    filename = "AnalysisResults_HL_289080.root";
    taskname = "dielectron_TPChadrejorTOFreq_minpt200_dca3d";
    suffix = "_LHC23_Thin_pp_13.6TeV_TPChadrejorTOFreq_minpt200_HL_289080";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dptdphi(filename, taskname, "uls", suffix);
    #draw_dptdphi(filename, taskname, "lspp", suffix);
    #draw_dptdphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);


    #filename = "AnalysisResults_HL_289636.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_reldiffPin02_minpt400";
    #suffix = "_FT0Coccupancy0_10000_TPCshared07_reldiffPin02_HL_289636";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    ##draw_dptdphi(filename, taskname, "uls", suffix);
    ##draw_dptdphi(filename, taskname, "lspp", suffix);
    ##draw_dptdphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.05, 3.2, suffix);


    #filename = "AnalysisResults_HL_289971.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_AA";
    #suffix = "_FT0Coccupancy0_10000_TPCshared07_HL_289971_AA";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    ###draw_dphi(filename, taskname, -0.5, -0.1, suffix);

    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_CC";
    #suffix = "_FT0Coccupancy0_10000_TPCshared07_HL_289971_CC";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    ###draw_dphi(filename, taskname, -0.5, -0.1, suffix);

    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_AC";
    #suffix = "_FT0Coccupancy0_10000_TPCshared07_HL_289971_AC";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    ###draw_dphi(filename, taskname, -0.5, -0.1, suffix);

    filename = "AnalysisResults_HL_289636.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_FT0Coccupancy0_10000_TPCshared07_HL_289636";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.05, 3.2, suffix);

    taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400";
    suffix = "_FT0Coccupancy10000_20000_TPCshared07_HL_289636";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.05, 3.2, suffix);

    taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400";
    suffix = "_FT0Coccupancy10000_20000_TPCshared07_HL_289636";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.05, 3.2, suffix);

    #taskname = "dielectron_occupancy50000_999999_TPChadrejorTOFreq_minpt400";
    #suffix = "_FT0Coccupancy50000_999999_TPCshared07_HL_289636";
    ##draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    ##draw_meephiv(filename, taskname, "lsmm", suffix);
    ##draw_mee(filename, taskname, 3.05, 3.2, suffix);

    filename = "AnalysisResults_HL_284904.root";
    taskname = "dielectron_3050_TPChadrejorTOFreq_TPCshared03_minpt400";
    suffix = "_TPChadrejorTOFreq_TPCshared03_HL_284904";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.05, 3.2, suffix);

    filename = "AnalysisResults_HL_290986.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_0090_ftococcupancy0_10000_TPChadrejorTOFreq_HL_290986";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dpt(filename, taskname, 0.1, 0.5, suffix);

    filename = "AnalysisResults_HL_290986.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_wophiv";
    suffix = "_0090_ftococcupancy0_10000_TPChadrejorTOFreq_wophiv_HL_290986";
    #draw_dpt(filename, taskname, 0.1, 0.5, suffix);

#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    ##draw_mee_100(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#
#    taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400_wophiv";
#    suffix = "_0090_ftococcupancy10000_20000_TPChadrejorTOFreq_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#
#    taskname = "dielectron_occupancy20000_30000_TPChadrejorTOFreq_minpt400_wophiv";
#    suffix = "_0090_ftococcupancy20000_30000_TPChadrejorTOFreq_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#
#    taskname = "dielectron_occupancy30000_40000_TPChadrejorTOFreq_minpt400_wophiv";
#    suffix = "_0090_ftococcupancy30000_40000_TPChadrejorTOFreq_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#
#    taskname = "dielectron_occupancy40000_50000_TPChadrejorTOFreq_minpt400_wophiv";
#    suffix = "_0090_ftococcupancy40000_50000_TPChadrejorTOFreq_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);
#
#    taskname = "dielectron_occupancy50000_999999_TPChadrejorTOFreq_minpt400_wophiv";
#    suffix = "_0090_ftococcupancy50000_999999_TPChadrejorTOFreq_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);


#    filename = "AnalysisResults_HL_291277.root";
#    taskname = "dielectron_TPChadrejorTOFreq_minpt200_wophiv";
#    suffix = "_pp_13.6TeV_LHC22o_TPChadrejorTOFreq_minpt200_wophiv_HL_290986";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.08, 3.2, suffix);
#    draw_mee_100(filename, taskname, 3.08, 3.2, suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    #filename = "AnalysisResults_HL_291572.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_TPChadrejorTOFreq_minpt400_wophiv_HL_291572";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);

    filename = "AnalysisResults_HL_294001_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_wophiv_HL_294001";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    ###draw_mee_100(filename, taskname, 3.11, 3.2, suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_minpt400_wophiv_HL_294001";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_TPChadrejorTOFreq_minpt400_wophiv_HL_294001";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    filename = "AnalysisResults_HL_294383_tmp.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_294383";
    ##draw_meephiv(filename, taskname, "uls", suffix);
    ##draw_meephiv(filename, taskname, "lspp", suffix);
    ##draw_meephiv(filename, taskname, "lsmm", suffix);
    ##draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400_wo_detadphi";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_294383";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);


    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_TPChadrejorTOFreq_minpt400_HL_294383";
    ##draw_meephiv(filename, taskname, "uls", suffix);
    ##draw_meephiv(filename, taskname, "lspp", suffix);
    ##draw_meephiv(filename, taskname, "lsmm", suffix);
    ##draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);


    filename = "AnalysisResults_HL_294383_tmp_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_294383";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dpt(filename, taskname, 0.1, 0.5, suffix);

    filename = "AnalysisResults_HL_294420_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_294420";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_dptdphi(filename, taskname, "uls", suffix);
    #draw_dptdphi(filename, taskname, "lspp", suffix);
    #draw_dptdphi(filename, taskname, "lsmm", suffix);


    filename = "AnalysisResults_HL_294420_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_HL_294420";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    #filename = "AnalysisResults_HL_294526_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_HL_294526";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    #filename = "AnalysisResults_HL_295259.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta005";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_deta005_HL_295259";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_deta010_HL_295259";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

#    filename = "AnalysisResults_HL_295155.root";
#    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta005";
#    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_deta005_HL_295155";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.11, 3.2, suffix);
#    draw_detadphi_geom(filename, taskname, "uls", suffix);
#    draw_detadphi_geom(filename, taskname, "lspp", suffix);
#    draw_detadphi_geom(filename, taskname, "lsmm", suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);

#    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta010";
#    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_minpt400_deta010_HL_295155";
#    draw_meephiv(filename, taskname, "uls", suffix);
#    draw_meephiv(filename, taskname, "lspp", suffix);
#    draw_meephiv(filename, taskname, "lsmm", suffix);
#    draw_mee(filename, taskname, 3.11, 3.2, suffix);
#    draw_detadphi_geom(filename, taskname, "uls", suffix);
#    draw_detadphi_geom(filename, taskname, "lspp", suffix);
#    draw_detadphi_geom(filename, taskname, "lsmm", suffix);
#    draw_detadphi(filename, taskname, "uls", suffix);
#    draw_detadphi(filename, taskname, "lspp", suffix);
#    draw_detadphi(filename, taskname, "lsmm", suffix);
#    draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    #filename = "AnalysisResults_HL_296423.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_occupancy0_10000_TPChadrejorTOFreq_minpt400_HL_296423";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    #taskname = "dielectron_occupancy50000_999999_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_occupancy50000_999999_TPChadrejorTOFreq_minpt400_HL_296423";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_detadphi_geom(filename, taskname, "uls", suffix);
    #draw_detadphi_geom(filename, taskname, "lspp", suffix);
    #draw_detadphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    filename = "AnalysisResults_HL_298365.root"; #without dz-rdphi cut
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_occupancy0_10000_TPChadrejorTOFreq_minpt400_HL_298365";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    filename = "AnalysisResults_HL_299066.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299066";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299066";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);


    filename = "AnalysisResults_HL_299065.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_occupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299065";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    filename = "AnalysisResults_HL_299220_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299220";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.1, +0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, +0.1, +0.5, suffix);

    #draw_dphi_yield(filename, taskname, "uls", 0.1, 0.5, suffix);
    #draw_dphi_yield(filename, taskname, "lspp", 0.1, 0.5, suffix);
    #draw_dphi_yield(filename, taskname, "lsmm", 0.1, 0.5, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299220";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);

    filename = "AnalysisResults_HL_299686_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_id20308"; #itsibany
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_299686";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.1, +0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, +0.1, +0.5, suffix);

    #filename = "AnalysisResults_HL_302828.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_302828";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_meephiv_prop(filename, taskname, "uls", suffix);
    #draw_meephiv_prop(filename, taskname, "lspp", suffix);
    #draw_meephiv_prop(filename, taskname, "lsmm", suffix);
    #draw_mee(filename, taskname, 3.11, 3.2, suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.5, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.5, -0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, -0.1, +0.1, suffix);
    #draw_rfactor_vs_dphi(filename, taskname, +0.1, +0.5, suffix);

    filename = "AnalysisResults_HL_302665.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_x50"; #without prefilter, x = 50, 83, 120, 200
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_x50_HL_302665";
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_x83"; #without prefilter, x = 50, 83, 120, 200
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_x83_HL_302665";
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_x120"; #without prefilter, x = 50, 83, 120, 200
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_x120_HL_302665";
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);


    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_x200"; #without prefilter, x = 50, 83, 120, 200
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_x200_HL_302665";
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);


    filename = "AnalysisResults_HL_305480.root";
    taskname = "dielectron_TPChadrejorTOFreq_minpt200"; #without prefilter
    suffix = "_pp_13.6TeV_LHC23_pass4_Thin_0_999_TPChadrejorTOFreq_minpt200_HL_305480";
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);

    filename = "AnalysisResults_HL_306301.root";
    taskname = "dielectron_TPChadrejorTOFreq_minpt200"; #without prefilter
    #taskname = "dielectron_TPChadrejorTOFreq_minpt200_wpf_mee"; #without prefilter
    suffix = "_pp_13.6TeV_LHC23_pass4_Thin_0_999_TPChadrejorTOFreq_minpt200_HL_306301";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    
    filename = "AnalysisResults_HL_306764.root";
    #taskname = "dielectron_TPChadrejorTOFreq_minpt200"; #without prefilter
    taskname = "dielectron_TPChadrejorTOFreq_minpt200_wpf_mee"; #without prefilter
    suffix = "_pp_13.6TeV_LHC22o_minbias_pass7_0_999_TPChadrejorTOFreq_minpt200_HL_306764";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    filename = "AnalysisResults_HL_306762.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS wide
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_306762";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    filename = "AnalysisResults_HL_307487.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_307487";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_307487";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);


    filename = "AnalysisResults_HL_307895.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_307895";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);


    filename = "AnalysisResults_HL_308878.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_308878";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    filename = "AnalysisResults_HL_309527_merged_centrality.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_309527";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);
    #draw_dphi(filename, taskname, 0.2, 999, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_309527";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);

    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_309527";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);
    #draw_dphi(filename, taskname, 0.2, 999, suffix);


    filename = "AnalysisResults_HL_310124.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_310124";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);

    filename = "AnalysisResults_HL_310124.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wopf_detadphi"; #without prefilter LS and ULS, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wopf_detadphi_HL_310124";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);

    filename = "AnalysisResults_HL_311554_merged_centrality.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_311554";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_311554";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    ##draw_dphi(filename, taskname, 0.1, 0.8, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_311554";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    ##draw_dphi(filename, taskname, 0.1, 0.8, suffix);

    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_311554";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dzrdphi_geom(filename, taskname, "uls", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lspp", suffix);
    #draw_dzrdphi_geom(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    filename = "AnalysisResults_HL_312076.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_312076";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);
    #draw_dphi(filename, taskname, 0.2, 999, suffix);
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_reldiffPin015"; #without prefilter LS and ULS, without deta dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_reldiffPin015_HL_312076";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);
    #draw_dphi(filename, taskname, 0.2, 999, suffix);

    filename = "AnalysisResults_HL_313263_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS, with deta > 0.04, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_313263";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS, with deta > 0.04, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_313263";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    filename = "AnalysisResults_HL_314544_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS, without detadphi, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_314544";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS without detadphi cut,with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_314544";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee_opangle(filename, taskname, "uls", suffix);
    #draw_mee_opangle(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400"; #without prefilter LS and ULS without detadphi cut,with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_314544";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_mee_opangle(filename, taskname, "uls", suffix);
    #draw_mee_opangle(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);


    filename = "AnalysisResults_HL_315520_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter, without deta-dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_315520";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    filename = "AnalysisResults_HL_316686_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter, with deta-dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_316686";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter, with deta-dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_316686";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400"; #with prefilter, with deta-dphi cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_316686";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);


    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wpf_onlyLS"; #with prefilter, with deta-dphi cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wpf_onlyLS_HL_316686";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    #filename = "AnalysisResults_HL_319300_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter 0.03, with deta-dphi cut 0.03, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_319300";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 999, suffix);

    #filename = "AnalysisResults_HL_319559_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter, without deta-dphi cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_319559";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel"; #without prefilter, without deta-dphi cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel_HL_319559";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi"; #without prefilter, without deta-dphi cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_319559";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel"; #without prefilter, without deta-dphi cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel_HL_319559";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    filename = "AnalysisResults_HL_319736_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.04, with deta > 0.04 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_319736";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.04, with detai > 0.04 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_319736";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.04, with deta > 0.04 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_319736";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.04, without deta > 0.04 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_319736";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    ###taskname = "dielectron_1090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.04, with deta > 0.04 cut, with ndiff_bc = 594
    ###suffix = "_PbPb_5.36TeV_LHC23_pass4_1090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_319736";
    ###draw_meephiv(filename, taskname, "uls", suffix);
    ###draw_meephiv(filename, taskname, "lspp", suffix);
    ###draw_meephiv(filename, taskname, "lsmm", suffix);
    ###draw_detadphi(filename, taskname, "uls", suffix);
    ###draw_detadphi(filename, taskname, "lspp", suffix);
    ###draw_detadphi(filename, taskname, "lsmm", suffix);
    ###draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    ###draw_dphi(filename, taskname, 0.8, 999, suffix);

    ###taskname = "dielectron_1090_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.04, without deta > 0.04 cut, with ndiff_bc = 594
    ###suffix = "_PbPb_5.36TeV_LHC23_pass4_1090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_319736";
    ###draw_meephiv(filename, taskname, "uls", suffix);
    ###draw_meephiv(filename, taskname, "lspp", suffix);
    ###draw_meephiv(filename, taskname, "lsmm", suffix);
    ###draw_detadphi(filename, taskname, "uls", suffix);
    ###draw_detadphi(filename, taskname, "lspp", suffix);
    ###draw_detadphi(filename, taskname, "lsmm", suffix);
    ###draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    ###draw_dphi(filename, taskname, 0.8, 999, suffix);

    taskname = "dielectron_5090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.04, with deta > 0.04 cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_5090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_319736";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #filename = "AnalysisResults_HL_320049_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.03, with deta > 0.03 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_320049";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.03, with detai > 0.03 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_320049";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.03, with deta > 0.03 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_320049";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.03, without deta > 0.03 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_320049";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #filename = "AnalysisResults_HL_320161_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.02, with deta > 0.02 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_320161";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.02, with detai > 0.02 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_320161";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.02, with deta > 0.02 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_320161";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.02, without deta > 0.02 cut, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_320161";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    filename = "AnalysisResults_HL_321612.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400"; #without prefilter, without deta cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_321612";
    #draw_meephiv(filename, taskname, "uls", suffix);
    #draw_meephiv(filename, taskname, "lspp", suffix);
    #draw_meephiv(filename, taskname, "lsmm", suffix);
    #draw_detadphi(filename, taskname, "uls", suffix);
    #draw_detadphi(filename, taskname, "lspp", suffix);
    #draw_detadphi(filename, taskname, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);


    filename = "AnalysisResults_HL_328076.root";
    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400"; #without prefilter, without deta cut, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC24as_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_328076";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 0.8, 999, suffix);

    filename = "AnalysisResults_HL_328095_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi"; #with prefilter deta > 0.04, withdeta cut > 0.04, with ndiff_bc = 594
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_328095";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel"; #with prefilter deta > 0.04, withdeta cut > 0.04, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel_HL_328095";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #filename = "AnalysisResults_HL_328470_merged_centrality.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #with prefilter deta > 0.04, withdeta cut > 0.04, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_328470";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel"; #with prefilter deta > 0.04, withdeta cut > 0.04, with ndiff_bc = 594
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_328470";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    filename = "AnalysisResults_HL_329829_merged_centrality.root";
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi";
    suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_HL_329829";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel_HL_329829";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_HL_329829";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_minpt400_TightEvSel_HL_329829";
    #draw_meephiv(filename,  taskname, 0, 90, "uls" , suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lspp", suffix);
    #draw_meephiv(filename,  taskname, 0, 90, "lsmm", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "uls", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lspp", suffix);
    #draw_detadphi(filename, taskname, 0, 90, "lsmm", suffix);
    #draw_dphi(filename, taskname, 0.1, 0.8, suffix);
    #draw_dphi(filename, taskname, 1.0, 1.6, suffix);


#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
