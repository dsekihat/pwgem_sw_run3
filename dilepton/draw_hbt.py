import datetime
import math
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TColor, TF1
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
gStyle.SetPalette(55);
#gStyle.SetOptStat(0);
#gStyle.SetOptTitle(0);

#__________________________________________________________
def draw_hbt_3d(filename, taskname, qmin, qmax, ktid, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();

    h3cf = rootdire.Get("h3cf_kt{0:d}".format(ktid));
    bin_x1 = h3cf.GetXaxis().FindBin(qmin + 1e-3);
    bin_x2 = h3cf.GetXaxis().FindBin(qmax - 1e-3);
    bin_y1 = h3cf.GetYaxis().FindBin(qmin + 1e-3);
    bin_y2 = h3cf.GetYaxis().FindBin(qmax - 1e-3);
    bin_z1 = h3cf.GetZaxis().FindBin(qmin + 1e-3);
    bin_z2 = h3cf.GetZaxis().FindBin(qmax - 1e-3);

    dbinx = float(bin_x2 - bin_x1 + 1);
    dbiny = float(bin_y2 - bin_y1 + 1);
    dbinz = float(bin_z2 - bin_z1 + 1);
    #print(dbinx, dbiny, dbinz);

    h1out  = h3cf.ProjectionX("h1out" , bin_y1, bin_y2, bin_z1, bin_z2, "");
    h1side = h3cf.ProjectionY("h1side", bin_x1, bin_x2, bin_z1, bin_z2, "");
    h1long = h3cf.ProjectionZ("h1long", bin_x1, bin_x2, bin_y1, bin_y2, "");

    h1out .Scale(1.0/dbiny/dbinz);
    h1side.Scale(1.0/dbinx/dbinz);
    h1long.Scale(1.0/dbinx/dbiny);

    h1out .SetDirectory(0);
    h1side.SetDirectory(0);
    h1long.SetDirectory(0);
    ROOT.SetOwnership(h1out , False);
    ROOT.SetOwnership(h1side, False);
    ROOT.SetOwnership(h1long, False);

    make_common_style(h1out , 20, 1.2, kBlack, 2, 0);
    make_common_style(h1side, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1long, 20, 1.2, kBlack, 2, 0);

    ymin = 0.8;
    ymax = 1.3;

    c1 = TCanvas("c1", "c1", 0, 0, 1500, 500);
    c1.Divide(3,1,0,0);
    p1 = c1.cd(1);
    p1.SetPad(0,0,0.36,1);
    p1.SetMargin(0.12, 0.0, 0.12, 0.03);
    p1.SetTicks(1, 1);
    frame1 = p1.DrawFrame(-0.29, ymin, +0.29, ymax);
    frame1.GetXaxis().SetTitle("#it{q}_{out} (GeV#it{c})");
    frame1.GetYaxis().SetTitle("C(#it{q})");
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.2);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    ROOT.SetOwnership(frame1, False);
    ROOT.SetOwnership(p1, False);
    h1out.Draw("E0same");

    p2 = c1.cd(2);
    p2.SetPad(0.36,0,0.68,1);
    p2.SetMargin(0.0,0.0,0.12,0.03);
    p2.SetTicks(1, 1);
    frame2 = p2.DrawFrame(-0.29, ymin, +0.29, ymax);
    frame2.GetXaxis().SetTitle("#it{q}_{side} (GeV#it{c})");
    frame2.GetYaxis().SetTitle("C(#it{q})");
    frame2.GetXaxis().SetTitleSize(0.05);
    frame2.GetYaxis().SetTitleSize(0.05);
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(1.2);
    frame2.GetXaxis().SetLabelSize(0.05);
    frame2.GetYaxis().SetLabelSize(0.05);
    ROOT.SetOwnership(frame2, False);
    ROOT.SetOwnership(p2, False);
    h1side.Draw("E0same");


    p3 = c1.cd(3);
    p3.SetPad(0.68,0,1,1);
    p3.SetMargin(0.0,0.02,0.12,0.03);
    p3.SetTicks(1, 1);
    frame3 = p3.DrawFrame(-0.29, ymin, +0.29, ymax);
    frame3.GetXaxis().SetTitle("#it{q}_{long} (GeV#it{c})");
    frame3.GetYaxis().SetTitle("C(#it{q})");
    frame3.GetXaxis().SetTitleSize(0.05);
    frame3.GetYaxis().SetTitleSize(0.05);
    frame3.GetXaxis().SetTitleOffset(1.0);
    frame3.GetYaxis().SetTitleOffset(1.2);
    frame3.GetXaxis().SetLabelSize(0.05);
    frame3.GetYaxis().SetLabelSize(0.05);
    ROOT.SetOwnership(frame3, False);
    ROOT.SetOwnership(p3, False);
    h1long.Draw("E0same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1, False);
    c1.SaveAs("{0}_gghbt_cf_3d_kT{1:d}_q{2:3.2f}_{3:3.2f}{4}.eps".format(date, ktid, qmin, qmax, suffix));
    c1.SaveAs("{0}_gghbt_cf_3d_kT{1:d}_q{2:3.2f}_{3:3.2f}{4}.pdf".format(date, ktid, qmin, qmax, suffix));
    c1.SaveAs("{0}_gghbt_cf_3d_kT{1:d}_q{2:3.2f}_{3:3.2f}{4}.png".format(date, ktid, qmin, qmax, suffix));
    rootfile.Close();

#__________________________________________________________
def draw_detadphi_photon(filename, taskname, suffix):
    rootfile = TFile(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    rootdire.Get("Pair/same").ls();

    h2same = rootdire.Get("Pair/same/hDeltaEtaDeltaPhi_Photon").Clone("h2same");
    h2mix  = rootdire.Get("Pair/mix/hDeltaEtaDeltaPhi_Photon").Clone("h2mix");

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
    h2ratio.SetMinimum(0.95);
    h2ratio.SetMaximum(1.05);
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

    #min_dphi = 0.2;
    #min_deta = 0.06;
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
    c1.SaveAs("{0}_detadphi_photon{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_detadphi_photon{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_detadphi_photon{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_detadphi(filename, taskname, suffix):
    rootfile = TFile(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    rootdire.Get("Pair/same").ls();

    h2same = rootdire.Get("Pair/same/hDeltaEtaDeltaPhi").Clone("h2same");
    h2mix  = rootdire.Get("Pair/mix/hDeltaEtaDeltaPhi").Clone("h2mix");

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

    frame1 = c1.DrawFrame(-math.pi, -1, +math.pi, +1);
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

    #min_dphi = 0.2;
    #min_deta = 0.06;
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
    c1.SaveAs("{0}_detadphi{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_detadphi{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_detadphi{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_drdz(filename, taskname, suffix):
    rootfile = TFile(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    rootdire.Get("Pair/same").ls();

    h2same = rootdire.Get("Pair/same/hDeltaRDeltaZ").Clone("h2same");
    h2mix  = rootdire.Get("Pair/mix/hDeltaRDeltaZ").Clone("h2mix");
    #h2same = rootdire_pair_same.Get("hsDeltaP").Projection(1,2).Clone("h2same");
    #h2mix  = rootdire_pair_mix .Get("hsDeltaP").Projection(1,2).Clone("h2mix");

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
    h2ratio.SetMinimum(0.95);
    h2ratio.SetMaximum(1.05);
    h2ratio.SetZTitle("same/mix");
    h2ratio.GetZaxis().SetTitleOffset(1.5);

    gStyle.SetPalette(55);
    #gStyle.SetPalette(87);
    #gStyle.SetPalette(104);
    ##gStyle.SetPalette(70);
    ##TColor.InvertPalette();
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.16, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetGridx(1);
    c1.SetGridy(1);
#    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0, 0, 50, 50);
    frame1.GetXaxis().SetTitle("|#Deltar| (cm)");
    frame1.GetYaxis().SetTitle("|#Deltaz| (cm)");
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

    min_dr = 15;
    min_dz = 20;
    f1up = TF1("f1up", "[1] * sqrt(1 - x*x/[0]/[0])", 0, +min_dr);
    f1up.SetNpx(min_dr * 1000);
    f1up.SetParameters(min_dr, min_dz);
    f1up.Draw("same");
    ROOT.SetOwnership(f1up, False);

    #f1down = TF1("f1down", "-[1] * sqrt(1 - x*x/[0]/[0])", -min_dphi, +min_dphi);
    #f1down.SetNpx(1000);
    #f1down.SetParameters(min_dphi, min_deta);
    #f1down.Draw("same");
    #ROOT.SetOwnership(f1down, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_drdz{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_drdz{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_drdz{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_test_hbt.root", "photon-hbt-pcmee_0010_mee0_140", -0.03, +0.03, 9, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_test_hbt.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 8, "");
    #draw_hbt_3d("photon_hbt_3d_pp_13.6TeV_LHC22o_pass7_HL_244379.root", "photon-hbt-pcmpcm", -0.03, +0.03, 0, "_pp_13.6TeV_LHC22o_pass7_pcmpcm");
    #draw_hbt_3d("photon_hbt_3d_pp_13.6TeV_LHC22o_pass7_HL_244379.root", "photon-hbt-pcmee_mee0_140", -0.03, +0.03, 0, "_pp_13.6TeV_LHC22o_pass7_pcmee");
    #draw_hbt_3d("photon_hbt_3d_pp_13.6TeV_LHC22o_pass7_HL_244379.root", "photon-hbt-pcmee_mee500_1100", -0.03, +0.03, 6, "_pp_13.6TeV_LHC22o_pass7_pcmee");

    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 1, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 2, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 3, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 4, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 5, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 6, "");

    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_3050", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_7090", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmee_0010", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmee_3050", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmee_7090", -0.03, +0.03, 0, "");

    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 0, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 1, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 2, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 3, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 4, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 5, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 6, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 7, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 8, "");
    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_0010", -0.03, +0.03, 9, "");

    #draw_hbt_3d("photon_hbt_3d_PbPb_5.36TeV_LHC23zzh_pass4_HL_263031.root", "photon-hbt-pcmpcm_3050", -0.03, +0.03, 4, "");

    #filename = "AnalysisResults_HL_309434.root";
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    ##suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_309434";
    ##draw_drdz(filename, taskname, suffix);
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_309434";
    #draw_detadphi(filename, taskname, suffix);
    ##taskname = "photon-hbt-pcmee_1d_0010";
    ##suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmee_HL_309434";
    ##draw_detadphi(filename, taskname, suffix);


    #filename = "AnalysisResults_HL_310198.root";
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_310198";
    ##draw_drdz(filename, taskname, suffix);
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_310198";
    #draw_detadphi(filename, taskname, suffix);
    ###taskname = "photon-hbt-pcmee_1d_0010";
    ###suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmee_HL_310198";
    ###draw_detadphi(filename, taskname, suffix);


    filename = "AnalysisResults_HL_311370.root";
    taskname = "photon-hbt-pcmpcm_1d_0010";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_311370";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);
    taskname = "photon-hbt-pcmpcm_1d_3050";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_311370";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);
    taskname = "photon-hbt-pcmpcm_1d_7090";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_311370";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);

    filename = "AnalysisResults_HL_315730.root";
    taskname = "photon-hbt-pcmpcm_1d_0010";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_315730";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);
    taskname = "photon-hbt-pcmpcm_1d_3050";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_315730";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);
    taskname = "photon-hbt-pcmpcm_1d_7090";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_315730";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);

    filename = "AnalysisResults_HL_316404.root";
    taskname = "photon-hbt-pcmpcm_1d_0010";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_316404";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);

    filename = "AnalysisResults_HL_316859.root";
    taskname = "photon-hbt-pcmpcm_1d_0010";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_316859";
    #draw_drdz(filename, taskname, suffix);
    #draw_detadphi(filename, taskname, suffix);


    filename = "AnalysisResults_HL_317565.root";
    taskname = "photon-hbt-pcmpcm_1d_0010";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_317565";
    #draw_detadphi_photon(filename, taskname, suffix);
    taskname = "photon-hbt-pcmpcm_1d_3050";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_317565";
    #draw_detadphi_photon(filename, taskname, suffix);


    #filename = "AnalysisResults_HL_317846.root";
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_317846";
    #draw_detadphi_photon(filename, taskname, suffix);

    #filename = "AnalysisResults_HL_318188.root";
    #taskname = "photon-hbt-pcmpcm_1d_0010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT100MeV_pcmpcm_HL_318188";
    #draw_detadphi_photon(filename, taskname, suffix);

    filename = "AnalysisResults_HL_336553.root";
    taskname = "photon-hbt-pcmpcm_1d_MB";
    suffix = "_pp_5.36TeV_LHC24ap_pass1_pT100MeV_pcmpcm_HL_336553";
    draw_detadphi_photon(filename, taskname, suffix);
    #draw_drdz(filename, taskname, suffix);

    taskname = "photon-hbt-pcmee_1d_MB";
    suffix = "_pp_5.36TeV_LHC24ap_pass1_pcmee_HL_336553";
    #draw_detadphi_photon(filename, taskname, suffix);

#__________________________________________________________
#__________________________________________________________
