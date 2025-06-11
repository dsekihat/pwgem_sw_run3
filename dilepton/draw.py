import datetime
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
from histo_manager import get_ratio, rebin_histogram, get_cumulative_histogram
from signal_extractor import get_significance, rebin_flow_thn
from dilepton_utils import extract_flow
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
def draw_mmumu_Rfactor(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_R = rootdire.Get("hs_R");

    bin0 = hs_R.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_R.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_R.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_R.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_R.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_R.GetAxis(2).SetRange(bin0, bin1);

    h1m_r = hs_R.Projection(0);

    h1m_r.SetDirectory(0);
    ROOT.SetOwnership(h1m_r, False);
    make_common_style(h1m_r, 20, 1.2, kBlack, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.18, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(0., 0.9, 12, 1.06);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line1_40 = TLine(0, 1, 12.0, 1);
    line1_40.SetLineColor(kBlack);
    line1_40.SetLineStyle(2);
    line1_40.SetLineWidth(2);
    line1_40.Draw("same");
    ROOT.SetOwnership(line1_40,False);
    h1m_r.Draw("E0same");

    #txt = TPaveText(0.45, 0.65, 0.65, 0.95, "NDC");
    txt = TPaveText(0.25, 0.2, 0.45, 0.5, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    minptmu = 0.2;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    if "pp" in suffix:
        #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("pp at #sqrt{#it{s}} = 5.36 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #p2 = TPad("p2", "p2", 0.35, 0.12, 0.92, 0.62);
    #p2.SetMargin(0.25, 0.03, 0.12, 0.02);
    #p2.SetTicks(1,1);
    #c1.cd();
    #p2.Draw("");
    #p2.cd();
    #frame2 = p2.DrawFrame(0.,0.9, 1.1, 1.02);
    #frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    #frame2.GetYaxis().SetTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
    #frame2.GetXaxis().SetTitleOffset(1.1);
    #frame2.GetYaxis().SetTitleOffset(2.3);
    #frame2.GetXaxis().SetTitleSize(0.05);
    #frame2.GetYaxis().SetTitleSize(0.05);
    #frame2.GetXaxis().SetLabelSize(0.05);
    #frame2.GetYaxis().SetLabelSize(0.05);
    #frame2.GetYaxis().SetMaxDigits(3);
    #frame2.GetXaxis().SetNdivisions(510);
    #frame2.GetYaxis().SetNdivisions(510);
    #line1_11 = TLine(0, 1, 1.1, 1);
    #line1_11.SetLineColor(kBlack);
    #line1_11.SetLineStyle(2);
    #line1_11.SetLineWidth(2);
    #line1_11.Draw("same");
    #ROOT.SetOwnership(line1_11, False);
    #ROOT.SetOwnership(p2,False);
    #ROOT.SetOwnership(frame2,False);


    h1m_r.Draw("E0same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mmumu_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();


#__________________________________________________________
def draw_mmumu_sbratio(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);

    h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 1e+3;
    y_range_min = 1e-1;

    #frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1 = c1.DrawFrame(0., y_range_min, 12, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.45, 0.65, 0.65, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    minptmu = 0.3;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt300" in suffix:
        minptmu = 0.3;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    if "pp" in suffix:
        #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("pp at #sqrt{#it{s}} = 5.36 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mmumu_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mmumu_significance(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);

    h1z = rootdire.Get("hZvtx");
    nev = h1z.GetEntries();
    h1norm = rootdire.Get("h1norm");
    if h1norm != None:
        nbin_tmp = h1norm.GetNbinsX();
        nev = h1norm.GetBinContent(nbin_tmp);
        print("h1norm is used");
    h1m_bkg.Scale(nev);
    h1m_sig.Scale(nev);

    h1significance = get_significance(h1m_sig, h1m_bkg);
    make_common_style(h1significance, 20, 1.2, kRed+1, 2, 0);
    h1significance.SetDirectory(0);
    ROOT.SetOwnership(h1significance, False);

    print(h1m_uls.GetBinContent(10));
    print(h1m_sig.GetBinContent(10));
    print(h1significance.GetBinContent(10));
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 1e+3;
    y_range_min = 1e-0;

    frame1 = c1.DrawFrame(0., y_range_min, 12, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/#sqrt{S+2B}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h1significance.Draw("E0same");

    txt = TPaveText(0.45, 0.65, 0.65, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    minptmu = 0.3;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt300" in suffix:
        minptmu = 0.3;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    if "pp" in suffix:
        #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("pp at #sqrt{#it{s}} = 5.36 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mmumu_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mmumu_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mmumu_uls_ls_1100MeV(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.03, 0.1, 0.05);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 1.1;
    y_range_min = max(h1m_sig.GetMinimum() * 0.5, 2e-4);
    if "standalone_muon" in suffix:
        y_range_max = h1m_uls.GetMaximum() * 2.0;
    y_range_max = h1m_uls.GetMaximum() * 10;

    frame1 = c1.DrawFrame(0.2, y_range_min, 1.1, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.45, 0.85, 0.6);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = ULS #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.45, 0.6, 0.65, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    minptmu = 0.2;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_1100MeV_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_1100MeV_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_1100MeV_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mmumu_uls_ls(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.165, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogy(1);

    h1m_sig.GetXaxis().SetRangeUser(0.21, 4);
    y_range_max = h1m_uls.GetMaximum() * 2;
    y_range_min = max(h1m_sig.GetMinimum() * 0.5, 2e-12);
    h1m_sig.GetXaxis().SetRangeUser(0.2, 12);

    frame1 = c1.DrawFrame(0.0, y_range_min, 4, y_range_max);
    #frame1 = c1.DrawFrame(0.0, y_range_min, 12, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.00);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    #leg = TLegend(0.65, 0.5, 0.85, 0.65);
    leg = TLegend(0.55, 0.15, 0.75, 0.30);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = ULS #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #txt = TPaveText(0.45, 0.65, 0.65, 0.95, "NDC");
    txt = TPaveText(0.2, 0.15, 0.4, 0.5, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    minptmu = 0.3;
    if "minpt200" in suffix:
        minptmu = 0.2;
    if "minpt300" in suffix:
        minptmu = 0.3;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    txt.AddText("ALICE Performance");
    if "pp" in suffix:
        #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("pp at #sqrt{#it{s}} = 5.36 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mmumu_uls_ls_pbpb_2panel(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    outfile_tmp = TFile("20241125_PbPb_5.36TeV_pass3_{0}.root".format(taskname), "RECREATE");
    outfile_tmp.WriteTObject(h1m_uls);
    outfile_tmp.WriteTObject(h1m_bkg);
    outfile_tmp.WriteTObject(h1m_sig);
    outfile_tmp.Close();

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.165, 0.025, 0.09, 0.01);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogy(1);

    h1m_sig.GetXaxis().SetRangeUser(0.22, 12);
    y_range_max = h1m_uls.GetMaximum() * 80;
    y_range_min = max(h1m_sig.GetMinimum() * 0.5, 2e-10);
    h1m_sig.GetXaxis().SetRangeUser(0.2, 12);

    frame1 = c1.DrawFrame(0.0, y_range_min, 12, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.00);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.DrawClone("E0same");
    ROOT.SetOwnership(h1m_sig, False);

    #leg = TLegend(0.65, 0.5, 0.85, 0.65);
    leg = TLegend(0.7, 0.38, 0.85, 0.5);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.AddEntry(h1m_uls, "#it{N}_{+#minus} = #it{S} + #it{B}", "P");
    leg.AddEntry(h1m_bkg, "#it{B} = 2#it{R} #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "#it{S} = #it{N}_{+#minus} #minus #it{B}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt1 = TPaveText(0.19, 0.8, 0.39, 0.95, "NDC");
    txt1.SetFillColor(kWhite);
    txt1.SetFillStyle(0);
    txt1.SetBorderSize(0);
    txt1.SetTextAlign(12);#middle,left
    txt1.SetTextFont(42);#helvetica
    txt1.SetTextSize(0.028);
    txt1.Draw();
    txt1.AddText("ALICE Performance");
    txt1.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt1.AddText("centrality {0:d}#minus{1:d}%".format(cen1, cen2));
    ROOT.SetOwnership(txt1, False);

    txt = TPaveText(0.19, 0.12, 0.39, 0.28, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.028);

    minptmu = 0.2;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MCH-MID";
    if "global_muon" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone_muon" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";
    txt.AddText(str_muon_type);

    #if "pp" in suffix:
    #    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #elif "PbPb" in suffix:
    #    txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #elif "pPb" in suffix:
    #    txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #else:
    #    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    #txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    h1m_sig.GetXaxis().SetRangeUser(0.5, 1.2);
    y_max2 = h1m_sig.GetMaximum() * 1.5;
    h1m_sig.GetXaxis().SetRangeUser(0, 12);

    p2 = TPad("p2", "p2", 0.5, 0.51, 0.95, 0.96);
    p2.SetMargin(0.22, 0.03, 0.12, 0.06);
    p2.SetTicks(1,1);
    c1.cd();
    p2.Draw("");
    p2.cd();
    frame2 = p2.DrawFrame(0.5, 0, 1.2, y_max2);
    frame2.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(1.8);
    frame2.GetXaxis().SetTitleSize(0.05);
    frame2.GetYaxis().SetTitleSize(0.05);
    frame2.GetXaxis().SetLabelSize(0.05);
    frame2.GetYaxis().SetLabelSize(0.05);
    frame2.GetYaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetNdivisions(505);
    frame2.GetXaxis().SetNdivisions(510);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    h1m_sig.GetXaxis().SetRangeUser(0.2, 1.2);
    h1m_sig.Draw("E0same");
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig, False);

    #f1 = TF1("f1", "crystalball(0) + crystalball(5) + [10]*exp(-x/[11]) + [12]", 0.2, 1.2);
    #f1.SetLineColor(kGreen+2);
    #f1.SetLineWidth(2);
    #f1.SetLineStyle(2);
    #f1.SetNpx(1000);
    #f1.SetParameters(1e-5, 0.782, 0.02, 0.1, 0.1, 1e-5, 1.019, 0.02, 0.1, 0.1);
    #f1.SetParameter(10, 1e-3);
    #f1.SetParameter(11, 10);
    #f1.SetParameter(12, 1e-3);
    #f1.SetParLimits(0, 0, 1e-4);
    #f1.SetParLimits(5, 0, 1e-4);
    #f1.SetParLimits(1, 0.7, 0.9);
    #f1.SetParLimits(6, 0.9, 1.1);
    #f1.SetParLimits(2, 1e-3, 0.1);
    #f1.SetParLimits(7, 1e-3, 0.1);
    #h1m_sig.Fit(f1, "SME", "", 0.5, 1.2);
    #ROOT.SetOwnership(f1, False);
    #print(f1.GetChisquare()/f1.GetNDF());

    #txt_omega = TPaveText(0.3, 0.78, 0.5, 0.9, "NDC");
    #txt_omega.SetFillColor(kWhite);
    #txt_omega.SetFillStyle(0);
    #txt_omega.SetBorderSize(0);
    #txt_omega.SetTextAlign(12);#middle,left
    #txt_omega.SetTextFont(42);#helvetica
    #txt_omega.SetTextSize(0.05);
    #txt_omega.AddText("#it{{m}}_{{#omega}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(1) * 1e+3, f1.GetParError(1) * 1e+3));
    #txt_omega.AddText("#sigma_{{#omega}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(2) * 1e+3, f1.GetParError(2) * 1e+3));
    #txt_omega.Draw();
    #ROOT.SetOwnership(txt_omega, False);

    #txt_phi = TPaveText(0.3, 0.64, 0.5, 0.76, "NDC");
    #txt_phi.SetFillColor(kWhite);
    #txt_phi.SetFillStyle(0);
    #txt_phi.SetBorderSize(0);
    #txt_phi.SetTextAlign(12);#middle,left
    #txt_phi.SetTextFont(42);#helvetica
    #txt_phi.SetTextSize(0.05);
    #txt_phi.AddText("#it{{m}}_{{#phi}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(6) * 1e+3, f1.GetParError(6) * 1e+3));
    #txt_phi.AddText("#sigma_{{#phi}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(7) * 1e+3, f1.GetParError(7) * 1e+3));
    #txt_phi.Draw();
    #ROOT.SetOwnership(txt_phi, False);


    ROOT.SetOwnership(p2,False);
    ROOT.SetOwnership(frame2,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mmumu_uls_ls_2panel(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    outfile_tmp = TFile("20241125_pp_13.6TeV_pass7_{0}.root".format(taskname), "RECREATE");
    outfile_tmp.WriteTObject(h1m_uls);
    outfile_tmp.WriteTObject(h1m_bkg);
    outfile_tmp.WriteTObject(h1m_sig);
    outfile_tmp.Close();


    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.165, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogy(1);

    h1m_sig.GetXaxis().SetRangeUser(0.21, 12);
    y_range_max = h1m_uls.GetMaximum() * 2;
    y_range_min = max(h1m_sig.GetMinimum() * 0.5, 2e-12);
    h1m_sig.GetXaxis().SetRangeUser(0.2, 12);

    #frame1 = c1.DrawFrame(0.0, y_range_min, 4, y_range_max);
    frame1 = c1.DrawFrame(0.0, y_range_min, 12, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.00);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.DrawClone("E0same");
    ROOT.SetOwnership(h1m_sig, False);

    #leg = TLegend(0.65, 0.5, 0.85, 0.65);
    leg = TLegend(0.7, 0.40, 0.85, 0.52);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.AddEntry(h1m_uls, "#it{N}_{+#minus}", "P");
    leg.AddEntry(h1m_bkg, "#it{B} = 2#it{R} #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "#it{S} = #it{N}_{+#minus} #minus #it{B}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #txt = TPaveText(0.45, 0.65, 0.65, 0.95, "NDC");
    txt = TPaveText(0.19, 0.12, 0.39, 0.42, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);

    minptmu = 0.3;
    if "minpt200" in suffix:
        minptmu = 0.2;
    elif "minpt300" in suffix:
        minptmu = 0.3;
    elif "minpt400" in suffix:
        minptmu = 0.4;
    elif "minpt500" in suffix:
        minptmu = 0.5;

    str_min_etamu = "#minus3.6";
    str_muon_type = "MFT-MCH-MID";
    if "global" in  suffix:
        str_muon_type = "MFT-MCH-MID";
        str_min_etamu = "#minus3.6";
    elif "standalone" in suffix:
        str_muon_type = "MCH-MID";
        str_min_etamu = "#minus4.0";

    txt.AddText("ALICE Performance");
    if "pp" in suffix:
        #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
        txt.AddText("pp at #sqrt{#it{s}} = 5.36 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    txt.AddText("#it{{p}}_{{T,#mu}} > {0:2.1f} GeV/#it{{c}}, {1} < #it{{#eta}}_{{#mu}} < #minus2.5".format(minptmu, str_min_etamu));
    txt.AddText("{0} < #it{{y}}_{{#mu#mu}} < #minus2.5".format(str_min_etamu));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{#mu#mu}}^{{XY}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.AddText(str_muon_type);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    h1m_sig.GetXaxis().SetRangeUser(0.5, 1.2);
    y_max2 = h1m_sig.GetMaximum() * 1.8;
    h1m_sig.GetXaxis().SetRangeUser(0, 12);

    p2 = TPad("p2", "p2", 0.5, 0.51, 0.95, 0.96);
    p2.SetMargin(0.22, 0.03, 0.12, 0.06);
    p2.SetTicks(1,1);
    c1.cd();
    p2.Draw("");
    p2.cd();
    frame2 = p2.DrawFrame(0.6, 0, 1.1, y_max2);
    #frame2 = p2.DrawFrame(0.5, 0, 1.2, y_max2);
    frame2.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{#mu#mu}} (GeV/#it{c}^{2})^{#minus1}");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(1.8);
    frame2.GetXaxis().SetTitleSize(0.05);
    frame2.GetYaxis().SetTitleSize(0.05);
    frame2.GetXaxis().SetLabelSize(0.05);
    frame2.GetYaxis().SetLabelSize(0.05);
    frame2.GetYaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetNdivisions(505);
    frame2.GetXaxis().SetNdivisions(510);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    #frame2.GetXaxis().SetRangeUser(0.1, 1.1);
    frame2.GetXaxis().SetNdivisions(505, False);
    h1m_sig.Draw("E0same");

    #f1 = TF1("f1", "crystalball(0) + crystalball(5) + [10]*exp(-x/[11]) + [12]", 0.2, 1.2);
    f1 = TF1("f1", "crystalball(0) + crystalball(5) + [10]+[11]*x", 0.2, 1.2);
    f1.SetLineColor(kGreen+2);
    f1.SetLineWidth(2);
    f1.SetLineStyle(2);
    f1.SetNpx(1000);
    f1.SetParameters(1e-5, 0.782, 0.02, 0.2, 0.2, 5e-6, 1.019, 0.02, 0.3, 0.3);
    f1.SetParameter(10, 2e-6);
    #f1.SetParameter(11, -0.1);
    #f1.SetParameter(12, 1e-3);
    f1.SetParLimits(0, 0, 1e-3);
    f1.SetParLimits(5, 0, 1e-3);
    f1.SetParLimits(1, 0.7, 0.9);
    f1.SetParLimits(6, 0.9, 1.1);
    f1.SetParLimits(2, 1e-3, 0.1);
    f1.SetParLimits(7, 1e-3, 0.1);
    h1m_sig.Fit(f1, "ISME", "", 0.6, 1.1);
    ROOT.SetOwnership(f1, False);
    print(f1.GetChisquare()/f1.GetNDF());

    txt_omega = TPaveText(0.3, 0.78, 0.5, 0.9, "NDC");
    txt_omega.SetFillColor(kWhite);
    txt_omega.SetFillStyle(0);
    txt_omega.SetBorderSize(0);
    txt_omega.SetTextAlign(12);#middle,left
    txt_omega.SetTextFont(42);#helvetica
    txt_omega.SetTextSize(0.05);
    txt_omega.AddText("#it{{m}}_{{#omega}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(1) * 1e+3, f1.GetParError(1) * 1e+3));
    txt_omega.AddText("#sigma_{{#omega}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(2) * 1e+3, f1.GetParError(2) * 1e+3));
    txt_omega.Draw();
    ROOT.SetOwnership(txt_omega, False);

    txt_phi = TPaveText(0.3, 0.64, 0.5, 0.76, "NDC");
    txt_phi.SetFillColor(kWhite);
    txt_phi.SetFillStyle(0);
    txt_phi.SetBorderSize(0);
    txt_phi.SetTextAlign(12);#middle,left
    txt_phi.SetTextFont(42);#helvetica
    txt_phi.SetTextSize(0.05);
    txt_phi.AddText("#it{{m}}_{{#phi}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(6) * 1e+3, f1.GetParError(6) * 1e+3));
    txt_phi.AddText("#sigma_{{#phi}} = {0:2.1f} #pm {1:2.1f} (MeV/#it{{c}}^{{2}})".format(f1.GetParameter(7) * 1e+3, f1.GetParError(7) * 1e+3));
    txt_phi.Draw();
    ROOT.SetOwnership(txt_phi, False);


    ROOT.SetOwnership(p2,False);
    ROOT.SetOwnership(frame2,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndmmumu_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_2panel.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_Rfactor(filename, taskname, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_R = rootdire.Get("hs_R");

    bin0 = hs_R.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_R.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_R.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_R.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_R.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_R.GetAxis(2).SetRange(bin0, bin1);

    h1m_r = hs_R.Projection(0);

    h1m_r.SetDirectory(0);
    ROOT.SetOwnership(h1m_r, False);
    make_common_style(h1m_r, 20, 1.2, kBlack, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.18, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(0., 0.9, 4, 1.06);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line1_40 = TLine(0, 1, 4.0, 1);
    line1_40.SetLineColor(kBlack);
    line1_40.SetLineStyle(2);
    line1_40.SetLineWidth(2);
    line1_40.Draw("same");
    ROOT.SetOwnership(line1_40,False);

    h1m_r.Draw("E0same");

    txt = TPaveText(0.25, 0.7, 0.6, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #p2 = TPad("p2", "p2", 0.35, 0.12, 0.92, 0.62);
    #p2.SetMargin(0.25, 0.03, 0.12, 0.02);
    #p2.SetTicks(1,1);
    #c1.cd();
    #p2.Draw("");
    #p2.cd();
    #frame2 = p2.DrawFrame(0.,0.9, 1.1, 1.02);
    #frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    #frame2.GetYaxis().SetTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
    #frame2.GetXaxis().SetTitleOffset(1.1);
    #frame2.GetYaxis().SetTitleOffset(2.3);
    #frame2.GetXaxis().SetTitleSize(0.05);
    #frame2.GetYaxis().SetTitleSize(0.05);
    #frame2.GetXaxis().SetLabelSize(0.05);
    #frame2.GetYaxis().SetLabelSize(0.05);
    #frame2.GetYaxis().SetMaxDigits(3);
    #frame2.GetXaxis().SetNdivisions(510);
    #frame2.GetYaxis().SetNdivisions(510);
    #line1_11 = TLine(0, 1, 1.1, 1);
    #line1_11.SetLineColor(kBlack);
    #line1_11.SetLineStyle(2);
    #line1_11.SetLineWidth(2);
    #line1_11.Draw("same");
    #ROOT.SetOwnership(line1_11, False);
    #h1m_r.Draw("E0same");
    #ROOT.SetOwnership(p2,False);
    #ROOT.SetOwnership(frame2,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_Rfactor_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();


#__________________________________________________________
def draw_mee_dcaee(filename, taskname, cen1, cen2, arr_m, ptmin, ptmax, suffix=""):

    print(filename);
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    h1m_uls_all_org = hs_uls.Projection(0);
    h1m_bkg_all_org = hs_bkg.Projection(0);
    h1m_sig_all_org = hs_sig.Projection(0);
    #h1m_uls_all.Scale(1, "width");
    #h1m_bkg_all.Scale(1, "width");
    #h1m_sig_all.Scale(1, "width");

    bin0 = hs_uls.GetAxis(2).FindBin(0 + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(1.0 - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_short_org = hs_uls.Projection(0);
    h1m_bkg_short_org = hs_bkg.Projection(0);
    h1m_sig_short_org = hs_sig.Projection(0);
    #h1m_uls_short.Scale(1, "width");
    #h1m_bkg_short.Scale(1, "width");
    #h1m_sig_short.Scale(1, "width");

    bin0 = hs_uls.GetAxis(2).FindBin(2 + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(99 - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_long_org = hs_uls.Projection(0);
    h1m_bkg_long_org = hs_bkg.Projection(0);
    h1m_sig_long_org = hs_sig.Projection(0);
    #h1m_uls_long.Scale(1, "width");
    #h1m_bkg_long.Scale(1, "width");
    #h1m_sig_long.Scale(1, "width");

    h1m_sig_all = rebin_histogram(h1m_sig_all_org, arr_m, False, False);
    h1m_sig_all.Scale(1, "width");
    h1m_sig_short = rebin_histogram(h1m_sig_short_org, arr_m, False, False);
    h1m_sig_short.Scale(1, "width");
    h1m_sig_long = rebin_histogram(h1m_sig_long_org, arr_m, False, False);
    h1m_sig_long.Scale(1, "width");

    #h1m_uls.SetDirectory(0);
    #h1m_bkg.SetDirectory(0);
    #h1m_sig.SetDirectory(0);
    #ROOT.SetOwnership(h1m_uls, False);
    #ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig_all, False);
    ROOT.SetOwnership(h1m_sig_short, False);
    ROOT.SetOwnership(h1m_sig_long, False);

    h1m_sig_all.SetDirectory(0);
    h1m_sig_short.SetDirectory(0);
    h1m_sig_long.SetDirectory(0);

    #make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    #make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig_all  , 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_sig_short, 21, 1.2, kRed+1, 2, 0);
    make_common_style(h1m_sig_long , 34, 1.5, kBlue+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.03, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = h1m_sig_all.GetMaximum() * 10;
    #y_range_min = max(h1m_sig_all.GetMinimum() * 0.1, 1e-8);
    y_range_max = 1e-2
    y_range_min = 1e-7;

    frame1 = c1.DrawFrame(0.14, y_range_min, 4.0, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    #h1m_uls.Draw("E0same");
    #h1m_bkg.Draw("E0same");

    h1m_sig_all.Draw("E0same");
    h1m_sig_long.Draw("E0same");
    h1m_sig_short.Draw("E0same");

    leg = TLegend(0.7, 0.75, 0.9, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    #leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    #leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig_all, "all DCA", "P");
    leg.AddEntry(h1m_sig_short, "DCA < 1.0 #sigma", "P");
    leg.AddEntry(h1m_sig_long, "DCA > 2.0 #sigma", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.75, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE WIP");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    #txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    #txt.AddText("#sigma = DCA resolution");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_dcaee_pt{1:2.1f}_{2:2.1f}GeV{3}.eps".format(date, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dndm_dcaee_pt{1:2.1f}_{2:2.1f}GeV{3}.pdf".format(date, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dndm_dcaee_pt{1:2.1f}_{2:2.1f}GeV{3}.png".format(date, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_v2_vs_ptee_old(filename, taskname, mmin, mmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    h2vn = rootdire.Get("h2vn_dca0");
    bin0 = h2vn.GetXaxis().FindBin(mmin + 1e-3);
    bin1 = h2vn.GetXaxis().FindBin(mmax - 1e-3);
    h1vn = h2vn.ProjectionY("h1vn", bin0, bin1, "");

    if mmax < 0.2:
        h1vn.SetBinContent(1, -1);
        h1vn.SetBinError(1, 0);

    h1vn.SetDirectory(0);
    ROOT.SetOwnership(h1vn, False);
    make_common_style(h1vn, 20, 1.5, kRed+1, 2, 0);

    #pion data
    filename_pion = "$HOME/analysis/EM/photon/HEPdata/identified_hadron_flow_PbPb_5.02TeV/charged_pion_vn_PbPb_5.02TeV.root";
    rootfile_pion = TFile.Open(filename_pion, "READ");
    h1v2_pi_3040_stat = rootfile_pion.Get("h1v2_pi_3040_stat");
    h1v2_pi_3040_syst = rootfile_pion.Get("h1v2_pi_3040_syst");
    h1v2_pi_4050_stat = rootfile_pion.Get("h1v2_pi_4050_stat");
    h1v2_pi_4050_syst = rootfile_pion.Get("h1v2_pi_4050_syst");

    h1v2_pi_3040_stat.SetDirectory(0);
    h1v2_pi_3040_syst.SetDirectory(0);
    h1v2_pi_4050_stat.SetDirectory(0);
    h1v2_pi_4050_syst.SetDirectory(0);
    ROOT.SetOwnership(h1v2_pi_3040_stat, False);
    ROOT.SetOwnership(h1v2_pi_3040_syst, False);
    ROOT.SetOwnership(h1v2_pi_4050_stat, False);
    ROOT.SetOwnership(h1v2_pi_4050_syst, False);

    make_common_style(h1v2_pi_3040_stat, 26, 1.5, kGray+2, 2, 0);
    make_common_style(h1v2_pi_3040_syst, 26, 1.5, kGray+2, 2, 0);
    make_common_style(h1v2_pi_4050_stat, 32, 1.5, kGray+2, 2, 0);
    make_common_style(h1v2_pi_4050_syst, 32, 1.5, kGray+2, 2, 0);

    #jpsi data
    filename_jpsi = "$HOME/analysis/EM/pwgem_sw_run3/dilepton/HEPData/jpsi_flow_PbPb_5.02TeV/jpsi_vn_PbPb_5.02TeV.root";
    rootfile_jpsi = TFile.Open(filename_jpsi, "READ");
    rootfile_jpsi.ls();
    g1v2_jpsi_3050_stat = rootfile_jpsi.Get("g1v2_jpsi_3050_stat");
    g1v2_jpsi_3050_syst = rootfile_jpsi.Get("g1v2_jpsi_3050_syst");
    ROOT.SetOwnership(g1v2_jpsi_3050_stat, False);
    ROOT.SetOwnership(g1v2_jpsi_3050_syst, False);
    make_common_style(g1v2_jpsi_3050_stat, 26, 1.5, kGray+2, 2, 0);
    make_common_style(g1v2_jpsi_3050_syst, 26, 1.5, kGray+2, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = +1;
    y_range_min = -0.1;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{v}_{2}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.1);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.15, 0.6, 0.35, 0.75);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn, "#it{v}_{2,ee}", "P");

    if mmax < 0.2:
        h1v2_pi_3040_stat.Draw("E0same");
        h1v2_pi_3040_syst.Draw("E2same");
        h1v2_pi_4050_stat.Draw("E0same");
        h1v2_pi_4050_syst.Draw("E2same");
        leg.AddEntry(h1v2_pi_3040_stat, "#it{v}_{2} #pi^{#pm}, 30-40%, Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV", "P");
        leg.AddEntry(h1v2_pi_4050_stat, "#it{v}_{2} #pi^{#pm}, 40-50%, Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV", "P");

    if 2.7 < mmin and mmax < 3.3:
        g1v2_jpsi_3050_stat.Draw("PZ");
        g1v2_jpsi_3050_syst.Draw("P2Z");
        leg.AddEntry(g1v2_jpsi_3050_stat, "#it{v}_{2} J/#psi, 30-50%, Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV", "P");

    h1vn.Draw("E0same");

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.75, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
    rootfile_pion.Close();
    rootfile_jpsi.Close();
#__________________________________________________________
def draw_v2_vs_dcaee(filename, taskname, cen1, cen2, arr_dca, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    print(arr_dca);

    hs_ry_uls = rootdire.Get("hs_uls"); #raw yield
    hs_ry_bkg = rootdire.Get("hs_bkg"); #raw yield
    hs_ry_sig = rootdire.Get("hs_sig"); #raw yield
    ROOT.SetOwnership(hs_ry_uls, False);
    ROOT.SetOwnership(hs_ry_bkg, False);
    ROOT.SetOwnership(hs_ry_sig, False);

    hs_vn_uls = rootdire.Get("hs_vn_tot"); #vn
    hs_vn_bkg = rootdire.Get("hs_vn_bkg"); #vn
    hs_vn_sig = rootdire.Get("hs_vn_sig"); #vn
    ROOT.SetOwnership(hs_vn_uls, False);
    ROOT.SetOwnership(hs_vn_bkg, False);
    ROOT.SetOwnership(hs_vn_sig, False);

    h1vn_uls = TH1D("h1vn_uls", "vn", len(arr_dca)-1, arr_dca);
    h1vn_uls.SetName("h1vn_sig");
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls, False);
    ROOT.SetOwnership(h1vn_bkg, False);
    ROOT.SetOwnership(h1vn_sig, False);

    for idca in range(0, len(arr_dca)-1):
        dcamin = arr_dca[idca];
        dcamax = arr_dca[idca+1];
        [vn_uls, vn_uls_err] = rebin_flow_thn(hs_ry_uls, hs_vn_uls, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_uls.SetBinContent(idca+1, vn_uls);
        h1vn_uls.SetBinError(idca+1, vn_uls_err);

        [vn_bkg, vn_bkg_err] = rebin_flow_thn(hs_ry_bkg, hs_vn_bkg, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_bkg.SetBinContent(idca+1, vn_bkg);
        h1vn_bkg.SetBinError(idca+1, vn_bkg_err);

        [vn_sig, vn_sig_err] = rebin_flow_thn(hs_ry_sig, hs_vn_sig, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_sig.SetBinContent(idca+1, vn_sig);
        h1vn_sig.SetBinError(idca+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = +0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 10, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    h1vn_uls.Draw("E0same");
    h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_dcaee_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_v2_vs_dcaee_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_v2_vs_dcaee_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_v2_vs_ptee(filename, taskname, cen1, cen2, arr_pt, mmin, mmax, dcamin, dcamax, suffix="", draw_pion=False, draw_jpsi=False):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    print(arr_pt);

    hs_ry_uls = rootdire.Get("hs_uls"); #raw yield
    hs_ry_bkg = rootdire.Get("hs_bkg"); #raw yield
    hs_ry_sig = rootdire.Get("hs_sig"); #raw yield
    ROOT.SetOwnership(hs_ry_uls, False);
    ROOT.SetOwnership(hs_ry_bkg, False);
    ROOT.SetOwnership(hs_ry_sig, False);

    hs_vn_uls = rootdire.Get("hs_vn_tot"); #vn
    hs_vn_bkg = rootdire.Get("hs_vn_bkg"); #vn
    hs_vn_sig = rootdire.Get("hs_vn_sig"); #vn
    ROOT.SetOwnership(hs_vn_uls, False);
    ROOT.SetOwnership(hs_vn_bkg, False);
    ROOT.SetOwnership(hs_vn_sig, False);

    h1vn_uls = TH1D("h1vn_uls", "vn", len(arr_pt)-1, arr_pt);
    h1vn_uls.SetName("h1vn_sig");
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls, False);
    ROOT.SetOwnership(h1vn_bkg, False);
    ROOT.SetOwnership(h1vn_sig, False);

    for ipt in range(0, len(arr_pt)-1):
        ptmin = arr_pt[ipt];
        ptmax = arr_pt[ipt+1];
        [vn_uls, vn_uls_err] = rebin_flow_thn(hs_ry_uls, hs_vn_uls, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_uls.SetBinContent(ipt+1, vn_uls);
        h1vn_uls.SetBinError(ipt+1, vn_uls_err);

        [vn_bkg, vn_bkg_err] = rebin_flow_thn(hs_ry_bkg, hs_vn_bkg, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_bkg.SetBinContent(ipt+1, vn_bkg);
        h1vn_bkg.SetBinError(ipt+1, vn_bkg_err);

        [vn_sig, vn_sig_err] = rebin_flow_thn(hs_ry_sig, hs_vn_sig, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_sig.SetBinContent(ipt+1, vn_sig);
        h1vn_sig.SetBinError(ipt+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = +0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 10, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    if draw_pion:
        #pion data
        filename_pion = "$HOME/analysis/EM/photon/HEPdata/identified_hadron_flow_PbPb_5.02TeV/charged_pion_vn_PbPb_5.02TeV.root";
        rootfile_pion = TFile.Open(filename_pion, "READ");
        h1v2_pi_3040_stat = rootfile_pion.Get("h1v2_pi_3040_stat");
        h1v2_pi_3040_syst = rootfile_pion.Get("h1v2_pi_3040_syst");
        h1v2_pi_4050_stat = rootfile_pion.Get("h1v2_pi_4050_stat");
        h1v2_pi_4050_syst = rootfile_pion.Get("h1v2_pi_4050_syst");

        h1v2_pi_3040_stat.SetDirectory(0);
        h1v2_pi_3040_syst.SetDirectory(0);
        h1v2_pi_4050_stat.SetDirectory(0);
        h1v2_pi_4050_syst.SetDirectory(0);
        ROOT.SetOwnership(h1v2_pi_3040_stat, False);
        ROOT.SetOwnership(h1v2_pi_3040_syst, False);
        ROOT.SetOwnership(h1v2_pi_4050_stat, False);
        ROOT.SetOwnership(h1v2_pi_4050_syst, False);

        make_common_style(h1v2_pi_3040_stat, 26, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_3040_syst, 26, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_4050_stat, 32, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_4050_syst, 32, 1.5, kGray+2, 2, 0);
        h1v2_pi_3040_stat.Draw("E0same");
        h1v2_pi_3040_syst.Draw("E2same");
        h1v2_pi_4050_stat.Draw("E0same");
        h1v2_pi_4050_syst.Draw("E2same");

        leg_pion = TLegend(0.15, 0.12, 0.35, 0.26);
        leg_pion.SetBorderSize(0);
        leg_pion.SetFillColor(kWhite);
        leg_pion.SetFillStyle(0);
        leg_pion.SetTextSize(0.035);
        leg_pion.SetHeader("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (JHEP 09 (2018) 006, 2018)");
        leg_pion.AddEntry(h1v2_pi_3040_stat, "#it{v}_{2} #pi^{#pm}, 30-40%", "P");
        leg_pion.AddEntry(h1v2_pi_4050_stat, "#it{v}_{2} #pi^{#pm}, 40-50%", "P");
        leg_pion.Draw("");
        ROOT.SetOwnership(leg_pion, False);

    if draw_jpsi:
        #jpsi data
        filename_jpsi = "$HOME/analysis/EM/pwgem_sw_run3/dilepton/HEPData/jpsi_flow_PbPb_5.02TeV/jpsi_vn_PbPb_5.02TeV.root";
        rootfile_jpsi = TFile.Open(filename_jpsi, "READ");
        rootfile_jpsi.ls();
        g1v2_jpsi2ee_3050_stat = rootfile_jpsi.Get("g1v2_jpsi2ee_3050_stat");
        g1v2_jpsi2ee_3050_syst = rootfile_jpsi.Get("g1v2_jpsi2ee_3050_syst");
        ROOT.SetOwnership(g1v2_jpsi2ee_3050_stat, False);
        ROOT.SetOwnership(g1v2_jpsi2ee_3050_syst, False);
        make_common_style(g1v2_jpsi2ee_3050_stat, 26, 1.5, kGray+2, 2, 0);
        make_common_style(g1v2_jpsi2ee_3050_syst, 26, 1.5, kGray+2, 2, 0);
        g1v2_jpsi2mumu_3050_stat = rootfile_jpsi.Get("g1v2_jpsi2mumu_3050_stat");
        g1v2_jpsi2mumu_3050_syst = rootfile_jpsi.Get("g1v2_jpsi2mumu_3050_syst");
        ROOT.SetOwnership(g1v2_jpsi2mumu_3050_stat, False);
        ROOT.SetOwnership(g1v2_jpsi2mumu_3050_syst, False);
        make_common_style(g1v2_jpsi2mumu_3050_stat, 32, 1.5, kGray+2, 2, 0);
        make_common_style(g1v2_jpsi2mumu_3050_syst, 32, 1.5, kGray+2, 2, 0);

        g1v2_jpsi2ee_3050_stat  .Draw("PZ");
        g1v2_jpsi2ee_3050_syst  .Draw("P2");
        g1v2_jpsi2mumu_3050_stat.Draw("PZ");
        g1v2_jpsi2mumu_3050_syst.Draw("P2");

        leg_jpsi = TLegend(0.15, 0.12, 0.35, 0.26);
        leg_jpsi.SetBorderSize(0);
        leg_jpsi.SetFillColor(kWhite);
        leg_jpsi.SetFillStyle(0);
        leg_jpsi.SetTextSize(0.035);
        leg_jpsi.SetHeader("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (JHEP 10 (2020) 141, 2020)");
        leg_jpsi.AddEntry(g1v2_jpsi2ee_3050_stat  , "#it{v}_{2} J/#psi, |#it{y}| < 0.9, 30-50%", "P");
        leg_jpsi.AddEntry(g1v2_jpsi2mumu_3050_stat, "#it{v}_{2} J/#psi, -4.0 < #it{y} < -2.5, 30-50%", "P");
        leg_jpsi.Draw("");
        ROOT.SetOwnership(leg_jpsi, False);

    h1vn_uls.Draw("E0same");
    h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_mee{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_v2_vs_mee(filename, taskname, cen1, cen2, arr_m, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    print(arr_m);

    hs_ry_uls = rootdire.Get("hs_uls"); #raw yield
    hs_ry_bkg = rootdire.Get("hs_bkg"); #raw yield
    hs_ry_sig = rootdire.Get("hs_sig"); #raw yield
    ROOT.SetOwnership(hs_ry_uls, False);
    ROOT.SetOwnership(hs_ry_bkg, False);
    ROOT.SetOwnership(hs_ry_sig, False);

    hs_vn_uls = rootdire.Get("hs_vn_tot"); #vn
    hs_vn_bkg = rootdire.Get("hs_vn_bkg"); #vn
    hs_vn_sig = rootdire.Get("hs_vn_sig"); #vn
    ROOT.SetOwnership(hs_vn_uls, False);
    ROOT.SetOwnership(hs_vn_bkg, False);
    ROOT.SetOwnership(hs_vn_sig, False);

    h1vn_uls = TH1D("h1vn_uls", "vn", len(arr_m)-1, arr_m);
    h1vn_uls.SetName("h1vn_sig");
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls, False);
    ROOT.SetOwnership(h1vn_bkg, False);
    ROOT.SetOwnership(h1vn_sig, False);

    for im in range(0, len(arr_m)-1):
        mmin = arr_m[im];
        mmax = arr_m[im+1];
        [vn_uls, vn_uls_err] = rebin_flow_thn(hs_ry_uls, hs_vn_uls, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_uls.SetBinContent(im+1, vn_uls);
        h1vn_uls.SetBinError(im+1, vn_uls_err);

        [vn_bkg, vn_bkg_err] = rebin_flow_thn(hs_ry_bkg, hs_vn_bkg, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_bkg.SetBinContent(im+1, vn_bkg);
        h1vn_bkg.SetBinError(im+1, vn_bkg_err);

        [vn_sig, vn_sig_err] = rebin_flow_thn(hs_ry_sig, hs_vn_sig, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        h1vn_sig.SetBinContent(im+1, vn_sig);
        h1vn_sig.SetBinError(im+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = 0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 4, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    h1vn_uls.Draw("E0same");
    h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_v2_vs_ptee_new_multiple_dca(filename, taskname, cen1, cen2, arr_pt, mmin, mmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();

    hs_uls_same  = rootdire.Get("Pair/same/uls/hs");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hs");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hs");
    hs_uls_mix  = rootdire.Get("Pair/mix/uls/hs");
    hs_lspp_mix = rootdire.Get("Pair/mix/lspp/hs");
    hs_lsmm_mix = rootdire.Get("Pair/mix/lsmm/hs");

    arr_dca_min = np.array([ 0, 0,  2], dtype=float);
    arr_dca_max = np.array([20, 1, 20], dtype=float);
    colors = [kBlack, kRed+1, kBlue+1];

    list_h1 = [];
    for idca in range(0, len(arr_dca_min)):
        h1vn = TH1D("h1vn_dca{0:d}".format(idca), "h1vn_dca{0:d}".format(idca), len(arr_pt) -1, arr_pt);
        h1vn.Sumw2();
        dcamin = arr_dca_min[idca];
        dcamax = arr_dca_max[idca];
        for ipt in range(0, len(arr_pt)-1):
            ptmin = arr_pt[ipt];
            ptmax = arr_pt[ipt+1];
            [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err] = extract_flow(hs_uls_same, hs_lspp_same, hs_lsmm_same, hs_uls_mix, hs_lspp_mix, hs_lsmm_mix, 1.0, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
            h1vn.SetBinContent(ipt+1, vn_sig);
            h1vn.SetBinError(ipt+1, vn_sig_err);
        make_common_style(h1vn, 20, 1.2, colors[idca], 2, 0);
        h1vn.SetDirectory(0);
        ROOT.SetOwnership(h1vn ,False);
        list_h1.append(h1vn);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = 0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 10, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    leg = TLegend(0.58, 0.75, 0.75, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.038);

    for ih in range(0, len(list_h1)):
        list_h1[ih].Draw("E1same");
        dcamin = arr_dca_min[ih];
        dcamax = arr_dca_max[ih];
        leg.AddEntry(list_h1[ih], "{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax), "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.038);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV.eps".format(date, mmin, mmax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV.pdf".format(date, mmin, mmax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV.png".format(date, mmin, mmax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_v2_vs_dcaee_new(filename, taskname, cen1, cen2, arr_dca, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    print(arr_dca);

    hs_uls_same  = rootdire.Get("Pair/same/uls/hs");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hs");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hs");
    hs_uls_mix  = rootdire.Get("Pair/mix/uls/hs");
    hs_lspp_mix = rootdire.Get("Pair/mix/lspp/hs");
    hs_lsmm_mix = rootdire.Get("Pair/mix/lsmm/hs");

    h1vn_uls = TH1D("h1vn_uls", "h1vn_uls", len(arr_dca) -1, arr_dca);
    h1vn_uls.Sumw2();
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    for idca in range(0, len(arr_dca)-1):
        dcamin = arr_dca[idca];
        dcamax = arr_dca[idca+1];
        [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err] = extract_flow(hs_uls_same, hs_lspp_same, hs_lsmm_same, hs_uls_mix, hs_lspp_mix, hs_lsmm_mix, 1.0, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        #print(vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err);
        h1vn_uls.SetBinContent(idca+1, vn_tot);
        h1vn_uls.SetBinError(idca+1, vn_tot_err);
        h1vn_bkg.SetBinContent(idca+1, vn_bkg);
        h1vn_bkg.SetBinError(idca+1, vn_bkg_err);
        h1vn_sig.SetBinContent(idca+1, vn_sig);
        h1vn_sig.SetBinError(idca+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);
    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls ,False);
    ROOT.SetOwnership(h1vn_bkg ,False);
    ROOT.SetOwnership(h1vn_sig ,False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = 0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 10, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    h1vn_uls.Draw("E0same");
    h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_dcaee_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}sigma{5}_new.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_v2_vs_dcaee_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}sigma{5}_new.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_v2_vs_dcaee_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}sigma{5}_new.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_v2_vs_ptee_new(filename, taskname, cen1, cen2, arr_pt, mmin, mmax, dcamin, dcamax, suffix="", draw_pion=False, draw_jpsi=False, draw_cocktail=False):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    print(arr_m);

    hs_uls_same  = rootdire.Get("Pair/same/uls/hs");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hs");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hs");
    hs_uls_mix  = rootdire.Get("Pair/mix/uls/hs");
    hs_lspp_mix = rootdire.Get("Pair/mix/lspp/hs");
    hs_lsmm_mix = rootdire.Get("Pair/mix/lsmm/hs");

    h1vn_uls = TH1D("h1vn_uls", "h1vn_uls", len(arr_pt) -1, arr_pt);
    h1vn_uls.Sumw2();
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    for ipt in range(0, len(arr_pt)-1):
        ptmin = arr_pt[ipt];
        ptmax = arr_pt[ipt+1];
        [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err] = extract_flow(hs_uls_same, hs_lspp_same, hs_lsmm_same, hs_uls_mix, hs_lspp_mix, hs_lsmm_mix, 1.0, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        #print(vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err);
        h1vn_uls.SetBinContent(ipt+1, vn_tot);
        h1vn_uls.SetBinError(ipt+1, vn_tot_err);
        h1vn_bkg.SetBinContent(ipt+1, vn_bkg);
        h1vn_bkg.SetBinError(ipt+1, vn_bkg_err);
        h1vn_sig.SetBinContent(ipt+1, vn_sig);
        h1vn_sig.SetBinError(ipt+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);
    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls ,False);
    ROOT.SetOwnership(h1vn_bkg ,False);
    ROOT.SetOwnership(h1vn_sig ,False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = 0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 10, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    if draw_pion:
        #pion data
        filename_pion = "$HOME/analysis/EM/photon/HEPdata/identified_hadron_flow_PbPb_5.02TeV/charged_pion_vn_PbPb_5.02TeV.root";
        rootfile_pion = TFile.Open(filename_pion, "READ");
        h1v2_pi_3040_stat = rootfile_pion.Get("h1v2_pi_3040_stat");
        h1v2_pi_3040_syst = rootfile_pion.Get("h1v2_pi_3040_syst");
        h1v2_pi_4050_stat = rootfile_pion.Get("h1v2_pi_4050_stat");
        h1v2_pi_4050_syst = rootfile_pion.Get("h1v2_pi_4050_syst");

        h1v2_pi_3040_stat.SetDirectory(0);
        h1v2_pi_3040_syst.SetDirectory(0);
        h1v2_pi_4050_stat.SetDirectory(0);
        h1v2_pi_4050_syst.SetDirectory(0);
        ROOT.SetOwnership(h1v2_pi_3040_stat, False);
        ROOT.SetOwnership(h1v2_pi_3040_syst, False);
        ROOT.SetOwnership(h1v2_pi_4050_stat, False);
        ROOT.SetOwnership(h1v2_pi_4050_syst, False);

        make_common_style(h1v2_pi_3040_stat, 26, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_3040_syst, 26, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_4050_stat, 32, 1.5, kGray+2, 2, 0);
        make_common_style(h1v2_pi_4050_syst, 32, 1.5, kGray+2, 2, 0);
        h1v2_pi_3040_stat.Draw("E0same");
        h1v2_pi_3040_syst.Draw("E2same");
        h1v2_pi_4050_stat.Draw("E0same");
        h1v2_pi_4050_syst.Draw("E2same");

        leg_pion = TLegend(0.15, 0.12, 0.35, 0.26);
        leg_pion.SetBorderSize(0);
        leg_pion.SetFillColor(kWhite);
        leg_pion.SetFillStyle(0);
        leg_pion.SetTextSize(0.035);
        leg_pion.SetHeader("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (JHEP 09 (2018) 006, 2018)");
        leg_pion.AddEntry(h1v2_pi_3040_stat, "#it{v}_{2} #pi^{#pm}, 30-40%", "P");
        leg_pion.AddEntry(h1v2_pi_4050_stat, "#it{v}_{2} #pi^{#pm}, 40-50%", "P");
        leg_pion.Draw("");
        ROOT.SetOwnership(leg_pion, False);

    if draw_jpsi:
        #jpsi data
        filename_jpsi = "$HOME/analysis/EM/pwgem_sw_run3/dilepton/HEPData/jpsi_flow_PbPb_5.02TeV/jpsi_vn_PbPb_5.02TeV.root";
        rootfile_jpsi = TFile.Open(filename_jpsi, "READ");
        rootfile_jpsi.ls();
        g1v2_jpsi2ee_3050_stat = rootfile_jpsi.Get("g1v2_jpsi2ee_3050_stat");
        g1v2_jpsi2ee_3050_syst = rootfile_jpsi.Get("g1v2_jpsi2ee_3050_syst");
        ROOT.SetOwnership(g1v2_jpsi2ee_3050_stat, False);
        ROOT.SetOwnership(g1v2_jpsi2ee_3050_syst, False);
        make_common_style(g1v2_jpsi2ee_3050_stat, 26, 1.5, kGray+2, 2, 0);
        make_common_style(g1v2_jpsi2ee_3050_syst, 26, 1.5, kGray+2, 2, 0);
        g1v2_jpsi2mumu_3050_stat = rootfile_jpsi.Get("g1v2_jpsi2mumu_3050_stat");
        g1v2_jpsi2mumu_3050_syst = rootfile_jpsi.Get("g1v2_jpsi2mumu_3050_syst");
        ROOT.SetOwnership(g1v2_jpsi2mumu_3050_stat, False);
        ROOT.SetOwnership(g1v2_jpsi2mumu_3050_syst, False);
        make_common_style(g1v2_jpsi2mumu_3050_stat, 32, 1.5, kGray+2, 2, 0);
        make_common_style(g1v2_jpsi2mumu_3050_syst, 32, 1.5, kGray+2, 2, 0);

        g1v2_jpsi2ee_3050_stat  .Draw("PZ");
        g1v2_jpsi2ee_3050_syst  .Draw("P2");
        g1v2_jpsi2mumu_3050_stat.Draw("PZ");
        g1v2_jpsi2mumu_3050_syst.Draw("P2");

        leg_jpsi = TLegend(0.15, 0.12, 0.35, 0.26);
        leg_jpsi.SetBorderSize(0);
        leg_jpsi.SetFillColor(kWhite);
        leg_jpsi.SetFillStyle(0);
        leg_jpsi.SetTextSize(0.035);
        leg_jpsi.SetHeader("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (JHEP 10 (2020) 141, 2020)");
        leg_jpsi.AddEntry(g1v2_jpsi2ee_3050_stat  , "#it{v}_{2} J/#psi, |#it{y}| < 0.9, 30-50%", "P");
        leg_jpsi.AddEntry(g1v2_jpsi2mumu_3050_stat, "#it{v}_{2} J/#psi, -4.0 < #it{y} < -2.5, 30-50%", "P");
        leg_jpsi.Draw("");
        ROOT.SetOwnership(leg_jpsi, False);

    if draw_cocktail:
        filename_cocktail = "lmee_lf_cocktail_PbPb_5.36TeV_3050_minpt400_maxeta08.root";
        if "pt200" in suffix  or "pT200" in suffix:
            filename_cocktail = "lmee_lf_cocktail_PbPb_5.36TeV_3050_minpt200_maxeta08.root";

        rootfile_cocktail = TFile.Open(filename_cocktail, "READ");
        hs_sum = rootfile_cocktail.Get("hs_mee_ptee_v2ee_sum");
        bin1 = hs_sum.GetAxis(0).FindBin(mmin + 1e-3);
        bin2 = hs_sum.GetAxis(0).FindBin(mmax - 1e-3);
        hs_sum.GetAxis(0).SetRange(bin1, bin2);
        h2sum = hs_sum.Projection(2, 1);
        h1v2 = h2sum.ProfileX("h1v2_prf").ProjectionX("h1v2");
        h1v2.SetDirectory(0);
        ROOT.SetOwnership(h1v2, False);
        h1v2.GetXaxis().SetRangeUser(arr_pt[0], arr_pt[-1]);
        make_common_style(h1v2, 26, 1.5, kGray+2, 2, 0);
        h1v2.Draw("chist,same");

        leg_cocktail = TLegend(0.15, 0.12, 0.35, 0.26);
        leg_cocktail.SetBorderSize(0);
        leg_cocktail.SetFillColor(kWhite);
        leg_cocktail.SetFillStyle(0);
        leg_cocktail.SetTextSize(0.035);
        #leg_cocktail.SetHeader("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
        leg_cocktail.AddEntry(h1v2, "LF cocktail without syst. unc.", "L");
        leg_cocktail.Draw("");
        ROOT.SetOwnership(leg_cocktail, False);

    #h1vn_uls.Draw("E0same");
    #h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    #leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    #leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_ptee_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_v2_vs_mee_new(filename, taskname, cen1, cen2, arr_m, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    print(arr_m);

    hs_uls_same  = rootdire.Get("Pair/same/uls/hs");
    hs_lspp_same = rootdire.Get("Pair/same/lspp/hs");
    hs_lsmm_same = rootdire.Get("Pair/same/lsmm/hs");
    hs_uls_mix  = rootdire.Get("Pair/mix/uls/hs");
    hs_lspp_mix = rootdire.Get("Pair/mix/lspp/hs");
    hs_lsmm_mix = rootdire.Get("Pair/mix/lsmm/hs");

    h1vn_uls = TH1D("h1vn_uls", "h1vn_uls", len(arr_m) -1, arr_m);
    h1vn_uls.Sumw2();
    h1vn_bkg = h1vn_uls.Clone("h1vn_bkg");
    h1vn_sig = h1vn_uls.Clone("h1vn_sig");

    for im in range(0, len(arr_m)-1):
        mmin = arr_m[im];
        mmax = arr_m[im+1];
        [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err] = extract_flow(hs_uls_same, hs_lspp_same, hs_lsmm_same, hs_uls_mix, hs_lspp_mix, hs_lsmm_mix, 1.0, mmin, mmax, ptmin, ptmax, dcamin, dcamax);
        #print(vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err);
        h1vn_uls.SetBinContent(im+1, vn_tot);
        h1vn_uls.SetBinError(im+1, vn_tot_err);
        h1vn_bkg.SetBinContent(im+1, vn_bkg);
        h1vn_bkg.SetBinError(im+1, vn_bkg_err);
        h1vn_sig.SetBinContent(im+1, vn_sig);
        h1vn_sig.SetBinError(im+1, vn_sig_err);

    make_common_style(h1vn_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1vn_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1vn_sig, 20, 1.2, kRed+1, 2, 0);
    h1vn_uls.SetDirectory(0);
    h1vn_bkg.SetDirectory(0);
    h1vn_sig.SetDirectory(0);
    ROOT.SetOwnership(h1vn_uls ,False);
    ROOT.SetOwnership(h1vn_bkg ,False);
    ROOT.SetOwnership(h1vn_sig ,False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    y_range_max = 0.6;
    y_range_min = -0.2;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#it{v}_{2,ee} {SP, |#Delta#it{#eta}| > 1.3}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    line0 = TLine(0, 0, 4, 0);
    line0.SetLineColor(kGray+2);
    line0.SetLineStyle(2);
    line0.SetLineWidth(2);
    line0.Draw("same");
    ROOT.SetOwnership(line0,False);

    h1vn_uls.Draw("E0same");
    h1vn_bkg.Draw("E0same");
    h1vn_sig.Draw("E0same");

    leg = TLegend(0.8, 0.7, 0.95, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1vn_uls, "#it{v}_{2,ee}^{S+B}", "P");
    leg.AddEntry(h1vn_bkg, "#it{v}_{2,ee}^{B}", "P");
    leg.AddEntry(h1vn_sig, "#it{v}_{2,ee}^{S}", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}%".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_v2_vs_mee_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}_new.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dcaee_test(filename, taskname, cen1, cen2, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_sig.GetAxis(0).FindBin(1.1 + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(2.7 - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);
    bin0 = hs_sig.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_sig.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    h1m_sig1127 = hs_sig.Projection(2);
    nint1127 = 0.0;
    for i in range(0, h1m_sig1127.GetNbinsX()):
        nint1127 += h1m_sig1127.GetBinContent(i+1);
    h1m_sig1127.Scale(1, "width");
    h1m_sig1127.Scale(1/nint1127);
    h1m_sig1127.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig1127, False);

    bin0 = hs_sig.GetAxis(0).FindBin(2.7 + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(3.2 - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);
    h1m_sig2732 = hs_sig.Projection(2);
    nint2732 = 0.0;
    for i in range(0, h1m_sig2732.GetNbinsX()):
        nint2732 += h1m_sig2732.GetBinContent(i+1);
    h1m_sig2732.Scale(1, "width");
    h1m_sig2732.Scale(1/nint2732);
    h1m_sig2732.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig2732, False);


    bin0 = hs_sig.GetAxis(0).FindBin(1.1 + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(3.2 - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);
    h1m_sig1132 = hs_sig.Projection(2);
    nint1132 = 0.0;
    for i in range(0, h1m_sig1132.GetNbinsX()):
        nint1132 += h1m_sig1132.GetBinContent(i+1);
    h1m_sig1132.Scale(1, "width");
    h1m_sig1132.Scale(1/nint1132);
    h1m_sig1132.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig1132, False);


    bin0 = hs_sig.GetAxis(0).FindBin(1.1 + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(2.7 - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);
    bin0 = hs_sig.GetAxis(1).FindBin(0.1 + 1e-3);
    bin1 = hs_sig.GetAxis(1).FindBin(10.0 - 1e-3);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    h1m_sig1127_np = hs_sig.Projection(2);
    nint1127_np = 0.0;
    for i in range(0, h1m_sig1127_np.GetNbinsX()):
        nint1127_np += h1m_sig1127_np.GetBinContent(i+1);
    h1m_sig1127_np.Scale(1, "width");
    h1m_sig1127_np.Scale(1/nint1127_np);
    h1m_sig1127_np.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig1127_np, False);

    make_common_style(h1m_sig1127   , 20, 1.4, kRed+1, 2, 0);
    make_common_style(h1m_sig2732   , 25, 1.4, kRed+1, 2, 0);
    make_common_style(h1m_sig1132   , 21, 1.4, kRed+1, 2, 0);
    make_common_style(h1m_sig1127_np, 20, 1.4, kBlue+1, 2, 0);

    h1m_sig1127   .SetName("h1dca_sig1127_p");
    h1m_sig2732   .SetName("h1dca_sig2732_p");
    h1m_sig1132   .SetName("h1dca_sig1132_p");
    h1m_sig1127_np.SetName("h1dca_sig1127_np");
    h1m_sig1127   .SetTitle("pair DCA template for prompt");
    h1m_sig2732   .SetTitle("pair DCA template for prompt");
    h1m_sig1132   .SetTitle("pair DCA template for prompt");
    h1m_sig1127_np.SetTitle("pair DCA template for nonprompt");
    h1m_sig1127   .SetXTitle("DCA_{ee}^{3D} (#sigma)");
    h1m_sig2732   .SetXTitle("DCA_{ee}^{3D} (#sigma)");
    h1m_sig1132   .SetXTitle("DCA_{ee}^{3D} (#sigma)");
    h1m_sig1127_np.SetXTitle("DCA_{ee}^{3D} (#sigma)");
    h1m_sig1127   .SetYTitle("#frac{1}{#it{N}_{ee}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    h1m_sig2732   .SetYTitle("#frac{1}{#it{N}_{ee}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    h1m_sig1132   .SetYTitle("#frac{1}{#it{N}_{ee}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    h1m_sig1127_np.SetYTitle("#frac{1}{#it{N}_{ee}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");



    outfile = TFile("20241016_dca_template_gammagamma2ee_7090.root", "RECREATE");
    outfile.WriteTObject(h1m_sig1127   );
    outfile.WriteTObject(h1m_sig2732   );
    outfile.WriteTObject(h1m_sig1132   );
    outfile.WriteTObject(h1m_sig1127_np);
    outfile.Clone();

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.15, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    y_range_max = h1m_sig1127.GetMaximum() * 2.0;
    y_range_min = -0.0;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ee}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_sig1127.Draw("E0same");
    h1m_sig2732.Draw("E0same");
    h1m_sig1127_np.Draw("E0same");

    leg = TLegend(0.15, 0.65, 0.35, 0.8);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.AddEntry(h1m_sig1127   , "{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}, #it{{p}}_{{T,ee}} < {2:2.1f} GeV/#it{{c}}, #gamma#gamma #rightarrow  e^{{+}}e^{{#minus}}".format(1.1, 2.7, 0.1), "P");
    leg.AddEntry(h1m_sig2732   , "{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}, #it{{p}}_{{T,ee}} < {2:2.1f} GeV/#it{{c}}, #gamma#gamma #rightarrow J/#psi #rightarrow  e^{{+}}e^{{#minus}}".format(2.7, 3.2, 0.1), "P");
    leg.AddEntry(h1m_sig1127_np, "{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}, #it{{p}}_{{T,ee}} > {2:2.1f} GeV/#it{{c}}, c#bar{{c}} #rightarrow e^{{+}}e^{{#minus}}".format(1.1, 2.7, 0.1), "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #h1ratio = h1m_sig1127.Clone("h1ratio");
    #h1ratio.Divide(h1m_sig1127, h1m_sig2732, 1., 1., "G");
    #h1ratio.Draw("E0");
    #h1ratio.Fit("pol0", "SME", "", 0, 2);
    #h1ratio.SetDirectory(0);
    #ROOT.SetOwnership(h1ratio, False);

    txt = TPaveText(0.2, 0.8, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    #txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    #txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dnddca_test_pt{1:2.1f}_{2:2.1f}GeV{3}.eps".format(date, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_test_pt{1:2.1f}_{2:2.1f}GeV{3}.pdf".format(date, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_test_pt{1:2.1f}_{2:2.1f}GeV{3}.png".format(date, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dcaee_uls_ls(filename, taskname, cen1, cen2, arr_dcaee, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
    hs_uls.GetAxis(0).SetRange(bin0, bin1);
    hs_bkg.GetAxis(0).SetRange(bin0, bin1);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(2);
    h1m_bkg_org = hs_bkg.Projection(2);
    h1m_sig_org = hs_sig.Projection(2);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_dcaee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_dcaee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_dcaee, False, False);

    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.17, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 50;
    #y_range_max = h1m_uls.GetMaximum() * 200;
    y_range_min = max(h1m_sig.GetMinimum() * 0.2, 2e-9);
    #y_range_max = 2e+0;
    #y_range_min = 2e-7;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.75, 0.85, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "#it{N}_{+#minus}", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = #it{N}_{+#minus} #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dnddca_uls_ls_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_uls_ls_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_uls_ls_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_ptee_uls_ls_upc(filename, taskname, cen1, cen2, mmin, mmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
    hs_uls.GetAxis(0).SetRange(bin0, bin1);
    hs_bkg.GetAxis(0).SetRange(bin0, bin1);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(1);
    h1m_bkg = hs_bkg.Projection(1);
    h1m_sig = hs_sig.Projection(1);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.05);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 1.5;
    #y_range_max = 15e-6;
    y_range_min = 0;

    frame1 = c1.DrawFrame(0., y_range_min, 1, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,ee}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetMaxDigits(2);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.7, 0.85, 0.85);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = ULS #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.65, 0.5, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndpt_lowpt_uls_ls_m{1:3.2f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_lowpt_uls_ls_m{1:3.2f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_lowpt_uls_ls_m{1:3.2f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_ptee_uls_ls(filename, taskname, cen1, cen2, arr_pt, mmin, mmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
    hs_uls.GetAxis(0).SetRange(bin0, bin1);
    hs_bkg.GetAxis(0).SetRange(bin0, bin1);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(1);
    h1m_bkg_org = hs_bkg.Projection(1);
    h1m_sig_org = hs_sig.Projection(1);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_pt, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_pt, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_pt, False, False);

    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.17, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 200;
    y_range_min = max(h1m_sig.GetMinimum() * 0.2, 2e-9);
    y_range_max = 2e+0;
    y_range_min = 2e-9;

    frame1 = c1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,ee}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.75, 0.85, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "#it{N}_{+#minus}", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = #it{N}_{+#minus} #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndpt_uls_ls_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_uls_ls_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_uls_ls_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_uls_ls_1100(filename, taskname, cen1, cen2, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    h1m_uls.GetXaxis().SetRangeUser(0, 1.1);
    h1m_sig.GetXaxis().SetRangeUser(0, 1.1);
    y_range_max = h1m_uls.GetMaximum() * 200;
    y_range_min = max(h1m_sig.GetMinimum() * 0.2, 2e-8);
    h1m_uls.GetXaxis().SetRange(0, 0);
    h1m_sig.GetXaxis().SetRange(0, 0);

    frame1 = c1.DrawFrame(0., y_range_min, 1.1, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.75, 0.85, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = ULS #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("centrality FT0C : 30#minus50 %");
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma_upto1100MeV_{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma_upto1100MeV_{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma_upto1100MeV_{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_uls_ls(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get(taskname);
    rootdire.ls();

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 200;
    y_range_min = max(h1m_sig.GetMinimum() * 0.2, 2e-7);
    if "pp" in suffix:
        y_range_max = 2e+0; #for pp
        y_range_min = 2e-7; #for pp
    elif "PbPb" in suffix:
        y_range_max = 2e+1; #for PbPb
        y_range_min = 2e-6; #for PbPb

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_bkg.Draw("E0same");
    h1m_uls.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.75, 0.85, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "#it{N}_{+#minus}", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = #it{N}_{+#minus} #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_signal_fraction_cumulative(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1m_uls_cumulative = get_cumulative_histogram(h1m_uls, False, False);
    h1m_bkg_cumulative = get_cumulative_histogram(h1m_bkg, False, False);
    h1m_sig_cumulative = get_cumulative_histogram(h1m_sig, False, False);

    h1sb = get_ratio(h1m_sig_cumulative, h1m_uls_cumulative, "B");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 1.05;
    y_range_min = 0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("integral range from 0 in #it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("cumulative S/(S+B)");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_sig_frac_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sig_frac_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sig_frac_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_sbratio_cumulative(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1m_bkg_cumulative = get_cumulative_histogram(h1m_bkg, False, False);
    h1m_sig_cumulative = get_cumulative_histogram(h1m_sig, False, False);

    h1sb = get_ratio(h1m_sig_cumulative, h1m_bkg_cumulative, "G");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 2e+4;
    y_range_min = 2e-2;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("integral range from 0 in #it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("cumulative S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_sbratio_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_cumulative_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_signal_fraction(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1sb = get_ratio(h1m_sig, h1m_uls, "B");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 1.05;
    y_range_min = 0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/(S+B)");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_sig_frac_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sig_frac_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sig_frac_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_sbratio(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);

    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 2e+4;
    y_range_min = 2e-5;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_significance(filename, taskname, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    h1z = rootdire.Get("hZvtx");
    nev = h1z.GetEntries();

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls_org = hs_uls.Projection(0);
    h1m_bkg_org = hs_bkg.Projection(0);
    h1m_sig_org = hs_sig.Projection(0);

    h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
    h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
    h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

    h1m_bkg.Scale(nev);
    h1m_sig.Scale(nev);

    h1significance = get_significance(h1m_sig, h1m_bkg);
    make_common_style(h1significance, 20, 1.2, kRed+1, 2, 0);
    h1significance.SetDirectory(0);
    ROOT.SetOwnership(h1significance, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.12, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 2e+3;
    y_range_min = 1e+0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("significance = S/#sqrt{S + 2B}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1significance.Draw("E0same");

    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE WIP");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_dcaee_sbratio_multiple(filename, tasknames, cen1, cen2, arr_dcaee, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kBlack, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e+4;
    y_range_min = 2e-5;
    if "pp" in suffix:
        y_range_max = 0.5;
        y_range_min = 0.001;
    elif "PbPb" in suffix:
        y_range_max = 2e+4;
        y_range_min = 2e-5;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.12, 0.02, 0., 0.02);
    p1.SetTicks(1,1);
    #p1.SetLogy(1);
    frame1 = p1.DrawFrame(arr_dcaee[0], y_range_min, arr_dcaee[-1], y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.2);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    ROOT.SetOwnership(p1, False);

    #leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg = TLegend(0.6, 0.85, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
        hs_uls.GetAxis(0).SetRange(bin0, bin1);
        hs_bkg.GetAxis(0).SetRange(bin0, bin1);
        hs_sig.GetAxis(0).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        h1m_uls_org = hs_uls.Projection(2);
        h1m_bkg_org = hs_bkg.Projection(2);
        h1m_sig_org = hs_sig.Projection(2);

        h1m_uls = rebin_histogram(h1m_uls_org, arr_dcaee, False, False);
        h1m_bkg = rebin_histogram(h1m_bkg_org, arr_dcaee, False, False);
        h1m_sig = rebin_histogram(h1m_sig_org, arr_dcaee, False, False);

        h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
        h1sb.SetDirectory(0);
        ROOT.SetOwnership(h1sb, False);
        make_common_style(h1sb, 20, 1.2, colors[it], 2, 0);
        h1sb.Draw("E0same");
        list_h1.append(h1sb);
        y_range_min = h1sb.GetMinimum();
        y_range_max = h1sb.GetMaximum();

        taskname_tmp = taskname.replace("dielectron_occupancy0_10000_TPChadrejorTOFreq_", "")
        if it == 0:
            taskname_tmp = "with TTCA";
        if it == 1:
            taskname_tmp = "without TTCA";
        if it == 0:
            taskname_tmp = "with TTCA, with weighting";
        leg.AddEntry(h1sb, taskname_tmp, "P");

    frame1.GetYaxis().SetRangeUser(y_range_min * 0.1, y_range_max * 2);
    txt = TPaveText(0.15, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}% ".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1., 0.35);
    p2.SetMargin(0.12, 0.02, 0.21, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(arr_dcaee[0], 0.9, arr_dcaee[-1], 1.55);
    frame2.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame2.GetYaxis().SetTitle("#frac{w}{wo}");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2, False);
    ROOT.SetOwnership(p2, False);
    line1 = TLine(arr_dcaee[0], 1, arr_dcaee[-1], 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1, False);

    for it in range(1, ntask):
        h1ratio = list_h1[it].Clone("h1ratio_{0}".format(it));
        h1ratio.Reset();
        h1ratio.Divide(list_h1[it], list_h1[0], 1., 1., "B");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);
        h1ratio.Draw("E0same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_dcaee_sbratio_multiple_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_dcaee_sbratio_multiple_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_dcaee_sbratio_multiple_mee{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_sbratio_multiple(filename, tasknames, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e+4;
    y_range_min = 2e-5;
    if "pp" in suffix:
        y_range_max = 2e+4;
        y_range_min = 2e-2;
    elif "PbPb" in suffix:
        y_range_max = 2e+4;
        y_range_min = 2e-5;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.12, 0.02, 0., 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.2);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    ROOT.SetOwnership(p1, False);

    #leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg = TLegend(0.6, 0.85, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls_org = hs_uls.Projection(0);
        h1m_bkg_org = hs_bkg.Projection(0);
        h1m_sig_org = hs_sig.Projection(0);

        h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
        h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
        h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

        h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
        h1sb.SetDirectory(0);
        ROOT.SetOwnership(h1sb, False);
        make_common_style(h1sb, 20, 1.2, colors[it], 2, 0);
        h1sb.Draw("E0same");
        list_h1.append(h1sb);

        taskname_tmp = taskname.replace("dielectron_occupancy0_10000_TPChadrejorTOFreq_", "")
        #if taskname == "dielectron-qc":
        #    taskname_tmp = "default";
        #if it == 0:
        #    taskname_tmp = "wo ITS cluster size cut";
        #elif it == 1:
        #    taskname_tmp = "w ITS cluster size cut";
        if it == 0:
            taskname_tmp = "without prefilter";
        if it == 1:
            taskname_tmp = "with #pi^{0} prefilter";
        leg.AddEntry(h1sb, taskname_tmp, "P");


    txt = TPaveText(0.15, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}% ".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1., 0.35);
    p2.SetMargin(0.12, 0.02, 0.21, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.9, 4, 1.55);
    frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("#frac{w}{wo}");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2, False);
    ROOT.SetOwnership(p2, False);
    line1 = TLine(0, 1, 4, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1, False);

    for it in range(1, ntask):
        h1ratio = list_h1[it].Clone("h1ratio_{0}".format(it));
        h1ratio.Reset();
        h1ratio.Divide(list_h1[it], list_h1[0], 1., 1., "B");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);
        h1ratio.Draw("E0same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_significance_multiple(filename, tasknames, cen1, cen2, arr_mee, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e+3;
    y_range_min = 8e-1;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.12, 0.02, 0., 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("significance = S/#sqrt{S + 2B}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.2);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1, False);
    ROOT.SetOwnership(p1, False);

    #leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg = TLegend(0.6, 0.8, 0.8, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");
        h1z = rootdire.Get("hZvtx");
        nev = h1z.GetEntries();

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls_org = hs_uls.Projection(0);
        h1m_bkg_org = hs_bkg.Projection(0);
        h1m_sig_org = hs_sig.Projection(0);

        h1m_uls = rebin_histogram(h1m_uls_org, arr_mee, False, False);
        h1m_bkg = rebin_histogram(h1m_bkg_org, arr_mee, False, False);
        h1m_sig = rebin_histogram(h1m_sig_org, arr_mee, False, False);

        h1m_bkg.Scale(nev);
        h1m_sig.Scale(nev);

        h1significance = get_significance(h1m_sig, h1m_bkg);
        h1significance.SetDirectory(0);
        ROOT.SetOwnership(h1significance, False);
        make_common_style(h1significance, 20, 1.2, colors[it], 2, 0);

        h1significance.Draw("E0same");
        list_h1.append(h1significance);

        #taskname_tmp = taskname.replace("dielectron-qc_", "")
        taskname_tmp = taskname.replace("dielectron_occupancy0_10000_TPChadrejorTOFreq_", "")
        #if taskname == "dielectron-qc":
        #    taskname_tmp = "default";
        #if it == 0:
        #    taskname_tmp = "wo ITS cluster size cut";
        #elif it == 1:
        #    taskname_tmp = "w ITS cluster size cut";
        leg.AddEntry(h1significance, taskname_tmp, "L");

    txt = TPaveText(0.15, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV, {0}#minus{1}% ".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1., 0.35);
    p2.SetMargin(0.12, 0.02, 0.2, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.8, 4, 1.55);
    frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("#frac{w}{wo}");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2, False);
    ROOT.SetOwnership(p2, False);
    line1 = TLine(0, 1, 4, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1, False);

    for it in range(1, ntask):
        h1ratio = list_h1[it].Clone("h1ratio_{0}".format(it));
        h1ratio.Reset();
        h1ratio.Divide(list_h1[it], list_h1[0], 1., 1., "B");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);
        h1ratio.Draw("E0same");


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_ttca_multiple_counting(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kGreen+2, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e-0;
    y_range_min = 7e-7;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.14, 0.02, 0.0, 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.55, 0.8, 0.75, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");
        list_h1.append(h1m_sig);

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 10;
            y_range_min = h1m_sig.GetMinimum() * 0.25;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");
        h1m_sig.SetDirectory(0);
        ROOT.SetOwnership(h1m_sig, False);

        taskname_tmp = taskname.replace("dielectron_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        if it == 0:
            taskname_tmp = "with TTCA, with weighting";
            make_common_style(h1m_sig, 20, 1.2, kBlue+1, 2, 0);
        h1m_sig.Draw("E0same");
        if it == 1:
            taskname_tmp = "with TTCA";
            make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    txt = TPaveText(0.18, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.14, 0.02, 0.2, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.99, 4, 1.022);
    frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("multiple counting");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2,False);
    line1 = TLine(0, 1, 4.0, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1,False);

    h1ratio = list_h1[1].Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Divide(list_h1[1], list_h1[0], 1., 1. ,"B");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_ttca_multiple_counting_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_ttca_multiple_counting_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_ttca_multiple_counting_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_ttca_recovery(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kGreen+2, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e-0;
    y_range_min = 7e-7;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.14, 0.02, 0.0, 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.55, 0.8, 0.75, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");
        list_h1.append(h1m_sig);

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 10;
            y_range_min = h1m_sig.GetMinimum() * 0.25;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");
        h1m_sig.SetDirectory(0);
        ROOT.SetOwnership(h1m_sig, False);

        taskname_tmp = taskname.replace("dielectron_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        if it == 0:
            taskname_tmp = "without TTCA";
            make_common_style(h1m_sig, 24, 1.2, kGreen+2, 2, 0);
        if it == 1:
            taskname_tmp = "with TTCA, with weighting";
            make_common_style(h1m_sig, 20, 1.2, kBlue+1, 2, 0);

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    txt = TPaveText(0.18, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.14, 0.02, 0.2, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.95, 4, 1.22);
    frame2.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame2.GetYaxis().SetTitle("signal recovery");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2,False);
    line1 = TLine(0, 1, 4.0, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1,False);

    h1ratio = list_h1[1].Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Divide(list_h1[1], list_h1[0], 1., 1. ,"B");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_ttca_recovery_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_ttca_recovery_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_ttca_recovery_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
#__________________________________________________________
def draw_dcaee_ttca_multiple_counting(filename, tasknames, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kGreen+2, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e-4;
    y_range_min = 7e-9;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.14, 0.02, 0.0, 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee} (#sigma)");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.55, 0.8, 0.75, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
        hs_uls.GetAxis(0).SetRange(bin0, bin1);
        hs_bkg.GetAxis(0).SetRange(bin0, bin1);
        hs_sig.GetAxis(0).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(2);
        h1m_bkg = hs_bkg.Projection(2);
        h1m_sig = hs_sig.Projection(2);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");
        list_h1.append(h1m_sig);

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 10;
            y_range_min = h1m_sig.GetMinimum() * 0.25;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");
        h1m_sig.SetDirectory(0);
        ROOT.SetOwnership(h1m_sig, False);

        taskname_tmp = taskname.replace("dielectron_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        if it == 1:
            taskname_tmp = "with TTCA";
            make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);
        if it == 0:
            taskname_tmp = "with TTCA, with weighting";
            make_common_style(h1m_sig, 20, 1.2, kBlue+1, 2, 0);

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    frame1.GetYaxis().SetRangeUser(y_range_min, y_range_max);
    txt = TPaveText(0.18, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:3.2f} < m_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.14, 0.02, 0.2, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.99, 10, 1.055);
    frame2.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame2.GetYaxis().SetTitle("multiple counting");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2,False);
    line1 = TLine(0, 1, 10.0, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1,False);

    h1ratio = list_h1[1].Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Divide(list_h1[1], list_h1[0], 1., 1. ,"B");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dnddca_ttca_multiple_counting_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_ttca_multiple_counting_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_ttca_multiple_counting_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_dcaee_ttca_recovery(filename, tasknames, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kGreen+2, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    y_range_max = 2e-4;
    y_range_min = 7e-9;

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.14, 0.02, 0.0, 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., y_range_min, 10, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee} (#sigma)");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.55, 0.8, 0.75, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs_uls.GetAxis(0).FindBin(mmax - 1e-3);
        hs_uls.GetAxis(0).SetRange(bin0, bin1);
        hs_bkg.GetAxis(0).SetRange(bin0, bin1);
        hs_sig.GetAxis(0).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(2);
        h1m_bkg = hs_bkg.Projection(2);
        h1m_sig = hs_sig.Projection(2);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");
        list_h1.append(h1m_sig);

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 10;
            y_range_min = h1m_sig.GetMinimum() * 0.25;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");
        h1m_sig.SetDirectory(0);
        ROOT.SetOwnership(h1m_sig, False);

        taskname_tmp = taskname.replace("dielectron_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        if it == 0:
            taskname_tmp = "without TTCA";
            make_common_style(h1m_sig, 24, 1.2, kGreen+2, 2, 0);
        if it == 1:
            taskname_tmp = "with TTCA, with weighting";
            make_common_style(h1m_sig, 20, 1.2, kBlue+1, 2, 0);

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    frame1.GetYaxis().SetRangeUser(y_range_min, y_range_max);
    txt = TPaveText(0.18, 0.6, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:3.2f} < m_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.14, 0.02, 0.2, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.95, 10, 2.7);
    frame2.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame2.GetYaxis().SetTitle("signal recovery");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2,False);
    line1 = TLine(0, 1, 10.0, 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1,False);

    h1ratio = list_h1[1].Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Divide(list_h1[1], list_h1[0], 1., 1. ,"B");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dnddca_ttca_recovery_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_ttca_recovery_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_ttca_recovery_m{1:3.2f}_{2:3.2f}_pt{3:2.1f}_{4:2.1f}{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_multiple(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    #colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    colors = [kGreen+2, kRed+1, kBlue+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = 0;
    #y_range_min = 0;
    y_range_max = 2e-2;
    y_range_min = 2e-8;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    #leg = TLegend(0.6, 0.85, 0.8, 0.95);
    leg = TLegend(0.2, 0.6, 0.4, 0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg,False);

    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 10;
            y_range_min = h1m_sig.GetMinimum() * 0.25;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");
        h1m_sig.SetDirectory(0);
        ROOT.SetOwnership(h1m_sig, False);

        taskname_tmp = taskname.replace("dielectron_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        #if it == 0:
        #    taskname_tmp = "without prefilter";
        #if it == 1:
        #    taskname_tmp = "with #pi^{0} prefilter";
        if it == 0:
            taskname_tmp = "with TTCA";
        if it == 1:
            taskname_tmp = "without TTCA";
            make_common_style(h1m_sig, 24, 1.2, colors[it], 2, 0);
        if it == 2:
            taskname_tmp = "with TTCA, with weighting";

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    #y_range_max = 2e-2;
    #y_range_min = 2e-8;
    frame1.GetYaxis().SetRangeUser(max(y_range_min, 2e-8), y_range_max);
    txt = TPaveText(0.2, 0.7, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_ptee_upc(filename, taskname, cen1, cen2, mmin, mmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_sig.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(mmax - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_sig.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_sig.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);
    h1m_sig = hs_sig.Projection(1);
    h1m_sig.Scale(1, "width");

    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.18, 0.02, 0.1, 0.05);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    y_range_max = h1m_sig.GetMaximum() * 5;
    y_range_min = 0.0;

    frame1 = c1.DrawFrame(0., y_range_min, 1, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,ee}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.1);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h1m_sig.Draw("E0same");

    txt = TPaveText(0.22, 0.62, 0.5, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE WIP");
    txt.AddText("{0}#minus{1}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_dcaee(filename, taskname, cen1, cen2, mmin, mmax, ptmin, ptmax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_sig.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(mmax - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_sig.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_sig.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);
    h1m_sig = hs_sig.Projection(2);
    h1m_sig.Scale(1, "width");

    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.18, 0.02, 0.1, 0.05);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    y_range_max = h1m_sig.GetMaximum() * 5;
    y_range_min = 0.0;

    frame1 = c1.DrawFrame(0., y_range_min, 5, y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{ee}^{3D} (#sigma)");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{dDCA_{ee}^{3D}} (#sigma)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.1);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h1m_sig.Draw("E0same");

    txt = TPaveText(0.22, 0.62, 0.5, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE WIP");
    txt.AddText("{0}#minus{1}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    filename = "mee_ptee_dcaee_pp_13.6TeV_LHC22o_pass6.root";
    taskname = "dielectron-qc";
    ptmin = 0;
    ptmax = 10;
    dcamin = 0;
    dcamax = 10;
    suffix = "_pp_13.6TeV_LHC22o_pass6";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_224982_0090_occupancy.root";
    #tasknames = [
    #    "dielectron-qc",
    #    "dielectron-qc_occupancy0_5000",
    #    "dielectron-qc_occupancy5001_99999",
    #];
    #draw_mee_sbratio_multiple(filename, tasknames, 0.0, 10.0, 0.0, 100.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0.0, 10.0, 0.0, 100.0, suffix);
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 100.0, suffix);

    ## for peripheral PbPb
    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_5070.root";
    #taskname = "dielectron-qc_5070_TOFreq_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_TOFreq_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0,  0.5, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0,  0.5, suffix);

    filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_occupancy.root";
    taskname = "dielectron-qc_0090_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_0090_occupancy0_1000_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy0_1000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_0090_occupancy0_5000_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy0_5000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_0090_occupancy5000_99999_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy5000_99999_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_occupancy_1bigbin.root";
    #taskname = "dielectron-qc_0090_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_newsel8";
    #draw_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_0090_occupancy0_1000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy0_1000_newsel8";
    #draw_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_0090_occupancy0_5000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy0_5000_newsel8";
    #draw_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_0090_occupancy5000_99999_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_0090_occupancy5000_99999_newsel8";
    #draw_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);

    taskname = "dielectron-qc_5090_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5090_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_5090_occupancy0_1000_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5090_occupancy0_1000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_5090_occupancy0_5000_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5090_occupancy0_5000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 1.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 1.1, 2.5, 0.0, 1.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    taskname = "dielectron-qc_5090_occupancy5000_99999_newsel8";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5090_occupancy5000_99999_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 90, 2.7, 3.2, 0.0, 10.0, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_occupancy.root";
    #taskname = "dielectron-qc_5070_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_5070_occupancy0_1000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_occupancy0_1000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_5070_occupancy0_5000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_occupancy0_5000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_5070_occupancy5000_99999_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_occupancy5000_99999_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_occupancy.root";
    #taskname = "dielectron-qc_7090_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_7090_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_7090_occupancy0_1000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_7090_occupancy0_1000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_7090_occupancy0_5000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_7090_occupancy0_5000_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron-qc_7090_occupancy5000_99999_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_7090_occupancy5000_99999_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10.0, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_occupancy.root";
    #taskname = "dielectron-qc_5090_occupancy0_5000_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5090_occupancy0_5000_newsel8";
    #draw_dcaee(filename, taskname, 50, 90, 1.1, 2.5, 0.0,  0.1, suffix);
    #draw_dcaee(filename, taskname, 50, 90, 1.1, 2.5, 0.1,  10, suffix);
    #draw_dcaee(filename, taskname, 50, 90, 2.7, 3.2, 0.0,  0.1, suffix);
    #draw_dcaee(filename, taskname, 50, 90, 2.7, 3.2, 0.1,  10, suffix);

    #filename = "mee_ptee_dcaee_pp_13.6TeV_LHC22o_pass6_HL_232106.root";
    ##taskname = "dimuon-qc_global";
    ##suffix = "_pp13.6TeV_LHC22o_pass6_global_muon";
    ##draw_mmumu_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon-qc_standalone";
    #suffix = "_pp13.6TeV_LHC22o_pass6_standalone_muon";
    #draw_mmumu_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);

    #filename = "mee_ptee_dcaee_pp_13.6TeV_LHC22o_pass6_HL_233027_big_ptmumu.root";
    #filename = "mee_ptee_dcaee_pp_13.6TeV_LHC22o_pass6_HL_233027.root";
    #taskname = "dimuon-qc_global_minpt200";
    #suffix = "_pp13.6TeV_LHC22o_pass6_global_muon_minpt200";
    #draw_mmumu_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 999, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 0.2, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 0.1, 0.0, 999.0, suffix);
    #taskname = "dimuon-qc_global_minpt400";
    #suffix = "_pp13.6TeV_LHC22o_pass6_global_muon_minpt400";
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon-qc_global_minpt500";
    #suffix = "_pp13.6TeV_LHC22o_pass6_global_muon_minpt500";
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);

    #taskname = "dimuon-qc_standalone_minpt200";
    #suffix = "_pp13.6TeV_LHC22o_pass6_standalone_muon_minpt200";
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_Rfactor(filename, taskname, 0.0, 10.0, 0.0, 999, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 0.2, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 0.1, 0.0, 999.0, suffix);
    #taskname = "dimuon-qc_standalone_minpt400";
    #suffix = "_pp13.6TeV_LHC22o_pass6_standalone_muon_minpt400";
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon-qc_standalone_minpt500";
    #suffix = "_pp13.6TeV_LHC22o_pass6_standalone_muon_minpt500";
    #draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0.0, 10.0, 0.0, 999.0, suffix);


#    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_HL_233073_big_ptmumu.root";
#    filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_HL_233073.root";
#    taskname = "dimuon-qc_0090_standalone_minpt200";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_standalone_muon_minpt200";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    taskname = "dimuon-qc_0090_standalone_minpt400";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_standalone_muon_minpt400";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    taskname = "dimuon-qc_0090_standalone_minpt500";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_standalone_muon_minpt500";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#
#    taskname = "dimuon-qc_0090_global_minpt200";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_global_muon_minpt200";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    taskname = "dimuon-qc_0090_global_minpt400";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_global_muon_minpt400";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    taskname = "dimuon-qc_0090_global_minpt500";
#    suffix = "_PbPb5.36TeV_LHC23zzh_pass3_0090_global_muon_minpt500";
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    draw_mmumu_uls_ls_1100MeV(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 0.0, 1.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 1.0, 2.0, 0.0, 999, suffix);
#    #draw_mmumu_Rfactor(filename, taskname, 0, 90, 2.0, 3.0, 0.0, 999, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_uls_ls(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 0.0, 10.0, 0.0, 999.0, suffix);
#    #draw_mmumu_sbratio(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);

    #filename = "dimuon_PbPb_5.36TeV_LHC23zzh_pass3_HL_243092_all.root";
    #taskname = "dimuon_standalone_0010";
    #suffix = "_PbPb5.36TeV_2023_pass3_0010_standalone_muon_minpt400";
    #draw_mmumu_uls_ls_pbpb_2panel(filename, taskname, 0, 10, 2.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_standalone_0090";
    #suffix = "_PbPb5.36TeV_2023_pass3_0090_standalone_muon_minpt400";
    #draw_mmumu_uls_ls_pbpb_2panel(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);

    #taskname = "dimuon_global_0010";
    #suffix = "_PbPb5.36TeV_2023_pass3_0010_global_muon_minpt400";
    #draw_mmumu_uls_ls_pbpb_2panel(filename, taskname, 0, 10, 2.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_global_0090";
    #suffix = "_PbPb5.36TeV_2023_pass3_0090_global_muon_minpt400";
    #draw_mmumu_uls_ls_pbpb_2panel(filename, taskname, 0, 90, 2.0, 10.0, 0.0, 999.0, suffix);


    #taskname = "dimuon_global_0010";
    #suffix = "_PbPb5.36TeV_2023_pass3_0010_global_muon_minpt400";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 2.0, 10.0, 0.0, 999.0, suffix);

    #filename = "dimuon_pp_13.6TeV_LHC22o_pass6_HL_242847.root";
    #taskname = "dimuon_standalone_id14533";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass6_standalone_muon_minpt200";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_global_id14533";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass6_global_muon_minpt200";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);

    filename = "dimuon_pp_13.6TeV_LHC22o_pass7_HL_249691.root";
    #taskname = "dimuon_standalone";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass7_standalone_muon_minpt200_HL_249691";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_global";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass7_global_muon_minpt200_HL_249691";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);

    #taskname = "dimuon_standalone_wWeighting";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass7_standalone_muon_minpt200_wWeighting_HL_249691";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_global_wWeighting";
    #suffix = "_pp13.6TeV_2022_LHC22o_pass7_global_muon_minpt200_wWeighting_HL_249691";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 999.0, suffix);


    #filename = "dimuon_PbPb_5.36TeV_LHC23zzh_pass4_HL_246171.root";
    #taskname = "dimuon_standalone_0010";
    #suffix = "_PbPb5.36TeV_2023_pass4_0010_standalone_muon_minpt400";
    #draw_mmumu_uls_ls(filename, taskname, 0, 10, 2.0, 10.0, 0.0, 999.0, suffix);
    #taskname = "dimuon_global_0010";
    #suffix = "_PbPb5.36TeV_2023_pass4_0010_global_muon_minpt400";
    #draw_mmumu_uls_ls(filename, taskname, 0, 10, 2.0, 10.0, 0.0, 999.0, suffix);

#    filename = "dielectron_pp_13.6TeV_LHC22o_pass7_HL_253469.root";
#    taskname = "dielectron_TPChadrejorTOFreq_minpt200";
#    suffix = "_kibanb";
#    draw_mee_dcaee(filename, taskname, 0.0, 10, suffix);

    #filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_258611.root";
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_3050_pT400MeV_tmp";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 1.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 2.0, 10.0, 0.0, 10.0, suffix);

    #draw_v2_vs_mee(filename, taskname, 2.0, 10.0, 0.0, 10.0, suffix);
    #draw_v2_vs_ptee(filename, taskname, 0.0, 0.14, 0.0, 10.0, suffix);
    #draw_v2_vs_ptee(filename, taskname, 2.8, 3.2, 0.0, 10.0, suffix);
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt800";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_3050_pT800MeV";
    #draw_v2_vs_ptee(filename, taskname, 2.8, 3.2, 0.0, 10.0, suffix);

    #taskname = "dielectron_v2_3050_TOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_3050_TOFreq_pT400MeV";
    #draw_v2_vs_mee(filename, taskname, 2.0, 10.0, 0.0, 10.0, suffix);
    #draw_v2_vs_ptee(filename, taskname, 0.0, 0.14, 0.0, 10.0, suffix);
    #draw_v2_vs_ptee(filename, taskname, 2.8, 3.2, 0.0, 10.0, suffix);

#    filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_261772.root";
#    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
#    suffix = "_PbPb_5.36TeV_LHC23zzh_pass4_3050_pT400MeV_HL_261772";
#    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
#    #draw_mee_uls_ls(filename, taskname, 1.0, 10.0, 0.0, 10.0, suffix);
#    #draw_mee_uls_ls(filename, taskname, 2.0, 10.0, 0.0, 10.0, suffix);
#    draw_v2_vs_mee( filename, taskname, 2.0, 10.0, 0.0, 10.0, suffix);
#    draw_v2_vs_ptee(filename, taskname, 0.0, 0.14, 0.0, 10.0, suffix);
#    draw_v2_vs_ptee(filename, taskname, 2.8,  3.2, 0.0, 10.0, suffix);

    #arr_m = np.array([0, 0.04, 0.14, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2, 2.4, 2.7, 3.0, 3.2, 4], dtype=float);
    arr_m = np.array([0, 0.14, 0.5, 1.1, 2.0, 2.7, 3.2, 4], dtype=float);
    filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_274632.root";
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23zzh_pass4_3050_pT400MeV_HL_274632";
    #draw_mee_uls_ls(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 30, 50, 0.0, 10.0, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 30, 50, 4.0, 6.0, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 30, 50, 6.0, 10.0, 2.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.0, 4.0, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.0,0.14, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.5, 1.1, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 1.1, 2.7, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_mee_dcaee(filename, taskname, 30, 50, 0.0, 10, suffix);

    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 10, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 10, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 10, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 2.0, 10, suffix);
    arr_pt = np.array([0, 1., 2, 3, 4, 6, 10], dtype=float);
    #arr_pt_pion = np.array([1., 2, 3, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 10, suffix, True, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.04, 0.0, 10, suffix, True, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 0.0, 10, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 2.0, 10, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 10, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 2.0, 10, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 10, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 2.0, 10, suffix);

    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 10, suffix, False, True);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 1 , suffix, False, True);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 2.0, 10, suffix, False, True);
    ##arr_dca = np.array([0, 0.5, 1, 1.5, 2, 4, 6, 10], dtype=float);
    #arr_dca = np.array([0, 1, 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.04, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.14, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.5,  1.1, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.5,  1.1, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 0.0, 10, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_277601.root";
    taskname = "dielectron_7090";
    suffix = "_PbPb_5.36TeV_LHC23pass4_7090_pT400MeV_HL_277601";
    #draw_mee_uls_ls(filename, taskname, 70, 90, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 70, 90, 0.0, 0.1, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls_upc(filename, taskname, 70, 90, 1.1, 3.2, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls_upc(filename, taskname, 70, 90, 0.8, 1.1, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls_upc(filename, taskname, 70, 90, 1.1, 2.7, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls_upc(filename, taskname, 70, 90, 3.2, 4.0, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 0.0, 0.14, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 1.1, 2.7, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 0.8, 1.1, 0.0, 0.1, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 1.1, 2.7, 0.0, 0.1, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 0.1, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 0.8, 1.1, 0.1, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 1.1, 2.7, 0.1, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 2.7, 3.2, 0.1, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, 1.1, 3.2, 0.0, 0.1, suffix);
    #draw_dcaee_test(filename, taskname, 70, 90, 0.0, 0.1, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_279716.root";
    taskname = "dielectron_3050";
    suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_HL_279716";

    tasknames = [
        # "dielectron_0010_TOFreq_ITSreq_minpt400"
        #,"dielectron_0010_TOFreq_minpt400"
        "dielectron_0010_TPChadrejorTOFreq_minpt400"
        ,"dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400"
    ];
    #draw_mee_uls_ls(filename, tasknames[0], 0, 10, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls(filename, tasknames[1], 0, 10, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_wITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[0], 0, 10, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[1], 0, 10, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_wITS_HL_279716");
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_HL_279716";
    tasknames = [
        # "dielectron_3050_TOFreq_ITSreq_minpt400"
        #,"dielectron_3050_TOFreq_minpt400"
        "dielectron_3050_TPChadrejorTOFreq_minpt400"
        ,"dielectron_3050_TPChadrejorTOFreq_ITSreq_minpt400"
    ];
    #draw_mee_uls_ls(filename, tasknames[0], 30, 50, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls(filename, tasknames[1], 30, 50, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_wITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[0], 30, 50, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[1], 30, 50, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_wITS_HL_279716");
    #draw_mee_sbratio_multiple(filename, tasknames, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);

    suffix = "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_HL_279716";
    tasknames = [
        # "dielectron_7090_TOFreq_ITSreq_minpt400"
        #,"dielectron_7090_TOFreq_minpt400"
        "dielectron_7090_TPChadrejorTOFreq_minpt400"
        ,"dielectron_7090_TPChadrejorTOFreq_ITSreq_minpt400"
    ];
    #draw_mee_uls_ls(filename, tasknames[0], 70, 90, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls(filename, tasknames[1], 70, 90, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_wITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[0], 70, 90, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_woITS_HL_279716");
    #draw_mee_uls_ls_1100(filename, tasknames[1], 70, 90, 0.0, 10.0, 0.0, 20.0, "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_wITS_HL_279716");
    #draw_mee_sbratio_multiple(filename, tasknames, 70, 90, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 70, 90, 0.0, 10.0, 0.0, 20.0, suffix);

    
    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_275984.root";
    #taskname = "dielectron_occupancy0_1000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_occupancy0_1000_HL_275984";
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.2, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_275984.root";
    #taskname = "dielectron_occupancy0_1000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_occupancy0_1000_HL_275984";
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.2, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_280731.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_HL_280731";
    ##draw_mee_uls_ls(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.2, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_HL_280731_tpcshared04";
    ##draw_mee_uls_ls(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.2, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_pT400MeV_HL_280731";
    ##draw_mee_uls_ls(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.2, 0.0, 20.0, suffix);
    #taskname = "dielectron_1030_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_pT400MeV_HL_280731_tpcshared04";
    ##draw_mee_uls_ls(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.2, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_HL_280731";
    ##draw_mee_uls_ls(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.2, 0.0, 20.0, suffix);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_HL_280731_tpcshared04";
    ##draw_mee_uls_ls(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.2, 0.0, 20.0, suffix);

    #taskname = "dielectron_5070_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_5070_pT400MeV_HL_280731";
    ##draw_mee_uls_ls(filename, taskname, 50, 70, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 50, 70, 0.1, 0.2, 0.0, 20.0, suffix);
    #taskname = "dielectron_5070_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_5070_pT400MeV_HL_280731_tpcshared04";
    ##draw_mee_uls_ls(filename, taskname, 50, 70, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 50, 70, 0.1, 0.2, 0.0, 20.0, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_ITSreq_minpt400_id18445";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_HL_280731";
    ##draw_mee_uls_ls(filename, taskname, 70, 90, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 70, 90, 0.1, 0.2, 0.0, 20.0, suffix);
    #taskname = "dielectron_7090_TPChadrejorTOFreq_ITSreq_minpt400_id18524";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_7090_pT400MeV_HL_280731_tpcshared04";
    ##draw_mee_uls_ls(filename, taskname, 70, 90, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 70, 90, 0.1, 0.2, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_281428.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_pT400MeV_HL_281428";
    ##draw_mee_uls_ls(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    ##draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.2, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_pT400MeV_HL_281428";
    ##draw_mee_uls_ls( filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    ##draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.2, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.2, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);



    #filename = "dielectron_pp_13.6TeV_LHC22o_pass7_HL_280891.root";
    #taskname = "dielectron_TPChadrejorTOFreq_minpt200";
    #suffix = "_pp_13.6TeV_LHC23_apass4_pT200MeV_HL_280891";
    #draw_mee_uls_ls( filename, taskname, 0, 999, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 20.0, suffix);



    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_282384.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_ITSreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_ITSreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);



    #taskname = "dielectron_1030_TPChadrejorTOFreq_ITSreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_TPChadrejorTOFreq_ITSreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_TPChadrejorTOFreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_ITSreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_TPChadrejorTOFreq_ITSreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_1030_TPChadrejorTOFreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 10, 30, 0.0, 10.0, 0.0, 20.0, suffix);



    #taskname = "dielectron_3050_TPChadrejorTOFreq_ITSreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_TPChadrejorTOFreq_ITSreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_TPCshared03_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_TPChadrejorTOFreq_TPCshared03_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_ITSreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_TPChadrejorTOFreq_ITSreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_3050_TPChadrejorTOFreq_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_282709.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared02_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_TPCshared02_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_ITSreq_TPCshared02_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_ITSreq_TPCshared02_pT400MeV_HL_281428";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_281955.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";# without shared03, without deta-dphi
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_pT400MeV_HL_281955";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_TPCshared03_minpt400"; # with shared03, without deta-dphi
    #suffix = "_PbPb_5.36TeV_LHC23_apass4_0010_TPChadrejorTOFreq_TPCshared03_pT400MeV_HL_281955";
    #draw_mee_uls_ls( filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, 0.1, 0.14, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_284153.root";
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_TPCshared0_pT400MeV_HL_284153";

    ##arr_m = np.array([0, 0.14, 0.5, 1.1, 2.0, 2.7, 3.2, 4], dtype=float);
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 2.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 2.0, 20, suffix);
    #arr_pt = np.array([0, 1., 2, 3, 4, 6, 10], dtype=float);
    #arr_pt_pion = np.array([1., 2, 3, 4, 6, 8, 10], dtype=float);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 0.5, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 2.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 2.0, 20, suffix);

    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 20, suffix, False, True);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 1 , suffix, False, True);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 2.0, 20, suffix, False, True);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 4, 6, 10], dtype=float);
    #arr_dca = np.array([0, 1, 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.04, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.14, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.5,  1.1, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.5,  1.1, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 0.0, 10, suffix);

    filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_284904.root";
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_TPCshared0_pT400MeV_HL_284904";

    ##arr_m = np.array([0, 0.14, 0.5, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 2.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 2.0, 20, suffix);

    #arr_pt = np.array([0, 1., 2, 3, 4, 6, 10], dtype=float);
    #arr_pt_pion = np.array([1., 2, 3, 4, 6, 10], dtype=float);
    #arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 1.1, 2.7, 2.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 20, suffix, False, True);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 0.0, 1 , suffix, False, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt, 2.7, 3.2, 2.0, 20, suffix, False, False);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.14, 2.7, 0.0, 20, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.14, 2.7, 0.0, 1, suffix);
    #draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.14, 2.7, 2.0, 20, suffix);

    ##arr_dca = np.array([0, 0.5, 1, 1.5, 2, 4, 6, 10], dtype=float);
    #arr_dca = np.array([0, 1, 2, 6, 10], dtype=float);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.14, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.04, 0.14, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.14,  1.1, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.14,  1.1, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 1.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 1.0, 10, suffix);


    #filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_284904_narrow.root";
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_TPCshared0_pT400MeV_HL_284904_narrow";
    #arr_m = np.array([0, 0.14, 0.5, 1.1, 2.7, 3.2, 4], dtype=float);
    ##arr_m = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);


    #filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_284904_ptintegrated.root";
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_TPCshared0_pT400MeV_HL_284904_ptintegrated";
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_dca = np.array([0, 1, 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 0.0, 10, suffix);

    filename = "dielectron_flow_PbPb_5.36TeV_LHC23zzh_pass4_HL_285692.root";
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_TPCshared0_pT400MeV_HL_285692";

    filename = "AnalysisResults_HL_286293.root";
    ##arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_m_hf = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    #arr_m_hf2 = np.array([0, 0.14, 4], dtype=float);
    ##arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);

    arr_pt_pion = np.array([1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10], dtype=float);
    arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 1, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.0, 0.14, 2.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.04, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.04, 0.14, 0.0, 1, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.04, 0.14, 2.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.0, 0.04, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.04, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.1, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_jpsi, 2.7, 3.2, 0.0, 20, suffix, False, True);

    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new_multiple_dca(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, suffix);

    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 4, 6, 10], dtype=float);
    #arr_dca = np.array([0, 1, 2, 6, 10], dtype=float);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 0.14, 1.1, 0.0, 10, suffix);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 0.14, 1.1, 1.0, 10, suffix);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 1.1, 2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 1.1, 2.7, 1.0, 10, suffix);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 0.14, 2.7, 0.0, 10, suffix);
    #draw_v2_vs_dcaee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_dca, 0.14, 2.7, 1.0, 10, suffix);


#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 20, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 1, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 2.0, 20, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 20, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 0.0, 1, suffix);
#    draw_v2_vs_mee(filename, taskname, 30, 50, arr_m, 2.0, 10.0, 2.0, 20, suffix);
#    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m_hf, 0.0, 10.0, 2.0, 20, suffix + "_HF");
#    #draw_v2_vs_mee(filename, taskname, 30, 50, arr_m_hf2, 0.0, 10.0, 2.0, 20, suffix + "_HF2");
#
#    arr_pt = np.array([0, 1., 2, 3, 4, 6, 10], dtype=float);
#    arr_pt_pion = np.array([1., 2, 3, 4, 6, 10], dtype=float);
#    arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
#    arr_pt_imr = np.array([0, 1., 2, 4, 6, 10], dtype=float);
#    arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_imr, 1.1, 2.7, 0.0, 20, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_imr, 1.1, 2.7, 0.0, 1, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_imr, 1.1, 2.7, 2.0, 20, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_jpsi, 2.7, 3.2, 0.0, 20, suffix, False, True);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_jpsi, 2.7, 3.2, 0.0, 1 , suffix, False, False);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_jpsi, 2.7, 3.2, 2.0, 20, suffix, False, False);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 2.7, 0.0, 20, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 2.7, 0.0, 1, suffix);
#    draw_v2_vs_ptee(filename, taskname, 30, 50, arr_pt_vm, 0.14, 2.7, 2.0, 20, suffix);
#
#    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 4, 6, 10], dtype=float);
#    arr_dca = np.array([0, 1, 2, 6, 10], dtype=float);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.0, 0.14, 1.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.04, 0.14, 1.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.14,  1.1, 0.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 0.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 0.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 0.14,  1.1, 1.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 1.1,  2.7, 1.0, 10, suffix);
#    draw_v2_vs_dcaee(filename, taskname, 30, 50, arr_dca, 2.7,  3.2, 1.0, 10, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_286293_yield.root";
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_286293";
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_uls_ls( filename, taskname, 30, 50, 0.0, 10.0, 2.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.0, 0.14, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 0.14, 0.2, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, 2.7, 3.2, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.14, 0.5, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.14, 1.1, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 0.5, 1.1, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 1.1, 2.7, 0.0, 10, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, 2.7, 3.2, 0.0, 10, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_284904.root";
    taskname = "dielectron_7090_TPChadrejorTOFreq_TPCshared03_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_TPChadrejorTOFreq_TPCshared03_pT400MeV_HL_284904";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.5, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_uls_ls( filename, taskname, 70, 90, arr_m, 0.0, 0.1, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 70, 90, 1.1, 2.7, 0.0, 20.0, suffix);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 1.1, 2.7, 0.0, 0.1, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 2.7, 3.2, 0.0, 0.1, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 1.1, 2.7, 0.2, 10, suffix);

    #filename = "dielectron_pp_13.6TeV_LHC22o_pass7_HL_275995.root";
    #taskname = "dielectron_TPChadrejorTOFreq_minpt200";
    #suffix = "_pp_13.6TeV_LHC22_LHC22o_apass7_pT200MeV_HL_275995";
    #arr_m = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
    #0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8,
    #0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1,
    #1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    #2.8, 2.85, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_288465.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_occupancy0_10000_TPChadrejorTOFreq_TPCshared07_pT400MeV_HL_288465";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.5, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 00, 90, arr_pt, 0.06, 0.2, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 00, 90, arr_pt, 2.7, 3.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    ##draw_dcaee_uls_ls(filename, taskname, 00, 90, arr_dca, 0.06, 0.2, 2.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 00, 90, arr_dca, 2.7, 3.2, 0.0, 10.0, suffix);


    filename = "AnalysisResults_HL_286547.root";
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    arr_m_hf = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    arr_m_hf2 = np.array([0, 0.14, 4], dtype=float);
    arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_286547";
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
    ##arr_m_wide = np.array([0, 0.14, 1.1, 2.7, 3.2, 8, 12], dtype=float);
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m_wide, 0.0, 10.0, 0.0, 20, suffix);
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m_wide, 0.0, 10.0, 0.0, 1, suffix);
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_m_wide, 0.0, 10.0, 2.0, 20, suffix);

    #arr_pt_pion = np.array([0.8, 1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10], dtype=float);
    #arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    ##draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.5, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt400", 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT200MeV_HL_286547";
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    ##draw_v2_vs_mee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);
    #arr_pt_pion = np.array([0.4, 0.6, 0.8, 1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10], dtype=float);
    #arr_pt_vm = np.array([0.4, 1., 2, 4, 6, 10], dtype=float);
    ##draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.5, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, "dielectron_v2_3050_TPChadrejorTOFreq_minpt200", 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_289636.root";
    tasknames = [
        "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400"
       ,"dielectron_occupancy0_10000_TPChadrejorTOFreq_chi2TOF1_minpt400"
       #,"dielectron_occupancy0_10000_TPChadrejorTOFreq_reldiffPin02_minpt400"
       #,"dielectron_occupancy0_10000_TPChadrejorTOFreq_chi2TOF1_reldiffPin02_minpt400"
    ];
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_HL_286547";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.5, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 100.0, suffix);
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_reldiffPin02_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_reldiffPin02_pT400MeV_HL_289636";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_290986.root";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);

    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.14, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.14, 0.4, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.08, 0.14, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.08, 0.14, 3.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.0, 0.14, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.0, 0.14, 3.0, 10.0, suffix);

    #taskname = "dielectron_occupancy10000_20000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy10000_20000_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_occupancy20000_30000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy20000_30000_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_occupancy30000_40000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy30000_40000_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_occupancy40000_50000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy40000_50000_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_occupancy50000_999999_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy50000_999999_TPChadrejorTOFreq_pT400MeV_HL_290986";
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_wophiv";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wophiv_HL_290986";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.14, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.14, 0.4, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.08, 0.14, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.08, 0.14, 3.0, 10.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_294001.root";
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_TPChadrejorTOFreq_pT400MeV_wophiv_HL_294001";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, arr_pt, 0.10, 0.14, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, arr_pt, 0.14, 0.4, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 00, 10, arr_dca, 0.1, 0.2, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 00, 10, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_294420.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_TPChadrejorTOFreq_pT400MeV_HL_294420";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_294526.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_HL_294526";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_295259.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta005";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_deta005_HL_295259";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_295259.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_deta010_HL_295259";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400_deta010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_deta010_HL_295259";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 40, 50, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400_deta005";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_deta005_HL_295259";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_295155.root";
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta005";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_deta005_HL_295155";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_deta010";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_TPChadrejorTOFreq_pT400MeV_deta010_HL_295155";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_296716.root"; #FT0C occupancy
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_296716";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_296716_flow.root"; #FT0C occupancy
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_296716";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt200";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT200MeV_HL_296716";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);


    #filename = "AnalysisResults_HL_296716.root"; #FT0C occupancy
    #arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_m_hf = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    #arr_m_hf2 = np.array([0, 0.14, 4], dtype=float);
    #arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_296716";
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);

    #arr_pt_pion = np.array([0.8, 1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10], dtype=float);
    #arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);


    filename = "AnalysisResults_HL_296817.root"; #track occupancy
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_296817";
    arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    arr_m_hf = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    arr_m_hf2 = np.array([0, 0.14, 4], dtype=float);
    arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);

    arr_pt_pion = np.array([1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10], dtype=float);
    arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, True, False);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_299066.root";
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299066";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299066";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 10, 30, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 10, 30, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299066";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299066";
    #arr_m = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 70, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_299065.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299065";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 2.4, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 3.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_298365.root";
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299065";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename, taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_299220.root";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 30, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 30, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);

    #taskname = "dielectron_5070_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_5070_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 50, 70, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 50, 70, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 50, 70, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 70, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);

    #taskname = "dielectron_1090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #draw_mee_uls_ls(filename  , taskname, 10, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_299686.root";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 30, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 30, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_5070_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_5070_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 50, 70, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 50, 70, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 50, 70, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 70, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_1090_TPChadrejorTOFreq_minpt400_id14815";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_299686";
    #draw_mee_uls_ls(filename  , taskname, 10, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    filename = "AnalysisResults_HL_299220.root"; #ft0c occupancy
    taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_299220";
    #arr_m = np.array([0, 0.14, 0.5, 1.1, 2.0, 2.7, 3.2, 4], dtype=float);
    arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    arr_m_hf = np.array([0, 0.14, 2.7, 3.2, 4], dtype=float);
    arr_m_hf2 = np.array([0, 0.14, 4], dtype=float);
    arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 1, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 2.0, 20, suffix);

    ##arr_pt_pion = np.array([1., 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 10], dtype=float);
    #arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    ##arr_pt_vm = np.array([1., 4, 10], dtype=float);
    #arr_pt_imr = np.array([0, 1., 2, 4, 6, 10], dtype=float);
    ###draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_pion, 0.0, 0.14, 0.0, 20, suffix, False, False, True);
    ##draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix, False, False, False);
    ###draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 1.1, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_imr, 1.1, 2.7, 0.0, 20, suffix, False, False, False);
    ##draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_imr, 1.1, 2.7, 2.0, 20, suffix, False, False, False);
    ###draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_jpsi, 2.7, 3.2, 0.0, 20, suffix, False, True, False);
    ###draw_v2_vs_ptee_new_multiple_dca(filename, taskname, 30, 50, arr_pt_vm, 0.5, 1.1, suffix);

    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.5, 2.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 20, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 0.0, 1.0, suffix);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_vm, 0.14, 0.3, 2.0, 20, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_302828.root";
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_302828";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_302828";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    filename = "dielectron_pp_13.6TeV_LHC23_pass4_thin_HL_289080.root";
    taskname = "dielectron_TPChadrejorTOFreq_minpt200_dca3d";
    suffix = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_HL_289080";
    arr_m = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
    0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.8, 2.85, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0], dtype=float);
    #draw_mee_uls_ls( filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_uls_ls( filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_uls_ls( filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 0.5, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 0.5, suffix);
    #draw_mee_sbratio_cumulative(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction_cumulative(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction(filename, taskname, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_306762.root";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_306762";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0,10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_306762";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30,50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio_cumulative(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_signal_fraction_cumulative(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_306762";
    #draw_mee_uls_ls(filename  , taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 70, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);


    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_307487.root";
    arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_307487";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_307487";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    filename = "dielectron_pp_13.6TeV_LHC23_pass4_thin_HL_307922.root";
    tasknames = [
        "dielectron_TPChadrejorTOFreq_minpt200",
        "dielectron_TPChadrejorTOFreq_minpt200_wpf_mee",
    ];
    arr_m = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
    0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.91, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.08, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.8, 2.85, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.9, 4.0], dtype=float);
    suffix = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_HL_307922";
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 100.0, suffix);
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 100.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 100.0, suffix);

    #draw_mee_sbratio_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 1.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_308878.root";
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_308878";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt200";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT200MeV_HL_308878";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_308878";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_309527.root";
    arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 2.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.14, 0.5, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.5, 1.1, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 1.1, 2.7, 0.0, 10.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 10, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 10, 30, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 10, 30, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #taskname = "dielectron_7090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_7090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_309527";
    #draw_mee_uls_ls(filename  , taskname, 70, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 70, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 70, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_310124.root"; #deta > 0.08
    #arr_m = np.array([0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_310124";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_310124";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);


    filename = "dielectron_pp_13.6TeV_LHC23_pass4_thin_HL_312813.root";
    tasknames = [
        "dielectron_TPChadrejorTOFreq_minpt200",
        "dielectron_TPChadrejorTOFreq_minpt200_wpf_mee",
        #"dielectron_TPChadrejorTOFreq_minpt200_itsibany",
        #"dielectron_TPChadrejorTOFreq_minpt200_itsibany_wpf_mee",
    ];
    arr_m = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
    0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.91, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.08, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.8, 2.85, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 4.0], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    suffix = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_HL_312813";
    suffix0 = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_HL_312813";
    suffix1 = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_wpf_mee_HL_312813";
    #draw_mee_uls_ls(filename  , tasknames[0], 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix0);
    #draw_mee_uls_ls(filename  , tasknames[1], 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix1);
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_dcaee_sbratio_multiple(filename, tasknames, 0, 999, arr_dca, 0.14, 0.3, 0.0, 10.0, suffix);
    #draw_dcaee_sbratio_multiple(filename, tasknames, 0, 999, arr_dca, 0.14, 0.5, 0.0, 10.0, suffix);
    #draw_dcaee_sbratio_multiple(filename, tasknames, 0, 999, arr_dca, 0.5, 1.1, 0.0, 10.0, suffix);
    #draw_dcaee_sbratio_multiple(filename, tasknames, 0, 999, arr_dca, 1.1, 2.7, 0.0, 10.0, suffix);
    #draw_dcaee_sbratio_multiple(filename, tasknames, 0, 999, arr_dca, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_sbratio_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 2.0, 20.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0, 999, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_dcaee(filename, tasknames[0], 0, 90, arr_m, 0.0, 10, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_312076.root"; #without prefilter, without deta-dphi cut
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_312076";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400_reldiffPin015";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_reldiffPin015_HL_312076";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_313263.root"; #without prefilter, with deta > 0.04
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_313263";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_316686.root"; #with prefilter, with deta > 0.03
    arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_313263";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    ##draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    ##draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_316686";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    #draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_1030_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_1030_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_316686";
    #draw_mee_uls_ls(filename  , taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    ##draw_ptee_uls_ls(filename , taskname, 10, 30, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    ##draw_dcaee_uls_ls(filename, taskname, 10, 30, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 10, 30, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_316686";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    ##draw_ptee_uls_ls(filename , taskname, 30, 50, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    ##draw_dcaee_uls_ls(filename, taskname, 30, 50, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_319300.root"; #with prefilter, with deta > 0.03
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_319300";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 90, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    ##draw_dcaee_uls_ls(filename, taskname, 0, 90, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_319300";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_ptee_uls_ls(filename , taskname, 0, 10, arr_pt, 0.10, 0.2, 0.0, 20.0, suffix);
    ##draw_dcaee_uls_ls(filename, taskname, 0, 10, arr_dca, 0.1, 0.2, 2.0, 10.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_319559.root"; #without prefilter, without detadphi cut
    #arr_m = np.array([0.0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_HL_319559";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_HL_319559";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_TightEvSel_HL_319559";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_TightEvSel_HL_319559";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_319736.root"; #with prefilter deta > 0.04, with deta > 0.04
    #arr_m = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_319736";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_319736";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_TightEvSel_HL_319736";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_TightEvSel_HL_319736";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_319736";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_320049.root"; #with prefilter deta > 0.03, with deta > 0.03
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta003_HL_320049";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta003_HL_320049";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta003_TightEvSel_HL_320049";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta003_TightEvSel_HL_320049";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta003_HL_320049";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #filename = "AnalysisResults_HL_320321.root"; #with prefilter deta > 0.02, with deta  > 0.02 cut
    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT400MeV_HL_320321";
    #arr_m = np.array([0, 0.14, 0.5, 1.1, 2.0, 2.7, 3.2, 4], dtype=float);
    ##arr_m = np.array([0, 0.14, 1.1, 2.7, 3.2, 4], dtype=float);
    #arr_pt_jpsi = np.array([0., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20, suffix);
    #draw_v2_vs_mee_new(filename, taskname, 30, 50, arr_m, 1.0, 10.0, 0.0, 20, suffix);

    #arr_pt_pion = np.array([1., 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10], dtype=float);
    #arr_pt_vm = np.array([1., 2, 4, 6, 10], dtype=float);
    #arr_pt_imr = np.array([0, 1., 2, 4, 6, 10], dtype=float);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_pion, 0.02, 0.14, 0.0, 20, suffix, False, False, True);

    #taskname = "dielectron_v2_3050_TPChadrejorTOFreq_minpt200";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_TPChadrejorTOFreq_pT200MeV_HL_320321";
    #arr_pt_pion = np.array([0.4, 0.6, 0.8, 1., 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10], dtype=float);
    #draw_v2_vs_ptee_new(filename, taskname, 30, 50, arr_pt_pion, 0.02, 0.14, 0.0, 20, suffix, False, False, True);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_320161.root"; #with prefilter deta > 0.02, with deta > 0.02
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta002_HL_320161";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta002_HL_320161";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta002_TightEvSel_HL_320161";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta002_TightEvSel_HL_320161";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta002_HL_320161";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_320499.root"; #with prefilter deta > 0.04, with deta > 0.04, ITS cluster size < 3
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    ##arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_320499";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #taskname = "dielectron_0010_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0010_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_320499";
    #draw_mee_uls_ls(filename  , taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 10, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #arr_m = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #taskname = "dielectron_3050_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_3050_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_pf_deta004_HL_320499";
    #draw_mee_uls_ls(filename  , taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 30, 50, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_321612.root"; #without prefilter, without deta
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);

    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_321612";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_328076.root"; #without prefilter, without deta
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    #taskname = "dielectron_occupancy0_10000_TPChadrejorTOFreq_minpt400";
    #suffix = "_PbPb_5.36TeV_LHC24as_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_328076";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_329829_330348.root"; #LHC24ar + LHC24as
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);


    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi"; #prefilter deta > 0.04
    #suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_HL_329289_330348";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_TightEvSel_HL_329289_330348";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    arr_m_tmp = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.6, 2.8, 2.9, 3.00, 3.10, 3.2, 4.0], dtype=float);
    taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #prefilter deta > 0.04
    suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_329289_330348";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m_tmp, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m_tmp, 0.0, 10.0, 0.0, 1.0, suffix);
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m_tmp, 0.0, 10.0, 2.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_dcaee(filename, taskname, 0, 90, arr_m_tmp, 0.0, 10, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC24_pass1_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_TightEvSel_HL_329289_330348";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    #filename = "dielectron_PbPb_5.36TeV_LHC23zzh_pass4_HL_328095_328470.root"; #LHC23
    #arr_m = np.array([0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.7, 2.8, 2.9, 3.00, 3.10, 3.2, 3.5, 4.0], dtype=float);
    #arr_pt = np.array([1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20], dtype=float);


    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi"; #prefilter deta > 0.04
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_HL_328095_328470";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_wo_detadphi_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_wo_detadphi_TightEvSel_HL_328095_328470";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400"; #prefilter deta > 0.04
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_HL_328095_328470";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);

    #taskname = "dielectron_0090_TPChadrejorTOFreq_minpt400_TightEvSel";
    #suffix = "_PbPb_5.36TeV_LHC23_pass4_0090_ft0coccupancy0_10000_TPChadrejorTOFreq_pT400MeV_TightEvSel_HL_328095_328470";
    #draw_mee_uls_ls(filename  , taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_sbratio(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_mee_significance(filename, taskname, 0, 90, arr_m, 0.0, 10.0, 0.0, 20.0, suffix);


    filename = "dielectron_pp_13.6TeV_LHC23_pass4_thin_HL_349170.root"; # for default, woTTCA, wWeighting
    tasknames = [
        "dielectron_TPChadrejorTOFreq_minpt200_woTTCA",
        "dielectron_TPChadrejorTOFreq_minpt200",
        "dielectron_TPChadrejorTOFreq_minpt200_wWeighting",
    ];
    arr_m = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
    0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.91, 0.92, 0.94, 0.95, 0.96, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.08, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.8, 2.85, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0], dtype=float);
    arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    suffix = "_pp_13.6TeV_LHC23_LHC23_Thin_apass4_pT200MeV_HL_349170";
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 20.0, suffix);
    tasknames = [
        "dielectron_TPChadrejorTOFreq_minpt200_woTTCA",
        "dielectron_TPChadrejorTOFreq_minpt200_wWeighting",
    ];
    #draw_mee_ttca_recovery(filename, tasknames, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_dcaee_ttca_recovery(filename, tasknames, 0.5, 1.1, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_recovery(filename, tasknames, 1.1, 2.7, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_recovery(filename, tasknames, 2.7, 3.2, 0.0, 10.0, suffix);
    tasknames = [
        "dielectron_TPChadrejorTOFreq_minpt200_wWeighting",
        "dielectron_TPChadrejorTOFreq_minpt200",
    ];
    #draw_mee_ttca_multiple_counting(filename, tasknames, 0.0, 10.0, 0.0, 20.0, suffix);
    #draw_dcaee_ttca_multiple_counting(filename, tasknames, 0.0, 0.14, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_multiple_counting(filename, tasknames, 0.14, 0.5, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_multiple_counting(filename, tasknames, 0.5, 1.1, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_multiple_counting(filename, tasknames, 1.1, 2.7, 0.0, 10.0, suffix);
    #draw_dcaee_ttca_multiple_counting(filename, tasknames, 2.7, 3.2, 0.0, 10.0, suffix);

    #filename = "dimuon_pp_5.36TeV_LHC24ap_pass1_HL_385409_385423.root";
    #taskname = "dimuon_global_tight_minpt300";
    #suffix = "_pp_5.36TeV_LHC24apaq_pass1_pT300MeV_HL_385409_385423";
    ##draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    ##draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    ##draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    ##draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    ##draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_significance(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);

    #taskname = "dimuon_global";
    #suffix = "_pp_5.36TeV_LHC24apaq_global_HL_385409_385423";
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_significance(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);


    #filename = "dimuon_pp_5.36TeV_LHC24ap_pass1_HL_389401.root";
    #taskname = "dimuon_global_tight_minpt300";
    #suffix = "_pp_13.6TeV_LHC22o_pass7_pT300MeV_HL_389401";
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_significance(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);


    filename = "dimuon_pp_13.60TeV_LHC23_Thin_pass4_HL_394826.root";
    taskname = "dimuon";
    suffix = "_pp_13.6TeV_LHC23_Thin_pass4_pT300MeV_global_HL_398826";
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    #draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_uls_ls_2panel(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    #draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    #draw_mmumu_significance(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);


    filename = "dimuon_pp_5.36TeV_LHC24ap_aq_pass1_HL_414591_414810.root";
    taskname = "dimuon";
    suffix = "_pp_5.36TeV_LHC24ap_aq_pass1_pT300MeV_global_HL_414591_414810";
    draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    draw_mmumu_uls_ls(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    draw_mmumu_uls_ls_2panel(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);
    draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 2.0, suffix);
    draw_mmumu_sbratio(filename, taskname, 0, 999, 0.0, 10.0, 4.0, 999.0, suffix);
    draw_mmumu_significance(filename, taskname, 0, 999, 0.0, 10.0, 0.0, 999.0, suffix);

#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
