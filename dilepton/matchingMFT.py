import os
import datetime
import sys
import shutil
import math
import numpy as np
import sys
sys.path.append("../common/");
from painter import make_common_style
from histo_manager import rebin_histogram
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross

#__________________________________________________________
def compare_correct_match(filename, tasknames, suffix=""):
    print("reading...", filename);
    print(tasknames);

    arr_pt = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);

    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kCyan+1, kMagenta+1];
    rootfile = TFile.Open(filename, "READ");

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.14, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(0., 0, 10, 1.1);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#mu} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{correct}{correct + fake}");
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

    leg = TLegend(0.5, 0.15, 0.7, 0.4);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);

    for idx, taskname in enumerate(tasknames):
        print(idx, taskname);
        rootdire_correct = rootfile.Get(taskname + "/MFTMCHMID/primary/correct");
        rootdire_wrong = rootfile.Get(taskname + "/MFTMCHMID/primary/wrong");
        #rootdire_correct.ls();
        #rootdire_wrong.ls();

        h1pt_correct_org = rootdire_correct.Get("hPt");
        h1pt_correct_org.SetName("h1Pt_correct_org_{0:d}".format(idx));
        h1pt_wrong_org = rootdire_wrong.Get("hPt");
        h1pt_wrong_org.SetName("h1Pt_wrong_org_{0:d}".format(idx));

        h1pt_correct_org.Sumw2();
        h1pt_wrong_org.Sumw2();

        h1pt_correct = rebin_histogram(h1pt_correct_org, arr_pt, False, False);
        h1pt_wrong = rebin_histogram(h1pt_wrong_org, arr_pt, False, False);
        h1pt_correct.SetName("h1Pt_correct_{0:d}".format(idx));
        h1pt_wrong.SetName("h1Pt_wrong_{0:d}".format(idx));

        h1pt_all = h1pt_correct.Clone("h1pt_all_{0:d}".format(idx));
        h1pt_all.Add(h1pt_wrong, 1.0);

        h1purity = h1pt_correct.Clone("h1purity_{0:d}".format(idx));
        h1purity.Sumw2();
        h1purity.Reset();
        h1purity.Divide(h1pt_correct, h1pt_all, 1, 1, "B");
        make_common_style(h1purity, 20, 1.2, colors[idx], 2, 0);
        h1purity.Draw("E0,same");
        ROOT.SetOwnership(h1purity, False);
        h1purity.SetDirectory(0);
        if taskname == "matching-mft":
            leg.AddEntry(h1purity, "No cut", "P");
        else:
            leg.AddEntry(h1purity, taskname.replace("matching-mft_", ""), "P");

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15, 0.15, 0.4, 0.40, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#minus3.6 < #it{#eta}_{#mu} < #minus2.5");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_global_muon_matching_purity{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_global_muon_matching_purity{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_global_muon_matching_purity{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
def compare_match_efficiency(filename, tasknames, suffix=""):
    print("reading...", filename);
    print(tasknames);

    arr_pt = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);

    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kCyan+1, kMagenta+1];
    rootfile = TFile.Open(filename, "READ");

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    p1 = c1.cd(1);
    p1.SetMargin(0.15, 0.03, 0., 0.03);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., 5e-7, 10, 1e-2);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#mu} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#mu}} (GeV/#it{c})^{#minus1}");
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

    leg = TLegend(0.6, 0.6, 0.8, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);

    list_h1 = [];
    for idx, taskname in enumerate(tasknames):
        print(idx, taskname);
        h1z = rootfile.Get(taskname + "/Event/hZvtx");
        nev = h1z.GetEntries();
        print("nev = {0:e}".format(nev));
        rootdire_correct = rootfile.Get(taskname + "/MFTMCHMID/primary/correct");

        h1pt_org = rootdire_correct.Get("hPt");
        h1pt_org.SetName("h1Pt_correct_org_{0:d}".format(idx));
        h1pt_org.Sumw2();
        h1pt = rebin_histogram(h1pt_org, arr_pt, False, False);
        h1pt.SetName("h1Pt_{0:d}".format(idx));
        make_common_style(h1pt, 20, 1.2, colors[idx], 2, 0);
        h1pt.Scale(1, "width");
        h1pt.Scale(1/nev);
        h1pt.Draw("E0,same");
        ROOT.SetOwnership(h1pt, False);
        h1pt.SetDirectory(0);
        if taskname == "matching-mft":
            leg.AddEntry(h1pt, "No cut", "P");
        else:
            leg.AddEntry(h1pt, taskname.replace("matching-mft_", ""), "P");
        list_h1.append(h1pt);

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.3, 0.7, 0.5, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#minus3.6 < #it{#eta}_{#mu} < #minus2.5");
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    p2 = c1.cd(2);
    p2.SetMargin(0.15, 0.03, 0.24, 0.00);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.5, 10, 1.05);
    frame2.GetXaxis().SetTitle("#it{p}_{T,#mu} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("ratio to No cut");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    ROOT.SetOwnership(frame2,False);

    h1pt_nocut = list_h1[0].Clone("h1pt_nocut");
    for idx in range(1, len(tasknames)):
        taskname = tasknames[idx];
        h1ratio = list_h1[idx].Clone("h1ratio_{0:d}".format(idx));
        h1ratio.Reset();
        h1ratio.Sumw2();
        h1ratio.Divide(list_h1[idx], list_h1[0], 1., 1., "B");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_global_muon_efficiency{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_global_muon_efficiency{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_global_muon_efficiency{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_388590.root";
    filename = "AnalysisResults_HL_389454.root";
    suffix = "_LHC23k4g_HL_388590";
    tasknames = [
        "matching-mft",
        "matching-mft_matchingchi240",
        "matching-mft_Rabs276",
        "matching-mft_chi240",
        #"matching-mft_dcaxy01",
    ];

    compare_correct_match(filename, tasknames, suffix);
    compare_match_efficiency(filename, tasknames, suffix);
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
