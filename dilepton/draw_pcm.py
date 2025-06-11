import numpy as np
import datetime
import math
import sys
sys.path.append("../common/");
import ROOT
import ctypes
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kGray
from painter import make_common_style
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#_____________________________________________________________________
def draw_run2_run3(filename_run2, filename_run3, taskname, suffix=""):
    print("reading...", filename_run2, filename_run3);

    rootfile_run2 = TFile.Open(filename_run2, "READ");
    h1rxy_run2 = rootfile_run2.Get("histoRData");
    make_common_style(h1rxy_run2, 21, 1.2, kGray+2, 2, 0);
    h1rxy_run2.SetDirectory(0);
    ROOT.SetOwnership(h1rxy_run2, False);
    bin74 = h1rxy_run2.FindBin(74.0 + 1e-3);
    bin86 = h1rxy_run2.FindBin(86.0 - 1e-3);
    ng_74_86_run2 = h1rxy_run2.Integral(bin74, bin86, "");

    #h1rxy_run2.GetXaxis().SetRangeUser(74, 86);
    #ng_74_86_run2 = h1rxy_run2.GetMaximum();
    #h1rxy_run2.GetXaxis().SetRangeUser(0, 90);


    rootfile_run3 = TFile.Open(filename_run3, "READ");
    rootdir = rootfile_run3.Get(taskname);
    rootdir.ls();
    list_ev = rootdir.Get("Event/after");
    list_ev.ls();
    list_v0 = rootdir.Get("V0");
    list_v0.ls();

    h1nch = list_ev.Get("hMultNTracksPV");
    nch = 0;
    for i in range(0, h1nch.GetNbinsX()):
        nch += h1nch.GetBinContent(i+1) * h1nch.GetBinCenter(i+1);
    print("nch = {0:e}".format(nch));

    h2 = list_v0.Get("hMassGamma");
    h1rxy_run3 = h2.ProjectionX("h1rxy_run3");
    h1rxy_run3.SetDirectory(0);
    ROOT.SetOwnership(h1rxy_run3, False);
    h1rxy_run3.Scale(1, "width");
    h1rxy_run3.Scale(1/nch);
    make_common_style(h1rxy_run3, 20, 1.2, kRed+1, 2, 0);

    bin74 = h1rxy_run3.FindBin(74.0 + 1e-3);
    bin86 = h1rxy_run3.FindBin(86.0 - 1e-3);
    ng_74_86_run3 = h1rxy_run3.Integral(bin74, bin86, "");
    #h1rxy_run3.GetXaxis().SetRangeUser(74, 86);
    #ng_74_86_run3 = h1rxy_run3.GetMaximum();
    #h1rxy_run3.GetXaxis().SetRangeUser(0, 90);

    sf = ng_74_86_run2 / ng_74_86_run3;
    print("sf = ", sf);
    h1rxy_run3.Scale(sf);
 
    c1 = TCanvas("c0","c0",0,0,850,850);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);
    c1.SetMargin(0.18, 0.03, 0.10, 0.05);

    frame1 = c1.DrawFrame(0, 0, 90, 2.5e-3);
    frame1.GetXaxis().SetTitle("R_{xy} (cm)");
    frame1.GetYaxis().SetTitle("#frac{1}{N_{ch}} #frac{dN_{#gamma #rightarrow e^{+}e^{#minus}}}{dR_{xy}} (cm)^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);
    frame1.GetYaxis().SetMaxDigits(3);
    h1rxy_run2.Draw("E0H,same");
    h1rxy_run3.Draw("E0H,same");

    leg = TLegend(0.2, 0.82, 0.4, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1rxy_run2 , "pp at #sqrt{s} = 5.02 TeV, JINST 18 (2023) P11032", "P");
    #leg.AddEntry(h1rxy_run3 , "pp at #sqrt{s} = 13.6 TeV", "P");
    leg.AddEntry(h1rxy_run3 , "pp at #sqrt{{s}} = 13.6 TeV #times {0:3.2f}".format(sf), "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #txt = TPaveText(0.22, 0.85, 0.4, 0.92,"NDC");
    #txt.SetFillColor(kWhite);
    #txt.SetFillStyle(0);
    #txt.SetBorderSize(0);
    #txt.SetTextAlign(12);#middle,left
    #txt.SetTextFont(42);#helvetica bold
    #txt.SetTextSize(0.04);
    #txt.AddText("0#minus90%, Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    #txt.Draw();
    #ROOT.SetOwnership(txt,False);

    txt_pcm = TPaveText(0.2, 0.75, 0.4, 0.82, "NDC");
    txt_pcm.SetFillColor(kWhite);
    txt_pcm.SetFillStyle(0);
    txt_pcm.SetBorderSize(0);
    txt_pcm.SetTextAlign(12);#middle,left
    txt_pcm.SetTextFont(42);#helvetica
    txt_pcm.SetTextSize(0.035);
    txt_pcm.AddText("p_{T,#gamma} > 0.1 GeV/c, |#eta_{#gamma}| < 0.8");
    txt_pcm.Draw();
    ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_PhotonConversionPointRxy_Run2_Run3{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_PhotonConversionPointRxy_Run2_Run3{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_PhotonConversionPointRxy_Run2_Run3{1}.png".format(date, suffix));

    rootfile_run2.Close();
    rootfile_run3.Close();

#_____________________________________________________________________
#_____________________________________________________________________
def draw_rxy(filename, taskname, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get(taskname);
    rootdir.ls();
    list_ev = rootdir.Get("Event/after");
    list_ev.ls();
    list_v0 = rootdir.Get("V0");
    list_v0.ls();

    h1ev = list_ev.Get("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h1nch = list_ev.Get("hMultNTracksPV");
    nch = h1nch.GetMean();
    print("nch = {0:e}".format(nch));

    h2 = list_v0.Get("hMassGamma");
    h1 = h2.ProjectionX("h1");
    ROOT.SetOwnership(h1, False);
    h1.SetDirectory(0);
    h1.Scale(1, "width");
    h1.Scale(1/nch);
    h1.Scale(1/nev);
    make_common_style(h1, 20, 1.2, kBlack, 2, 0);
 
    c1 = TCanvas("c0","c0",0,0,850,850);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);
    c1.SetMargin(0.19, 0.03, 0.10, 0.05);

    frame1 = c1.DrawFrame(0, 0, 90, 4e-4);
    frame1.GetXaxis().SetTitle("r_{xy} (cm)");
    #frame1.GetYaxis().SetTitle("#frac{1}{N_{ch}} #frac{dN_{#gamma #rightarrow e^{+}e^{#minus}}}{dR_{xy}} (cm)^{#minus1}");
    frame1.GetYaxis().SetTitle("#frac{dN_{#gamma #rightarrow e^{+}e^{#minus}}}{dr_{xy}} (a.u.)");
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.2);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);
    frame1.GetYaxis().SetMaxDigits(3);
    h1.Draw("E0H,same");

    #leg = TLegend(0.6,0.7,0.8,0.9);
    #leg.SetBorderSize(0);
    #leg.SetFillColor(kWhite);
    #leg.SetTextSize(0.035);
    #leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    #leg.AddEntry(h1same , "same","LP");
    #leg.AddEntry(h1bkg  , "bkg (mixed event)","LP");
    #leg.AddEntry(h1sig  , "signal","LP");
    #leg.Draw("");
    #ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.22, 0.85, 0.4, 0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.04);
    txt.AddText("0#minus90%, Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

#    txt_pcm = TPaveText(0.4,0.92,0.7,1.0,"NDC");
#    txt_pcm.SetFillColor(kWhite);
#    txt_pcm.SetFillStyle(0);
#    txt_pcm.SetBorderSize(0);
#    txt_pcm.SetTextAlign(12);#middle,left
#    txt_pcm.SetTextFont(42);#helvetica
#    txt_pcm.SetTextSize(0.035);
#    txt_pcm.AddText("Photon Conversion Method (#gamma #rightarrow e^{+}e^{#minus})");
#    txt_pcm.AddText("|#it{#eta}_{#gamma}| < 0.8");
#    txt_pcm.Draw();
#    ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_conversion_point_rxy{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_conversion_point_rxy{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_conversion_point_rxy{1}.png".format(date, suffix));

    rootfile.Close();

#_____________________________________________________________________
#_____________________________________________________________________
def draw_xy(filename, taskname, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get(taskname);
    rootdir.ls();
    list_ev = rootdir.Get("Event/after");
    list_v0 = rootdir.Get("V0");
    list_v0.Print();

    h1ev = list_ev.Get("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h2 = list_v0.Get("hGammaRxy");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    h2.GetZaxis().SetTitle("N_{#gamma}/N_{ev}");
    h2.GetZaxis().SetTitleSize(0.04);
    h2.GetZaxis().SetTitleOffset(0.7);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,850,850);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.12,0.10,0.09);

    frame1 = c1.DrawFrame(-100,-100,100,100);
    frame1.GetXaxis().SetTitle("conversion point X (cm)");
    frame1.GetYaxis().SetTitle("conversion point Y (cm)");
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);
    h2.Scale(1/nev);
    print(h2.GetMinimum());
    print(h2.GetMaximum());
    h2.SetMinimum(1e-8);
    h2.SetMaximum(9e-4);
    h2.Draw("colz, same");

    #leg = TLegend(0.6,0.7,0.8,0.9);
    #leg.SetBorderSize(0);
    #leg.SetFillColor(kWhite);
    #leg.SetTextSize(0.035);
    #leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    #leg.AddEntry(h1same , "same","LP");
    #leg.AddEntry(h1bkg  , "bkg (mixed event)","LP");
    #leg.AddEntry(h1sig  , "signal","LP");
    #leg.Draw("");
    #ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.10,0.92,0.4,1.0,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.035);
    txt.AddText("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("centrality FT0C 0-90 %");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

#    txt_pcm = TPaveText(0.4,0.92,0.7,1.0,"NDC");
#    txt_pcm.SetFillColor(kWhite);
#    txt_pcm.SetFillStyle(0);
#    txt_pcm.SetBorderSize(0);
#    txt_pcm.SetTextAlign(12);#middle,left
#    txt_pcm.SetTextFont(42);#helvetica
#    txt_pcm.SetTextSize(0.035);
#    txt_pcm.AddText("Photon Conversion Method (#gamma #rightarrow e^{+}e^{#minus})");
#    txt_pcm.AddText("|#it{#eta}_{#gamma}| < 0.8");
#    txt_pcm.Draw();
#    ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_conversion_point_xy{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_conversion_point_xy{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_conversion_point_xy{1}.png".format(date, suffix));

    rootfile.Close();

#_____________________________________________________________________
#_____________________________________________________________________
if __name__ == "__main__":
    #draw_xy("AnalysisResults_HL_246294.root", "pcm-qc_occupancy0_1000", "_PbPb_2023_apass4");
    #draw_rxy("AnalysisResults_HL_246294.root", "pcm-qc_occupancy0_1000", "_PbPb_2023_apass4");
    #draw_run2_run3("MaterialRun2_paper.root", "AnalysisResults_HL_375352.root", "pcm-qc", "_absolute");
    #draw_run2_run3("MaterialRun2_paper.root", "AnalysisResults_HL_375352.root", "pcm-qc", "_scaled");
    #draw_run2_run3("MaterialRun2_paper.root", "AnalysisResults_HL_375352.root", "pcm-qc_woITSonly", "_woITSonly_absolute");
    draw_run2_run3("MaterialRun2_paper.root", "AnalysisResults_HL_375352.root", "pcm-qc_woITSonly", "_woITSonly_scaled");

