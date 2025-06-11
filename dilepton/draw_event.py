import datetime
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TPaletteAxis, TLine
from ROOT import gStyle, gROOT, gSystem, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
from histo_manager import get_ratio
from signal_extractor import get_significance
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);


#__________________________________________________________
def draw_occupancy_1d(filename, taskname, period, apass_number, suffix):
    rootfile = TFile.Open(filename, "REAAD");
    rootdire_task = rootfile.Get(taskname);
    rootdire_ev = rootdire_task.Get("Event/before");
    rootdire_ev.ls();

    h1z = rootdire_ev.Get("hZvtx");
    nev = h1z.GetEntries();

    h2 = rootdire_ev.Get("hMultFT0CvsOccupancy");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h1 = h2.ProjectionY("h1");
    h1.SetDirectory(0);
    ROOT.SetOwnership(h1, False);
    make_common_style(h1, 20, 1.0, kBlack, 2, 0);
    h1.RebinX(10);
    h1.Scale(1/nev);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.10, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    frame1 = c1.DrawFrame(0., 1e-8, 2e+4,  1);
    frame1.GetXaxis().SetTitle("occupancy");
    frame1.GetYaxis().SetTitle("1/N_{ev}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetXaxis().SetMaxDigits(3);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    h1.Draw("E0same");

    txt = TPaveText(0.5, 0.77, 0.8, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.AddText(suffix.replace("_", "", 1));
    txt.AddText("N_{{ev}} = {0:3.1f} M events".format(nev/1e+6));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_occupancy_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_occupancy_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_occupancy_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_ncontrib_centrality(filename, taskname, period, apass_number, suffix):
    rootfile = TFile.Open(filename, "REAAD");
    rootdire_task = rootfile.Get(taskname);
    rootdire_ev = rootdire_task.Get("Event/after");
    #rootdire_ev = rootdire_task.Get("Event/before");
    rootdire_ev.ls();

    h1z = rootdire_ev.Get("hZvtx");
    nev = h1z.GetEntries();

    h2 = rootdire_ev.Get("hCentFT0CvsMultNTracksPV");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.11, 0.1, 0.05);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0., 0, 100,  5000);
    frame1.GetXaxis().SetTitle("centrality FT0C (%)");
    frame1.GetYaxis().SetTitle("N_{contributors} to PV");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetXaxis().SetMaxDigits(3);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    h2.SetContour(1000);
    h2.Draw("colz,same");

    h2.SetMinimum(1);
    h2.SetMaximum(1e+7);

#    gPad.Update();
#    palette = h2.GetListOfFunctions().FindObject("palette");
#    #print(palette.GetX1NDC());
#    #print(palette.GetX2NDC());
#    #print(palette.GetY1NDC());
#    #print(palette.GetY2NDC());
#    #palette.SetX1NDC(0.9);
#    #palette.SetX2NDC(0.99);
#    palette.SetY1NDC(0.16);
#    #palette.SetY2NDC(0.95);

    txt = TPaveText(0.5, 0.77, 0.8, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.AddText(suffix.replace("_", "", 1));
    txt.AddText("N_{{ev}} = {0:3.1f} M events".format(nev/1e+6));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_NPV_eta08_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_eta08_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_eta08_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_ncontrib_ft0c(filename, taskname, period, apass_number, suffix):
    rootfile = TFile.Open(filename, "REAAD");
    rootdire_task = rootfile.Get(taskname);
    rootdire_ev = rootdire_task.Get("Event/after");
    #rootdire_ev = rootdire_task.Get("Event/before");
    rootdire_ev.ls();

    h1z = rootdire_ev.Get("hZvtx");
    nev = h1z.GetEntries();

    h2 = rootdire_ev.Get("hMultFT0CvsMultNTracksPV");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.11, 0.1, 0.05);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0., 0, 60000,  5000);
    frame1.GetXaxis().SetTitle("multiplicity FT0C");
    frame1.GetYaxis().SetTitle("N_{contributors} to PV");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetXaxis().SetMaxDigits(3);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    h2.SetContour(1000);
    h2.Draw("colz,same");

    h2.SetMinimum(1);
    h2.SetMaximum(1e+7);

    gPad.Update();
    palette = h2.GetListOfFunctions().FindObject("palette");
    #print(palette.GetX1NDC());
    #print(palette.GetX2NDC());
    #print(palette.GetY1NDC());
    #print(palette.GetY2NDC());
    #palette.SetX1NDC(0.9);
    #palette.SetX2NDC(0.99);
    palette.SetY1NDC(0.16);
    #palette.SetY2NDC(0.95);

    txt = TPaveText(0.15, 0.77, 0.56, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.AddText(suffix.replace("_", "", 1));
    txt.AddText("N_{{ev}} = {0:3.1f} M events".format(nev/1e+6));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_NPV_eta08_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_eta08_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_eta08_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_ncontrib_1d(filename, tasknames, period, apass_number, suffix, arr_occ_min, arr_occ_max):
    rootfile = TFile.Open(filename, "REAAD");

    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.13, 0.08, 0.0, 0.03);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0., 2e-6, 4000, 1e-0);
    frame1.GetXaxis().SetTitle("N_{contributors} to PV");
    frame1.GetYaxis().SetTitle("1/N_{ev}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.1);
    frame1.GetXaxis().SetTitleSize(0.06);
    frame1.GetYaxis().SetTitleSize(0.06);
    frame1.GetXaxis().SetLabelSize(0.06);
    frame1.GetYaxis().SetLabelSize(0.06);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    #frame1.GetXaxis().SetMaxDigits(3);
    #frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(p1,False);

    leg = TLegend(0.2, 0.2, 0.4, 0.5);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);

    list_h1mult = [];
    for idx, taskname in enumerate(tasknames):
        rootdire_task = rootfile.Get(taskname);
        rootdire_ev = rootdire_task.Get("Event/after");
        #rootdire_ev.ls();

        h2 = rootdire_ev.Get("hMultFT0CvsMultNTracksPV");
        h1mult = h2.ProjectionY("h1mult_{0:d}".format(idx));
        nev = h1mult.GetEntries();
        h1mult.RebinX(10);
        h1mult.Scale(1/nev);
        h1mult.SetDirectory(0);
        ROOT.SetOwnership(h1mult, False);
        h1mult.Draw("E0same");
        make_common_style(h1mult, 20, 1.0, colors[idx], 1, 0);
        list_h1mult.append(h1mult);
        leg.AddEntry(h1mult, "{0:d} #leq occupancy < {1:d}".format(arr_occ_min[idx], arr_occ_max[idx]), "P");

    txt = TPaveText(0.5, 0.77, 0.8, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.13, 0.08, 0.25, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.6, 4000, 1.25);
    frame2.GetXaxis().SetTitle("N_{contributors} to PV");
    frame2.GetYaxis().SetTitle("ratio to 0 #leq occ. < 1000");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.5);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetXaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetMaxDigits(3);
    #frame2.GetYaxis().SetNdivisions(505);
    ROOT.SetOwnership(frame2,False);
    ROOT.SetOwnership(p2,False);

    line = TLine(0, 1, 4000, 1);
    line.SetLineColor(kBlack);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw("");
    ROOT.SetOwnership(line, False);

    for i in range(1, len(list_h1mult)):
        h1 = list_h1mult[i].Clone("h1mult_{0:d}".format(i));
        h1.Divide(list_h1mult[0]);
        h1.Draw("E0same");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_NPV_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_NPV_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, apass_number, suffix, arr_occ_min, arr_occ_max):
    rootfile = TFile.Open(filename, "REAAD");

    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.12, 0.08, 0.0, 0.07);
    p1.SetTicks(1,1);
    frame1 = p1.DrawFrame(0., -100, 60000, 4000);
    frame1.GetXaxis().SetTitle("mult. FT0C");
    frame1.GetYaxis().SetTitle("<N_{contributors} to PV>");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.06);
    frame1.GetYaxis().SetTitleSize(0.06);
    frame1.GetXaxis().SetLabelSize(0.06);
    frame1.GetYaxis().SetLabelSize(0.06);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    #frame1.GetXaxis().SetMaxDigits(3);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(p1,False);

    leg = TLegend(0.45, 0.05, 0.65, 0.35);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);

    list_h1mult = [];
    for idx, taskname in enumerate(tasknames):
        rootdire_task = rootfile.Get(taskname);
        rootdire_ev = rootdire_task.Get("Event/after");
        #rootdire_ev.ls();

        h2 = rootdire_ev.Get("hMultFT0CvsMultNTracksPV");
        h1mult = h2.ProfileX("h1prf_mult_{0:d}".format(idx)).ProjectionX("h1prj_mult_{0:d}".format(idx));
        h1mult.SetDirectory(0);
        ROOT.SetOwnership(h1mult, False);
        h1mult.Draw("E0same");
        make_common_style(h1mult, 20, 1.0, colors[idx], 1, 0);
        list_h1mult.append(h1mult);
        leg.AddEntry(h1mult, "{0:d} #leq occupancy < {1:d}".format(arr_occ_min[idx], arr_occ_max[idx]), "P");

    txt = TPaveText(0.2, 0.75, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.12, 0.08, 0.25, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.9, 60000, 1.04);
    frame2.GetXaxis().SetTitle("mult. FT0C");
    frame2.GetYaxis().SetTitle("ratio to 0 #leq occ. < 1000");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetXaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetNdivisions(505);
    ROOT.SetOwnership(frame2,False);
    ROOT.SetOwnership(p2,False);

    line = TLine(0, 1, 60000, 1);
    line.SetLineColor(kBlack);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw("");
    ROOT.SetOwnership(line, False);

    for i in range(1, len(list_h1mult)):
        h1 = list_h1mult[i].Clone("h1mult_{0:d}".format(i));
        h1.Divide(list_h1mult[0]);
        h1.Draw("E0same");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_profile_NPV_vs_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_profile_NPV_vs_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_profile_NPV_vs_multFT0C_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, apass_number, suffix, arr_occ_min, arr_occ_max):
    rootfile = TFile.Open(filename, "REAAD");

    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];
    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1, 2, 0, 0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.12, 0.03, 0.0, 0.07);
    p1.SetTicks(1,1);
    frame1 = p1.DrawFrame(0., -100, 100, 4000);
    frame1.GetXaxis().SetTitle("centrality FT0C (%d)");
    frame1.GetYaxis().SetTitle("<N_{contributors} to PV>");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.06);
    frame1.GetYaxis().SetTitleSize(0.06);
    frame1.GetXaxis().SetLabelSize(0.06);
    frame1.GetYaxis().SetLabelSize(0.06);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    #frame1.GetXaxis().SetMaxDigits(3);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(p1,False);

    leg = TLegend(0.45, 0.4, 0.65, 0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);

    list_h1mult = [];
    for idx, taskname in enumerate(tasknames):
        rootdire_task = rootfile.Get(taskname);
        rootdire_ev = rootdire_task.Get("Event/after");
        #rootdire_ev.ls();

        h2 = rootdire_ev.Get("hCentFT0CvsMultNTracksPV");
        h1mult = h2.ProfileX("h1prf_mult_{0:d}".format(idx)).ProjectionX("h1prj_mult_{0:d}".format(idx));
        h1mult.SetDirectory(0);
        ROOT.SetOwnership(h1mult, False);
        h1mult.Draw("E0same");
        make_common_style(h1mult, 20, 1.0, colors[idx], 1, 0);
        list_h1mult.append(h1mult);
        leg.AddEntry(h1mult, "{0:d} #leq occupancy < {1:d}".format(arr_occ_min[idx], arr_occ_max[idx]), "P");

    txt = TPaveText(0.2, 0.75, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(1001);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("{0}_{1}".format(period, apass_number));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.12, 0.03, 0.25, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0.9, 100, 1.04);
    frame2.GetXaxis().SetTitle("centrality FT0C (%)");
    frame2.GetYaxis().SetTitle("ratio to 0 #leq occ. < 1000");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetXaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetMaxDigits(3);
    frame2.GetYaxis().SetNdivisions(505);
    ROOT.SetOwnership(frame2,False);
    ROOT.SetOwnership(p2,False);

    line = TLine(0, 1, 100, 1);
    line.SetLineColor(kBlack);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw("");
    ROOT.SetOwnership(line, False);

    for i in range(1, len(list_h1mult)):
        h1 = list_h1mult[i].Clone("h1mult_{0:d}".format(i));
        h1.Divide(list_h1mult[0]);
        h1.Draw("E0same");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_profile_NPV_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.eps".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_profile_NPV_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.pdf".format(date, period, apass_number, suffix));
    c1.SaveAs("{0}_profile_NPV_vs_centFT0C_PbPb_5.36TeV_{1}_{2}{3}.png".format(date, period, apass_number, suffix));

    rootfile.Close();
#__________________________________________________________
if __name__ == "__main__":
    filename = "AnalysisResults_HL_246294.root";
    period = "LHC23_PbPb";
    pass_number = "apass4";

    taskname = "single-electron-qc_occupancy0_1000_TOFreq_minpt400";
    suffix = "_nocut";
    #draw_occupancy_1d(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy0_1000_TOFreq_minpt400";
    suffix = "_occupancy0_1000";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy1000_2000_TOFreq_minpt400";
    suffix = "_occupancy1000_2000";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy2000_3000_TOFreq_minpt400";
    suffix = "_occupancy2000_3000";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy3000_4000_TOFreq_minpt400";
    suffix = "_occupancy3000_4000";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy4000_5000_TOFreq_minpt400";
    suffix = "_occupancy4000_5000";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    taskname = "single-electron-qc_occupancy5000_99999_TOFreq_minpt400";
    suffix = "_occupancy5000_99999";
    #draw_ncontrib_ft0c(filename, taskname, period, pass_number, suffix);
    #draw_ncontrib_centrality(filename, taskname, period, pass_number, suffix);

    tasknames = [
        "single-electron-qc_occupancy0_1000_TOFreq_minpt400",
        "single-electron-qc_occupancy1000_2000_TOFreq_minpt400",
        "single-electron-qc_occupancy2000_3000_TOFreq_minpt400",
        "single-electron-qc_occupancy3000_4000_TOFreq_minpt400",
        "single-electron-qc_occupancy4000_5000_TOFreq_minpt400",
        "single-electron-qc_occupancy5000_99999_TOFreq_minpt400",
    ];

    arr_occ_min = np.array([   0, 1000, 2000, 3000, 4000, 5000], dtype=int);
    arr_occ_max = np.array([1000, 2000, 3000, 4000, 5000, 99999], dtype=int);
    #draw_ncontrib_1d(filename, tasknames, period, pass_number, "", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "", arr_occ_min, arr_occ_max);

    filename = "AnalysisResults_HL_246294_544116.root";
    #draw_occupancy_1d(filename, taskname, period, pass_number, "_HL_246294_544116");
    #draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_246294_544116", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_246294_544116", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_246294_544116", arr_occ_min, arr_occ_max);

    filename = "AnalysisResults_HL_246294_544123.root";
    #draw_occupancy_1d(filename, taskname, period, pass_number, "_HL_246294_544123");
    #draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_246294_544123", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_246294_544123", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_246294_544123", arr_occ_min, arr_occ_max);


    tasknames = [
        "single-electron-qc_occupancy0_1000_TOFreq_minpt400_id15336",
        "single-electron-qc_occupancy1000_2000_TOFreq_minpt400_id15336",
        "single-electron-qc_occupancy2000_3000_TOFreq_minpt400_id15336",
        "single-electron-qc_occupancy3000_4000_TOFreq_minpt400_id15336",
        "single-electron-qc_occupancy4000_5000_TOFreq_minpt400_id15336",
        "single-electron-qc_occupancy5000_99999_TOFreq_minpt400_id15336",
    ];

    #filename = "AnalysisResults_HL_250911.root";
    #draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_250911", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_250911", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_250911", arr_occ_min, arr_occ_max);

    tasknames = [
        "single-electron-qc_occupancy0_1000_TOFreq_minpt400_NoCollInTRS_id15700",
        "single-electron-qc_occupancy1000_2000_TOFreq_minpt400_NoCollInTRS_id15700",
        "single-electron-qc_occupancy2000_3000_TOFreq_minpt400_NoCollInTRS_id15700",
        "single-electron-qc_occupancy3000_4000_TOFreq_minpt400_NoCollInTRS_id15700",
        "single-electron-qc_occupancy4000_5000_TOFreq_minpt400_NoCollInTRS_id15700",
        "single-electron-qc_occupancy5000_99999_TOFreq_minpt40_NoCollInTRS0_id15700",
    ];

    #filename = "AnalysisResults_HL_250911.root";
    #draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_250911_NoCollInTRS", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_250911_NoCollInTRS", arr_occ_min, arr_occ_max);
    #draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_250911_NoCollInTRS", arr_occ_min, arr_occ_max);

    filename = "AnalysisResults_HL_250911_544116.root";
    draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_250911_544116_NoCollInTRS", arr_occ_min, arr_occ_max);
    draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_250911_544116_NoCollInTRS", arr_occ_min, arr_occ_max);
    draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_250911_544116_NoCollInTRS", arr_occ_min, arr_occ_max);

    filename = "AnalysisResults_HL_250911_544123.root";
    draw_ncontrib_1d(filename, tasknames, period, pass_number, "_HL_250911_544123_NoCollInTRS", arr_occ_min, arr_occ_max);
    draw_profile_ncontrib_vs_multft0c(filename, tasknames, period, pass_number, "_HL_250911_544123_NoCollInTRS", arr_occ_min, arr_occ_max);
    draw_profile_ncontrib_vs_centft0c(filename, tasknames, period, pass_number, "_HL_250911_544123_NoCollInTRS", arr_occ_min, arr_occ_max);

