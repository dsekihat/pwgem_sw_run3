import datetime
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kViolet, kAzure
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
from histo_manager import get_ratio
from signal_extractor import get_significance
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
def draw_itsob_cluster_size_slice(fliename, taskname, cen1, cen2, suffix):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire_track = rootdire.Get("Track");
    rootdire_track.ls();

    arr_sizemax = np.array([15, 4, 3.5, 3, 2.5, 2, 1.5], dtype=float);
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1, kViolet+1, kAzure+7];

    rootdire_event = rootdire.Get("Event/after");
    #rootdire_event.ls();
    h1z = rootdire_event.Get("hZvtx");
    nev = h1z.GetEntries();

    h2pos = rootdire_track.Get("positive/hMeanClusterSizeITSob");
    h2neg = rootdire_track.Get("negative/hMeanClusterSizeITSob");
    h2pos.Sumw2();
    h2neg.Sumw2();

    h2 = h2pos.Clone("h2");
    h2.Reset();
    h2.Add(h2pos, 1.);
    h2.Add(h2neg, 1.);
    h2.SetContour(1000);
    h2.Scale(1/nev);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1,2,0,0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.1, 0.02, 0.0, 0.06);
    p1.SetTicks(1,1);
    #p1.SetLogy(1);

    y_max = 0;
    y_min = -0.001;
    frame1 = p1.DrawFrame(0, y_min, 2, 0.04);
    frame1.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame1.GetYaxis().SetTitle("N_{track}/N_{ev}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.5, 0.4, 0.7, 0.8);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    for i in range(0, len(arr_sizemax)):
        sizemax = arr_sizemax[i];
        bin1 = h2.GetYaxis().FindBin(0 + 1e-3);
        bin2 = h2.GetYaxis().FindBin(sizemax - 1e-3);
        h1 = h2.ProjectionX("h1_{0:d}".format(i), bin1, bin2, "");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);
        h1.Draw("E0,same");
        make_common_style(h1, 20, 1.0, colors[i], 2, 0);
        leg.AddEntry(h1, "<ITSob cluster size> #times cos(#lambda) < {0:2.1f}".format(sizemax), "P");
        list_h1.append(h1);
        if i == 0:
            y_max = h1.GetMaximum();
    leg.Draw("");
    frame1.GetYaxis().SetRangeUser(y_max * -0.05, y_max * 1.2);

    txt = TPaveText(0.15, 0.8, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.1, 0.02, 0.22, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0, 0.2, 2, 1.05);
    frame2.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame2.GetYaxis().SetTitle("ratio to no cut");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.5);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    for i in range(1, len(arr_sizemax)):
        h1tmp = list_h1[i];
        h1ratio = h1tmp.Clone("h1tmp_{0:d}".format(i));
        h1ratio.Divide(h1tmp, list_h1[0], 1., 1., "B");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_itsob_cluster_size_slice{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_itsob_cluster_size_slice{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_itsob_cluster_size_slice{1}.png".format(date, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_its_cluster_size_slice(fliename, taskname, cen1, cen2, suffix):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire_track = rootdire.Get("Track");
    rootdire_track.ls();

    arr_sizemax = np.array([15, 4, 3.5, 3, 2.5], dtype=float);
    #arr_sizemax = np.array([15, 4, 3.5, 3.2, 3], dtype=float);
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1, kViolet+1, kAzure+7];

    rootdire_event = rootdire.Get("Event/after");
    #rootdire_event.ls();
    h1z = rootdire_event.Get("hZvtx");
    nev = h1z.GetEntries();

    #h2pos = rootdire_track.Get("positive/hMeanClusterSizeITS");
    #h2neg = rootdire_track.Get("negative/hMeanClusterSizeITS");
    #h2pos.Sumw2();
    #h2neg.Sumw2();
    #h2 = h2pos.Clone("h2");
    #h2.Reset();
    #h2.Add(h2pos, 1.);
    #h2.Add(h2neg, 1.);

    h2 = rootdire_track.Get("hMeanClusterSizeITS");
    h2.Sumw2();

    h2.SetContour(1000);
    h2.Scale(1/nev);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1,2,0,0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.1, 0.02, 0.0, 0.06);
    p1.SetTicks(1,1);
    #p1.SetLogy(1);

    y_max = 0;
    y_min = -0.001;
    frame1 = p1.DrawFrame(0, y_min, 2, 0.04);
    frame1.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame1.GetYaxis().SetTitle("N_{track}/N_{ev}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.5, 0.4, 0.7, 0.8);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    for i in range(0, len(arr_sizemax)):
        sizemax = arr_sizemax[i];
        bin1 = h2.GetYaxis().FindBin(0 + 1e-3);
        bin2 = h2.GetYaxis().FindBin(sizemax - 1e-3);
        h1 = h2.ProjectionX("h1_{0:d}".format(i), bin1, bin2, "");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);
        h1.Draw("E0,same");
        make_common_style(h1, 20, 1.0, colors[i], 2, 0);
        leg.AddEntry(h1, "<ITS cluster size> #times cos(#lambda) < {0:2.1f}".format(sizemax), "P");
        list_h1.append(h1);
        if i == 0:
            y_max = h1.GetMaximum();
    leg.Draw("");
    frame1.GetYaxis().SetRangeUser(y_max * -0.05, y_max * 1.2);

    txt = TPaveText(0.15, 0.8, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.1, 0.02, 0.22, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0, 0.2, 2, 1.05);
    frame2.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame2.GetYaxis().SetTitle("ratio to no cut");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.5);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    for i in range(1, len(arr_sizemax)):
        h1tmp = list_h1[i];
        h1ratio = h1tmp.Clone("h1tmp_{0:d}".format(i));
        h1ratio.Divide(h1tmp, list_h1[0], 1., 1., "B");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_its_cluster_size_slice{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_its_cluster_size_slice{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_its_cluster_size_slice{1}.png".format(date, suffix));

    rootfile.Close();

#__________________________________________________________
#__________________________________________________________
def draw_tofchi2_slice(fliename, taskname, cen1, cen2, suffix):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire_track = rootdire.Get("Track");
    rootdire_track.ls();

    #arr_chi2max = np.array([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], dtype=float);
    arr_chi2max = np.array([10, 5, 3, 2, 1, 0.5, 0.1], dtype=float);
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1, kViolet+1, kAzure+7];

    rootdire_event = rootdire.Get("Event/after");
    rootdire_event.ls();
    h1z = rootdire_event.Get("hZvtx");
    nev = h1z.GetEntries();

    h2pos = rootdire_track.Get("positive/hChi2TOF");
    h2neg = rootdire_track.Get("negative/hChi2TOF");
    h2pos.Sumw2();
    h2neg.Sumw2();

    h2 = h2pos.Clone("h2");
    h2.Reset();
    h2.Add(h2pos, 1.);
    h2.Add(h2neg, 1.);
    h2.SetContour(1000);
    h2.Scale(1/nev);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1,2,0,0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.1, 0.02, 0.0, 0.06);
    p1.SetTicks(1,1);
    #p1.SetLogy(1);

    y_max = 0;
    y_min = -0.001;
    frame1 = p1.DrawFrame(0, y_min, 2, 0.03);
    frame1.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame1.GetYaxis().SetTitle("N_{track}/N_{ev}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.7, 0.4, 0.9, 0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    ROOT.SetOwnership(leg,False);

    list_h1 = [];
    for i in range(0, len(arr_chi2max)):
        chi2max = arr_chi2max[i];
        bin1 = h2.GetYaxis().FindBin(0 + 1e-3);
        bin2 = h2.GetYaxis().FindBin(chi2max - 1e-3);
        h1 = h2.ProjectionX("h1_{0:d}".format(i), bin1, bin2, "");
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);
        h1.Draw("E0,same");
        make_common_style(h1, 20, 1.0, colors[i], 2, 0);
        leg.AddEntry(h1, "TOF #chi^{{2}} < {0:2.1f}".format(chi2max), "P");
        list_h1.append(h1);
        if i == 0:
            y_max = h1.GetMaximum();
            print(y_max);

    frame1.GetYaxis().SetRangeUser(y_max * -0.05, y_max * 1.2);
    leg.Draw("");

    txt = TPaveText(0.15, 0.8, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    #txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1, 0.35);
    p2.SetMargin(0.1, 0.02, 0.22, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0, 0.5, 2, 1.05);
    frame2.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame2.GetYaxis().SetTitle("ratio to no cut");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.5);
    frame2.GetXaxis().SetTitleSize(0.09);
    frame2.GetYaxis().SetTitleSize(0.09);
    frame2.GetXaxis().SetLabelSize(0.09);
    frame2.GetYaxis().SetLabelSize(0.09);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    for i in range(1, len(arr_chi2max)):
        h1tmp = list_h1[i];
        h1ratio = h1tmp.Clone("h1tmp_{0:d}".format(i));
        h1ratio.Divide(h1tmp, list_h1[0], 1., 1., "B");
        h1ratio.Draw("E0,same");
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_tofchi2_slice{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_tofchi2_slice{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_tofchi2_slice{1}.png".format(date, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_tofchi2(fliename, taskname, cen1, cen2, suffix):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire_track = rootdire.Get("Track");
    rootdire_track.ls();

    rootdire_event = rootdire.Get("Event/after");
    rootdire_event.ls();
    h1z = rootdire_event.Get("hZvtx");
    nev = h1z.GetEntries();

    h2pos = rootdire_track.Get("positive/hChi2TOF");
    h2neg = rootdire_track.Get("negative/hChi2TOF");
    h2pos.Sumw2();
    h2neg.Sumw2();

    h2 = h2pos.Clone("h2");
    h2.Reset();
    h2.Add(h2pos, 1.);
    h2.Add(h2neg, 1.);
    h2.SetContour(1000);
    h2.Scale(1/nev);
    h2.SetMinimum(1e-8);
    h2.SetMaximum(1e-2);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.12, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0, 0, 10, 10);
    frame1.GetXaxis().SetTitle("p_{pv} (GeV/c)");
    frame1.GetYaxis().SetTitle("TOF #chi^{2}");
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
    h2.Draw("colz,same");

    txt = TPaveText(0.15, 0.85, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_tofchi2{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_tofchi2{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_tofchi2{1}.png".format(date, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_nsigma_tpc_el(fliename, taskname, cen1, cen2, suffix):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire_track = rootdire.Get("Track");
    rootdire_track.ls();

    h2pos = rootdire_track.Get("positive/hTPCNsigmaEl");
    h2neg = rootdire_track.Get("negative/hTPCNsigmaEl");
    h2pos.Sumw2();
    h2neg.Sumw2();

    h2 = h2pos.Clone("h2");
    h2.Reset();
    h2.Add(h2pos, 1.);
    h2.Add(h2neg, 1.);
    h2.SetContour(1000);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.12, 0.12, 0.05);
    c1.SetTicks(1,1);
    c1.SetLogx(1);

    frame1 = c1.DrawFrame(0.1, -5, 10, +5);
    frame1.GetXaxis().SetTitle("p_{TPC inner wall} (GeV/c)");
    frame1.GetYaxis().SetTitle("n #sigma_{e}^{TPC}");
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(1.1);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    txt = TPaveText(0.15, 0.8, 0.4, 0.9, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_nsigma_tpc_el{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_nsigma_tpc_el{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_nsigma_tpc_el{1}.png".format(date, suffix));

    rootfile.Close();

#__________________________________________________________
if __name__ == "__main__":
    filename = "AnalysisResults_HL_281428.root";
    taskname = "single-electron-qc_0010_TPChadrejorTOFreq_TPCshared03_minpt400";
    suffix = "_LHC23_PbPb_apass4_0010";
    #draw_nsigma_tpc_el(filename, taskname, 0, 10, suffix);

    #filename = "AnalysisResults_HL_318815.root";
    #taskname = "single-electron-qc_0010_TOFreq_minpt400";
    #suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_0010_TOFreq_minpt400";
    ##draw_tofchi2(filename, taskname, 0, 10, suffix);
    ##draw_tofchi2_slice(filename, taskname, 0, 10, suffix);
    #draw_its_cluster_size_slice(filename, taskname, 0, 10, suffix);

    #taskname = "single-electron-qc_3050_TOFreq_minpt400";
    #suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_3050_TOFreq_minpt400";
    ##draw_tofchi2_slice(filename, taskname, 30, 50, suffix);
    #draw_its_cluster_size_slice(filename, taskname, 30, 50, suffix);

    #taskname = "single-electron-qc_7090_TOFreq_minpt400";
    #suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_7090_TOFreq_minpt400";
    ##draw_tofchi2_slice(filename, taskname, 70, 90, suffix);
    #draw_its_cluster_size_slice(filename, taskname, 70, 90, suffix);

    filename = "AnalysisResults_HL_319429_544116.root";
    taskname = "event-qc_0010";
    suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_0010_TOFreq_minpt400_HL_319429_544116";
    draw_its_cluster_size_slice(filename, taskname, 0, 10, suffix);
    taskname = "event-qc_3050";
    suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_3050_TOFreq_minpt400_HL_319429_544116";
    draw_its_cluster_size_slice(filename, taskname, 30, 50, suffix);
    taskname = "event-qc_7090";
    suffix = "_LHC23_PbPb_apass4_ft0cocc0_10000_7090_TOFreq_minpt400_HL_319429_544116";
    draw_its_cluster_size_slice(filename, taskname, 70, 90, suffix);



