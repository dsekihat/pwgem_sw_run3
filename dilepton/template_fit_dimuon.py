import datetime
import sys
import math
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText, TLine, TF1, TH1D, TGraphErrors
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray, kViolet
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
from histo_manager import get_ratio, rebin_histogram, convert_dn2n
from signal_extractor import get_significance, rebin_flow_thn
from true_shape_fitter import TrueShapeFitter
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
def template_fit_dca(filename, taskname, cen1, cen2, arr_dca, mmin, mmax, ptmin, ptmax, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_sig.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(mmax - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_sig.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_sig.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    h1m_sig_org = hs_sig.Projection(2);
    h1m_sig = rebin_histogram(h1m_sig_org , arr_dca, False, False);
    h1m_sig.Scale(1, "width");

    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig, False);
    make_common_style(h1m_sig, 20, 1.2, kBlack, 2, 0);

    list_h1dca_photon = [];
    list_h1dca_lf = [];
    list_h1dca_ccbar_uls = [];
    list_h1dca_bbbar_uls = [];
    list_h1dca_bbbar_ls = [];
    list_h1dca_promptjpsi = [];
    list_h1dca_nonpromptjpsi = [];

    rootfile_template = TFile.Open(filename_template, "READ");

    #for photon
    for ip in range(0, len(parnames_photon)):
        parname = parnames_photon[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_photon.append(h1);
    h1dca_photon_org = None;
    if len(list_h1dca_photon) != 0:
        h1dca_photon_org = list_h1dca_photon[0].Clone("h1dca_photon_org");
        for ip in range(1, len(parnames_photon)):
            h1dca_photon_org.Add(list_h1dca_photon[ip]);

    #for lf
    for ip in range(0, len(parnames_lf)):
        parname = parnames_lf[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_lf.append(h1);

    h1dca_lf_org = None;
    if len(list_h1dca_lf) != 0:
        h1dca_lf_org = list_h1dca_lf[0].Clone("h1dca_lf_org");
        for ip in range(1, len(parnames_lf)):
            h1dca_lf_org.Add(list_h1dca_lf[ip]);

    #for ccbar
    for ip in range(0, len(parnames_ccbar_uls)):
        parname = parnames_ccbar_uls[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_ccbar_uls.append(h1);

    h1dca_ccbar_org = None;
    if len(list_h1dca_ccbar_uls) != 0:
        h1dca_ccbar_org = list_h1dca_ccbar_uls[0].Clone("h1dca_ccbar_org");
        for ip in range(1, len(parnames_ccbar_uls)):
            h1dca_ccbar_org.Add(list_h1dca_ccbar_uls[ip], 1);

    #for bbbar
    for ip in range(0, len(parnames_bbbar_uls)):
        parname = parnames_bbbar_uls[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_bbbar_uls.append(h1);

    for ip in range(0, len(parnames_bbbar_ls)):
        parname = parnames_bbbar_ls[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_bbbar_ls.append(h1);

    h1dca_bbbar_org = None;
    if len(list_h1dca_bbbar_uls) != 0:
        h1dca_bbbar_org = list_h1dca_bbbar_uls[0].Clone("h1dca_bbbar_org");
        for ip in range(1, len(parnames_bbbar_uls)):
            h1dca_bbbar_org.Add(list_h1dca_bbbar_uls[ip], 1);
    if len(list_h1dca_bbbar_ls) != 0:
        for ip in range(0, len(parnames_bbbar_ls)):
            h1dca_bbbar_org.Add(list_h1dca_bbbar_ls[ip], -1);

    #for promptjpsi
    for ip in range(0, len(parnames_promptjpsi)):
        parname = parnames_promptjpsi[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_promptjpsi.append(h1);
    h1dca_promptjpsi_org = None;
    if len(list_h1dca_promptjpsi) != 0:
        h1dca_promptjpsi_org = list_h1dca_promptjpsi[0].Clone("h1dca_promptjpsi_org");

    #for nonpromptjpsi
    for ip in range(0, len(parnames_nonpromptjpsi)):
        parname = parnames_nonpromptjpsi[ip];
        hs = rootfile_template.Get("hs_template_{0}".format(parname));
        bin0 = hs.GetAxis(0).FindBin(mmin + 1e-3);
        bin1 = hs.GetAxis(0).FindBin(mmax - 1e-3);
        hs.GetAxis(0).SetRange(bin0, bin1);
        bin0 = hs.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs.GetAxis(1).FindBin(ptmax - 1e-3);
        hs.GetAxis(1).SetRange(bin0, bin1);
        h1 = hs.Projection(2);
        h1.Sumw2();
        h1.SetName("h1dca_{0}".format(parname));
        list_h1dca_nonpromptjpsi.append(h1);
    h1dca_nonpromptjpsi_org = None;
    if len(list_h1dca_nonpromptjpsi) != 0:
        h1dca_nonpromptjpsi_org = list_h1dca_nonpromptjpsi[0].Clone("h1dca_nonpromptjpsi_org");

    h1dca_photon        = None;
    h1dca_lf            = None;
    h1dca_ccbar         = None;
    h1dca_bbbar         = None;
    h1dca_promptjpsi    = None;
    h1dca_nonpromptjpsi = None;

    list_dca_templates = [];

    if len(list_h1dca_photon) != 0:
        h1dca_photon = rebin_histogram(h1dca_photon_org, arr_dca, False, False);
        h1dca_photon.Scale(1, "width");
        h1dca_photon.SetDirectory(0);
        ROOT.SetOwnership(h1dca_photon, False);

    if len(list_h1dca_lf) != 0:
        h1dca_lf            = rebin_histogram(h1dca_lf_org            , arr_dca, False, False);
        h1dca_lf           .Scale(1, "width");
        h1dca_lf           .SetDirectory(0);
        ROOT.SetOwnership(h1dca_lf           , False);
    if len(list_h1dca_ccbar_uls) != 0:
        h1dca_ccbar         = rebin_histogram(h1dca_ccbar_org         , arr_dca, False, False);
    if len(list_h1dca_bbbar_uls) != 0:
        h1dca_bbbar         = rebin_histogram(h1dca_bbbar_org         , arr_dca, False, False);

    if len(list_h1dca_promptjpsi) != 0:
        h1dca_promptjpsi    = rebin_histogram(h1dca_promptjpsi_org     , arr_dca, False, False);
        h1dca_promptjpsi   .Scale(1, "width");
        h1dca_promptjpsi   .SetDirectory(0);
        ROOT.SetOwnership(h1dca_promptjpsi   , False);
    if len(list_h1dca_nonpromptjpsi) != 0:
        h1dca_nonpromptjpsi = rebin_histogram(h1dca_nonpromptjpsi_org, arr_dca, False, False);
        h1dca_nonpromptjpsi .Scale(1, "width");
        ROOT.SetOwnership(h1dca_nonpromptjpsi, False);
        h1dca_nonpromptjpsi.SetDirectory(0);

    if len(list_h1dca_ccbar_uls) != 0:
        h1dca_ccbar        .Scale(1, "width");
        h1dca_ccbar        .SetDirectory(0);
        ROOT.SetOwnership(h1dca_ccbar        , False);
    if len(list_h1dca_bbbar_uls) != 0:
        h1dca_bbbar        .Scale(1, "width");
        h1dca_bbbar        .SetDirectory(0);
        ROOT.SetOwnership(h1dca_bbbar        , False);

    list_dca_templates = [];
    if len(list_h1dca_photon) != 0:
        list_dca_templates.append(h1dca_photon);
    if len(list_h1dca_lf) != 0:
        list_dca_templates.append(h1dca_lf);
    if len(list_h1dca_ccbar_uls) != 0:
        list_dca_templates.append(h1dca_ccbar);
    if len(list_h1dca_bbbar_uls) != 0:
        list_dca_templates.append(h1dca_bbbar);
    if len(list_h1dca_promptjpsi) != 0:
        list_dca_templates.append(h1dca_promptjpsi);
    if len(list_h1dca_nonpromptjpsi) != 0:
        list_dca_templates.append(h1dca_nonpromptjpsi);
    rootfile_template.Close();

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.Divide(1,2,0,0);

    y_range_max = h1m_sig.GetMaximum() * 100;
    y_range_min = max(h1m_sig.GetMinimum() * 0.1, 2e-9);

    str_dca = "3D";
    if "dca3d" in suffix:
        str_dca = "3D";
    elif "dcaxy" in suffix:
        str_dca = "XY";
    elif "dcaz" in suffix:
        str_dca = "Z";

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1, 1);
    p1.SetMargin(0.15, 0.02, 0.0, 0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(arr_dca[0], y_range_min, arr_dca[-1], y_range_max);
    frame1.GetXaxis().SetTitle("DCA_{{#mu#mu}}^{{{0}}} (#sigma)".format(str_dca));
    frame1.GetYaxis().SetTitle("#frac{{1}}{{#it{{N}}_{{ev}}}} #frac{{d#it{{N}}}}{{dDCA_{{#mu#mu}}^{{{0}}}}} (#sigma)^{{#minus1}}".format(str_dca));
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(p1,False);
    ROOT.SetOwnership(frame1,False);

    h1m_sig.Draw("E0same");

    tsf = TrueShapeFitter(list_dca_templates);
    f1dca = TF1("f1dca", tsf, arr_dca[0], arr_dca[-1], len(list_dca_templates));
    f1dca.SetNpx(1000);
    f1dca.SetLineColor(kBlack);
    f1dca.SetLineWidth(2);
    f1dca.SetParLimits(0, 0, 10);
    f1dca.SetParLimits(1, 0, 10);
    f1dca.SetParLimits(2, 0, 10);
    f1dca.SetParLimits(3, 0, 10);
    f1dca.SetParameter(0, 0.5);
    f1dca.SetParameter(1, 0.5);
    f1dca.SetParameter(2, 0.5);
    f1dca.SetParameter(3, 0.5);
    ROOT.SetOwnership(f1dca, False);
    fitter = h1m_sig.Fit(f1dca, "SMEI", "", arr_dca[0], arr_dca[-1]);

    arr_cl68 = np.array(fitter.GetConfidenceIntervals(0.68), dtype=float);
    print(arr_cl68);
    h1cl68 = h1m_sig.Clone("h1cl68");
    h1cl68.Reset();
    h1cl68.SetDirectory(0);
    ROOT.SetOwnership(h1cl68, False);
    make_common_style(h1cl68 , 20, 0.0, kGray+1, 2, 1001);
    for i in range(0, len(arr_cl68)):
        center = h1cl68.GetBinCenter(i+1);
        h1cl68.SetBinContent(i+1, f1dca.Eval(center));
        h1cl68.SetBinError(i+1, arr_cl68[i]);

    #h1dca_p.Scale(sf_prompt);
    #h1dca_np.Scale(sf_nonprompt);
    #make_common_style(h1dca_p , 20, 0.0, kRed+1 , 2, 0);
    #make_common_style(h1dca_np, 20, 0.0, kBlue+1, 2, 0);
    #h1dca_p .SetLineStyle(2);
    #h1dca_np.SetLineStyle(2);
    #h1dca_p .Draw("Hist,same");
    #h1dca_np.Draw("Hist,same");

    h1cl68.Draw("E2,same");
    h1m_sig.Draw("E0same");

    h1cl68.SetDirectory(0);
    ROOT.SetOwnership(h1cl68, False);

    leg = TLegend(0.55, 0.4, 0.75, 0.8);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    leg.AddEntry(h1m_sig, "Data", "P");
    leg.AddEntry(f1dca, "Total fit", "L");
    leg.AddEntry(h1cl68, "Fit uncertainty 68% C.L.", "F");

    sf = [];
    sf_err = [];
    for i in range(0, len(list_dca_templates)):
        sf     = f1dca.GetParameter(i);
        sf_err = f1dca.GetParError(i);
        #print(list_dca_templates[i], sf, sf_err);
        list_dca_templates[i].Scale(sf);
        if "photon" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kViolet+2, 2, 0, 1);
            leg.AddEntry(list_dca_templates[i], "#gamma #rightarrow #mu^{+}#mu^{#minus}", "L");
        elif "lf" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kGreen+2, 2, 0, 1);
            leg.AddEntry(list_dca_templates[i], "LF #rightarrow #mu^{+}#mu^{#minus}", "L");
        elif "ccbar" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kRed+1, 2, 0, 1);
            leg.AddEntry(list_dca_templates[i], "c#bar{c} #rightarrow #mu^{+}#mu^{#minus}", "L");
        elif "bbbar" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kMagenta+1, 2, 0, 1);
            leg.AddEntry(list_dca_templates[i], "b#bar{b} #rightarrow #mu^{+}#mu^{#minus}", "L");
        elif "nonpromptjpsi" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kYellow+1, 2, 0, 2);
            leg.AddEntry(list_dca_templates[i], "Nonprompt J/#psi #rightarrow #mu^{+}#mu^{#minus}", "L");
        elif "promptjpsi" in list_dca_templates[i].GetName():
            make_common_style(list_dca_templates[i], 20, 0.0,  kYellow+1, 2, 0, 1);
            leg.AddEntry(list_dca_templates[i], "Prompt J/#psi #rightarrow #mu^{+}#mu^{#minus}", "L");
        list_dca_templates[i].Draw("Hist,same");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    #n_p = n_p_org * sf_prompt;
    #n_np = n_np_org * sf_nonprompt;
    #frac = n_p/(n_p + n_np);

    #n_p_err = n_p_org * sf_prompt_err;
    #n_np_err = n_np_org * sf_nonprompt_err;
    #frac_err = math.sqrt(math.pow(n_np * n_p_err, 2) + math.pow(n_p * n_np_err,2))/math.pow(n_p + n_np, 2);

    txt_chi2ndf = TPaveText(0.7, 0.35, 0.9, 0.4, "NDC");
    txt_chi2ndf.SetFillColor(kWhite);
    txt_chi2ndf.SetFillStyle(0);
    txt_chi2ndf.SetBorderSize(0);
    txt_chi2ndf.SetTextAlign(12);#middle,left
    txt_chi2ndf.SetTextFont(42);#helvetica
    txt_chi2ndf.SetTextSize(0.05);
    txt_chi2ndf.AddText("#chi^{{2}}/ndf = {0:3.1f}/{1:d}".format(f1dca.GetChisquare(), f1dca.GetNDF()));
    #txt_chi2ndf.AddText("#frac{{N_{{Prompt}}}}{{N_{{Prompt}} + N_{{Nonprompt}}}} = {0:3.2f} #pm {1:3.2f}".format(frac, frac_err));
    txt_chi2ndf.Draw();
    ROOT.SetOwnership(txt_chi2ndf,False);


    txt = TPaveText(0.18, 0.6, 0.4, 0.94, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0}#minus{1}%, Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));

    txt.AddText("#it{p}_{T,#mu} > 0.3 GeV/#it{c}, #minus3.6 < #it{#eta}_{#mu} < #minus2.5");
    txt.AddText("#minus3.6 < #it{y}_{#mu#mu} < #minus2.5");
    txt.AddText("{0:3.2f} < #it{{m}}_{{#mu#mu}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,#mu#mu}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0., 1, 0.35);
    p2.SetMargin(0.15, 0.02, 0.23, 0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(arr_dca[0], 0.4, arr_dca[-1], 2.1);
    frame2.GetXaxis().SetTitle("DCA_{{#mu#mu}}^{{{0}}} (#sigma)".format(str_dca));
    frame2.GetYaxis().SetTitle("#frac{Data}{Fit}");
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.1);
    frame2.GetYaxis().SetTitleSize(0.1);
    frame2.GetXaxis().SetLabelSize(0.1);
    frame2.GetYaxis().SetLabelSize(0.1);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle();
    ROOT.SetOwnership(p2,False);
    ROOT.SetOwnership(frame2,False);

    h1ratio_cl68 = h1cl68.Clone("h1ratio_cl68");
    h1ratio_cl68.Divide(f1dca);
    h1ratio_cl68.Draw("E2same");
    h1ratio_cl68.SetDirectory(0);
    ROOT.SetOwnership(h1ratio_cl68, False);

    line1 = TLine(arr_dca[0], 1, arr_dca[-1], 1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(2);
    line1.SetLineWidth(2);
    line1.Draw("same");
    ROOT.SetOwnership(line1, False);

    h1ratio = h1m_sig.Clone("h1ratio");
    h1ratio.Divide(f1dca);
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_template_fit_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.eps".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_template_fit_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.pdf".format(date, mmin, mmax, ptmin, ptmax, suffix));
    c1.SaveAs("{0}_template_fit_raw_dnddca_m{1:3.2f}_{2:3.2f}GeV_pt{3:2.1f}_{4:2.1f}GeV{5}.png".format(date, mmin, mmax, ptmin, ptmax, suffix));

    rootfile.Close();
#__________________________________________________________
if __name__ == "__main__":
    #filename = "dimuon_pp_5.36TeV_LHC24ap_pass1_HL_389401.root";
    #taskname = "dimuon_global_tight_minpt300";
    #suffix = "_dimuon_pp_13.6TeV_LHC22o_pass7_pT300MeV_dcaxy_HL_389401";
    #filename_template = "output_template_pp_13.6TeV_dcaxy_LHC24f3c_fix_HL_389477.root";

    filename = "dimuon_pp_13.60TeV_LHC23_Thin_pass4_HL_394826.root";
    taskname = "dimuon";
    suffix = "_dimuon_pp_13.6TeV_LHC23_Thin_pass4_pT300MeV_dcaxy_HL_394826";
    filename_template = "output_template_pp_13.6TeV_dcaxy_LHC23k4g_HL_394914.root";

    #arr_dca = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    arr_dca = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);
    #arr_dca = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10], dtype=float);

    parnames_lf = ["eta", "etaprime", "rho", "omega", "phi"];
    parnames_photon = [];
    parnames_ccbar_uls = ["c2l_c2l"];
    parnames_bbbar_uls = ["b2l_b2l", "b2c2l_b2c2l", "b2c2l_b2l_sameb"];
    parnames_bbbar_ls = ["b2c2l_b2l_diffb"];
    parnames_promptjpsi = [];
    parnames_nonpromptjpsi = [];
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 0.21, 0.5, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 0.70, 0.84, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 0.94, 1.08, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);

    parnames_lf = [];
    parnames_photon = [];
    parnames_ccbar_uls = ["c2l_c2l"];
    parnames_bbbar_uls = ["b2l_b2l", "b2c2l_b2c2l", "b2c2l_b2l_sameb"];
    parnames_bbbar_ls = ["b2c2l_b2l_diffb"];
    parnames_promptjpsi = ["promptjpsi"];
    parnames_nonpromptjpsi = ["nonpromptjpsi"];
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 2.7, 3.2, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 2.7, 3.2, 2.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);

    parnames_lf = [];
    parnames_photon = [];
    parnames_ccbar_uls = ["c2l_c2l"];
    parnames_bbbar_uls = ["b2l_b2l", "b2c2l_b2c2l", "b2c2l_b2l_sameb"];
    parnames_bbbar_ls = ["b2c2l_b2l_diffb"];
    parnames_promptjpsi = [];
    parnames_nonpromptjpsi = [];
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 1.1, 2.7, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);

    parnames_lf = [];
    parnames_photon = [];
    parnames_ccbar_uls = ["c2l_c2l"];
    parnames_bbbar_uls = ["b2l_b2l", "b2c2l_b2c2l", "b2c2l_b2l_sameb"];
    parnames_bbbar_ls = ["b2c2l_b2l_diffb"];
    parnames_promptjpsi = [];
    parnames_nonpromptjpsi = [];
    template_fit_dca(filename, taskname, 0, 999, arr_dca, 4.0, 8.0, 0.0, 10, filename_template, parnames_photon, parnames_lf, parnames_ccbar_uls, parnames_bbbar_uls, parnames_bbbar_ls, parnames_promptjpsi, parnames_nonpromptjpsi, suffix);

#__________________________________________________________
#__________________________________________________________
#__________________________________________________________

