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
from histo_manager import rebin_histogram
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
def draw_cospa(filename, taskname, pairname, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();
    #roottask.Get("{0}".format(pairname)).ls();

    colors = [kBlack, kRed+1, kBlue+1];
    markers = [20, 24];

    parnames = ["all", "D0", "Dpm"];

    leg_strings = [
        "all eK",
        "eK from D^{0} and c.c.",
        "eK from D^{+} and c.c.",
        "eK from #Lambda_{c}^{+} and c.c.",
    ];
    
    list_h1eff = [];
    y_max = 100;
    str_mass = "";
    if "e_Kpm" in pairname:
        str_mass = "m_{e^{+}K^{#minus}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";

    for ipar, parname in enumerate(parnames):
        h1 = roottask.Get("{0}/{1}/hCosPA".format(pairname, parname));
        h1.Sumw2();
        make_common_style(h1, markers[0], 1.2, colors[ipar], 2, 0);
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);
        list_h1eff.append(h1);
        if ipar == 0:
            y_max = h1.GetMaximum() * 2;

    leg = TLegend(0.15, 0.55, 0.35, 0.70);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.03, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.9, 1, 1.0, y_max);
    frame1.GetXaxis().SetTitle("cosine of pointing angle");
    frame1.GetYaxis().SetTitle("counts");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ie in range(0, len(list_h1eff)):
        list_h1eff[ie].Draw("same,E0");
        leg.AddEntry(list_h1eff[ie], leg_strings[ie], "P");

    leg.Draw("");
    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("#it{{p}}_{{T,h_{{s}}}} > {0:2.1f} GeV/#it{{c}}, |#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.AddText("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(0.5));
    txt.AddText("|#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_cosPA{2}.eps".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_cosPA{2}.pdf".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_cosPA{2}.png".format(date, pairname, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_dca2legs(filename, taskname, pairname, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();
    #roottask.Get("{0}".format(pairname)).ls();

    colors = [kBlack, kRed+1, kBlue+1];
    markers = [20, 24];

    parnames = ["all", "D0", "Dpm"];

    leg_strings = [
        "all eK",
        "eK from D^{0} and c.c.",
        "eK from D^{+} and c.c.",
        "eK from #Lambda_{c}^{+} and c.c.",
    ];
    
    list_h1eff = [];
    y_max = 100;
    str_mass = "";
    if "e_Kpm" in pairname:
        str_mass = "m_{e^{+}K^{#minus}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";

    for ipar, parname in enumerate(parnames):
        h1 = roottask.Get("{0}/{1}/hDCA2Legs".format(pairname, parname));
        h1.Sumw2();
        make_common_style(h1, markers[0], 1.2, colors[ipar], 2, 0);
        h1.SetDirectory(0);
        ROOT.SetOwnership(h1, False);
        list_h1eff.append(h1);
        if ipar == 0:
            y_max = h1.GetMaximum() * 2;

    leg = TLegend(0.15, 0.55, 0.35, 0.70);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.03, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.0, 1, 0.5, y_max);
    frame1.GetXaxis().SetTitle("DCA between 2 legs (cm)");
    frame1.GetYaxis().SetTitle("counts");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ie in range(0, len(list_h1eff)):
        list_h1eff[ie].Draw("same,E0");
        leg.AddEntry(list_h1eff[ie], leg_strings[ie], "P");

    leg.Draw("");
    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("#it{{p}}_{{T,h_{{s}}}} > {0:2.1f} GeV/#it{{c}}, |#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.AddText("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(0.5));
    txt.AddText("|#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_distance_between_2legs{2}.eps".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_distance_between_2legs{2}.pdf".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_distance_between_2legs{2}.png".format(date, pairname, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_mass(filename, taskname, pairname, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();
    #roottask.Get("{0}".format(pairname)).ls();

    colors = [kBlack, kRed+1, kBlue+1];
    markers = [20, 24];

    parnames = ["all", "D0", "Dpm"];

    leg_strings = [
        "all eK",
        "eK from D^{0} and c.c.",
        "eK from D^{+} and c.c.",
        "eK from #Lambda_{c}^{+} and c.c.",
    ];
    
    list_h1eff = [];
    y_max = 100;
    str_mass = "";
    if "e_Kpm" in pairname:
        str_mass = "m_{e^{+}K^{#minus}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";
    elif "e_K0S" in pairname:
        str_mass = "m_{e^{+}#bar{K^{0}_{S}}} and c.c. (GeV/#it{c}^{2})";

    for ipar, parname in enumerate(parnames):
        h2 = roottask.Get("{0}/{1}/hMass".format(pairname, parname));
        h2.Sumw2();
        h1m  = h2.ProjectionX("h1m_{0}".format(parname));
        make_common_style(h1m, markers[0], 1.2, colors[ipar], 2, 0);
        h1m.SetDirectory(0);
        ROOT.SetOwnership(h1m, False);
        list_h1eff.append(h1m);
        if ipar == 0:
            y_max = h1m.GetMaximum() * 2;

    leg = TLegend(0.15, 0.55, 0.35, 0.70);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.5, 1, 2, y_max);
    frame1.GetXaxis().SetTitle(str_mass);
    frame1.GetYaxis().SetTitle("counts");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ie in range(0, len(list_h1eff)):
        list_h1eff[ie].Draw("same,E0");
        leg.AddEntry(list_h1eff[ie], leg_strings[ie], "P");

    leg.Draw("");
    txt = TPaveText(0.15, 0.70, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("#it{{p}}_{{T,h_{{s}}}} > {0:2.1f} GeV/#it{{c}}, |#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.AddText("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(0.5));
    txt.AddText("|#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_{1}_mass{2}.eps".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_mass{2}.pdf".format(date, pairname, suffix));
    c1.SaveAs("{0}_{1}_mass{2}.png".format(date, pairname, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_tagging_efficiency(filename, taskname, etamax_e, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    roottask.ls();

    arr_pte = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6,7,8,9,10], dtype=float);
    hs_gen_prompt_d0 = roottask.Get("Generated/D0/prompt/hs");
    hs_gen_prompt_d0.Sumw2();
    ngen = hs_gen_prompt_d0.GetEntries();
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1];
    markers = [20, 24];

    parnames = ["D0", "Dpm", "Ds", "Lc"];
    sourcetypes = ["prompt", "nonprompt"];

    leg_strings = [
        "e from D^{0} #rightarrow e^{+}#nu_{e}K^{#minus}, e^{+}#nu_{e}K^{#minus}#pi^{0}, e^{+}#nu_{e}#bar{K^{0}_{S}}#pi^{#minus} and c.c.",
        "e from D^{+} #rightarrow e^{+}#nu_{e}#bar{K^{0}_{S}}, e^{+}#nu_{e}K^{#minus}#pi^{+} and c.c.",
        "e from D^{+}_{s} #rightarrow e^{+}#nu_{e}#phi with #phi #rightarrow K^{+}K^{#minus} and c.c.",
        "e from #Lambda_{c}^{+} #rightarrow e^{+}#nu_{e}#Lambda and c.c.",
    ];
    
    list_h1eff = [];
    #for ipar, parname in enumerate(parnames):
    #    for isource, stype in enumerate(sourcetypes):
    #        print(parname, stype);
    #        hs_findable = roottask.Get("{0}/electron/{1}/findable/hs".format(parname, stype));
    #        hs_correct = roottask.Get("{0}/electron/{1}/correct/hs".format(parname, stype));
    #        h1pt_findable = hs_findable.Projection(0);
    #        h1pt_correct = hs_correct.Projection(0);
    #        h1eff = h1pt_correct.Clone("h1eff_{0}_{1}".format(stype, parname));
    #        h1eff.Reset();
    #        h1eff.Divide(h1pt_correct, h1pt_findable, 1., 1., "B");
    #        make_common_style(h1eff, markers[isource], 1.2, colors[ipar], 2, 0);
    #        h1eff.SetDirectory(0);
    #        ROOT.SetOwnership(h1eff, False);
    #        list_h1eff.append(h1eff);

    for ipar, parname in enumerate(parnames):
        hs_prompt_findable = roottask.Get("{0}/electron/prompt/findable/hs".format(parname));
        hs_prompt_correct  = roottask.Get("{0}/electron/prompt/correct/hs".format(parname));
        hs_nonprompt_findable = roottask.Get("{0}/electron/nonprompt/findable/hs".format(parname));
        hs_nonprompt_correct  = roottask.Get("{0}/electron/nonprompt/correct/hs".format(parname));

        hs_prompt_findable.Add(hs_nonprompt_findable, 1);
        hs_prompt_correct.Add(hs_nonprompt_correct, 1);

        h1pt_findable = hs_prompt_findable.Projection(0);
        h1pt_correct  = hs_prompt_correct.Projection(0);
        h1eff = h1pt_correct.Clone("h1eff_{0}".format(parname));
        h1eff.Reset();
        h1eff.Divide(h1pt_correct, h1pt_findable, 1., 1., "B");
        make_common_style(h1eff, markers[0], 1.2, colors[ipar], 2, 0);
        h1eff.SetDirectory(0);
        ROOT.SetOwnership(h1eff, False);
        list_h1eff.append(h1eff);

    leg = TLegend(0.15, 0.55, 0.35, 0.75);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.SetHeader("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(etamax_e));
    ROOT.SetOwnership(leg, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    frame1 = c1.DrawFrame(0.0, 0, 10, 1.0);
    frame1.GetXaxis().SetTitle("#it{p}_{T,e} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("tagging efficiency");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ie in range(0, len(list_h1eff)):
        list_h1eff[ie].Draw("same,E1");
        leg.AddEntry(list_h1eff[ie], leg_strings[ie], "P");

    leg.Draw("");
    txt = TPaveText(0.15, 0.75, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#it{{p}}_{{T,h_{{s}}}} > {0:2.1f} GeV/#it{{c}}, |#it{{#eta}}_{{h_{{s}}}}| < {1:2.1f}".format(0.1, 1.2));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_tagging_eff{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_tagging_eff{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_tagging_eff{1}.png".format(date, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_conditional_acceptance(filename, taskname, pt_min_k, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();

    arr_pte = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6,7,8,9,10], dtype=float);
    hs_gen_prompt_d0 = roottask.Get("Generated/D0/prompt/hs");
    hs_gen_prompt_d0.Sumw2();
    ngen = hs_gen_prompt_d0.GetEntries();
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1];
    markers = [20, 20, 21, 21, 24, 24];

    arr_eta_max_e = np.array([0.5, 0.5, 0.8, 0.8, 1.0, 4.0], dtype=float);
    arr_eta_max_k = np.array([0.8, 1.2, 0.8, 1.2, 4.0, 4.0], dtype=float);
    list_h1ratio = [];

    leg = TLegend(0.6, 0.15, 0.8, 0.45);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.SetHeader("#it{{p}}_{{T,K}} > {0:3.2f} GeV/#it{{c}}".format(pt_min_k));
    ROOT.SetOwnership(leg,False);

    for ieta in range(0, len(arr_eta_max_e)):
        hs_gen_prompt_d0.GetAxis(0).SetRange(0, 0); #reset
        hs_gen_prompt_d0.GetAxis(1).SetRange(0, 0); #reset
        hs_gen_prompt_d0.GetAxis(2).SetRange(0, 0); #reset
        hs_gen_prompt_d0.GetAxis(3).SetRange(0, 0); #reset

        eta_max_e = arr_eta_max_e[ieta];
        eta_max_k = arr_eta_max_k[ieta];

        bin0_e = hs_gen_prompt_d0.GetAxis(2).FindBin(-1 * eta_max_e + 1e-3);
        bin1_e = hs_gen_prompt_d0.GetAxis(2).FindBin(+1 * eta_max_e - 1e-3);
        hs_gen_prompt_d0.GetAxis(2).SetRange(bin0_e, bin1_e);
        h1pt_denominator_org = hs_gen_prompt_d0.Projection(0);
        h1pt_denominator = rebin_histogram(h1pt_denominator_org, arr_pte, False, False);

        bin0_k = hs_gen_prompt_d0.GetAxis(1).FindBin(pt_min_k + 1e-3);
        bin1_k = hs_gen_prompt_d0.GetAxis(1).FindBin(100      - 1e-3);
        hs_gen_prompt_d0.GetAxis(1).SetRange(bin0_k, bin1_k);
        bin0_k = hs_gen_prompt_d0.GetAxis(3).FindBin(-1 * eta_max_k + 1e-3);
        bin1_k = hs_gen_prompt_d0.GetAxis(3).FindBin(+1 * eta_max_k - 1e-3);
        hs_gen_prompt_d0.GetAxis(3).SetRange(bin0_k, bin1_k);
        h1pt_numerator_org = hs_gen_prompt_d0.Projection(0);
        h1pt_numerator_org = hs_gen_prompt_d0.Projection(0);
        h1pt_numerator = rebin_histogram(h1pt_numerator_org, arr_pte, False, False);

        h1ratio = h1pt_numerator.Clone("h1ratio_{0:d}".format(ieta));
        h1ratio.Reset();
        h1ratio.Divide(h1pt_numerator, h1pt_denominator, 1., 1., "B");
        make_common_style(h1ratio, markers[ieta], 1.2, colors[ieta], 2, 0);
        list_h1ratio.append(h1ratio);
        h1ratio.SetDirectory(0);
        ROOT.SetOwnership(h1ratio, False);
        leg.AddEntry(h1ratio, "|#it{{#eta}}_{{e}}| < {0:2.1f}, |#it{{#eta}}_{{K}}| < {1:2.1f}".format(eta_max_e, eta_max_k), "P");

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.11, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(0.0, 0, 10, 1);
    frame1.GetXaxis().SetTitle("#it{p}_{T,e} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("conditional acceptance");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ieta in range(0, len(arr_eta_max_e)):
        list_h1ratio[ieta].Draw("E0,same");

    leg.Draw("");
    txt = TPaveText(0.12, 0.75, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("D^{0} #rightarrow e^{+}#nu_{e}K^{#minus} and c.c.");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_conditional_acc_promptD0_minptK{1:3.2f}GeV{2}.eps".format(date, pt_min_k, suffix));
    c1.SaveAs("{0}_conditional_acc_promptD0_minptK{1:3.2f}GeV{2}.pdf".format(date, pt_min_k, suffix));
    c1.SaveAs("{0}_conditional_acc_promptD0_minptK{1:3.2f}GeV{2}.png".format(date, pt_min_k, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_eta_correlation(filename, taskname, ptmin_e, ptmin_k, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();

    hs_gen_prompt_d0 = roottask.Get("Generated/D0/prompt/hs");
    ngen = hs_gen_prompt_d0.GetEntries();

    bin0_e = hs_gen_prompt_d0.GetAxis(0).FindBin(ptmin_e + 1e-3);
    bin0_k = hs_gen_prompt_d0.GetAxis(1).FindBin(ptmin_k + 1e-3);
    hs_gen_prompt_d0.GetAxis(0).SetRange(bin0_e, 10000);
    hs_gen_prompt_d0.GetAxis(1).SetRange(bin0_k, 10000);

    h2 = hs_gen_prompt_d0.Projection(3, 2);
    h2.SetContour(1000);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.Scale(1/ngen);
    h2.GetZaxis().SetTitle("1/#it{N}_{D^{0}}");
    h2.GetZaxis().SetTitleOffset(1.3);
    h2.GetXaxis().SetTitle("#it{#eta}_{e}");
    h2.GetYaxis().SetTitle("#it{#eta}_{K}");
    h2.GetXaxis().SetTitleOffset(1.0);
    h2.GetYaxis().SetTitleOffset(0.9);
    h2.GetXaxis().SetTitleSize(0.04);
    h2.GetYaxis().SetTitleSize(0.04);
    h2.GetXaxis().SetLabelSize(0.04);
    h2.GetYaxis().SetLabelSize(0.04);
    h2.GetXaxis().SetLabelOffset(0.01);
    h2.GetYaxis().SetLabelOffset(0.01);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.08, 0.15, 0.1, 0.04);
    c1.SetTicks(1,1);

    frame1 = c1.DrawFrame(-5, -5, +5, +5);
    frame1.GetXaxis().SetTitle("#it{#eta}_{e}");
    frame1.GetYaxis().SetTitle("#it{#eta}_{K}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz");

    txt = TPaveText(0.1, 0.64, 0.4, 0.94, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("D^{0} #rightarrow e^{+}#nu_{e}K^{#minus} and c.c.");
    txt.AddText("#it{{p}}_{{T,e}} > {0:3.2f} GeV/#it{{c}}".format(ptmin_e));
    txt.AddText("#it{{p}}_{{T,K}} > {0:3.2f} GeV/#it{{c}}".format(ptmin_k));
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_eta_gen_ptminE{1:3.2f}_ptminK{2:3.2f}GeV_promptD0{3}.eps".format(date, ptmin_e, ptmin_k, suffix));
    c1.SaveAs("{0}_eta_gen_ptminE{1:3.2f}_ptminK{2:3.2f}GeV_promptD0{3}.pdf".format(date, ptmin_e, ptmin_k, suffix));
    c1.SaveAs("{0}_eta_gen_ptminE{1:3.2f}_ptminK{2:3.2f}GeV_promptD0{3}.png".format(date, ptmin_e, ptmin_k, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_pt_correlation(filename, taskname, etamax_e, etamax_k, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();

    hs_gen_prompt_d0 = roottask.Get("Generated/D0/prompt/hs");
    ngen = hs_gen_prompt_d0.GetEntries();

    bin0_e = hs_gen_prompt_d0.GetAxis(2).FindBin(-1 * etamax_e + 1e-3);
    bin0_k = hs_gen_prompt_d0.GetAxis(3).FindBin(-1 * etamax_k + 1e-3);
    bin1_e = hs_gen_prompt_d0.GetAxis(2).FindBin(1 * etamax_e - 1e-3);
    bin1_k = hs_gen_prompt_d0.GetAxis(3).FindBin(1 * etamax_k - 1e-3);

    hs_gen_prompt_d0.GetAxis(2).SetRange(bin0_e, bin1_e);
    hs_gen_prompt_d0.GetAxis(3).SetRange(bin0_k, bin1_k);

    h2 = hs_gen_prompt_d0.Projection(1,0);
    h2.SetContour(1000);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.Scale(1/ngen);
    h2.GetZaxis().SetTitle("1/#it{N}_{D^{0}}");
    h2.GetZaxis().SetTitleOffset(1.3);
    h2.GetXaxis().SetTitle("#it{p}_{T,e} (GeV/#it{c})");
    h2.GetYaxis().SetTitle("#it{p}_{T,K} (GeV/#it{c})");
    h2.GetXaxis().SetTitleOffset(1.1);
    h2.GetYaxis().SetTitleOffset(1.1);
    h2.GetXaxis().SetTitleSize(0.04);
    h2.GetYaxis().SetTitleSize(0.04);
    h2.GetXaxis().SetLabelSize(0.04);
    h2.GetYaxis().SetLabelSize(0.04);
    h2.GetXaxis().SetLabelOffset(0.01);
    h2.GetYaxis().SetLabelOffset(0.01);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.15, 0.1, 0.04);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(0, 0, 10, 10);
    frame1.GetXaxis().SetTitle("#it{p}_{T,e}");
    frame1.GetYaxis().SetTitle("#it{p}_{T,K}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz");

    txt = TPaveText(0.5, 0.7, 0.8, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("D^{0} #rightarrow e^{+}#nu_{e}K^{#minus} and c.c.");
    txt.AddText("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(etamax_e));
    txt.AddText("|#it{{#eta}}_{{K}}| < {0:2.1f}".format(etamax_k));
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pt_gen_etamaxE{1:2.1f}_etamaxnK{2:2.1f}_promptD0{3}.eps".format(date, etamax_e, etamax_k, suffix));
    c1.SaveAs("{0}_pt_gen_etamaxE{1:2.1f}_etamaxnK{2:2.1f}_promptD0{3}.pdf".format(date, etamax_e, etamax_k, suffix));
    c1.SaveAs("{0}_pt_gen_etamaxE{1:2.1f}_etamaxnK{2:2.1f}_promptD0{3}.png".format(date, etamax_e, etamax_k, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
def draw_pt_eta_correlation_kaon(filename, taskname, ptmin_e, etamax_e, suffix=""):
    rootfile = TFile.Open(filename, "READ");
    roottask = rootfile.Get(taskname);
    #roottask.ls();

    hs_gen_prompt_d0 = roottask.Get("Generated/D0/prompt/hs");
    ngen = hs_gen_prompt_d0.GetEntries();

    bin0_e = hs_gen_prompt_d0.GetAxis(0).FindBin(ptmin_e + 1e-3);
    bin1_e = hs_gen_prompt_d0.GetAxis(0).FindBin(100 - 1e-3);
    hs_gen_prompt_d0.GetAxis(0).SetRange(bin0_e, bin1_e);

    bin0_e = hs_gen_prompt_d0.GetAxis(2).FindBin(-1 * etamax_e + 1e-3);
    bin1_e = hs_gen_prompt_d0.GetAxis(2).FindBin(1 * etamax_e - 1e-3);
    hs_gen_prompt_d0.GetAxis(2).SetRange(bin0_e, bin1_e);

    h2 = hs_gen_prompt_d0.Projection(1,3);
    h2.SetContour(1000);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.Scale(1/ngen);
    h2.GetZaxis().SetTitle("1/#it{N}_{D^{0}}");
    h2.GetZaxis().SetTitleOffset(1.3);
    h2.GetXaxis().SetTitle("#it{#eta}_{K}");
    h2.GetYaxis().SetTitle("#it{p}_{T,K} (GeV/#it{c})");
    h2.GetXaxis().SetTitleOffset(1.1);
    h2.GetYaxis().SetTitleOffset(1.1);
    h2.GetXaxis().SetTitleSize(0.04);
    h2.GetYaxis().SetTitleSize(0.04);
    h2.GetXaxis().SetLabelSize(0.04);
    h2.GetYaxis().SetLabelSize(0.04);
    h2.GetXaxis().SetLabelOffset(0.01);
    h2.GetYaxis().SetLabelOffset(0.01);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.1, 0.15, 0.1, 0.04);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-5, 0, +5, 10);
    frame1.GetXaxis().SetTitle("#it{p}_{T,e}");
    frame1.GetYaxis().SetTitle("#it{p}_{T,K}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.0);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz");

    txt = TPaveText(0.12, 0.64, 0.4, 0.94, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");
    txt.AddText("LHC23k4g_small");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("D^{0} #rightarrow e^{+}#nu_{e}K^{#minus} and c.c.");
    txt.AddText("#it{{p}}_{{T,e}} > {0:3.2f} GeV/#it{{c}}".format(ptmin_e));
    txt.AddText("|#it{{#eta}}_{{e}}| < {0:2.1f}".format(etamax_e));
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pt_eta_ptminE{1:3.2f}_etamaxE{2:2.1f}_promptD0{3}.eps".format(date, ptmin_e, etamax_e, suffix));
    c1.SaveAs("{0}_pt_eta_ptminE{1:3.2f}_etamaxE{2:2.1f}_promptD0{3}.pdf".format(date, ptmin_e, etamax_e, suffix));
    c1.SaveAs("{0}_pt_eta_ptminE{1:3.2f}_etamaxE{2:2.1f}_promptD0{3}.png".format(date, ptmin_e, etamax_e, suffix));

    rootfile.Close();
#__________________________________________________________
if __name__ == "__main__":
    filename = "AnalysisResults_HL_409020.root";
    taskname = "tagging-hfe_all";
    #draw_pt_eta_correlation_kaon(filename, taskname, 0.2, 0.8, "_HL_409020");
    #draw_pt_eta_correlation_kaon(filename, taskname, 0.2, 0.5, "_HL_409020");
    #draw_pt_eta_correlation_kaon(filename, taskname, 0.4, 0.8, "_HL_409020");
    #draw_pt_eta_correlation_kaon(filename, taskname, 0.1, 4.0, "_HL_409020");
    #draw_pt_correlation(filename, taskname, 0.8, 0.8, "_HL_409020");
    #draw_eta_correlation(filename, taskname, 0.1, 0.1, "_HL_409020");
    #draw_eta_correlation(filename, taskname, 0.2, 0.2, "_HL_409020");
    #draw_conditional_acceptance(filename, taskname, 0.1, "_HL_409020");
    #draw_conditional_acceptance(filename, taskname, 0.2, "_HL_409020");
    #draw_conditional_acceptance(filename, taskname, 0.4, "_HL_409020");
    draw_tagging_efficiency(filename, taskname, 0.5, "_HL_409020");
    #draw_mass(filename, taskname, "e_Kpm", "_HL_409020");
    #draw_dca2legs(filename, taskname, "e_Kpm", "_HL_409020");
    #draw_cospa(filename, taskname, "e_Kpm", "_HL_409020");
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________

