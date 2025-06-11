import sys
sys.path.append("../common/");
import numpy as np
import datetime
import ROOT
from ROOT import TFile, TList, TCanvas, TPaveText, gStyle
from ROOT import kWhite, kBlack
from histo_manager import rebin_histogram_2d
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
def extract(filename, taskname, parname, outname, arr_pt):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    h2pt_org = rootdire.Get("Ptgen_RelDeltaPt");
    h2eta_org = rootdire.Get("Ptgen_DeltaEta");
    h2phi_pos_org = rootdire.Get("Ptgen_DeltaPhi_Pos");
    h2phi_neg_org = rootdire.Get("Ptgen_DeltaPhi_Neg");

    arr_rel_reso_pt = np.linspace(h2pt_org.GetYaxis().GetXmin(), h2pt_org.GetYaxis().GetXmax(), h2pt_org.GetYaxis().GetNbins() + 1);
    arr_reso_eta = np.linspace(h2eta_org.GetYaxis().GetXmin(), h2eta_org.GetYaxis().GetXmax(), h2eta_org.GetYaxis().GetNbins() + 1);
    arr_reso_phi_pos = np.linspace(h2phi_pos_org.GetYaxis().GetXmin(), h2phi_pos_org.GetYaxis().GetXmax(), h2phi_pos_org.GetYaxis().GetNbins() + 1);
    arr_reso_phi_neg = np.linspace(h2phi_neg_org.GetYaxis().GetXmin(), h2phi_neg_org.GetYaxis().GetXmax(), h2phi_neg_org.GetYaxis().GetNbins() + 1);
    #print(arr_rel_reso_pt);
    #print(arr_reso_eta);
    #print(arr_reso_phi_pos);
    #print(arr_reso_phi_neg);

    h2pt      = rebin_histogram_2d(h2pt_org     , arr_pt, arr_rel_reso_pt);
    h2eta     = rebin_histogram_2d(h2eta_org    , arr_pt, arr_reso_eta);
    h2phi_pos = rebin_histogram_2d(h2phi_pos_org, arr_pt, arr_reso_phi_pos);
    h2phi_neg = rebin_histogram_2d(h2phi_neg_org, arr_pt, arr_reso_phi_neg);

    h2pt     .SetName("Ptgen_RelDeltaPt");
    h2eta    .SetName("Ptgen_DeltaEta");
    h2phi_pos.SetName("Ptgen_DeltaPhi_Pos");
    h2phi_neg.SetName("Ptgen_DeltaPhi_Neg");

    outfile = TFile(outname, "RECREATE");
    outlist = TList();
    outlist.SetName("ccdb_object"); #Don't touch.
    outlist.SetOwner(True);
    outlist.Add(h2pt )
    outlist.Add(h2eta)
    outlist.Add(h2phi_pos)
    outlist.Add(h2phi_neg)
    outfile.WriteTObject(outlist);
    outfile.Close();
    outlist.Clear();

    rootfile.Close();
#__________________________________________________________
def extract_hs(filename, taskname, parname, outname, arr_pt):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    outfile = TFile(outname, "RECREATE");
    outlist = TList();
    outlist.SetName("ccdb_object"); #Don't touch.
    outlist.SetOwner(True);
    outlist.Add(hs_reso)
    outfile.WriteTObject(outlist);
    outfile.Close();
    outlist.Clear();

    rootfile.Close();
#__________________________________________________________
def check_ptgen_vs_etagen_deltaeta(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h2 = hs_reso.Projection(5, 1);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.SetContour(1000);

    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.12, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-1.5, -1, +1.5, +1);
    frame1.GetXaxis().SetTitle("#eta_{l}^{gen}");
    frame1.GetYaxis().SetTitle("#eta_{l}^{gen} #minus #eta_{l}^{rec}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    sign_str = "pos";
    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_etagen_deta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etagen_deta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etagen_deta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
def check_ptgen_vs_etagen_deltapt(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h2 = hs_reso.Projection(4, 1);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.SetContour(1000);

    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.12, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-1.5, -1, +1.5, +1);
    frame1.GetXaxis().SetTitle("#eta_{l}^{gen}");
    frame1.GetYaxis().SetTitle("(p_{T,l}^{gen} #minus p_{T,l}^{rec})/p_{T,l}^{gen}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    sign_str = "pos";
    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_etagen_dpt_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etagen_dpt_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etagen_dpt_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
#__________________________________________________________
def check_ptgen_etagen_vs_pteta(filename, taskname, parname, minptgen, maxptgen, minetagen, maxetagen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);

    bin1 = hs_reso.GetAxis(1).FindBin(minetagen + 1e-3);
    bin2 = hs_reso.GetAxis(1).FindBin(maxetagen - 1e-3);
    hs_reso.GetAxis(1).SetRange(bin1, bin2);

    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h2 = hs_reso.Projection(5, 4);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.SetContour(1000);

    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.12, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-1, -0.1, +1, +0.1);
    frame1.GetXaxis().SetTitle("(p_{T,l}^{gen} #minus p_{T,l}^{rec})/p_{T,l}^{gen}");
    frame1.GetYaxis().SetTitle("#eta_{l}^{gen} #minus #eta_{l}^{rec}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    sign_str = "pos";
    txt = TPaveText(0.16, 0.7, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    txt.AddText("{0:2.1f} < #eta_{{l}}^{{gen}} < {1:2.1f}".format(minetagen, maxetagen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_etagen{3:2.1f}_{4:2.1f}_{5}{6}.eps".format(date, minptgen, maxptgen, minetagen, maxetagen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_etagen{3:2.1f}_{4:2.1f}_{5}{6}.pdf".format(date, minptgen, maxptgen, minetagen, maxetagen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_etagen{3:2.1f}_{4:2.1f}_{5}{6}.png".format(date, minptgen, maxptgen, minetagen, maxetagen, sign_str, suffix));

#__________________________________________________________
def check_ptgen_vs_pteta(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h2 = hs_reso.Projection(5, 4);
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);
    h2.SetContour(1000);

    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.12, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-1, -0.1, +1, +0.1);
    frame1.GetXaxis().SetTitle("(p_{T,l}^{gen} #minus p_{T,l}^{rec})/p_{T,l}^{gen}");
    frame1.GetYaxis().SetTitle("#eta_{l}^{gen} #minus #eta_{l}^{rec}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    sign_str = "pos";
    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_pteta_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
def check_ptgen_vs_ptphi(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);
    hs_reso.Sumw2();

    h2 = hs_reso.Projection(6, 4);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    ROOT.SetOwnership(h2, False);
    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.12, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-1, -0.1, +1, +0.1);
    frame1.GetXaxis().SetTitle("(p_{T,l}^{gen} #minus p_{T,l}^{rec})/p_{T,l}^{gen}");
    frame1.GetYaxis().SetTitle("#varphi_{l}^{gen} #minus #varphi_{l}^{rec} (rad.)");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    sign_str = "pos";

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_ptphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_ptphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_ptphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
def check_ptgen_vs_etaphi(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h2 = hs_reso.Projection(6, 5);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    ROOT.SetOwnership(h2, False);
    h2.Scale(1/h2.GetEntries());
    h2.Scale(1, "width");
    max_entry = h2.GetMaximum();
    h2.SetMaximum(max_entry * 1.1);
    h2.SetMinimum(max_entry * 1e-6);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.13, 0.11, 0.10, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogz(1);

    frame1 = c1.DrawFrame(-0.1, -0.1, +0.1, +0.1);
    frame1.GetXaxis().SetTitle("#eta_{l}^{gen} #minus #eta_{l}^{rec}");
    frame1.GetYaxis().SetTitle("#varphi_{l}^{gen} #minus #varphi_{l}^{rec} (rad.)");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h2.Draw("colz,same");

    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    sign_str = "pos";

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_etaphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etaphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_etaphi_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
def check_ptgen_vs_ptetaphi(filename, taskname, parname, minptgen, maxptgen, sign, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("{0}/{1}".format(taskname, parname));
    rootdire.ls();

    hs_reso = rootdire.Get("hs_reso");
    bin1 = hs_reso.GetAxis(0).FindBin(minptgen + 1e-3);
    bin2 = hs_reso.GetAxis(0).FindBin(maxptgen - 1e-3);
    hs_reso.GetAxis(0).SetRange(bin1, bin2);
    bin1 = hs_reso.GetAxis(3).FindBin(sign + 1e-3);
    bin2 = hs_reso.GetAxis(3).FindBin(sign - 1e-3);
    hs_reso.GetAxis(3).SetRange(bin1, bin2);

    h3 = hs_reso.Projection(4, 5, 6);
    h3.SetDirectory(0);
    #h3.SetContour(1000);
    ROOT.SetOwnership(h3, False);
    h3.Scale(1/h3.GetEntries());
    h3.Scale(1, "width");

    h3.GetXaxis().SetTitleOffset(1.9);
    h3.GetYaxis().SetTitleOffset(2.5);
    h3.GetZaxis().SetTitleOffset(2.0);
    h3.GetXaxis().SetTitleSize(0.035);
    h3.GetYaxis().SetTitleSize(0.035);
    h3.GetZaxis().SetTitleSize(0.035);
    h3.GetXaxis().SetLabelSize(0.035);
    h3.GetYaxis().SetLabelSize(0.035);
    h3.GetZaxis().SetLabelSize(0.035);
    h3.GetXaxis().SetLabelOffset(0.01);
    h3.GetYaxis().SetLabelOffset(0.01);
    h3.GetZaxis().SetLabelOffset(0.01);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.14, 0.03, 0.07, 0.02);
    c1.SetTicks(1,1);
    c1.SetGrid(1,1);
    #c1.SetLogz(1);

    #frame1 = c1.DrawFrame(-0.1, -0.1, +0.1, +0.1);
    #frame1.GetXaxis().SetTitle("#eta_{l}^{gen} #minus #eta_{l}^{rec}");
    #frame1.GetYaxis().SetTitle("#varphi_{l}^{gen} #minus #varphi_{l}^{rec} (rad.)");
    #frame1.GetXaxis().SetTitleOffset(1.2);
    #frame1.GetYaxis().SetTitleOffset(1.8);
    #frame1.GetXaxis().SetTitleSize(0.035);
    #frame1.GetYaxis().SetTitleSize(0.035);
    #frame1.GetXaxis().SetLabelSize(0.035);
    #frame1.GetYaxis().SetLabelSize(0.035);
    #frame1.GetXaxis().SetLabelOffset(0.01);
    #frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(1);
    #ROOT.SetOwnership(frame1,False);
    h3.Draw("");

    txt = TPaveText(0.16, 0.72, 0.4, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE simulation");

    sign_str = "pos";

    if "pp" in suffix:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    elif "PbPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    elif "pPb" in suffix:
        txt.AddText("{0:d}#minus{1:d}% p#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    else:
        txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");

    if "LHC23k4g" in suffix:
        txt.AddText("LHC23k4g");

    txt.AddText("{0:3.2f} < p_{{T,l}}^{{gen}} < {1:3.2f} GeV/c".format(minptgen, maxptgen));
    if sign > 0:
        txt.AddText("positive leptons");
        sign_str = "pos";
    else:
        txt.AddText("negative leptons");
        sign_str = "neg";
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_reso_map_3d_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.eps".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_3d_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.pdf".format(date, minptgen, maxptgen, sign_str, suffix));
    c1.SaveAs("{0}_reso_map_3d_ptgen{1:3.2f}_{2:3.2f}GeV_{3}{4}.png".format(date, minptgen, maxptgen, sign_str, suffix));

#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_245050.root";
    #extract(filename, "Electron"      , "20240805_resolution_map_electron_LHC24b1b.root");
    #extract(filename, "StandaloneMuon", "20240805_resolution_map_sa_muon_LHC24b1b.root");
    #extract(filename, "GlobalMuon"    , "20240805_resolution_map_gl_muon_LHC24b1b.root");

    arr_pt = np.array([0.00, 0.05, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00], dtype=float);
    # arr_pt = np.array([0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00], dtype=float); # best binning, if there is a lot of statistics in MC.

    #filename = "AnalysisResults_reso_map_LHC24f3.root";
    #extract(filename, "Electron"      , "20240815_resolution_map_electron_LHC24f3.root", arr_pt);
    #extract(filename, "StandaloneMuon", "20240815_resolution_map_sa_muon_LHC24f3.root", arr_pt);
    #extract(filename, "GlobalMuon"    , "20240815_resolution_map_gl_muon_LHC24f3.root", arr_pt);

    #filename = "AnalysisResults_HL_280634.root";
    #extract(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20241025_resolution_map_electron_LHC23k4g.root", arr_pt);
    #extract(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20241025_resolution_map_sa_muon_LHC23k4g.root" , arr_pt);
    #extract(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20241025_resolution_map_gl_muon_LHC23k4g.root" , arr_pt);

    #filename = "AnalysisResults_HL_289802.root";
    #extract(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20241109_resolution_map_electron_LHC24g3.root", arr_pt);
    #extract(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20241109_resolution_map_sa_muon_LHC24g3.root" , arr_pt);
    #extract(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20241109_resolution_map_gl_muon_LHC24g3.root" , arr_pt);

    #filename = "AnalysisResults_HL_317206.root";
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20241225_resolution_map_electron_LHC24g3.root", arr_pt);
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20241225_resolution_map_sa_muon_LHC24g3.root" , arr_pt);
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20241225_resolution_map_gl_muon_LHC24g3.root" , arr_pt);

    #filename = "AnalysisResults_HL_318619.root"; # LHC24g3_small
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20241230_resolution_map_electron_LHC24g3.root", arr_pt);
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20241230_resolution_map_sa_muon_LHC24g3.root" , arr_pt);
    #extract_hs(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20241230_resolution_map_gl_muon_LHC24g3.root" , arr_pt);

    filename = "AnalysisResults_HL_317370.root"; # LHC23k4g
    suffix = "_pp_13.6TeV_LHC23k4g";
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, +1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6  , +1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1  , +1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 5.0, 10.0 , +1, suffix);

    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, -1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6  , -1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1  , -1, suffix);
    #check_ptgen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 5.0, 10.0 , -1, suffix);

    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, +1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6  , +1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1  , +1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 5.0, 10.0 , +1, suffix);

    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, -1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6  , -1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1  , -1, suffix);
    #check_ptgen_vs_ptphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 5.0, 10.0 , -1, suffix);

    #check_ptgen_vs_etaphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, +1, suffix);
    #check_ptgen_vs_etaphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6  , +1, suffix);
    #check_ptgen_vs_etaphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1  , +1, suffix);
    #check_ptgen_vs_etaphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 5.0, 10.0 , +1, suffix);

    #check_ptgen_vs_ptetaphi(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, +1, suffix);

    #check_ptgen_vs_etagen_deltapt(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, +1, suffix);
    #check_ptgen_vs_etagen_deltapt(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, +1, suffix);
    #check_ptgen_vs_etagen_deltapt(filename, "create-resolution-map_allgen_itsib1st", "Electron", 1.0, 1.1, +1, suffix);
    #check_ptgen_vs_etagen_deltaeta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, +1, suffix);
    #check_ptgen_vs_etagen_deltaeta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, +1, suffix);

    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, 0.0, 0.1, +1, suffix);
    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, 0.4, 0.5, +1, suffix);
    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.14, 0.16, 1.0, 1.1, +1, suffix);
    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, 0.0, 0.1, +1, suffix);
    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, 0.4, 0.5, +1, suffix);
    #check_ptgen_etagen_vs_pteta(filename, "create-resolution-map_allgen_itsib1st", "Electron", 0.5, 0.6, 1.0, 1.1, +1, suffix);

    filename = "AnalysisResults_HL_321449.root"; # LHC23k4g
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20250110_resolution_map_electron_LHC23k4g.root", arr_pt);
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20250110_resolution_map_sa_muon_LHC23k4g.root" , arr_pt);
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20250110_resolution_map_gl_muon_LHC23k4g.root" , arr_pt);

    filename = "AnalysisResults_HL_320757.root"; # LHC24f3
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "Electron"      , "20250110_resolution_map_electron_LHC24f3.root", arr_pt);
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "StandaloneMuon", "20250110_resolution_map_sa_muon_LHC24f3.root" , arr_pt);
    extract_hs(filename, "create-resolution-map_allgen_itsib1st", "GlobalMuon"    , "20250110_resolution_map_gl_muon_LHC24f3.root" , arr_pt);
#__________________________________________________________
