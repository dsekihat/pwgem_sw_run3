import os
import sys
import shutil
import math
import numpy as np
import yaml
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory, TH1D, TH2D, THnSparseD
from dilepton_utils import get_yield_1d
from histo_manager import rebin_histogram, get_bkg_subtracted, rebin_profile
import sys
sys.path.append("../common/");

def run_efficiency_mll_ptll_dcall(filename, tasks, arr_mll, arr_ptll, arr_dcall, system, energy, suffix):
    print(sys._getframe().f_code.co_name);
    print("filename = ", filename);
    print("arr_mll = ", arr_mll);
    print("arr_ptll = ", arr_ptll);
    print("arr_dcall = ", arr_dcall);
    __delta = 1e-3;

    rootfile = TFile.Open(filename, "READ");
    outname = "";
    axis_title_ll = "ll";

    if "dielectron" in tasks[0]["name"]:
        axis_title_ll = "ll";
        if "pp" in system:
            outname = "dielectron_mc_{0}_{1:3.1f}TeV{2}.root".format(system, energy, suffix)
        else:
            outname = "dielectron_mc_{0}_{1:3.2f}TeV{2}.root".format(system, energy, suffix)
    elif "dimuon" in tasks[0]["name"]:
        axis_title_ll = "#mu#mu";
        if "pp" in system:
            outname = "dimuon_{0}_{1:3.1f}TeV{2}.root".format(system, energy, suffix)
        else:
            outname = "dimuon_{0}_{1:3.2f}TeV{2}.root".format(system, energy, suffix)

    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    for itask in range(0, len(tasks)):
        taskname = tasks[itask]["name"];
        cent_min = tasks[itask]["cent_min"];
        cent_max = tasks[itask]["cent_min"];
        print("analyzing ... ", taskname, cent_min, cent_max); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        rootdir_gen = rootdir_task.Get("Generated");
        rootdir_event.ls();
        rootdir_pair .ls();
        rootdir_gen  .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();

        parnames_org = [
            "sm/Pi0/hs",
            "sm/Eta/hs",
            "sm/EtaPrime/hs",
            "sm/Rho/hs",
            "sm/Omega/hs",
            "sm/Phi/hs",
            "sm/PromptJPsi/hs",
            "sm/NonPromptJPsi/hs",
            "sm/PromptPsi2S/hs",
            "sm/NonPromptPsi2S/hs",
            "ccbar/c2l_c2l/hadron_hadron/hs",
            "bbbar/b2l_b2l/hadron_hadron/hs",
            "bbbar/b2c2l_b2c2l/hadron_hadron/hs",
            "bbbar/b2c2l_b2l_sameb/hadron_hadron/hs",
            "bbbar/b2c2l_b2l_diffb/hadron_hadron/hs",
        ];

        parnames = [
            "pi0",
            "eta",
            "etaprime",
            "rho",
            "omega",
            "phi",
            "promptjpsi",
            "nonpromptjpsi",
            "promptpsi2s",
            "nonpromptpsi2s",
            "c2l_c2l",
            "b2l_b2l",
            "b2c2l_b2c2l",
            "b2c2l_b2l_sameb",
            "b2c2l_b2l_diffb",
        ];

        nbin = np.array([len(arr_mll)-1, len(arr_ptll)-1, len(arr_dcall)-1], dtype=np.int32);
        xmin = np.array([arr_mll[0], arr_ptll[0], arr_dcall[0]], dtype=np.float64);
        xmax = np.array([arr_mll[-1], arr_ptll[-1], arr_dcall[-1]], dtype=np.float64);
        hs_tmp = THnSparseD("hs_tmp", "#it{N}_{ll}^{all}/#it{N}_{ev};#it{m}_{ll} (GeV/#it{c}^{2});#it{p}_{T,ll} (GeV/#it{c});DCA_{ll} (#sigma);", 3, nbin, xmin, xmax);
        hs_tmp.Sumw2();
        hs_tmp.SetBinEdges(0, arr_mll);
        hs_tmp.SetBinEdges(1, arr_ptll);
        hs_tmp.SetBinEdges(2, arr_dcall);

        #rec. part
        for ipar, parname_org in enumerate(parnames_org):
            print(ipar, parname_org);
            hs_rec = hs_tmp.Clone("hs_rec_{0}".format(parnames[ipar]));
            hs_rec.SetTitle("#it{N}_{ll}^{rec}/#it{N}_{ev}");
            hs_rec_org = rootdir_pair.Get(parname_org);
            hs_rec_org.Sumw2();

            for idca in range(0, len(arr_dcall)-1):
                dca_min = arr_dcall[idca];
                dca_max = arr_dcall[idca+1];
                dca_center = (dca_min + dca_max)/2;
                bin_dca1 = hs_rec_org.GetAxis(10).FindBin(dca_min + __delta);
                bin_dca2 = hs_rec_org.GetAxis(10).FindBin(dca_max - __delta);
                hs_rec_org.GetAxis(10).SetRange(bin_dca1, bin_dca2);
    
                for ipt in range(0, len(arr_ptll)-1):
                    pt_min = arr_ptll[ipt];
                    pt_max = arr_ptll[ipt+1];
                    pt_center = (pt_min + pt_max)/2;
                    bin_pt1 = hs_rec_org.GetAxis(1).FindBin(pt_min + __delta);
                    bin_pt2 = hs_rec_org.GetAxis(1).FindBin(pt_max - __delta);
                    hs_rec_org.GetAxis(1).SetRange(bin_pt1, bin_pt2);
    
                    h1m_rec_org = hs_rec_org.Projection(0);
                    h1m_rec_org.SetName("h1m_rec_{0}_pt{1}_dca{2}_org".format(parnames[ipar], ipt, idca));
                    h1m_rec_org.SetTitle("#it{{N}}^{{gen}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                    h1m_rec_org.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                    h1m_rec_org.SetYTitle("counts per bin");
    
                    h1m_rec_rebin  = rebin_histogram(h1m_rec_org , arr_mll, False, False);
                    h1m_rec_rebin .SetName("h1m_rec_{0}_pt{1}_dca{2}".format(parnames[ipar], ipt, idca));
                    h1m_rec_rebin .SetTitle("#it{{N}}^{{gen}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
    
                    for im in range(0, len(arr_mll)-1):
                        m_min = arr_mll[im];
                        m_max = arr_mll[im+1];
                        m_center = (m_min + m_max)/2;
                        global_bin_id = hs_rec.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                        hs_rec.SetBinContent(global_bin_id, h1m_rec_rebin.GetBinContent(im+1));
                        hs_rec.SetBinError(global_bin_id, h1m_rec_rebin.GetBinError(im+1));
            hs_rec.Scale(1/nev);
            outdir.WriteTObject(hs_rec);

        hs_rec_pi0 = outdir.Get("hs_rec_pi0");
        hs_rec_eta = outdir.Get("hs_rec_eta");
        hs_rec_etaprime = outdir.Get("hs_rec_etaprime");
        hs_rec_rho = outdir.Get("hs_rec_rho");
        hs_rec_omega = outdir.Get("hs_rec_omega");
        hs_rec_phi = outdir.Get("hs_rec_phi");
        hs_rec_promptjpsi = outdir.Get("hs_rec_promptjpsi");
        hs_rec_nonpromptjpsi = outdir.Get("hs_rec_nonpromptjpsi");
        hs_rec_promptpsi2s = outdir.Get("hs_rec_promptpsi2s");
        hs_rec_nonpromptpsi2s = outdir.Get("hs_rec_nonpromptpsi2s");

        hs_rec_sm = hs_rec_pi0.Clone("hs_rec_sm");
        hs_rec_sm.Add(hs_rec_eta, 1.);
        hs_rec_sm.Add(hs_rec_etaprime, 1.);
        hs_rec_sm.Add(hs_rec_rho, 1.);
        hs_rec_sm.Add(hs_rec_omega, 1.);
        hs_rec_sm.Add(hs_rec_phi, 1.);
        hs_rec_sm.Add(hs_rec_promptjpsi, 1.);
        hs_rec_sm.Add(hs_rec_nonpromptjpsi, 1.);
        hs_rec_sm.Add(hs_rec_promptpsi2s, 1.);
        hs_rec_sm.Add(hs_rec_nonpromptpsi2s, 1.);
        outdir.WriteTObject(hs_rec_sm);
        h1m_rec_sm = hs_rec_sm.Projection(0);
        h1m_rec_sm.SetName("h1m_rec_sm");
        outdir.WriteTObject(h1m_rec_sm);

        hs_rec_c2l_c2l = outdir.Get("hs_rec_c2l_c2l");
        hs_rec_ccbar = hs_rec_c2l_c2l.Clone("hs_rec_ccbar");
        outdir.WriteTObject(hs_rec_ccbar);
        h1m_rec_ccbar = hs_rec_ccbar.Projection(0);
        h1m_rec_ccbar.SetName("h1m_rec_ccbar");
        outdir.WriteTObject(h1m_rec_ccbar);

        hs_rec_b2l_b2l = outdir.Get("hs_rec_b2l_b2l");
        hs_rec_b2c2l_b2c2l = outdir.Get("hs_rec_b2c2l_b2c2l");
        hs_rec_b2c2l_b2l_sameb = outdir.Get("hs_rec_b2c2l_b2l_sameb");
        hs_rec_b2c2l_b2l_diffb = outdir.Get("hs_rec_b2c2l_b2l_diffb"); #LS
        hs_rec_bbbar = hs_rec_b2l_b2l.Clone("hs_rec_bbbar");
        hs_rec_bbbar.Add(hs_rec_b2c2l_b2c2l, 1.);
        hs_rec_bbbar.Add(hs_rec_b2c2l_b2l_sameb, 1.);
        hs_rec_bbbar.Add(hs_rec_b2c2l_b2l_diffb, -1.);
        outdir.WriteTObject(hs_rec_bbbar);
        h1m_rec_bbbar = hs_rec_bbbar.Projection(0);
        h1m_rec_bbbar.SetName("h1m_rec_bbbar");
        outdir.WriteTObject(h1m_rec_bbbar);

        nbin = np.array([len(arr_mll)-1, len(arr_ptll)-1], dtype=np.int32);
        xmin = np.array([arr_mll[0], arr_ptll[0]], dtype=np.float64);
        xmax = np.array([arr_mll[-1], arr_ptll[-1]], dtype=np.float64);
        hs_tmp_gen = THnSparseD("hs_tmp_gen", "#it{N}_{ll}^{all}/#it{N}_{ev};#it{m}_{ll} (GeV/#it{c}^{2});#it{p}_{T,ll} (GeV/#it{c});", 2, nbin, xmin, xmax);
        hs_tmp_gen.Sumw2();
        hs_tmp_gen.SetBinEdges(0, arr_mll);
        hs_tmp_gen.SetBinEdges(1, arr_ptll);
        #gen. part
        for ipar, parname_org in enumerate(parnames_org):
            print(ipar, parname_org);
            hs_gen = hs_tmp_gen.Clone("hs_gen_{0}".format(parnames[ipar]));
            hs_gen.SetTitle("#it{N}_{ll}^{gen}/#it{N}_{ev}");
            hs_gen_org = rootdir_gen.Get(parname_org);
            hs_gen_org.Sumw2();

            for ipt in range(0, len(arr_ptll)-1):
                pt_min = arr_ptll[ipt];
                pt_max = arr_ptll[ipt+1];
                pt_center = (pt_min + pt_max)/2;
                bin_pt1 = hs_gen_org.GetAxis(1).FindBin(pt_min + __delta);
                bin_pt2 = hs_gen_org.GetAxis(1).FindBin(pt_max - __delta);
                hs_gen_org.GetAxis(1).SetRange(bin_pt1, bin_pt2);
    
                h1m_gen_org = hs_gen_org.Projection(0);
                h1m_gen_org.SetName("h1m_gen_{0}_pt{1}_org".format(parnames[ipar], ipt));
                h1m_gen_org.SetTitle("#it{{N}}^{{gen}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c".format(pt_min, pt_max));
                h1m_gen_org.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_gen_org.SetYTitle("counts per bin");
    
                h1m_gen_rebin  = rebin_histogram(h1m_gen_org , arr_mll, False, False);
                h1m_gen_rebin .SetName("h1m_gen_{0}_pt{1}".format(parnames[ipar], ipt));
                h1m_gen_rebin .SetTitle("#it{{N}}^{{gen}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c".format(pt_min, pt_max));
    
                for im in range(0, len(arr_mll)-1):
                    m_min = arr_mll[im];
                    m_max = arr_mll[im+1];
                    m_center = (m_min + m_max)/2;
                    global_bin_id = hs_gen.GetBin(np.array([m_center, pt_center], dtype=np.float64), True);
                    hs_gen.SetBinContent(global_bin_id, h1m_gen_rebin.GetBinContent(im+1));
                    hs_gen.SetBinError(global_bin_id, h1m_gen_rebin.GetBinError(im+1));
            hs_gen.Scale(1/nev);
            outdir.WriteTObject(hs_gen);

        hs_gen_pi0 = outdir.Get("hs_gen_pi0");
        hs_gen_eta = outdir.Get("hs_gen_eta");
        hs_gen_etaprime = outdir.Get("hs_gen_etaprime");
        hs_gen_rho = outdir.Get("hs_gen_rho");
        hs_gen_omega = outdir.Get("hs_gen_omega");
        hs_gen_phi = outdir.Get("hs_gen_phi");
        hs_gen_promptjpsi = outdir.Get("hs_gen_promptjpsi");
        hs_gen_nonpromptjpsi = outdir.Get("hs_gen_nonpromptjpsi");
        hs_gen_promptpsi2s = outdir.Get("hs_gen_promptpsi2s");
        hs_gen_nonpromptpsi2s = outdir.Get("hs_gen_nonpromptpsi2s");

        hs_gen_sm = hs_gen_pi0.Clone("hs_gen_sm");
        hs_gen_sm.Add(hs_gen_eta, 1.);
        hs_gen_sm.Add(hs_gen_etaprime, 1.);
        hs_gen_sm.Add(hs_gen_rho, 1.);
        hs_gen_sm.Add(hs_gen_omega, 1.);
        hs_gen_sm.Add(hs_gen_phi, 1.);
        hs_gen_sm.Add(hs_gen_promptjpsi, 1.);
        hs_gen_sm.Add(hs_gen_nonpromptjpsi, 1.);
        hs_gen_sm.Add(hs_gen_promptpsi2s, 1.);
        hs_gen_sm.Add(hs_gen_nonpromptpsi2s, 1.);
        outdir.WriteTObject(hs_gen_sm);
        h1m_gen_sm = hs_gen_sm.Projection(0);
        h1m_gen_sm.SetName("h1m_gen_sm");
        outdir.WriteTObject(h1m_gen_sm);

        hs_gen_c2l_c2l = outdir.Get("hs_gen_c2l_c2l");
        hs_gen_ccbar = hs_gen_c2l_c2l.Clone("hs_gen_ccbar");
        outdir.WriteTObject(hs_gen_ccbar);
        h1m_gen_ccbar = hs_gen_ccbar.Projection(0);
        h1m_gen_ccbar.SetName("h1m_gen_ccbar");
        outdir.WriteTObject(h1m_gen_ccbar);

        hs_gen_b2l_b2l = outdir.Get("hs_gen_b2l_b2l");
        hs_gen_b2c2l_b2c2l = outdir.Get("hs_gen_b2c2l_b2c2l");
        hs_gen_b2c2l_b2l_sameb = outdir.Get("hs_gen_b2c2l_b2l_sameb");
        hs_gen_b2c2l_b2l_diffb = outdir.Get("hs_gen_b2c2l_b2l_diffb"); #LS
        hs_gen_bbbar = hs_gen_b2l_b2l.Clone("hs_gen_bbbar");
        hs_gen_bbbar.Add(hs_gen_b2c2l_b2c2l, 1.);
        hs_gen_bbbar.Add(hs_gen_b2c2l_b2l_sameb, 1.);
        hs_gen_bbbar.Add(hs_gen_b2c2l_b2l_diffb, -1.);
        outdir.WriteTObject(hs_gen_bbbar);
        h1m_gen_bbbar = hs_gen_bbbar.Projection(0);
        h1m_gen_bbbar.SetName("h1m_gen_bbbar");
        outdir.WriteTObject(h1m_gen_bbbar);

        h1m_eff_sm = h1m_rec_sm.Clone("h1m_eff_sm");
        h1m_eff_sm.Reset();
        h1m_eff_sm.Divide(h1m_rec_sm, h1m_gen_sm, 1., 1., "B");
        h1m_eff_sm.SetYTitle("pair rec. efficiency");
        outdir.WriteTObject(h1m_eff_sm);

        h1m_eff_ccbar = h1m_rec_ccbar.Clone("h1m_eff_ccbar");
        h1m_eff_ccbar.Reset();
        h1m_eff_ccbar.Divide(h1m_rec_ccbar, h1m_gen_ccbar, 1., 1., "B");
        h1m_eff_ccbar.SetYTitle("pair rec. efficiency");
        outdir.WriteTObject(h1m_eff_ccbar);

        h1m_eff_bbbar = h1m_rec_bbbar.Clone("h1m_eff_bbbar");
        h1m_eff_bbbar.Reset();
        h1m_eff_bbbar.Divide(h1m_rec_bbbar, h1m_gen_bbbar, 1., 1., "B");
        h1m_eff_bbbar.SetYTitle("pair rec. efficiency");
        outdir.WriteTObject(h1m_eff_bbbar);

    outfile.Close();
    rootfile.Close();

