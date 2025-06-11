import os
import sys
import shutil
import math
import numpy as np
import yaml
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory, TH1F, TH2F, THnSparseD
from dilepton_utils import get_yield_1d, get_flow
import sys
sys.path.append("../common/");
from histo_manager import rebin_profile
from signal_extractor import get_Rn
from dilepton_efficiency import run_efficiency_mll_ptll_dcall

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

config = "";
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
print(config)
system = config["common"]["system"];
energy = config["common"]["energy"]; #float
suffix = args.suffix;

period = "";
is_mc = False;
if "mc" in args.type:
    is_mc = True;
    period = config["period_mc"];
else:
    is_mc = False;
    period = config["period_data"];
pass_number = config["pass_number"];

#__________________________________________________________________
def set_range_00(hs, list_axis_ids):
    for ia in list_axis_ids:
        hs.GetAxis(ia).SetRange(0, 0);

#__________________________________________________________________
def run_mll_ptll_dcall(filename, tasks, arr_mll, arr_ptll, arr_dcall):
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
            outname = "dielectron_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
        else:
            outname = "dielectron_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    elif "dimuon" in tasks[0]["name"]:
        axis_title_ll = "#mu#mu";
        if "pp" in system:
            outname = "dimuon_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
        else:
            outname = "dimuon_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)

    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    for itask in range(0, len(tasks)):
        taskname = tasks[itask]["name"];
        cent_min = tasks[itask]["cent_min"];
        cent_max = tasks[itask]["cent_max"];
        print("analyzing ... ", taskname, cent_min, cent_max); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        #rootdir_event.ls();
        #rootdir_pair .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();
        if nev < 0.5:
            return None;

        h1norm = rootdir_task.Get("Event/norm/hCollisionCounter");
        if h1norm != None:
            nbin_tmp = h1norm.GetNbinsX();
            nev = h1norm.GetBinContent(nbin_tmp);
            print("nev from norm = {0:e}".format(nev));
            h1norm.SetName("h1norm");
            outdir.WriteTObject(h1norm);

        hs_uls_same  = rootdir_pair.Get("same/uls/hs") .Clone("hs_uls_same")
        hs_lspp_same = rootdir_pair.Get("same/lspp/hs").Clone("hs_lspp_same")
        hs_lsmm_same = rootdir_pair.Get("same/lsmm/hs").Clone("hs_lsmm_same")
        hs_uls_mix   = rootdir_pair.Get("mix/uls/hs")  .Clone("hs_uls_mix")
        hs_lspp_mix  = rootdir_pair.Get("mix/lspp/hs") .Clone("hs_lspp_mix")
        hs_lsmm_mix  = rootdir_pair.Get("mix/lsmm/hs") .Clone("hs_lsmm_mix")
        hs_uls_same .Sumw2();
        hs_lspp_same.Sumw2();
        hs_lsmm_same.Sumw2();
        hs_uls_mix  .Sumw2();
        hs_lspp_mix .Sumw2();
        hs_lsmm_mix .Sumw2();

        #outdir.WriteTObject(hs_uls_same);
        #outdir.WriteTObject(hs_lspp_same);
        #outdir.WriteTObject(hs_lsmm_same);
        #outdir.WriteTObject(hs_uls_mix);
        #outdir.WriteTObject(hs_lspp_mix);
        #outdir.WriteTObject(hs_lsmm_mix);

        set_range_00(hs_uls_same , [0, 1, 2]);
        set_range_00(hs_lspp_same, [0, 1, 2]);
        set_range_00(hs_lsmm_same, [0, 1, 2]);
        set_range_00(hs_uls_mix  , [0, 1, 2]);
        set_range_00(hs_lspp_mix , [0, 1, 2]);
        set_range_00(hs_lsmm_mix , [0, 1, 2]);

        __ndim = 3;
        nbin = np.array([len(arr_mll)-1, len(arr_ptll)-1, len(arr_dcall)-1], dtype=np.int32);
        xmin = np.array([arr_mll[0], arr_ptll[0], arr_dcall[0]], dtype=np.float64);
        xmax = np.array([arr_mll[-1], arr_ptll[-1], arr_dcall[-1]], dtype=np.float64);
        hs_uls = THnSparseD("hs_uls", "#it{N}_{ll}^{all}/#it{N}_{ev};#it{m}_{ll} (GeV/#it{c}^{2});#it{p}_{T,ll} (GeV/#it{c});DCA_{ll}^{3D} (#sigma);", __ndim, nbin, xmin, xmax);
        hs_uls.Sumw2();
        hs_uls.SetBinEdges(0, arr_mll);
        hs_uls.SetBinEdges(1, arr_ptll);
        hs_uls.SetBinEdges(2, arr_dcall);
        hs_R = hs_uls.Clone("hs_R");
        hs_R.SetTitle("R factor");
        hs_bkg = hs_uls.Clone("hs_bkg");
        hs_bkg.SetTitle("#it{N}_{ll}^{bkg}/#it{N}_{ev}");
        hs_sig = hs_uls.Clone("hs_sig");
        hs_sig.SetTitle("#it{N}_{ll}^{signal}/#it{N}_{ev}");

        for idca in range(0, len(arr_dcall)-1):
            dca_min = arr_dcall[idca];
            dca_max = arr_dcall[idca+1];
            dca_center = (dca_min + dca_max)/2;
            bin_dca1 = hs_uls_same.GetAxis(2).FindBin(dca_min + __delta);
            bin_dca2 = hs_uls_same.GetAxis(2).FindBin(dca_max - __delta);
            hs_uls_same .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_uls_mix  .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);

            h2R = TH2F("h2R_dca{0}".format(idca), "R factor in {0:2.1f} < DCA_{{ll}}^{{3D}} < {1:2.1f} #sigma;#it{{m}}_{{ll}} (GeV/#it{{c}}^{{2}});#it{{p}}_{{T,ll}} (GeV/#it{{c}})".format(dca_min, dca_max) , len(arr_mll)-1, arr_mll, len(arr_ptll)-1, arr_ptll);
            h2R.Sumw2();
            h2R.SetContour(100);

            for ipt in range(0, len(arr_ptll)-1):
                pt_min = arr_ptll[ipt];
                pt_max = arr_ptll[ipt+1];
                pt_center = (pt_min + pt_max)/2;
                bin_pt1 = hs_uls_same.GetAxis(1).FindBin(pt_min + __delta);
                bin_pt2 = hs_uls_same.GetAxis(1).FindBin(pt_max - __delta);
                hs_uls_same .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_uls_mix  .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);

                h1m_uls_same  = hs_uls_same.Projection(0);
                h1m_lspp_same = hs_lspp_same.Projection(0);
                h1m_lsmm_same = hs_lsmm_same.Projection(0);
                h1m_uls_mix   = hs_uls_mix.Projection(0);
                h1m_lspp_mix  = hs_lspp_mix.Projection(0);
                h1m_lsmm_mix  = hs_lsmm_mix.Projection(0);

                h1m_uls_same .SetName("h1m_uls_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_same.SetName("h1m_lspp_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_same.SetName("h1m_lsmm_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_uls_mix  .SetName("h1m_uls_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_mix .SetName("h1m_lspp_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_mix .SetName("h1m_lsmm_mix_pt{0}_dca{1}".format(ipt, idca));

                h1m_uls_same .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1m_uls_same .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lspp_same.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lsmm_same.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_uls_mix  .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lspp_mix .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lsmm_mix .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");

                h1m_uls_same .SetYTitle("counts per bin");
                h1m_lspp_same.SetYTitle("counts per bin");
                h1m_lsmm_same.SetYTitle("counts per bin");
                h1m_uls_mix  .SetYTitle("counts per bin");
                h1m_lspp_mix .SetYTitle("counts per bin");
                h1m_lsmm_mix .SetYTitle("counts per bin");

                [h1m_uls_same_rebin, h1m_lspp_same_rebin, h1m_lsmm_same_rebin, h1m_uls_mix_rebin, h1m_lspp_mix_rebin, h1m_lsmm_mix_rebin, h1R, h1bkg, h1sig] = get_yield_1d(h1m_uls_same, h1m_lspp_same, h1m_lsmm_same, h1m_uls_mix, h1m_lspp_mix, h1m_lsmm_mix, arr_mll);

                h1m_uls_same_rebin .SetName("h1m_uls_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_same_rebin.SetName("h1m_lspp_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_same_rebin.SetName("h1m_lsmm_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_uls_mix_rebin  .SetName("h1m_uls_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_mix_rebin .SetName("h1m_lspp_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_mix_rebin .SetName("h1m_lsmm_mix_pt{0}_dca{1}_rebin".format(ipt, idca));

                h1m_uls_same_rebin .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same_rebin.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same_rebin.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix_rebin  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix_rebin .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix_rebin .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetName("h1R_pt{0}_dca{1}"  .format(ipt, idca));
                h1bkg.SetName("h1bkg_pt{0}_dca{1}".format(ipt, idca));
                h1sig.SetName("h1sig_pt{0}_dca{1}".format(ipt, idca));
                h1R  .SetTitle("R factor in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1bkg.SetTitle("#it{{N}}_{{ll}}^{{bkg}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1sig.SetTitle("#it{{N}}_{{ll}}^{{signal}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetYTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
                h1bkg.SetYTitle("counts per bin");
                h1sig.SetYTitle("counts per bin");

                #outdir.WriteTObject(h1m_uls_same_rebin );
                #outdir.WriteTObject(h1m_lspp_same_rebin);
                #outdir.WriteTObject(h1m_lsmm_same_rebin);
                #outdir.WriteTObject(h1m_uls_mix_rebin  );
                #outdir.WriteTObject(h1m_lspp_mix_rebin );
                #outdir.WriteTObject(h1m_lsmm_mix_rebin );
                #outdir.WriteTObject(h1R);
                #outdir.WriteTObject(h1bkg);
                #outdir.WriteTObject(h1sig);

                for im in range(0, len(arr_mll)-1):
                    m_min = arr_mll[im];
                    m_max = arr_mll[im+1];
                    m_center = (m_min + m_max)/2;
                    global_bin_id = hs_uls.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_R  .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_bkg.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_sig.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    #print(m_center, pt_center, dca_center, global_bin_id);
                    hs_uls.SetBinContent(global_bin_id, h1m_uls_same_rebin.GetBinContent(im+1));
                    hs_uls.SetBinError(global_bin_id, h1m_uls_same_rebin.GetBinError(im+1));
                    hs_R  .SetBinContent(global_bin_id, h1R.GetBinContent(im+1));
                    hs_R  .SetBinError(global_bin_id, h1R.GetBinError(im+1));
                    hs_bkg.SetBinContent(global_bin_id, h1bkg.GetBinContent(im+1));
                    hs_bkg.SetBinError(global_bin_id, h1bkg.GetBinError(im+1));
                    hs_sig.SetBinContent(global_bin_id, h1sig.GetBinContent(im+1));
                    hs_sig.SetBinError(global_bin_id, h1sig.GetBinError(im+1));
                    h2R.SetBinContent(im+1, ipt+1, h1R.GetBinContent(im+1));
                    h2R.SetBinError(im+1, ipt+1, h1R.GetBinError(im+1));
            outdir.WriteTObject(h2R);
        hs_uls.Scale(1/nev);
        hs_bkg.Scale(1/nev);
        hs_sig.Scale(1/nev);
        outdir.WriteTObject(hs_uls);
        outdir.WriteTObject(hs_R  );
        outdir.WriteTObject(hs_bkg);
        outdir.WriteTObject(hs_sig);

    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
def run_flow(filename, tasks, arr_mll, arr_ptll, arr_dcall):
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
            outname = "dielectron_flow_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
        else:
            outname = "dielectron_flow_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    elif "dimuon" in tasks[0]["name"]:
        axis_title_ll = "#mu#mu";
        if "pp" in system:
            outname = "dimuon_flow_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
        else:
            outname = "dimuon_flow_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)

    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    for itask in range(0, len(tasks)):
        taskname = tasks[itask]["name"];
        cent_min = tasks[itask]["cent_min"];
        cent_max = tasks[itask]["cent_max"];
        print("analyzing ... ", taskname, cent_min, cent_max); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        #rootdir_pair_same = rootdir_task.Get("Pair/same/uls").ls();
        #rootdir_pair_mix = rootdir_task.Get("Pair/mix/uls").ls();
        #rootdir_event.ls();
        #rootdir_pair .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();
        if nev < 0.5:
            return None;

        hPrfQ2FT0MQ2BPos_CentFT0C = rootdir_event.Get("hPrfQ2FT0MQ2BPos_CentFT0C");
        hPrfQ2FT0MQ2BNeg_CentFT0C = rootdir_event.Get("hPrfQ2FT0MQ2BNeg_CentFT0C");
        hPrfQ2BPosQ2BNeg_CentFT0C = rootdir_event.Get("hPrfQ2BPosQ2BNeg_CentFT0C");
        outdir.WriteTObject(hPrfQ2FT0MQ2BPos_CentFT0C);
        outdir.WriteTObject(hPrfQ2FT0MQ2BNeg_CentFT0C);
        outdir.WriteTObject(hPrfQ2BPosQ2BNeg_CentFT0C);
        h1R2 = get_Rn(hPrfQ2FT0MQ2BPos_CentFT0C, hPrfQ2FT0MQ2BNeg_CentFT0C, hPrfQ2BPosQ2BNeg_CentFT0C);
        h1R2.SetName("h1R2");
        outdir.WriteTObject(h1R2);

        hPrfQ2FT0MQ2BPos_CentFT0C_clone = hPrfQ2FT0MQ2BPos_CentFT0C.Clone("hPrfQ2FT0MQ2BPos_CentFT0C_clone");
        hPrfQ2FT0MQ2BNeg_CentFT0C_clone = hPrfQ2FT0MQ2BNeg_CentFT0C.Clone("hPrfQ2FT0MQ2BNeg_CentFT0C_clone");
        hPrfQ2BPosQ2BNeg_CentFT0C_clone = hPrfQ2BPosQ2BNeg_CentFT0C.Clone("hPrfQ2BPosQ2BNeg_CentFT0C_clone");
        hPrfQ2FT0MQ2BPos_CentFT0C_rebin = hPrfQ2FT0MQ2BPos_CentFT0C_clone.Rebin(1, "hPrfQ2FT0MQ2BPos_CentFT0C_rebin", np.array([cent_min, cent_max], dtype=float));
        hPrfQ2FT0MQ2BNeg_CentFT0C_rebin = hPrfQ2FT0MQ2BNeg_CentFT0C_clone.Rebin(1, "hPrfQ2FT0MQ2BNeg_CentFT0C_rebin", np.array([cent_min, cent_max], dtype=float));
        hPrfQ2BPosQ2BNeg_CentFT0C_rebin = hPrfQ2BPosQ2BNeg_CentFT0C_clone.Rebin(1, "hPrfQ2BPosQ2BNeg_CentFT0C_rebin", np.array([cent_min, cent_max], dtype=float));
        outdir.WriteTObject(hPrfQ2FT0MQ2BPos_CentFT0C_rebin);
        outdir.WriteTObject(hPrfQ2FT0MQ2BNeg_CentFT0C_rebin);
        outdir.WriteTObject(hPrfQ2BPosQ2BNeg_CentFT0C_rebin);
        h1R2_rebin = get_Rn(hPrfQ2FT0MQ2BPos_CentFT0C_rebin, hPrfQ2FT0MQ2BNeg_CentFT0C_rebin, hPrfQ2BPosQ2BNeg_CentFT0C_rebin);
        h1R2_rebin.SetName("h1R2_rebin");
        outdir.WriteTObject(h1R2_rebin);
        R2 = h1R2_rebin.GetBinContent(1);
        R2_err = h1R2_rebin.GetBinError(1);
        print("R2 = ", R2);

        hs_uls_same  = rootdir_pair.Get("same/uls/hs") .Clone("hs_uls_same"); #for yield
        hs_lspp_same = rootdir_pair.Get("same/lspp/hs").Clone("hs_lspp_same"); #for yield
        hs_lsmm_same = rootdir_pair.Get("same/lsmm/hs").Clone("hs_lsmm_same"); #for yield
        hs_uls_mix   = rootdir_pair.Get("mix/uls/hs")  .Clone("hs_uls_mix"); #for yield
        hs_lspp_mix  = rootdir_pair.Get("mix/lspp/hs") .Clone("hs_lspp_mix"); #for yield
        hs_lsmm_mix  = rootdir_pair.Get("mix/lsmm/hs") .Clone("hs_lsmm_mix"); #for yield
        hs_uls_same .Sumw2();
        hs_lspp_same.Sumw2();
        hs_lsmm_same.Sumw2();
        hs_uls_mix  .Sumw2();
        hs_lspp_mix .Sumw2();
        hs_lsmm_mix .Sumw2();

        #outdir.WriteTObject(hs_uls_same);
        #outdir.WriteTObject(hs_lspp_same);
        #outdir.WriteTObject(hs_lsmm_same);
        #outdir.WriteTObject(hs_uls_mix);
        #outdir.WriteTObject(hs_lspp_mix);
        #outdir.WriteTObject(hs_lsmm_mix);

        set_range_00(hs_uls_same , [0, 1, 2, 3]);
        set_range_00(hs_lspp_same, [0, 1, 2, 3]);
        set_range_00(hs_lsmm_same, [0, 1, 2, 3]);
        set_range_00(hs_uls_mix  , [0, 1, 2]);
        set_range_00(hs_lspp_mix , [0, 1, 2]);
        set_range_00(hs_lsmm_mix , [0, 1, 2]);

        __ndim = 3;
        nbin = np.array([len(arr_mll)-1, len(arr_ptll)-1, len(arr_dcall)-1], dtype=np.int32);
        xmin = np.array([arr_mll[0], arr_ptll[0], arr_dcall[0]], dtype=np.float64);
        xmax = np.array([arr_mll[-1], arr_ptll[-1], arr_dcall[-1]], dtype=np.float64);
        hs_uls = THnSparseD("hs_uls", "#it{N}_{ll}^{all}/#it{N}_{ev};#it{m}_{ll} (GeV/#it{c}^{2});#it{p}_{T,ll} (GeV/#it{c});DCA_{ll}^{3D} (#sigma);", __ndim, nbin, xmin, xmax);
        hs_uls.Sumw2();
        hs_uls.SetBinEdges(0, arr_mll);
        hs_uls.SetBinEdges(1, arr_ptll);
        hs_uls.SetBinEdges(2, arr_dcall);
        hs_R = hs_uls.Clone("hs_R");
        hs_R.SetTitle("R factor");
        hs_bkg = hs_uls.Clone("hs_bkg");
        hs_bkg.SetTitle("#it{N}_{ll}^{bkg}/#it{N}_{ev}");
        hs_sig = hs_uls.Clone("hs_sig");
        hs_sig.SetTitle("#it{N}_{ll}^{signal}/#it{N}_{ev}");
        hs_vn_sig = hs_uls.Clone("hs_vn_sig");
        hs_vn_sig.SetTitle("v_{n}^{S}");
        hs_vn_bkg = hs_uls.Clone("hs_vn_bkg");
        hs_vn_bkg.SetTitle("v_{n}^{B}");
        hs_vn_tot = hs_uls.Clone("hs_vn_tot");
        hs_vn_tot.SetTitle("v_{n}^{S+B}");

        for idca in range(0, len(arr_dcall)-1):
            dca_min = arr_dcall[idca];
            dca_max = arr_dcall[idca+1];
            dca_center = (dca_min + dca_max)/2;

            bin_dca1 = hs_uls_same.GetAxis(2).FindBin(dca_min + __delta);
            bin_dca2 = hs_uls_same.GetAxis(2).FindBin(dca_max - __delta);
            hs_uls_same .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_uls_mix  .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);

            h2R = TH2F("h2R_dca{0}".format(idca), "R factor in {0:2.1f} < DCA_{{ll}}^{{3D}} < {1:2.1f} #sigma;#it{{m}}_{{ll}} (GeV/#it{{c}}^{{2}});#it{{p}}_{{T,ll}} (GeV/#it{{c}})".format(dca_min, dca_max) , len(arr_mll)-1, arr_mll, len(arr_ptll)-1, arr_ptll);
            h2R.Sumw2();
            h2R.SetContour(100);

            for ipt in range(0, len(arr_ptll)-1):
                pt_min = arr_ptll[ipt];
                pt_max = arr_ptll[ipt+1];
                pt_center = (pt_min + pt_max)/2;

                bin_pt1 = hs_uls_same.GetAxis(1).FindBin(pt_min + __delta);
                bin_pt2 = hs_uls_same.GetAxis(1).FindBin(pt_max - __delta);
                hs_uls_same .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_uls_mix  .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);

                h1m_uls_same  = hs_uls_same.Projection(0);
                h1m_lspp_same = hs_lspp_same.Projection(0);
                h1m_lsmm_same = hs_lsmm_same.Projection(0);
                h1m_uls_mix   = hs_uls_mix.Projection(0);
                h1m_lspp_mix  = hs_lspp_mix.Projection(0);
                h1m_lsmm_mix  = hs_lsmm_mix.Projection(0);

                h1m_uls_same .SetName("h1m_uls_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_same.SetName("h1m_lspp_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_same.SetName("h1m_lsmm_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_uls_mix  .SetName("h1m_uls_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_mix .SetName("h1m_lspp_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_mix .SetName("h1m_lsmm_mix_pt{0}_dca{1}".format(ipt, idca));

                h1m_uls_same .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1m_uls_same .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lspp_same.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lsmm_same.SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_uls_mix  .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lspp_mix .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");
                h1m_lsmm_mix .SetXTitle("#it{m}_{ll} (GeV/#it{c}^{2})");

                h1m_uls_same .SetYTitle("counts per bin");
                h1m_lspp_same.SetYTitle("counts per bin");
                h1m_lsmm_same.SetYTitle("counts per bin");
                h1m_uls_mix  .SetYTitle("counts per bin");
                h1m_lspp_mix .SetYTitle("counts per bin");
                h1m_lsmm_mix .SetYTitle("counts per bin");

                #print(h1m_uls_same );
                #print(h1m_lspp_same);
                #print(h1m_lsmm_same);
                #print(h1m_uls_mix  );
                #print(h1m_lspp_mix );
                #print(h1m_lsmm_mix );

                [h1m_uls_same_rebin, h1m_lspp_same_rebin, h1m_lsmm_same_rebin, h1m_uls_mix_rebin, h1m_lspp_mix_rebin, h1m_lsmm_mix_rebin, h1R, h1bkg, h1sig] = get_yield_1d(h1m_uls_same, h1m_lspp_same, h1m_lsmm_same, h1m_uls_mix, h1m_lspp_mix, h1m_lsmm_mix, arr_mll);
                h1m_uls_same_rebin .SetName("h1m_uls_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_same_rebin.SetName("h1m_lspp_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_same_rebin.SetName("h1m_lsmm_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_uls_mix_rebin  .SetName("h1m_uls_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_mix_rebin .SetName("h1m_lspp_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_mix_rebin .SetName("h1m_lsmm_mix_pt{0}_dca{1}_rebin".format(ipt, idca));

                h1m_uls_same_rebin .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same_rebin.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same_rebin.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix_rebin  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix_rebin .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix_rebin .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetName("h1R_pt{0}_dca{1}"  .format(ipt, idca));
                h1bkg.SetName("h1bkg_pt{0}_dca{1}".format(ipt, idca));
                h1sig.SetName("h1sig_pt{0}_dca{1}".format(ipt, idca));
                h1R  .SetTitle("R factor in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1bkg.SetTitle("#it{{N}}_{{ll}}^{{bkg}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1sig.SetTitle("#it{{N}}_{{ll}}^{{signal}} in {0:2.1f} < #it{{p}}_{{T,ll}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ll}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetYTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
                h1bkg.SetYTitle("counts per bin");
                h1sig.SetYTitle("counts per bin");

                h1frac_sig = h1sig.Clone("h1farc_sig_pt{0}_dca{1}".format(ipt, idca));
                h1frac_sig.Reset();
                h1frac_sig.Divide(h1sig, h1m_uls_same_rebin, 1., 1., "B");
                outdir.WriteTObject(h1frac_sig);

                for im in range(0, len(arr_mll)-1):
                    m_min = arr_mll[im];
                    m_max = arr_mll[im+1];
                    m_center = (m_min + m_max)/2;

                    bin_m1 = hs_uls_same.GetAxis(0).FindBin(m_min + __delta);
                    bin_m2 = hs_uls_same.GetAxis(0).FindBin(m_max - __delta);
                    hs_uls_same .GetAxis(0).SetRange(bin_m1, bin_m2);
                    hs_lspp_same.GetAxis(0).SetRange(bin_m1, bin_m2);
                    hs_lsmm_same.GetAxis(0).SetRange(bin_m1, bin_m2);

                    h1uq_uls_same  = hs_uls_same .Projection(3);
                    h1uq_lspp_same = hs_lspp_same.Projection(3);
                    h1uq_lsmm_same = hs_lsmm_same.Projection(3);
                    h1uq_uls_same .SetName("h1uq_uls_same_m{0}_pt{1}_dca{2}" .format(im, ipt, idca));
                    h1uq_lspp_same.SetName("h1uq_lspp_same_m{0}_pt{1}_dca{2}".format(im, ipt, idca));
                    h1uq_lsmm_same.SetName("h1uq_lsmm_same_m{0}_pt{1}_dca{2}".format(im, ipt, idca));
                    h1uq_uls_same .RebinX(2);
                    h1uq_lspp_same.RebinX(2);
                    h1uq_lsmm_same.RebinX(2);

                    frac_sig     = h1frac_sig.GetBinContent(im+1);
                    frac_sig_err = h1frac_sig.GetBinError(im+1);
                    R = h1R.GetBinContent(im+1);
                    R_err = h1R.GetBinError(im+1);
                    [vn_tot, vn_tot_err, vn_bkg, vn_bkg_err, vn_sig, vn_sig_err, h1uq_bkg] = get_flow(h1uq_uls_same, h1uq_lspp_same, h1uq_lsmm_same, frac_sig, frac_sig_err, R, R_err, R2);

                    h1uq_bkg.SetName("h1uq_bkg_m{0}_pt{1}_dca{2}".format(im, ipt, idca));
                    #h1uq_sig.SetName("h1uq_sig_m{0}_pt{1}_dca{2}".format(im, ipt, idca));

                    h1uq_uls_same.SetTitle("uQ in {0:3.2f} < #it{{m}}_{{ll}} < {1:3.2f} GeV/c^{{2}}, {2:2.1f} < #it{{p}}_{{T,ll}} < {3:2.1f} GeV/c, {4:2.1f} < DCA_{{ll}}^{{3D}} < {5:2.1f} #sigma".format(m_min, m_max, pt_min, pt_max, dca_min, dca_max));
                    h1uq_bkg     .SetTitle("uQ in {0:3.2f} < #it{{m}}_{{ll}} < {1:3.2f} GeV/c^{{2}}, {2:2.1f} < #it{{p}}_{{T,ll}} < {3:2.1f} GeV/c, {4:2.1f} < DCA_{{ll}}^{{3D}} < {5:2.1f} #sigma".format(m_min, m_max, pt_min, pt_max, dca_min, dca_max));
                    #h1uq_sig     .SetTitle("uQ in {0:3.2f} < #it{{m}}_{{ll}} < {1:3.2f} GeV/c^{{2}}, {2:2.1f} < #it{{p}}_{{T,ll}} < {3:2.1f} GeV/c, {4:2.1f} < DCA_{{ll}}^{{3D}} < {5:2.1f} #sigma".format(m_min, m_max, pt_min, pt_max, dca_min, dca_max));

                    #outdir.WriteTObject(h1uq_uls_same);
                    #outdir.WriteTObject(h1uq_bkg);
                    #outdir.WriteTObject(h1uq_sig);

                    global_bin_id = hs_uls   .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_R     .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_bkg   .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_sig   .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_vn_tot.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_vn_bkg.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_vn_sig.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);

                    #print(m_center, pt_center, dca_center, global_bin_id);
                    #print("vn_tot: ", h1uq_uls_same.Integral(), vn_tot, vn_tot_err);
                    #print("vn_bkg: ", h1uq_bkg.Integral(), vn_bkg, vn_bkg_err);
                    hs_uls.SetBinContent(global_bin_id, h1m_uls_same_rebin.GetBinContent(im+1));
                    hs_uls.SetBinError(global_bin_id, h1m_uls_same_rebin.GetBinError(im+1));
                    hs_R  .SetBinContent(global_bin_id, h1R.GetBinContent(im+1));
                    hs_R  .SetBinError(global_bin_id, h1R.GetBinError(im+1));
                    hs_bkg.SetBinContent(global_bin_id, h1bkg.GetBinContent(im+1));
                    hs_bkg.SetBinError(global_bin_id, h1bkg.GetBinError(im+1));
                    hs_sig.SetBinContent(global_bin_id, h1sig.GetBinContent(im+1));
                    hs_sig.SetBinError(global_bin_id, h1sig.GetBinError(im+1));
                    h2R.SetBinContent(im+1, ipt+1, h1R.GetBinContent(im+1));
                    h2R.SetBinError(im+1, ipt+1, h1R.GetBinError(im+1));
                    hs_vn_tot.SetBinContent(global_bin_id, vn_tot);
                    hs_vn_tot.SetBinError(global_bin_id  , vn_tot_err);
                    hs_vn_bkg.SetBinContent(global_bin_id, vn_bkg);
                    hs_vn_bkg.SetBinError(global_bin_id  , vn_bkg_err);
                    hs_vn_sig.SetBinContent(global_bin_id, vn_sig);
                    hs_vn_sig.SetBinError(global_bin_id  , vn_sig_err);
                    set_range_00(hs_uls_same , [0]);
                    set_range_00(hs_lspp_same, [0]);
                    set_range_00(hs_lsmm_same, [0]);
                    set_range_00(hs_uls_mix  , [0]);
                    set_range_00(hs_lspp_mix , [0]);
                    set_range_00(hs_lsmm_mix , [0]);

            outdir.WriteTObject(h2R);
        hs_uls.Scale(1/nev);
        hs_bkg.Scale(1/nev);
        hs_sig.Scale(1/nev);
        outdir.WriteTObject(hs_uls);
        outdir.WriteTObject(hs_R  );
        outdir.WriteTObject(hs_bkg);
        outdir.WriteTObject(hs_sig);
        outdir.WriteTObject(hs_vn_tot);
        outdir.WriteTObject(hs_vn_bkg);
        outdir.WriteTObject(hs_vn_sig);

    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
if __name__ == "__main__":
    filename = args.input;

    arr_mll = np.array(config["common"]["mll_bin"], dtype=float);
    arr_ptll = np.array(config["common"]["ptll_bin"], dtype=float);
    arr_dcall = np.array(config["common"]["dcall_bin"], dtype=float);

    if is_mc :
        tasks = config["mc"]["tasks"];
        print(tasks);
        if config["common"]["do_mll_ptll_dcall"]:
            run_efficiency_mll_ptll_dcall(filename, tasks, arr_mll, arr_ptll, arr_dcall, system, energy, suffix);
    else:
        tasks = config["data"]["tasks"];
        #print(tasks);
        if config["common"]["do_mll_ptll_dcall"]:
            run_mll_ptll_dcall(filename, tasks, arr_mll, arr_ptll, arr_dcall);
        if config["common"]["do_flow"]:
            run_flow(filename, tasks, arr_mll, arr_ptll, arr_dcall);

