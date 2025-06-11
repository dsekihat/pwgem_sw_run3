#!/bin/bash


python3 run_dilepton.py -i AnalysisResults_HL_414591_414810.root -t data -c dimuon/configs/config_pp_5.36TeV.yml -s "_HL_414591_414810"; # for dimuon LHC24ap + aq pass1
#python3 run_dilepton.py -i AnalysisResults_HL_414810.root -t data -c dimuon/configs/config_pp_5.36TeV.yml -s "_HL_414810"; # for dimuon LHC24aq pass1
#python3 run_dilepton.py -i AnalysisResults_HL_414591.root -t data -c dimuon/configs/config_pp_5.36TeV.yml -s "_HL_414591"; # for dimuon LHC24ap pass1
#python3 run_dilepton.py -i AnalysisResults_HL_405933.root -t data -c dielectron/configs/config_pp_5.36TeV.yml -s "_HL_405933"; # for dielectron LHC24aq pass1
#python3 run_dilepton.py -i AnalysisResults_HL_405917.root -t data -c dielectron/configs/config_pp_5.36TeV.yml -s "_HL_405917"; # for dielectron LHC24ap pass1
#python3 run_dilepton.py -i AnalysisResults_HL_395023.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_395023"; # for dimuon LHC23_Thin pass4 with TTCA weighting
#python3 run_dilepton.py -i AnalysisResults_HL_394896_394927.root -t data -c dimuon/configs/config_pp_5.36TeV.yml -s "_HL_394896_394827"; # for dimuon LHC24ap,aq pass1
#python3 run_dilepton.py -i AnalysisResults_HL_394826.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_394826"; # for dimuon LHC23_Thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_389401.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_389401"; # for dimuon LHC22o MB pass7
#python3 run_dilepton.py -i AnalysisResults_HL_385409_385423.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_385409_385423"; # for dimuon LHC24ap,aq pass1
#python3 run_dilepton.py -i AnalysisResults_HL_349170.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_349170"; # for dielectron LHC23 thin pass4, woTTCA, wWeighting, default
#python3 run_hbt.py -i AnalysisResults_HL_337703.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_337703"; #with deta-dphi cut > 0.15
#python3 run_hbt.py -i AnalysisResults_HL_336553.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_336553"; #without deta-dphi cut
#python3 run_dilepton.py -i AnalysisResults_HL_328095_328470_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_328095_330348"; #2023 PbPb
#python3 run_dilepton.py -i AnalysisResults_HL_329829_330348_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_329829_330348"; #2024 PbPb, LHC24ar + LHC24as
#python3 run_dilepton.py -i AnalysisResults_HL_330348_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_330348"; #2024 PbPb, LHC24as
#python3 run_dilepton.py -i AnalysisResults_HL_329829_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_329829"; #2024 PbPb, LHC24ar

#python3 run_dilepton.py -i AnalysisResults_HL_328076.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_328076";
#python3 run_dilepton.py -i AnalysisResults_HL_321612.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_321612";
#python3 run_dilepton.py -i AnalysisResults_HL_321577.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_321577";
#python3 run_dilepton.py -i AnalysisResults_HL_320499_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_320499"; # with prefilter deta > 0.04, with deta > 0.04, ITS cluster size < 3
#python3 run_dilepton.py -i AnalysisResults_HL_320321.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_320321";
#python3 run_dilepton.py -i AnalysisResults_HL_320161_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_320161"; # with prefilter deta > 0.02, with deta > 0.02
#python3 run_dilepton.py -i AnalysisResults_HL_320049_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_320049"; # with prefilter deta > 0.03, with deta > 0.03
#python3 run_dilepton.py -i AnalysisResults_HL_319736_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_319736"; # with prefilter deta > 0.04, with deta > 0.04
#python3 run_dilepton.py -i AnalysisResults_HL_319559_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_319559"; # without prefilter, without detadphi
#python3 run_dilepton.py -i AnalysisResults_HL_319300_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_319300"; # prefilter deta > 0.03
#python3 run_hbt.py -i AnalysisResults_HL_318188.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_318188";
#python3 run_hbt.py -i AnalysisResults_HL_317846.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_317846";
#python3 run_hbt.py -i AnalysisResults_HL_316859.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_316859";
#python3 run_dilepton.py -i AnalysisResults_HL_316686_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_316686"; # prefilter deta > 0.03
#python3 run_hbt.py -i AnalysisResults_HL_316404.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_316404";
#python3 run_hbt.py -i AnalysisResults_HL_315730.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_315730";
#python3 run_hbt.py -i AnalysisResults_HL_311370.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_311370";
#python3 run_dilepton.py -i AnalysisResults_HL_313263_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_313263"; # deta 0.04 cut
#python3 run_dilepton.py -i AnalysisResults_HL_312076.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_312076"; #rel diff pin < 0.15
#python3 run_dilepton.py -i AnalysisResults_HL_312813.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_312813"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_310124.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_310124"; # deta 0.08, dphi > 1e+10 with prefilter
#python3 run_hbt.py -i AnalysisResults_HL_310198.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_310198";
#python3 run_hbt.py -i AnalysisResults_HL_309434.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_309434";
#python3 run_dilepton.py -i AnalysisResults_HL_309527_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_309527"; # deta 0.04, dphi > 1e+10 with prefilter
#python3 run_dilepton.py -i AnalysisResults_HL_308878.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_308878"; # deta 0.04, dphi > 0.2 with prefilter
#python3 run_dilepton.py -i AnalysisResults_HL_307922.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_307922"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_307487.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_307487"; # centrality bin for mixing is every 1%
#python3 run_dilepton.py -i AnalysisResults_HL_306762.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_306762"; # centrality bin for mixing is every 1%
#python3 run_dilepton.py -i AnalysisResults_HL_306301.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_306301"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_302828.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_302828"; # PIDTPC, -1 < n sigma e TPC < 2
#python3 run_dilepton.py -i AnalysisResults_HL_299686_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_299686"; # PIDTPC, -0.5 < n sigma e TPC < 2
#python3 run_dilepton.py -i AnalysisResults_HL_299220_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_299220"; # PIDTPC, -1 < n sigma e TPC < 2
#python3 run_dilepton.py -i AnalysisResults_HL_298365.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_298365"; # PIDTPC, -1 < n sigma e TPC < 2
#python3 run_dilepton.py -i AnalysisResults_HL_299065.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_299065";
#python3 run_dilepton.py -i AnalysisResults_HL_299066.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_299066";
#python3 run_dilepton.py -i AnalysisResults_HL_296817.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_296817_flow";
#python3 run_dilepton.py -i AnalysisResults_HL_296716.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_296716_flow";
#python3 run_dilepton.py -i AnalysisResults_HL_296817.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_296817";
#python3 run_dilepton.py -i AnalysisResults_HL_296716.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_296716";
#python3 run_dilepton.py -i AnalysisResults_HL_295155.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_295155";
#python3 run_dilepton.py -i AnalysisResults_HL_295259.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_295259";
#python3 run_dilepton.py -i AnalysisResults_HL_294526_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_294526";
#python3 run_dilepton.py -i AnalysisResults_HL_294420_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_294420";
#python3 run_dilepton.py -i AnalysisResults_HL_294001.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_294001";
#python3 run_dilepton.py -i AnalysisResults_HL_289080.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_289080"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_291104.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_291104";
#python3 run_dilepton.py -i AnalysisResults_HL_290986.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_290986";
#python3 run_dilepton.py -i AnalysisResults_HL_289636.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_289636";
#python3 run_dilepton.py -i AnalysisResults_HL_286547.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_286547";
#python3 run_dilepton.py -i AnalysisResults_HL_288465.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_288465";
#python3 run_dilepton.py -i AnalysisResults_HL_275995.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_275995"; # for dielectron pass7
#python3 run_dilepton.py -i AnalysisResults_HL_284904.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_284904";
#python3 run_dilepton.py -i AnalysisResults_HL_286293.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_286293";
#python3 run_dilepton.py -i AnalysisResults_HL_285248.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_285248";
#python3 run_dilepton.py -i AnalysisResults_HL_285692.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_285692";
#python3 run_dilepton.py -i AnalysisResults_HL_284904.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_284904_narrow_m";
#python3 run_dilepton.py -i AnalysisResults_HL_284904.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_284904_ptintegrated";
#python3 run_dilepton.py -i AnalysisResults_HL_284904.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_284904_narrow";
#python3 run_dilepton.py -i AnalysisResults_HL_284153.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_284153";
#python3 run_dilepton.py -i AnalysisResults_HL_282709.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_282709";
#python3 run_dilepton.py -i AnalysisResults_HL_282608.root -t mc -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_282608"; # for dielectron LHC24k4g
#python3 run_dilepton.py -i AnalysisResults_HL_282384.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_282384";
#python3 run_dilepton.py -i AnalysisResults_HL_281025.root -t mc -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_281025"; # for dielectron LHC24k4g
#python3 run_dilepton.py -i AnalysisResults_HL_281955.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_281955";
#python3 run_dilepton.py -i AnalysisResults_HL_281906.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_281906"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_281428.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_281428";
#python3 run_dilepton.py -i AnalysisResults_HL_280891.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_280891"; # for dielectron LHC23 thin pass4
#python3 run_dilepton.py -i AnalysisResults_HL_280731.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_280731";
#python3 run_dilepton.py -i AnalysisResults_test_detadphi.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_test_detadphi";
#python3 run_dilepton.py -i AnalysisResults_HL_280225_merged_centrality.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_280225_merged_centrality";
#python3 run_dilepton.py -i AnalysisResults_HL_280225.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_280225";
#python3 run_dilepton.py -i AnalysisResults_HL_279716.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_279716";
#python3 run_dilepton.py -i AnalysisResults_HL_275984.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_275984";
#python3 run_dilepton.py -i AnalysisResults_HL_277601.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_277601";
#python3 run_dilepton.py -i AnalysisResults_HL_277106.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_277106";
#python3 run_dilepton.py -i AnalysisResults_HL_276138.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_276138";
#python3 run_dilepton.py -i AnalysisResults_HL_274137.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_HL_274137"; #first attempt of deta-dphi cut
#python3 run_dilepton.py -i AnalysisResults_HL_274632.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_274632";
#python3 run_dilepton.py -i AnalysisResults_HL_273587.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_273587";
#python3 run_hbt.py -i AnalysisResults_HL_266490.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_266490";
#python3 run_hbt.py -i AnalysisResults_HL_265475.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_265475";
#python3 run_hbt.py -i AnalysisResults_HL_263811.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_263811";
#python3 run_dilepton.py -i AnalysisResults_HL_243092_merged_centrality_global.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_243092_all"; # for dimuon full stat. for HP24
#python3 run_dilepton.py -i AnalysisResults_HL_243092_merged_centrality_standalone.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_243092"; # for dimuon full stat. for HP24
#python3 run_hbt.py -i AnalysisResults_HL_263299.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_263299";
#python3 run_hbt.py -i AnalysisResults_HL_263031.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_263031";
#python3 run_dilepton.py -i AnalysisResults_HL_262447.root -t data -c dielectron/configs/config_PbPb_5.36TeV_6080_UPC.yml -s "_HL_262447";
#python3 run_dilepton.py -i AnalysisResults_HL_261772.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_261772";
#python3 run_dilepton.py -i AnalysisResults_HL_259881.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_259881";
#python3 run_dilepton.py -i AnalysisResults_HL_254915.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090_jj.yml -s "_HL_254915_woMixOcc";
#python3 run_dilepton.py -i AnalysisResults_HL_253953.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_253953";
#python3 run_dilepton.py -i AnalysisResults_HL_253788.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_253788";
#python3 run_dilepton.py -i AnalysisResults_HL_253469.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_253469"; # for dielectron pass7
#python3 run_hbt.py -i AnalysisResults_HL_251833.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_251833";
#python3 run_dilepton.py -i AnalysisResults_test.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090_jj.yml -s "_test_R";
#python3 run_dilepton.py -i AnalysisResults_HL_250594.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090_jj.yml -s "_HL_250594_jj";
#python3 run_dilepton.py -i AnalysisResults_HL_251678.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_251678";
#python3 run_hbt.py -i AnalysisResults_HL_250608.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_250608";
#python3 run_dilepton.py -i AnalysisResults_HL_249691.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_249691"; # for dielectron pass7
#python3 run_dilepton.py -i AnalysisResults_HL_249691.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_249691"; # for dimuon pass7
#python3 run_hbt.py -i AnalysisResults_HL_250160.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_250160";
#python3 run_dilepton.py -i AnalysisResults_tmp_v2.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_tmp_v2";
#python3 run_dilepton.py -i AnalysisResults_tmp_forR.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_tmp_forR"; # for dielectron pass7
#python3 run_dilepton.py -i AnalysisResults_HL_248963.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_248963";
#python3 run_dilepton.py -i AnalysisResults_HL_248396.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_248396"; # for dielectron pass7
#python3 run_dilepton.py -i AnalysisResults_HL_248396.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_248396"; # for dimuon pass7
#python3 run_hbt.py -i AnalysisResults_HL_248545.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_248545";
#python3 run_dilepton.py -i AnalysisResults_HL_247038.root -t mc -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_247038"; # for dielectron LHC24f3
#python3 run_hbt.py -i AnalysisResults_HL_247995.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_247995";
#python3 run_dilepton.py -i AnalysisResults_HL_246171.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_246171";
#python3 run_hbt.py -i AnalysisResults_HL_247241.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_HL_247241";
#python3 run_hbt.py -i AnalysisResults_HL_247198.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_247198";
#python3 run_dilepton.py -i AnalysisResults_HL_247038_LHC22o_minbias.root -t mc -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_247038"; # for dielectron LHC24f3
#python3 run_dilepton.py -i AnalysisResults_HL_245648.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_HL_245648"; # for dielectron pass7
#python3 run_dilepton.py -i AnalysisResults_HL_245648.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_245648"; # for dimuon pass7
#python3 run_dilepton.py -i AnalysisResults_HL_246171.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_246171"; # for dimuon pass4
#python3 run_dilepton.py -i AnalysisResults_HL_245689.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_245689";
#python3 run_dilepton.py -i AnalysisResults_HL_245160.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_245160";
#python3 run_hbt.py -i AnalysisResults_HL_244379.root -t data -c hbt/configs/config_pp_13.6TeV_HBT.yml -s "_HL_244379";
#python3 run_hbt.py -i AnalysisResults_test_hbt.root -t data -c hbt/configs/config_PbPb_5.36TeV_0090_HBT.yml -s "_test_hbt";
#python3 run_dilepton.py -i AnalysisResults_HL_244958.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_244958";
#python3 run_dilepton.py -i AnalysisResults_test_v2.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_test_v2";
#python3 run_dilepton.py -i AnalysisResults_HL_243492.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "_HL_243492";
#python3 run_dilepton.py -i AnalysisResults_HL_242847.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_242847"; # for dimuon LHC22o pass6
#python3 run_dilepton.py -i AnalysisResults_HL_243092.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_243092"; # for dimuon full stat.
#python3 run_dilepton.py -i AnalysisResults_HL_242210_242211.root -t data -c dielectron/configs/config_PbPb_5.36TeV_3050.yml -s "";
#python3 run_lmee.py -i AnalysisResults_HL_233073.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_233073_big_ptmumu";
#python3 run_lmee.py -i AnalysisResults_HL_233073.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_233073";
#python3 run_lmee.py -i AnalysisResults_HL_233027.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_233027_big_ptmumu";
#python3 run_lmee.py -i AnalysisResults_HL_233027.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_233027_tmp";
#python3 run_lmee.py -i AnalysisResults_HL_232981.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_232981";
#python3 run_lmee.py -i AnalysisResults_HL_232159.root -t data -c dimuon/configs/config_PbPb_5.36TeV.yml -s "_HL_232159";
#python3 run_lmee.py -i AnalysisResults_HL_232390.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_232390";
#python3 run_lmee.py -i AnalysisResults_HL_232106.root -t data -c dimuon/configs/config_pp_13.6TeV.yml -s "_HL_232106";
#python3 run_lmee.py -i AnalysisResults_HL_PbPb_LHC23zzh_pass3_20240617_merged_centrality_occupancy.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_occupancy_1bigbin";
#python3 run_lmee.py -i AnalysisResults_HL_224982.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_224982_0090_occupancy";
#python3 run_lmee.py -i AnalysisResults_HL_224921.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_224921_noDCA_forR";
#python3 run_lmee.py -i AnalysisResults_HL_224874.root -t data -c dielectron/configs/config_PbPb_5.36TeV_0090.yml -s "_224874_0090_TOFreq";
#python3 run_lmee.py -i AnalysisResults_HL_224398.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_224398_noDCA";
#python3 run_lmee.py -i AnalysisResults_HL_222576.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "_noDCA";
#python3 run_lmee.py -i AnalysisResults_HL_223383.root -t data -c dielectron/configs/config_PbPb_5.36TeV_5070.yml -s "_5070";
#python3 run_lmee.py -i AnalysisResults_HL_223383.root -t data -c dielectron/configs/config_PbPb_5.36TeV_7090.yml -s "_7090";
#python3 run_lmee.py -i AnalysisResults_HL_5090.root -t data -c dielectron/configs/config_PbPb_5.36TeV_5090.yml -s "_5090";
