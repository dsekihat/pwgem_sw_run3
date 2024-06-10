#!/bin/bash

#python run_lmee.py -i AnalysisResults_HL_222576.root -t data -c dielectron/configs/config_pp_13.6TeV.yml -s "";
python run_lmee.py -i AnalysisResults_test.root -t data -c dielectron/configs/config_PbPb_5.36TeV_5070.yml -s "_5070";
python run_lmee.py -i AnalysisResults_test.root -t data -c dielectron/configs/config_PbPb_5.36TeV_7090.yml -s "_7090";
python run_lmee.py -i AnalysisResults_HL_5090.root -t data -c dielectron/configs/config_PbPb_5.36TeV_5090.yml -s "_5090";
