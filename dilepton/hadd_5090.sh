#!/bin/bash

rootcp -r AnalysisResults_test.root:dielectron-qc_5070 AnalysisResults_HL_5070.root:dielectron-qc_5090
rootcp -r AnalysisResults_test.root:dielectron-qc_5070_TOFreq AnalysisResults_HL_5070.root:dielectron-qc_5090_TOFreq

rootcp -r AnalysisResults_test.root:dielectron-qc_7090 AnalysisResults_HL_7090.root:dielectron-qc_5090
rootcp -r AnalysisResults_test.root:dielectron-qc_7090_TOFreq AnalysisResults_HL_7090.root:dielectron-qc_5090_TOFreq

hadd -f AnalysisResults_HL_5090.root AnalysisResults_HL_5070.root AnalysisResults_HL_7090.root;

