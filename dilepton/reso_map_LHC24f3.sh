#!/bin/bash



ls AnalysisResults_526641.root AnalysisResults_526964.root AnalysisResults_527041.root AnalysisResults_527057.root AnalysisResults_527109.root AnalysisResults_527240.root AnalysisResults_527850.root AnalysisResults_527871.root AnalysisResults_527895.root AnalysisResults_527899.root AnalysisResults_528292.root AnalysisResults_528531.root | xargs hadd -f AnalysisResults_reso_map_LHC24f3.root
