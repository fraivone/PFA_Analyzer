import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import ROOT
import math
import numpy as np
import os.path
import sys

from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.DAQ_Chamber_Mapping import *


"""
Script to compare the final AVG THR applied to the front-end
Taking into account both the general THR ARM DAC (set per VFAT)
and the single channel TRIM DAC (set per channel)
Typical execution (from top folder) python -m helperScript.Overall_VFAT_Threshold
"""

THR_ARM_folder = "/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/helperScript/THR_ARM_DAC/SBit100_Trimming/"
THR_ARM_old = "/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/helperScript/THR_ARM_DAC/SBit100"

## ROOT
ROOT.gStyle.SetOptTitle(0) #Don't print titles
ROOT.gStyle.SetLineScalePS(1)
c = ROOT.TCanvas("C1","C1")

data = pd.read_csv("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/helperScript/TRIM_ARM_DAC/TrimmingData.csv")
chamber_list = list( set( data["Chamber_ID"].tolist() ) )


AfterTrimming = ROOT.TH1F("Overall After Trimming","Overall After Trimming",255,0,255)
BeforeTrimming = ROOT.TH1F("Overall Before Trimming","Overall Before Trimming",255,0,255)
BeforeTrimming.SetLineColor(ROOT.kRed)

AfterTrimmingfC = ROOT.TH1F("Overall After Trimming","Overall After Trimming",255,0,31)
BeforeTrimmingfC = ROOT.TH1F("Overall Before Trimming","Overall Before Trimming",255,0,31)
BeforeTrimmingfC.SetLineColor(ROOT.kRed)

for chamber in sorted(chamber_list):
    print chamber
    for VFAT_N in range(24):
        # Trimming data
        pandas_column_with_trim_fc = data[( (data["Chamber_ID"] == chamber) & (data["VFATN"] == VFAT_N) )]["TRIM_fC"]
        # THR data
        try:
            THR_ARM_DAC = GetVFAT_THRDAC(THR_ARM_folder,chamber,VFAT_N)
        except:
            continue
        overall_thresholdsfc =  pandas_column_with_trim_fc + THR_ARM_DAC*THRDAC2fC
        
        overall_thresholds_dac =  pandas_column_with_trim_fc / THRDAC2fC + THR_ARM_DAC  ## Converting TRIM fC into DAC THR (conversion factors in PFA_Analyzer_Utils)

        AfterTrimming.Fill(overall_thresholds_dac.mean())
        AfterTrimmingfC.Fill(overall_thresholdsfc.mean())


for chamber in sorted(chamber_list):
    print chamber
    for VFAT_N in range(24):
        # THR data
        try:
            THR_ARM_DAC = GetVFAT_THRDAC(THR_ARM_old,chamber,VFAT_N)
        except:
            continue
        overall_thresholds = THR_ARM_DAC
        BeforeTrimming.Fill(overall_thresholds)
        BeforeTrimmingfC.Fill(overall_thresholds*THRDAC2fC)


BeforeTrimming.GetXaxis().SetTitle("Global THR in DAC Units")
BeforeTrimming.GetYaxis().SetTitle("Entries")
BeforeTrimming.SetStats(False)
BeforeTrimming.SetFillColor(ROOT.kRed)
BeforeTrimming.SetFillStyle(3003)
BeforeTrimming.SetTitle("w/o Trim AVG THR = "+str(round(BeforeTrimming.GetMean(),2))+" DAC Units")

AfterTrimming.SetTitle("w/ Trim AVG THR = "+str(round(AfterTrimming.GetMean(),2))+" DAC Units")
AfterTrimming.SetStats(False)
AfterTrimming.SetFillColor(ROOT.kBlue)
AfterTrimming.SetFillStyle(3003)

BeforeTrimming.Draw("HIST F")
AfterTrimming.Draw("SAME HIST F")

c.BuildLegend(0.5,0.5,0.95,0.9)
c.Modified()
c.Update()

c.SaveAs("Output/Overall_THR_perVFAT.pdf")
raw_input("Press any key to continue")


BeforeTrimmingfC.GetXaxis().SetTitle("Global THR in fC")
BeforeTrimmingfC.GetYaxis().SetTitle("Entries")
BeforeTrimmingfC.SetStats(False)
BeforeTrimmingfC.SetFillColor(ROOT.kRed)
BeforeTrimmingfC.SetFillStyle(3003)
BeforeTrimmingfC.SetTitle("w/o Trim AVG THR = "+str(round(BeforeTrimmingfC.GetMean(),2))+"fC")

AfterTrimmingfC.SetTitle("w/ Trim AVG THR = "+str(round(AfterTrimmingfC.GetMean(),2))+" fC")
AfterTrimmingfC.SetStats(False)
AfterTrimmingfC.SetFillColor(ROOT.kBlue)
AfterTrimmingfC.SetFillStyle(3003)

BeforeTrimmingfC.Draw("HIST F")
AfterTrimmingfC.Draw("SAME HIST F")

c.BuildLegend(0.5,0.5,0.95,0.9)
c.Modified()
c.Update()
c.SaveAs("Output/Overall_fC_perVFAT.pdf")

raw_input("Press any key to close")
