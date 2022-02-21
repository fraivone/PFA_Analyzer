import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import ROOT
import math
import numpy as np


data = pd.read_csv("./TRIM_ARM_DAC/TrimmingData.csv")

Plot2D = ROOT.TH2F("AVG TRIM DAC","AVG TRIM DAC",128,-0.5,127.5,24,-0.5,23.5)
Plot2D.GetXaxis().SetTitle("Channel Number")
Plot2D.GetYaxis().SetTitle("VFAT_N")
Plot2D.SetMinimum(-64)
Plot2D.SetMaximum(64)
Plot2D.SetStats(False)
c1 = ROOT.TCanvas("c1","c1")

VFATN = data["VFATN"].tolist()
CHANNEL = data["Channel"].tolist()
TRIM_DAC = data["TRIM_DAC"].tolist()
for vfat in range( 24 ):
    print vfat
    for channel in range(128):        
        mean_TRIM_DAC = data[((data["Channel"] == channel) & (data["VFATN"] == vfat))]["TRIM_DAC"].mean()
        Plot2D.SetBinContent(channel+1,vfat+1,mean_TRIM_DAC)
        

c1.cd()
Plot2D.Draw("COLZ")
Plot2D.SaveAs("t.png")
c1.Modified()
c1.Update()
raw_input()
