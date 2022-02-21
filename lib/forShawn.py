import ROOT
import csv
import os.path
from os import mkdir
import subprocess
import numpy as np
import math
import sys
import time
import argparse
import pandas as pd
from argparse import RawTextHelpFormatter
### Let's add some more from different folder
lib_folder = os.path.expandvars('$myLIB')
sys.path.insert(1, lib_folder)
try:
    from ROOT_Utils import *
    

except:
  print("ERROR:\n\tCan't find the package CMS_lumi and tdrstlye\n\tPlease verify that this file are placed in the path $myLIB/ROOT_Utils/ \n\tAdditionally keep in mind to export the environmental variable $myLIB\nEXITING...\n") 
  sys.exit(0)
try:
    from PFA_Analyzer_Utils import *
except:
    print ("ERROR:\n\tCan't find the package PFA_Analyzer_Utils\nEXITING...\n")
    sys.exit(0)


parser = argparse.ArgumentParser(
        description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Checks for chamber having GEMRecHits with BX!=0''',
        epilog="""Typical exectuion\n\t python forShawn.py  --dataset /eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/CRUZET/Run_343496/""",
        formatter_class=RawTextHelpFormatter
)
parser.add_argument('--dataset','-ds', type=str,help="Path to the folder containing the NTuples to be analyzed",required=True)
args = parser.parse_args()

ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gROOT.SetBatch(True)

start_time = time.time()
## Input data
folder = args.dataset
Run_Number = folder.split("/Run_")[1].replace("/","")
print "####\t Run ",Run_Number,"\t####"
files = files_in_folder(folder)
files = [f for f in files if ".root" in f]

chain = ROOT.TChain("muNtupleProducer/MuDPGTree")
for fl in files:
    chain.Add(fl)
MaxLS = int(chain.GetMaximum("event_lumiBlock"))
## Input data 

## Container for bad chambers and LS 
badChamberInLS = {}
for LS in range(1,MaxLS+1):
    badChamberInLS.setdefault(LS,[])
## Container for bad chambers and LS 

## OutputFile 
OutF = ROOT.TFile("./BXInfo_Run_"+Run_Number+".root","RECREATE")
## OutputFile 


## Info on Bins

# On X axis the bin N corresponds to LS N-1
# On Y axis the bin N corresponds to LS N-1 - 4000
binyz = 4001
binybot = 1
binytop = 8000

## Info on Bins

## Identify bad chambers (speeds up the execution)
histo = ROOT.TH2F("Temp","Temp",36,0.5,36.5,8000,-4000,4000)
chambersToBeChecked = []
for region in [-1,1]:
    for layer in [1,2]:
        selection = "gemRecHit_region == "+str(region)+" && gemRecHit_layer =="+str(layer)
        chain.Draw("gemRecHit_bx:gemRecHit_chamber >>Temp",selection, "goff")

        for ch in range(1,37):
            numberOfBXNonZero = histo.Integral(ch,ch,binyz+1,binytop) + histo.Integral(ch,ch,binybot,binyz-1)
            if numberOfBXNonZero!=0:
                print ReChLa2chamberName(region,ch,layer)+"\thas Delta BX!=0"
                chambersToBeChecked.append([region,ch,layer])

## Identify bad chambers (speeds up the execution)


## Operative Loop
for identifier in chambersToBeChecked:
    region = identifier[0]
    chamber = identifier[1]
    layer = identifier[2]
    

    chID = ReChLa2chamberName(region,chamber,layer)
    print "Processing \t"+chID


    histo = ROOT.TH2F("TempTitle","TempTitle",MaxLS,0.5,MaxLS+0.5,8000,-4000,4000)
    histo.SetStats(False)
    histo.GetXaxis().SetTitle("LumiSection")
    histo.GetYaxis().SetTitle("GEMRecHit #DeltaBX")

    
    selection = "gemRecHit_chamber == "+str(chamber)+" && gemRecHit_region == "+str(region)+" && gemRecHit_layer =="+str(layer)
    chain.Draw("gemRecHit_bx:event_lumiBlock >> TempTitle" ,selection, "goff")
    histo.SetName(chID)
    histo.SetTitle(chID)
    totalEntriesWithBXNonZero = histo.Integral(1,MaxLS+1,binyz+1,binytop) + histo.Integral(1,MaxLS+1,binybot,binyz-1) 
    if totalEntriesWithBXNonZero!= 0:
        for LS in range(1,MaxLS+1):
            numberOfBXNonZero = histo.Integral(LS,LS,binyz+1,binytop) + histo.Integral(LS,LS,binybot,binyz-1) 
            if numberOfBXNonZero != 0:
                badChamberInLS[LS].append(chID)
    ## Rebinning the Histo for sake of readability
    histo.RebinX(20)
    histo.RebinY(100)
    writeToTFile(OutF,histo,"BadChambers")
## Operative Loop

## Summary Plots
graph = ROOT.TGraph()
graph.SetName("SummaryTGraph")
graph.SetTitle("SummaryTGraph")
graph.SetMarkerStyle(7)
graph.SetMarkerColor(ROOT.kRed)
graph.GetXaxis().SetTitle("LumiSection")
graph.GetYaxis().SetTitle("Number of Ch w/ #DeltaBX!=0")


SummaryTH2 = {}
for region in [-1,1]:
    for layer in [1,2]:
        endcap = "M" if region == -1 else "P"
        key = endcap + "L" + str(layer)
        SummaryTH2[key] = ROOT.TH2F("BX!=0 for GE11-"+key,"BX!=0 for GE11-"+key,MaxLS,0.5,MaxLS+0.5,36,0.5,36.5)
        SummaryTH2[key].SetStats(False)
        SummaryTH2[key].GetXaxis().SetTitle("LumiSection")
        SummaryTH2[key].GetYaxis().SetTitle("Chamber Number")
        SummaryTH2[key].GetYaxis().SetNdivisions(36,0,0)
        SummaryTH2[key].GetYaxis().SetTickLength(0.005)
        SummaryTH2[key].GetYaxis().SetLabelSize(0.03)

for LS in range(1,MaxLS+1):
    graph.SetPoint(LS-1,LS,len(badChamberInLS[LS]))
    badChambersInLS = list(set(badChamberInLS[LS]))
    for chID in badChambersInLS:
        region,chamber,layer = chamberName2ReChLa(chID)[0],chamberName2ReChLa(chID)[1],chamberName2ReChLa(chID)[2]
        endcap = "M" if region == -1 else "P"
        key = endcap + "L" + str(layer)
        SummaryTH2[key].Fill(LS,chamber)

c1 = setUpCanvas("SummaryCanvas",1600,900)

c1.Divide(2,2)
c1.cd(1)
ROOT.gPad.SetGrid()
SummaryTH2['ML1'].Draw("COLZ")
c1.cd(2)
ROOT.gPad.SetGrid()
SummaryTH2['PL1'].Draw("COLZ")
c1.cd(3)
ROOT.gPad.SetGrid()
SummaryTH2['ML2'].Draw("COLZ")
c1.cd(4)
ROOT.gPad.SetGrid()
SummaryTH2['PL2'].Draw("COLZ")

## Adding text
c1.cd()
tempPad=ROOT.TPad("tempPad","tempPad",0,0,1,1)
tempPad.SetFillStyle(4000)
tempPad.Draw()
tempPad.cd()
s = "Run "+str(Run_Number)
t = ROOT.TLatex(0.44,0.48,s)
t.SetTextSize(0.04)
t.Draw()


c1.Modified()
c1.Update()

writeToTFile(OutF,graph)
writeToTFile(OutF,SummaryTH2['ML1'])
writeToTFile(OutF,SummaryTH2['ML2'])
writeToTFile(OutF,SummaryTH2['PL1'])
writeToTFile(OutF,SummaryTH2['PL2'])
writeToTFile(OutF,c1)
c1.SaveAs("./Plot/Run_"+Run_Number+"_BXSummaryCanvas.png")
## Summary Plots


print("--- %s seconds ---" % (time.time() - start_time))
