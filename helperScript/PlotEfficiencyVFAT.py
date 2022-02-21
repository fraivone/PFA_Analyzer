import ROOT
import argparse
from argparse import RawTextHelpFormatter
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.THR_Utils import *
import pandas as pd
import subprocess
import numpy as np


parser = argparse.ArgumentParser(
        description='''Scripts that prouces efficiency per VFAT/Chamber''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.PlotEfficiencyVFAT --input 690uA_MergedRuns_NoTrim_byVFAT""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--input', type=str , help="Tag of the run to be plotted",required=True)
parser.add_argument('--verbose', default=False, action='store_true',help="Verbose printing",required=False)
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)



args = parser.parse_args()
inputTag = args.input
inputFile = "./Output/PFA_Analyzer_Output/CSV/"+inputTag+"/MatchingSummary_glb_rdphi_byVFAT.csv"

outputFolder = "Output/PFA_Analyzer_Output/Plot/"+inputTag+"/"
subprocess.call(["mkdir", "-p", outputFolder])

OutF = ROOT.TFile(outputFolder+"Summary.root","RECREATE")
Text_Dict = {}
PositionDict = {}

### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)
if not args.verbose: ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=2001;") #suppressed everything less-than-or-equal-to kWarning
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetOptTitle(0) #Don't print titles
## ROOT Objects


c1 = setUpCanvas("Eff_byVFAT",2000,2000)                     
c1.SetLeftMargin(0.1)
c1.SetRightMargin(0.1)      


df = pd.read_csv(inputFile, sep=',')
EfficiencyDictVFAT = generateVFATDict(['glb_rdphi'])


for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    
    endcapTag = EndcapLayer2label(re,la)
    Text_Dict[endcapTag] = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 "+endcapTag, "brNDC")
    Text_Dict[endcapTag].SetBorderSize(0)
    Text_Dict[endcapTag].SetFillColor(ROOT.gStyle.GetTitleFillColor())
    Text_Dict[endcapTag].SetFillStyle(0)

    PositionDict[endcapTag] = {}


    for chamber in range(1,37):
        current_chamber_ID = ReChLa2chamberName(re,chamber,la)
        for VFATN in range(24):
            EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['num'] = df[ (df['chamberID']==current_chamber_ID) &(df['VFATN']==VFATN)]["matchedRecHit"].values[0]
            EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['den'] = df[ (df['chamberID']==current_chamber_ID) &(df['VFATN']==VFATN)]["propHit"].values[0]


            angle_deg = (chamber - 1)*10## degrees
            angle_rad = np.radians(angle_deg)

            PositionDict[endcapTag][current_chamber_ID] = ROOT.TPaveLabel(0.5+0.195*np.cos(angle_rad)-0.015, 0.52+0.195*np.sin(angle_rad)-0.015, 0.5+0.195*np.cos(angle_rad)+0.015,0.52+0.195*np.sin(angle_rad)+0.015,str(chamber), "brNDC")
            PositionDict[endcapTag][current_chamber_ID].SetTextAngle(angle_deg-90)
            PositionDict[endcapTag][current_chamber_ID].SetBorderSize(0)
            PositionDict[endcapTag][current_chamber_ID].SetFillColor(ROOT.gStyle.GetTitleFillColor())
            PositionDict[endcapTag][current_chamber_ID].SetFillStyle(0)
        
        writeToTFile(OutF,generate1DEfficiencyPlotbyVFAT(EfficiencyDictVFAT['glb_rdphi'],current_chamber_ID),endcapTag+"/")
        writeToTFile(OutF,generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT['glb_rdphi'],current_chamber_ID),endcapTag+"/")

for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    endcapTag = EndcapLayer2label(re,la)
    s = generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT["glb_rdphi"],endcapTag)
    s.Draw("COLZ TEXT")
    Text_Dict[endcapTag].Draw()
    for chamberID in PositionDict[endcapTag].keys():
        PositionDict[endcapTag][chamberID].Draw()

    c1.Modified()
    c1.Update()
    c1.SaveAs(outputFolder+endcapTag+".pdf")


print "Your ouput plots \t"+outputFolder
print "Your output ROOT \t",outputFolder+"Summary.root"
raw_input("Press any key to terminate...")