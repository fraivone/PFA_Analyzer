import ROOT
import argparse
from argparse import RawTextHelpFormatter
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.THR_Utils import *
import pandas as pd
import subprocess
import numpy as np
from lib.THR_Utils_v2 import *

parser = argparse.ArgumentParser(
        description='''Scripts that prouces efficiency per VFAT/Chamber''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.PlotEfficiencyVFAT --input 350107_Express --output 350107 --batch""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--input', type=str , help="Tag of the run to be plotted",required=True)
parser.add_argument('--verbose', default=False, action='store_true',help="Verbose printing",required=False)
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode. OFF by default",required=False)
parser.add_argument('--THR', default=False,action='store_true', help="Enables THR comparison plot on the secondary y axis",required=False)
parser.add_argument('--output', type=str , help="Output file name",required=False)



args = parser.parse_args()
enable_THR = args.THR
inputTag = args.input
inputFile = "./Output/PFA_Analyzer_Output/CSV/"+inputTag+"/MatchingSummary_glb_rdphi_byVFAT.csv"
run_number = GetRunNumber(inputTag)
output = args.output if args.output is not None else run_number
outputFolder = "/eos/user/f/fivone/www/run"+output+"/VFAT/"
CreatEOSFolder(outputFolder)

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


for re,la in [(-1,2),(1,1),(-1,1),(1,2)]: 
    
    endcapTag = EndcapLayer2label(re,la)
    CreatEOSFolder(outputFolder+"/"+endcapTag)
    print "Processing ", endcapTag
    Text_Dict[endcapTag] = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 "+endcapTag, "brNDC")
    Text_Dict[endcapTag].SetBorderSize(0)
    Text_Dict[endcapTag].SetFillColor(ROOT.gStyle.GetTitleFillColor())
    Text_Dict[endcapTag].SetFillStyle(0)

    PositionDict[endcapTag] = {}


    for chamber in range(1,37):
        current_chamber_ID = ReChLa2chamberName(re,chamber,la)
        current_output_folder = outputFolder+"/"+endcapTag + "/" + current_chamber_ID +"/"
        CreatEOSFolder(current_output_folder)

        if enable_THR:
            THR = []
            VFAT_N = []
        if args.verbose: print "Processing ", current_chamber_ID
        for VFATN in range(24):
            if current_chamber_ID in df['chamberID'].tolist():
                EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['num'] = df[ (df['chamberID']==current_chamber_ID) &(df['VFATN']==VFATN)]["matchedRecHit"].values[0]
                EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['den'] = df[ (df['chamberID']==current_chamber_ID) &(df['VFATN']==VFATN)]["propHit"].values[0]
            else: 
                EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['num'] = 0
                EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFATN]['den'] = 0

            if enable_THR:
                THR.append(GetVFAT_THRDAC(current_chamber_ID,VFATN))
                VFAT_N.append(VFATN)

            angle_deg = (chamber - 1)*10## degrees
            angle_rad = np.radians(angle_deg)

            PositionDict[endcapTag][current_chamber_ID] = ROOT.TPaveLabel(0.5+0.195*np.cos(angle_rad)-0.015, 0.52+0.195*np.sin(angle_rad)-0.015, 0.5+0.195*np.cos(angle_rad)+0.015,0.52+0.195*np.sin(angle_rad)+0.015,str(chamber), "brNDC")
            PositionDict[endcapTag][current_chamber_ID].SetTextAngle(angle_deg-90)
            PositionDict[endcapTag][current_chamber_ID].SetBorderSize(0)
            PositionDict[endcapTag][current_chamber_ID].SetFillColor(ROOT.gStyle.GetTitleFillColor())
            PositionDict[endcapTag][current_chamber_ID].SetFillStyle(0)
        
        oneD = generate1DEfficiencyPlotbyVFAT(EfficiencyDictVFAT['glb_rdphi'],current_chamber_ID)
        

        twoD = generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT['glb_rdphi'],current_chamber_ID)
        writeToTFile(OutF,oneD,endcapTag+"/")
        writeToTFile(OutF,twoD,endcapTag+"/")
        oneD.Draw("APE")
        
        if enable_THR:
            y_max = oneD.GetMaximum()
            THR = [float(i)*y_max/255 for i in THR]
            oneD_Thr = ROOT.TGraph(len(VFAT_N),np.asarray(VFAT_N,dtype=float),np.asarray(THR,dtype=float))
            oneD_Thr.SetTitle(current_chamber_ID+" THR DAC")
            oneD_Thr.SetName(current_chamber_ID+ " THR DAC")
            oneD_Thr.SetMarkerStyle(45)
            oneD_Thr.SetMarkerSize(3)
            oneD_Thr.SetMarkerColor(ROOT.kMagenta)

            secondary_axis = ROOT.TGaxis(23.5,0,23.5, y_max,0,255,505,"+LS")
            secondary_axis.SetLineColor(ROOT.kMagenta)
            secondary_axis.SetLabelColor(ROOT.kMagenta)
            secondary_axis.SetTitleColor(ROOT.kMagenta)
            secondary_axis.SetTitle("VFAT THR (DAC)")
            secondary_axis.SetTickLength(0.015)
            secondary_axis.SetTitleOffset(1.35)

            secondary_axis.Draw()
            oneD_Thr.Draw("SAME P")
        
        c1.Modified()
        c1.Update()
        c1.SaveAs(current_output_folder+"Eff1D_"+current_chamber_ID+".pdf")
        Convert2png(current_output_folder+"Eff1D_"+current_chamber_ID+".pdf")

        twoD.Draw("COLZ TEXT")
        c1.Modified()
        c1.Update()
        c1.SaveAs(current_output_folder+"Eff2D_"+current_chamber_ID+".pdf")
        Convert2png(current_output_folder+"Eff2D_"+current_chamber_ID+".pdf")
        print current_chamber_ID+" done..."

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
    Convert2png(outputFolder+endcapTag+".pdf")
    


print "Your ouput plots \t"+outputFolder
print "Your output ROOT \t",outputFolder+"NoCUT_Summary.root"