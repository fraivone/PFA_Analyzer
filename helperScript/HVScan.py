import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter

from ROOT_Utils import *
from PFA_Analyzer_Utils import *
import os
import array



parser = argparse.ArgumentParser(
        description='''Scripts that prouces HV Scan plotsstarting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.HVScan --inputs 360019_700 360019_600 --HV 700 600 """,
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--HV', type=int , help="HV values",required=True,nargs="*")
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)


args = parser.parse_args()
inputs = ["./Output/PFA_Analyzer_Output/CSV/"+i+"/" for i in args.inputs]
hv_list = [ k for k in args.HV ]

if len(hv_list) != len(inputs):
    print("Parsed inputs and labels are different in number...\nExiting ..")
    print(f"{hv_list}")
    print(f"{inputs}")
    sys.exit(0)    


QC8_THR = []
PRIMARY_XAXIS_MIN = 580
PRIMARY_XAXIS_MAX = 700
PRIMARY_YAXIS_MIN = 0
PRIMARY_YAXIS_MAX = 105


### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)

ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0) #Don't print titles
## ROOT Objects
TGraphError_Dict = {}
Text_Dict = {}
TGraph_THR_Dict = {}

for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    label = EndcapLayer2label(re,la)
    TGraphError_Dict[label] = {}

## ROOTStyle
color = [ROOT.TColor.GetColorPalette( int(255/(len(inputs))*(i))) for i in range(len(inputs))]
c2 = setUpCanvas("Comparison",1400,1400)
c2.SetGrid()

## Operative loop

for re in [-1,1]:
    for la in [1,2]:
        label = EndcapLayer2label(re,la)           
        for chamber in range(1,37):
            x,y,exl,exh,eyl,eyh = [],[],[],[],[],[]
            chamberID = getChamberName(re,chamber,la)
            for index,file_path in enumerate(inputs):
                HV = hv_list[index]
                matched,propagated,efficiency_CL68 = ChamberEfficiencyFromCSV(file_path,chamberID)
                if propagated == 0:
                    continue
                else:
                    center_eff = 100*float(matched)/float(propagated)

                x.append(HV)
                y.append(center_eff)
                eyh.append(efficiency_CL68[1]*100 - center_eff)
                eyl.append(center_eff - 100*efficiency_CL68[0])
                exl.append(0.5)
                exh.append(0.5)
            
            if len(x) == 0: continue
            TGraphError_Dict[label][chamberID] = ROOT.TGraphAsymmErrors(    len(x),
                                                                    array.array('d',x),
                                                                    array.array('d',y),
                                                                    array.array('d',exl),
                                                                    array.array('d',exh),
                                                                    array.array('d',eyl),
                                                                    array.array('d',eyh)
                                                                )
            TGraphError_Dict[label][chamberID].SetName(chamberID)
            TGraphError_Dict[label][chamberID].SetTitle(chamberID)
            TGraphError_Dict[label][chamberID].SetMarkerColor(ROOT.kBlack)
            TGraphError_Dict[label][chamberID].SetLineColor(ROOT.kBlack)
            TGraphError_Dict[label][chamberID].SetMarkerStyle(20)
            TGraphError_Dict[label][chamberID].SetFillColorAlpha(ROOT.kBlack,0.5)
            TGraphError_Dict[label][chamberID].SetMarkerSize(1)
            TGraphError_Dict[label][chamberID].SetLineWidth(2)
            TGraphError_Dict[label][chamberID].GetXaxis().SetLimits(PRIMARY_XAXIS_MIN,PRIMARY_XAXIS_MAX)
            TGraphError_Dict[label][chamberID].GetXaxis().SetLabelSize(0.03)
            TGraphError_Dict[label][chamberID].GetXaxis().SetTitle("Equivalent Divider Current (uA)")
            TGraphError_Dict[label][chamberID].SetMaximum(PRIMARY_YAXIS_MAX)
            TGraphError_Dict[label][chamberID].SetMinimum(PRIMARY_YAXIS_MIN)
            TGraphError_Dict[label][chamberID].GetYaxis().SetLabelSize(0.03)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTickLength(0.015)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTitleOffset(0.9)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTitle("Efficiency [%]")

            c2.cd()
            TGraphError_Dict[label][chamberID].Draw("APEC")
            c2.Modified()
            c2.Update()

            ## Output files
            folder_name = "/eos/user/f/fivone/www/P5_Operations/Run3/RunMerge/HVScan/"+label
            CreatEOSFolder(folder_name)
            outputPath = folder_name+"/"+chamberID +".pdf"
            c2.SaveAs(outputPath)
            Convert2png(outputPath)