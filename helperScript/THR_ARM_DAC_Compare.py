import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import ROOT
import math
import numpy as np
import os.path
import sys
from collections import OrderedDict
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.DAQ_Chamber_Mapping import *


parser = argparse.ArgumentParser(
        description='''Scripts that given that plos AVG THR ARM DAC per chamber''',
        epilog="""Typical execution (from top folder)\n\t python -m helperScript.THR_ARM_DAC_Compare --inputs SBit100 SBit100_Trimming SBit10000_Trimming""",
        formatter_class=RawTextHelpFormatter
)

parser.add_argument('--inputs', type=str,help="path to the top folder where THR ARM DAC are stored",required=True,nargs='*')
parser.add_argument('--labels', type=str , help="Label with which the folder should be listed in the legend (according to inputs order). If not provided, input names will be used",required=False,nargs='*')

args = parser.parse_args()

inputs = ["./helperScript/THR_ARM_DAC/"+i for i in args.inputs]
label_list = args.labels if args.labels is not None else args.inputs
n_inputs = len(inputs)
n_label = len(label_list)
if n_label != n_inputs:
    print "Parsed inputs and labels are different in number...\nExiting .."
    sys.exit(0)


color = [ROOT.TColor.GetColorPalette(255/n_inputs*(i)) for i in range(n_inputs)]

### ROOT style settings
ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0) #Don't print titles

## ROOT Objects
c2 = setUpCanvas("Comparison",1400,900)
c2.Divide(2,2)
Multigraph_Dict = OrderedDict()
TGraphError_Dict = {}
Text_Dict = {}
for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    endcap_layer = EndcapLayer2label(re,la)
    Multigraph_Dict[endcap_layer]=ROOT.TMultiGraph("A","A")
    TGraphError_Dict[endcap_layer] = {}
    Text_Dict[endcap_layer] = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 "+endcap_layer, "brNDC")


## ROOTStyle
color = [ROOT.TColor.GetColorPalette(255/(len(inputs))*(i)) for i in range(len(inputs))]
for key in Multigraph_Dict.keys():
    Multigraph_Dict[key].GetXaxis().SetNdivisions(80,1,1,True)
    Multigraph_Dict[key].GetXaxis().SetLimits(0.,36.5)
    Multigraph_Dict[key].GetXaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetXaxis().SetTitle("Chamber Number")
    Multigraph_Dict[key].SetMinimum(0)
    Multigraph_Dict[key].SetMaximum(150)
    Multigraph_Dict[key].GetYaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetYaxis().SetTitleOffset(0.9)
    Multigraph_Dict[key].GetYaxis().SetTitle("AVG THR ARM DAC")

    Text_Dict[key].SetBorderSize(0)
    Text_Dict[key].SetFillColor(ROOT.gStyle.GetTitleFillColor())
    Text_Dict[key].SetFillStyle(0)

for index,input_folder in enumerate(inputs):
    legend_name = label_list[index]
    for re,layer in [(-1,1),(1,1),(-1,2),(1,2)]:
        x_data = []
        y_data = []
        ey_data = []
        ex_data = []

        endcap_layer = EndcapLayer2label(re,layer)

        for ch in range(1,37):
            chamber_ID = ReChLa2chamberName(re,ch,layer)
            
            tupl = chamber_mapping[chamber_ID]
            crate = int(tupl[0])
            amc = int(tupl[1])
            OHLink = int(tupl[2])
 
            file_path = input_folder+"/crate%02d" % crate + "-amc%02d" % amc +".txt"
            data = pd.read_csv(file_path, sep=':')
            mean = data[data["OH/I"]==OHLink]["threshold/I"].mean()
            if math.isnan(mean) == False:
                x_data.append(ch)
                y_data.append(mean)
                ex_data.append(0)
                ey_data.append(0)
        
        n = len(x_data)
        x = np.asarray(x_data,dtype=float)
        y = np.asarray(y_data,dtype=float)
        ex = np.asarray(ex_data,dtype=float)
        ey = np.asarray(ey_data,dtype=float)
        print endcap_layer, legend_name, np.mean(y)

        TGraphError_Dict[endcap_layer][legend_name] = ROOT.TGraphErrors(n,x,y,ex,ey)
        TGraphError_Dict[endcap_layer][legend_name].SetMarkerColor(color[index])
        TGraphError_Dict[endcap_layer][legend_name].SetTitle(legend_name)
        TGraphError_Dict[endcap_layer][legend_name].SetName(legend_name)
        TGraphError_Dict[endcap_layer][legend_name].SetMarkerStyle(21)
        Multigraph_Dict[endcap_layer].Add(TGraphError_Dict[endcap_layer][legend_name])
        
for index,key in enumerate(Multigraph_Dict.keys()):
    c2.cd(index+1).SetGrid()
    Multigraph_Dict[key].Draw("0APE")
    c2.cd(index+1).BuildLegend().SetBorderSize(2)
    Text_Dict[key].Draw()

c2.Modified()
c2.Update()
c2.SaveAs("./Output/CompareTHR.pdf")
raw_input()
