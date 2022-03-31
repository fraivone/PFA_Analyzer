import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.THR_Utils_v2 import *
import os

parser = argparse.ArgumentParser(
        description='''Scripts that prouces efficiency summary plot of different runs, starting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.CompareCSV --inputs day26_344420_690uA  day27_344421_690uA  day28_344423_690uA day29_344518_690uA --output 690uA_NoTrim""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--output', type=str , help="Output file name",required=False)
parser.add_argument('--labels', type=str , help="Label with which the runs should be listed in the legend (according to inputs order). If not provided, input names will be used",required=False,nargs='*')
parser.add_argument('--THR', default=False,action='store_true', help="Enables THR comparison plot on the secondary y axis",required=False)
parser.add_argument('--Short', default=False,action='store_true', help="Enables plotting of the chambers in short",required=False)
parser.add_argument('--verbose', default=False, action='store_true',help="Verbose printing",required=False)
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)


args = parser.parse_args()
enable_THR = args.THR
inputs = ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi.csv" for i in args.inputs]
run_numbers = [GetRunNumber(i) for i in inputs]
output = args.output if args.output is not None else run_numbers[0]
label_list = args.labels if args.labels is not None else args.inputs
plot_short = args.Short
items_in_the_legend = 0

if len(label_list) != len(inputs):
    print "Parsed inputs and labels are different in number...\nExiting .."
    print label_list
    print inputs
    sys.exit(0)    


QC8_THR = ["GE11-M-31L1-S","GE11-M-31L2-S","GE11-P-07L1-S","GE11-P-07L2-S","GE11-P-11L1-S","GE11-P-11L2-S","GE11-P-25L1-S","GE11-P-25L2-S","GE11-P-27L1-S","GE11-P-27L2-S","GE11-P-29L1-S","GE11-P-29L2-S"]


if plot_short:
    items_in_the_legend+=2 

    chamberNames_withShort = []
    if args.verbose: print "Fetching latest GE11 short status"
    cmd = " wget 'https://docs.google.com/spreadsheets/d/1m_OqvUmCpz6ge8rljOpFvAVRUmRCXPSGtEvBphsajBo/gviz/tq?tqx=out:csv&sheet=List of chambers with HV anomalies' -O Output/GE11Short.csv"
    os.system(cmd)
    file1 = open('Output/GE11Short.csv', 'r')
    Lines = file1.readlines()
    for line in Lines: 
        line = line.replace("\"","")
        line = line.strip()
        line = line.split(",")
        if "GE" in line[1] and "disappeared" not in line[0] and "anomaly" not in line[0] and line[2].lstrip("-").isdigit():  #check if the short is there, documented in its layer and chamber
            if "-" in line[1]: region_short = -1
            else: region_short = 1
            line[1] = line[1].replace(" ","")
            chamber_short = int(line[1][-2:])
            layer_short = int(line[3])
    
            chamberNames_withShort.append(ReChLa2chamberName(region_short,chamber_short,layer_short))
            if args.verbose:
                print chamberNames_withShort[-1],"\t has a short"
### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)
if not args.verbose: ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=2001;") #suppressed everything less-than-or-equal-to kWarning
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0) #Don't print titles
## ROOT Objects
Multigraph_Dict = OrderedDict() 
MultiHist_Dict = OrderedDict() 
ChamberWQC8_THR = OrderedDict() 
TGraphError_Dict = {}
Text_Dict = {}
TGraph_THR_Dict = {}
for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    label = EndcapLayer2label(re,la)
    Multigraph_Dict[label]=ROOT.TMultiGraph()
    TGraphError_Dict[label] = {}
    Text_Dict[label] = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 "+label, "brNDC")
    if enable_THR:
        TGraph_THR_Dict[label] = ROOT.TGraph()
    if plot_short:
        MultiHist_Dict[label] = ROOT.TH1F("HasShort","HasShort",36,0.5,36.5)
        MultiHist_Dict[label].SetFillColorAlpha(ROOT.kGray+3,0.35)
        MultiHist_Dict[label].SetLineColor(ROOT.kGray+3)

        ChamberWQC8_THR[label] = ROOT.TH1F("QC8THR","QC8THR",36,0.5,36.5)
        ChamberWQC8_THR[label].SetFillColorAlpha(ROOT.kGreen+3,0.15)
        ChamberWQC8_THR[label].SetLineColor(ROOT.kGreen+3)

## ROOTStyle
color = [ROOT.TColor.GetColorPalette(255/(len(inputs))*(i)) for i in range(len(inputs))]
c2 = setUpCanvas("Comparison",1400,900)
c2.Divide(2,2)

for key in Multigraph_Dict.keys():
    Multigraph_Dict[key].GetXaxis().SetNdivisions(80,1,1,True)
    Multigraph_Dict[key].GetXaxis().SetLimits(0.,36.5)
    Multigraph_Dict[key].GetXaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetXaxis().SetTitle("Chamber Number")
    Multigraph_Dict[key].SetMaximum(1.1)
    Multigraph_Dict[key].SetMinimum(0.)
    Multigraph_Dict[key].GetYaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetYaxis().SetTickLength(0.015)
    Multigraph_Dict[key].GetYaxis().SetTitleOffset(0.9)
    Multigraph_Dict[key].GetYaxis().SetTitle("MIP Efficiency")

    Text_Dict[key].SetBorderSize(0)
    Text_Dict[key].SetFillColor(ROOT.gStyle.GetTitleFillColor())
    Text_Dict[key].SetFillStyle(0)

## Operative loop
for index,file_path in enumerate(inputs):
    df = pd.read_csv(file_path, sep=',')
    for re in [-1,1]:
        for la in [1,2]:
            
            Num_TH1F=ROOT.TH1F("Denominator","Denominator",36,0.5,36.5)
            Den_TH1F=ROOT.TH1F("Denominator","Denominator",36,0.5,36.5)
            label = EndcapLayer2label(re,la)            

            TGraphError_Dict[label][index] = ROOT.TGraphAsymmErrors()
            

            for chamber in range(1,37):
                chamberID = ReChLa2chamberName(re,chamber,la)
                matched = sum(df[ df['chamberID']==chamberID ]["matchedRecHit"].tolist())
                prop = sum(df[ df['chamberID']==chamberID ]["propHit"].tolist())
                if prop != 0:
                    Num_TH1F.SetBinContent(chamber,matched)
                    Den_TH1F.SetBinContent(chamber,prop)
                
                if plot_short and chamberID in chamberNames_withShort:
                    MultiHist_Dict[label].SetBinContent(chamber,1.1)
                if plot_short and chamberID in QC8_THR:
                    ChamberWQC8_THR[label].SetBinContent(chamber,1.1)

            TGraphError_Dict[label][index].Divide(Num_TH1F,Den_TH1F,"B")
            TGraphError_Dict[label][index].SetName(label_list[index])
            TGraphError_Dict[label][index].SetTitle(label_list[index])
            TGraphError_Dict[label][index].SetMarkerColor(color[index])
            TGraphError_Dict[label][index].SetMarkerStyle(20)
            TGraphError_Dict[label][index].SetLineColor(color[index])
            TGraphError_Dict[label][index].SetLineWidth(2)

            Multigraph_Dict[label].Add(TGraphError_Dict[label][index])

if enable_THR:
    items_in_the_legend += 1
    print "Adding a secondary y-axis for THR comparison... using data from"
    print "\t"+"THR_Utils_v2.py"

    secondary_axis = ROOT.TGaxis(36.5,0,36.5, 1.1,0,255,505,"+LS")
    secondary_axis.SetLineColor(ROOT.kMagenta)
    secondary_axis.SetLabelColor(ROOT.kMagenta)
    secondary_axis.SetTitleColor(ROOT.kMagenta)
    secondary_axis.SetTitle("Chamber Overall THR")
    secondary_axis.SetTickLength(0.015)
    for region,layer in [(-1,1),(1,1),(-1,2),(1,2)]:
        label = EndcapLayer2label(region,layer)
        print label
        for chamber in range(1,37):
            
            chamberID = ReChLa2chamberName(region,chamber,layer)
            old_thr = GetOverallChamberThreshold(chamberID)
            if old_thr == None : graph_point = 0
            else: graph_point = old_thr
            
            TGraph_THR_Dict[label].SetPoint(chamber-1,chamber,1.1*float(graph_point)/255)
        
        TGraph_THR_Dict[label].SetTitle("Chamber Overall THR")
        TGraph_THR_Dict[label].SetName("Chamber Overall THR")
        TGraph_THR_Dict[label].SetMarkerStyle(45)
        TGraph_THR_Dict[label].SetMarkerColor(ROOT.kMagenta)
        Multigraph_Dict[label].Add(TGraph_THR_Dict[label])

leg_list = {}
for index,key in enumerate(Multigraph_Dict.keys()):
    c2.cd(index+1).SetGrid()    
    Multigraph_Dict[key].Draw("0APE")
    
    leg_list[key] = ROOT.TLegend(0.1,0.1,0.3,0.2+0.03*(len(inputs)-1 + items_in_the_legend))
    leg_list[key].SetBorderSize(2)
    leg_list[key].SetFillStyle(0)
    for j,file_path in enumerate(inputs):
        leg_list[key].AddEntry(TGraphError_Dict[key][j],TGraphError_Dict[key][j].GetName(),"PE")

    if enable_THR: 
        secondary_axis.Draw()
        leg_list[key].AddEntry(TGraph_THR_Dict[key],TGraph_THR_Dict[key].GetName(),"P")
    if plot_short:
        MultiHist_Dict[key].Draw("SAME")
        ChamberWQC8_THR[key].Draw("same")
        leg_list[key].AddEntry(MultiHist_Dict[key],MultiHist_Dict[key].GetName(),"F")
        leg_list[key].AddEntry(ChamberWQC8_THR[key],ChamberWQC8_THR[key].GetName(),"F")
        
        
    leg_list[key].Draw()
    


    Text_Dict[key].Draw()


c2.Modified()
c2.Update()

folder_name = "/eos/user/f/fivone/www/run"+output
CreatEOSFolder(folder_name)

outputPath = folder_name+"/"+output+".pdf"
c2.SaveAs(outputPath)
Convert2png(outputPath)
print "Your ouput \t",outputPath
