import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from ROOT_Utils import *
from PFA_Analyzer_Utils import *
from THR_Utils_v2 import *
import CMS_lumi 
import os
import array



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


thr_folder_old = "/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/THR_Data/THR_ARM_DAC/SBit100_Trimming/"

args = parser.parse_args()
enable_THR = args.THR
inputs = [f"{OUTPUT_PATH}/PFA_Analyzer_Output/CSV/{i}/" for i in args.inputs]
run_numbers = [GetRunNumber(i) for i in inputs]
output = args.output if args.output is not None else run_numbers[0]
label_list = args.labels if args.labels is not None else args.inputs
plot_short = args.Short
items_in_the_legend = 0

if len(label_list) != len(inputs):
    print("Parsed inputs and labels are different in number...\nExiting ..")
    print(f"{label_list}")
    print(f"{inputs}")
    sys.exit(0)    


QC8_THR = []
PRIMARY_YAXIS_MIN = 0
PRIMARY_YAXIS_MAX = 105
SECONDARY_YAXIS_MIN = 0.
SECONDARY_YAXIS_MAX = 165.


if plot_short:
    items_in_the_legend+=2 

    chamberNames_withShort = []
    if args.verbose: print(f"Fetching latest GE11 short status")
    cmd = " wget 'https://docs.google.com/spreadsheets/d/1m_OqvUmCpz6ge8rljOpFvAVRUmRCXPSGtEvBphsajBo/gviz/tq?tqx=out:csv&sheet=List of chambers with HV anomalies' -O Output/GE11Short.csv"
    os.system(cmd)
    file1 = open(f'{OUTPUT_PATH}/GE11Short.csv', 'r')
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
            if line[3] == '':
                chamberNames_withShort.append(getChamberName(region_short,chamber_short,1))
                chamberNames_withShort.append(getChamberName(region_short,chamber_short,2))
            else:
                layer_short = int(line[3])
                chamberNames_withShort.append(getChamberName(region_short,chamber_short,layer_short))

            if args.verbose:
                print(f"{chamberNames_withShort[-1]}\t has a short")
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

## Output files
folder_name = "/eos/user/f/fivone/www/P5_Operations/Run3/RunMerge/"+output
print(f"Creating Folder")
CreatEOSFolder(folder_name)

for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    label = EndcapLayer2label(re,la)
    Multigraph_Dict[label]=ROOT.TMultiGraph()
    TGraphError_Dict[label] = {}
    
    Text_Dict[label] = "GE" + "+1/1 Layer"+str(la) if re > 0 else "GE" + "-1/1 Layer"+str(la)
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
color = [ROOT.TColor.GetColorPalette( int(255/(len(inputs))*(i))) for i in range(len(inputs))]
c2 = setUpCanvas("Comparison",1400,900)
c2.Divide(2,2)

for key in Multigraph_Dict.keys():
    Multigraph_Dict[key].GetXaxis().SetNdivisions(80,1,1,True)
    Multigraph_Dict[key].GetXaxis().SetLimits(0.,36.5)
    Multigraph_Dict[key].GetXaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetXaxis().SetTitle("Chamber Number")
    Multigraph_Dict[key].SetMaximum(PRIMARY_YAXIS_MAX)
    Multigraph_Dict[key].SetMinimum(PRIMARY_YAXIS_MIN)
    Multigraph_Dict[key].GetYaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetYaxis().SetTickLength(0.015)
    Multigraph_Dict[key].GetYaxis().SetTitleOffset(0.9)
    Multigraph_Dict[key].GetYaxis().SetTitle("Efficiency [%]")
## Operative loop
for index,file_path in enumerate(inputs):
    #df = pd.read_csv(file_path, sep=',')
    for re in [-1,1]:
        for la in [1,2]:
            
            label = EndcapLayer2label(re,la)           

            x,y,exl,exh,eyl,eyh = [],[],[],[],[],[]
            for chamber in range(1,37):
                chamberID = getChamberName(re,chamber,la)
                matched,propagated,efficiency_CL68 = ChamberEfficiencyFromCSV(file_path,chamberID)
                if propagated == 0:
                    continue

                else:
                    center_eff = 100*float(matched)/float(propagated)

                x.append(chamber)
                y.append(center_eff)
                eyh.append(efficiency_CL68[1]*100 - center_eff)
                eyl.append(center_eff - 100*efficiency_CL68[0])
                exl.append(0.5)
                exh.append(0.5)
                
                if plot_short and chamberID in chamberNames_withShort:
                    MultiHist_Dict[label].SetBinContent(chamber,PRIMARY_YAXIS_MAX)
                if plot_short and chamberID in QC8_THR:
                    ChamberWQC8_THR[label].SetBinContent(chamber,PRIMARY_YAXIS_MAX)
            TGraphError_Dict[label][index] = ROOT.TGraphAsymmErrors(    len(x),
                                                                        array.array('d',x),
                                                                        array.array('d',y),
                                                                        array.array('d',exl),
                                                                        array.array('d',exh),
                                                                        array.array('d',eyl),
                                                                        array.array('d',eyh)
                                                                    )
            TGraphError_Dict[label][index].SetName(label_list[index])
            TGraphError_Dict[label][index].SetTitle(label_list[index])
            TGraphError_Dict[label][index].SetMarkerColor(color[index])
            TGraphError_Dict[label][index].SetMarkerStyle(20)
            TGraphError_Dict[label][index].SetLineColor(color[index])
            TGraphError_Dict[label][index].SetFillColorAlpha(ROOT.kBlue,0.99)
            TGraphError_Dict[label][index].SetMarkerSize(1)
            TGraphError_Dict[label][index].SetLineWidth(2)

            Multigraph_Dict[label].Add(TGraphError_Dict[label][index])

if enable_THR:
    items_in_the_legend += 1
    print(f"Adding a secondary y-axis for THR comparison... using data from")
    print(f"\tTHR_Utils_v2.py")

    secondary_axis = ROOT.TGaxis(36.5,PRIMARY_YAXIS_MIN,36.5, PRIMARY_YAXIS_MAX,SECONDARY_YAXIS_MIN,SECONDARY_YAXIS_MAX,505,"+LS")
    secondary_axis.SetLineColor(ROOT.kMagenta)
    secondary_axis.SetLabelColor(ROOT.kMagenta)
    secondary_axis.SetTitleColor(ROOT.kMagenta)
    secondary_axis.SetTitle("DAC THR Delta(0Hz - 100Hz)")
    secondary_axis.SetTickLength(0.015)
    for region,layer in [(-1,1),(1,1),(-1,2),(1,2)]:
        label = EndcapLayer2label(region,layer)
        print(f"{label}")
        for chamber in range(1,37):
            
            chamberID = getChamberName(region,chamber,layer)
            thr = GetOverallChamberThreshold(chamberID,run=357329)
            graph_point = thr
            if graph_point is  None:
                graph_point = 10**4
            scaled_secondary_y_point = (PRIMARY_YAXIS_MAX-PRIMARY_YAXIS_MIN)*(float(graph_point)/(SECONDARY_YAXIS_MAX-SECONDARY_YAXIS_MIN))

            TGraph_THR_Dict[label].SetPoint(chamber-1,chamber,scaled_secondary_y_point)
        
        TGraph_THR_Dict[label].SetTitle("Chamber Overall THR")
        TGraph_THR_Dict[label].SetName("Chamber Overall THR")
        TGraph_THR_Dict[label].SetMarkerStyle(45)
        TGraph_THR_Dict[label].SetMarkerColor(ROOT.kMagenta)
        Multigraph_Dict[label].Add(TGraph_THR_Dict[label])

leg_list = {}

for index,key in enumerate(Multigraph_Dict.keys()):
    virtual_Pad = c2.cd(index+1)
    virtual_Pad.SetGrid()    
    Multigraph_Dict[key].Draw("0AEP")        
    CMS_lumi.CMS_lumi(virtual_Pad,  0,  0,"(13.6 TeV)")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)    
    
    latex.SetTextFont(61)
    latex.SetTextAlign(21) 
    latex.SetTextSize(0.045)    
    latex.DrawLatex(.5,0.909,Text_Dict[key])
    leg_list[key] = ROOT.TLegend()
    leg_list[key].SetBorderSize(2)
    leg_list[key].SetFillStyle(1)
    for j,file_path in enumerate(inputs):
        leg_list[key].AddEntry(TGraphError_Dict[key][j],TGraphError_Dict[key][j].GetTitle(),"PE")

    if enable_THR: 
        c2.cd(index+1)
        secondary_axis.Draw()
        leg_list[key].AddEntry(TGraph_THR_Dict[key],TGraph_THR_Dict[key].GetName(),"P")

    if plot_short:
        c2.cd(index+1)
        MultiHist_Dict[key].Draw("SAME")
        leg_list[key].AddEntry(MultiHist_Dict[key],"Chamber with short in GEM Foil HV sector","F")
        
    leg_list[key].Draw()

c2.Modified()
c2.Update()



outputPath = folder_name+"/Merged.pdf"
c2.SaveAs(outputPath)
Convert2png(outputPath)
print(f"Your ouput \t{outputPath}")
