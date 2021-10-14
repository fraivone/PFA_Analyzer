import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.THR_Utils import *

parser = argparse.ArgumentParser(
        description='''Scripts that prouces efficiency summary plot of different runs, starting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.CompareCSV --inputs day26_344420_690uA  day27_344421_690uA  day28_344423_690uA day29_344518_690uA --output 690uA_NoTrim""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--output', type=str , help="Output file name",required=False)
parser.add_argument('--labels', type=str , help="Label with which the runs should be listed in the legend (according to inputs order). If not provided, input names will be used",required=False,nargs='*')
parser.add_argument('--THR', default=False,action='store_true', help="Enables THR comparison plot on the secondary y axis",required=False)
parser.add_argument('--verbose', default=False, action='store_true',help="Verbose printing",required=False)
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)


thr_folder_old = "/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/THR_Data/THR_ARM_DAC/SBit100_Trimming/"

args = parser.parse_args()
enable_THR = args.THR
inputs = ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi.csv" for i in args.inputs]
output = args.output if args.output is not None else "CompareCSV"
label_list = args.labels if args.labels is not None else args.inputs

if len(label_list) != len(inputs):
    print "Parsed inputs and labels are different in number...\nExiting .."
    sys.exit(0)    

### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)
if not args.verbose: ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=2001;") #suppressed everything less-than-or-equal-to kWarning
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0) #Don't print titles
## ROOT Objects
Multigraph_Dict = OrderedDict() 
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

## ROOTStyle
color = [ROOT.TColor.GetColorPalette(255/(len(inputs))*(i)) for i in range(len(inputs))]
c2 = setUpCanvas("Comparison",1400,900)
c2.Divide(2,2)

for key in Multigraph_Dict.keys():
    Multigraph_Dict[key].GetXaxis().SetNdivisions(80,1,1,True)
    Multigraph_Dict[key].GetXaxis().SetLimits(0.,36.5)
    Multigraph_Dict[key].GetXaxis().SetLabelSize(0.03)
    Multigraph_Dict[key].GetXaxis().SetTitle("Chamber Number")
    Multigraph_Dict[key].SetMinimum(0)
    Multigraph_Dict[key].SetMaximum(1.1)
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

            TGraphError_Dict[label][index].Divide(Num_TH1F,Den_TH1F,"B")
            TGraphError_Dict[label][index].SetName(label_list[index])
            TGraphError_Dict[label][index].SetTitle(label_list[index])
            TGraphError_Dict[label][index].SetMarkerColor(color[index])
            TGraphError_Dict[label][index].SetMarkerStyle(20)
            TGraphError_Dict[label][index].SetLineColor(color[index])
            TGraphError_Dict[label][index].SetLineWidth(2)

            Multigraph_Dict[label].Add(TGraphError_Dict[label][index])

if enable_THR:
    print "Adding a secondary y-axis for THR comparison... using data from"
    print "\t"+thr_folder_old

    secondary_axis = ROOT.TGaxis(36.5,0,36.5, 1.1,0,255,505,"+LS")
    secondary_axis.SetLineColor(ROOT.kRed)
    secondary_axis.SetLabelColor(ROOT.kRed)
    secondary_axis.SetTitle("Chamber Overall THR")
    secondary_axis.SetTickLength(0.015)
    for region,layer in [(-1,1),(1,1),(-1,2),(1,2)]:
        label = EndcapLayer2label(region,layer)
        print label
        for chamber in range(1,37):
            
            chamberID = ReChLa2chamberName(region,chamber,layer)
            old_thr = GetOverallChamberThreshold(thr_folder_old,chamberID,verbose=args.verbose)
            if old_thr == None : graph_point = 0
            else: graph_point = old_thr
            
            TGraph_THR_Dict[label].SetPoint(chamber-1,chamber,1.1*float(graph_point)/255)
            print graph_point
        
        TGraph_THR_Dict[label].SetTitle("Chamber Overall THR")
        TGraph_THR_Dict[label].SetName("Chamber Overall THR")
        TGraph_THR_Dict[label].SetMarkerStyle(45)
        TGraph_THR_Dict[label].SetMarkerColor(ROOT.kRed)
        Multigraph_Dict[label].Add(TGraph_THR_Dict[label])


for index,key in enumerate(Multigraph_Dict.keys()):
    c2.cd(index+1).SetGrid()
    Multigraph_Dict[key].Draw("0APE")
    
    leg = c2.cd(index+1).BuildLegend(.1,.8,.4,1.)
    leg.SetBorderSize(2)
    if enable_THR: 
        secondary_axis.Draw()
        leg.GetListOfPrimitives()[-1].SetOption("P")

    Text_Dict[key].Draw()


c2.Modified()
c2.Update()
c2.SaveAs("./Output/CompareCSV/"+output+".pdf")
print "Your ouput \t","./Output/CompareCSV/"+output+".pdf"