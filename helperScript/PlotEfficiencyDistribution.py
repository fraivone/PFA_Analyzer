import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from ROOT_Utils import *
from PFA_Analyzer_Utils import *
import CMS_lumi 

ROOT.gStyle.SetOptTitle(0) #Don't print titles
parser = argparse.ArgumentParser(
        description='''Scripts that produces efficiency distribution different runs, starting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.PlotEfficiencyDistribution --inputs MergedRuns_700uA_100HZ_Trim  MergedRuns_690uA_100HZ_Trim  --output EfficiencyDistributuion""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--output', type=str , help="Output file name",required=True)
parser.add_argument('--labels', type=str , help="Label with which the runs should be listed in the legend (according to inputs order). If not provided, input names will be used",required=False,nargs='*')
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)

chamber_with_short =[] # ["GE11-P-07L2-S", "GE11-P-09L2-S", "GE11-P-10L2-L", "GE11-P-12L1-L", "GE11-P-14L1-L", "GE11-P-15L1-S", "GE11-P-15L2-S", "GE11-P-18L1-L", "GE11-P-24L2-L", "GE11-P-34L2-L", "GE11-P-36L1-L", "GE11-P-36L2-L", "GE11-M-05L2-S", "GE11-M-05L1-S", "GE11-M-06L1-L", "GE11-M-07L1-S", "GE11-M-10L2-L", "GE11-M-17L2-S", "GE11-M-19L2-S", "GE11-M-21L1-S", "GE11-M-23L1-S", "GE11-M-29L1-S", "GE11-M-30L1-L", "GE11-M-31L2-S", "GE11-M-33L1-S", "GE11-M-35L1-S"] 

args = parser.parse_args()
inputs = [f"{OUTPUT_PATH}/PFA_Analyzer_Output/CSV/"+i+"/" for i in args.inputs]
output = args.output
label_list = args.labels if args.labels is not None else args.inputs

if len(label_list) != len(inputs):
    print("Parsed inputs and labels are different in number...\nExiting ..")
    sys.exit(0)    

### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
#ROOT.gStyle.SetOptTitle(0) #Don't print titles

## ROOT Objects
TH1F_Container = [ i for i in range(len(inputs))]
Cumulative_Container = [ i for i in range(len(inputs))]
hs = ROOT.THStack("hs",output+"_EfficiencyDistribution")
hCumulative = ROOT.THStack("hCumulative",output+"_CumulativeDistribution")
hCumulative.SetMaximum(144)



## ROOTStyle
color = [ROOT.kBlue +2 , ROOT.kGreen + 1]
c2 = setUpCanvas("Comparison",1400,900)
c2.SetTopMargin(0.1)
c2.SetBottomMargin(0.1)



## Operative loop
for index_input_file,file_path in enumerate(inputs):

    TH1F_Container[index_input_file] = generateEffDistr(file_path)
    
    fa1 = ROOT.TF1("f11","1",0,1) ## dumb function to add a constant to TH1F
    Cumulative_Container[index_input_file] = -1*TH1F_Container[index_input_file].GetCumulative()
    Cumulative_Container[index_input_file].Add(fa1,-Cumulative_Container[index_input_file].GetMinimum())
    # Cumulative_Container[index_input_file].Scale(100./Cumulative_Container[index_input_file].GetMaximum())

    TH1F_Container[index_input_file].SetFillColorAlpha(color[index_input_file],0.3)
    TH1F_Container[index_input_file].SetLineWidth(4)
    TH1F_Container[index_input_file].SetLineColor(color[index_input_file])
    Cumulative_Container[index_input_file].SetLineColor(color[index_input_file])
    Cumulative_Container[index_input_file].SetLineWidth(4)

    TH1F_Container[index_input_file].SetName(label_list[index_input_file])
    TH1F_Container[index_input_file].SetTitle(label_list[index_input_file])
    Cumulative_Container[index_input_file].SetName(label_list[index_input_file])
    Cumulative_Container[index_input_file].SetTitle(label_list[index_input_file])

    hs.Add(TH1F_Container[index_input_file])
    hCumulative.Add(Cumulative_Container[index_input_file])
    print(f"{file_path}\t85%\t{Cumulative_Container[index_input_file].GetBinContent(86)} {TH1F_Container[index_input_file].GetEntries()}")




for h in [hs,hCumulative]:
    h.Draw("nostack")
    CMS_lumi.CMS_lumi(c2,  0,  0.3,"1.83 fb ^{-1} (13.6 TeV)")

    h.GetXaxis().SetTitle("Efficiency [%]")
    h.GetXaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitle("Number of Chambers")
    h.GetYaxis().SetTitleSize(0.04)
    if h.GetName() == "hCumulative": h.GetYaxis().SetTitle("%Ch w/ efficiency greater than X")
    leg = ROOT.TLegend(0.15,0.65,0.55,0.75)
    leg.AddEntry(h.GetHists()[0],"GE1/1 AVG Efficiency "+str(round(h.GetHists()[0].GetMean(),2))+"%","F")
    leg.Draw()
    c2.Modified()
    c2.Update()
    c2.SaveAs(f"{OUTPUT_PATH}/EfficiencyDistribution/"+h.GetTitle()+".pdf")
    print(f"Your ouput \t./Output/EfficiencyDistribution/{h.GetTitle()}.pdf")
