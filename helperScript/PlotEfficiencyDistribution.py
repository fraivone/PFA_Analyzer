import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *
from lib.THR_Utils import *

parser = argparse.ArgumentParser(
        description='''Scripts that produces efficiency distribution different runs, starting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.PlotEfficiencyDistribution --inputs MergedRuns_700uA_100HZ_Trim  MergedRuns_690uA_100HZ_Trim  --output EfficiencyDistributuion""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--output', type=str , help="Output file name",required=False)
parser.add_argument('--labels', type=str , help="Label with which the runs should be listed in the legend (according to inputs order). If not provided, input names will be used",required=False,nargs='*')
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)

chamber_with_short =[] # ["GE11-P-07L2-S", "GE11-P-09L2-S", "GE11-P-10L2-L", "GE11-P-12L1-L", "GE11-P-14L1-L", "GE11-P-15L1-S", "GE11-P-15L2-S", "GE11-P-18L1-L", "GE11-P-24L2-L", "GE11-P-34L2-L", "GE11-P-36L1-L", "GE11-P-36L2-L", "GE11-M-05L2-S", "GE11-M-05L1-S", "GE11-M-06L1-L", "GE11-M-07L1-S", "GE11-M-10L2-L", "GE11-M-17L2-S", "GE11-M-19L2-S", "GE11-M-21L1-S", "GE11-M-23L1-S", "GE11-M-29L1-S", "GE11-M-30L1-L", "GE11-M-31L2-S", "GE11-M-33L1-S", "GE11-M-35L1-S"] 

args = parser.parse_args()
inputs = ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi.csv" for i in args.inputs]
output = args.output
label_list = args.labels if args.labels is not None else args.inputs

if len(label_list) != len(inputs):
    print "Parsed inputs and labels are different in number...\nExiting .."
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
hCumulative.SetMaximum(110)



## ROOTStyle
color = [ROOT.TColor.GetColorPalette(255/(len(inputs))*(i)) for i in range(len(inputs))]
c2 = setUpCanvas("Comparison",1400,900)



## Operative loop
for index_input_file,file_path in enumerate(inputs):
    print index_input_file,file_path
    df = pd.read_csv(file_path, sep=',')
    dictionary = {}

    for index, row in df.iterrows():
        region = row['region']
        chamber = row['chamber']
        layer = row['layer']
        etaP = row['etaPartition']
        propHit = row['propHit']
        matchedRecHit = row['matchedRecHit']

        if ReChLa2chamberName(region,chamber,layer) in chamber_with_short:
            # print "Skipping ",ReChLa2chamberName(region,chamber,layer)," cause it has a short"
            continue
        
        etaPartitionID =  region*(100*chamber+10*layer+etaP)

        dictionary[etaPartitionID] = {1:{'den':propHit,'num':matchedRecHit}}

    TH1F_Container[index_input_file] = generateEfficiencyDistribution(dictionary)
    
    fa1 = ROOT.TF1("f11","1",0,1) ## dumb function to add a constant to TH1F
    Cumulative_Container[index_input_file] = -1*TH1F_Container[index_input_file].GetCumulative()
    Cumulative_Container[index_input_file].Add(fa1,-Cumulative_Container[index_input_file].GetMinimum())
    # Cumulative_Container[index_input_file].Scale(100./Cumulative_Container[index_input_file].GetMaximum())

    TH1F_Container[index_input_file].SetFillColorAlpha(color[index_input_file],.6)
    TH1F_Container[index_input_file].SetLineColor(color[index_input_file])
    Cumulative_Container[index_input_file].SetLineColor(color[index_input_file])
    Cumulative_Container[index_input_file].SetLineWidth(4)

    TH1F_Container[index_input_file].SetName(label_list[index_input_file])
    TH1F_Container[index_input_file].SetTitle(label_list[index_input_file])
    Cumulative_Container[index_input_file].SetName(label_list[index_input_file])
    Cumulative_Container[index_input_file].SetTitle(label_list[index_input_file])

    hs.Add(TH1F_Container[index_input_file])
    hCumulative.Add(Cumulative_Container[index_input_file])
    print file_path,"\t85%\t",Cumulative_Container[index_input_file].GetBinContent(86),TH1F_Container[index_input_file].GetEntries()




for h in [hs,hCumulative]:
    h.Draw("nostack")
    h.GetXaxis().SetTitle("Chamber MIP Efficiency")
    if h.GetName() == "hCumulative": h.GetYaxis().SetTitle("%Ch w/ efficiency greater than X")
    c2.BuildLegend()
    c2.Modified()
    c2.Update()
    c2.SaveAs("./Output/EfficiencyDistribution/"+h.GetTitle()+".pdf")
    print "Your ouput \t","./Output/EfficiencyDistribution/"+h.GetTitle()+".pdf"




raw_input("Press any key to terminate...")
