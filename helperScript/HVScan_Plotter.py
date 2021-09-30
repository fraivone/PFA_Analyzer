import ROOT
import argparse
import os.path
import sys
import pandas as pd
from argparse import RawTextHelpFormatter
import subprocess
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *

parser = argparse.ArgumentParser(
        description='''Scripts that produces the HV Scan plot for the selected list_of_chambers. Takes as input the output csv from PFA_Analyzer.py''',
        epilog="""Typical exectuion\n\t  python -m helperScript.HVScan_Plotter --inputs  700uA_100HzTHR_Trimmed 690uA_100HzTHR CRUZET_680_HVScan CRUZET_660_HVScan --HV 700 690 680 660""",
        formatter_class=RawTextHelpFormatter
        )


parser.add_argument('--inputs', type=str , help="file_path to the files to be compared",required=True,nargs='*')
parser.add_argument('--HV', type=int , help="correspective equivalent divider current of the inputs",required=True,nargs='*')
parser.add_argument('--GMM', default=False,action='store_true', help="Enables special plotting for GMM. False by default",required=False)


args = parser.parse_args()
inputs = inputs = ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi.csv" for i in args.inputs]
HV = args.HV

bin_width = 1/100.
HistogramMinHV = min(HV) - 10 + bin_width/2
HistogramMaxHV = max(HV) +  10  + bin_width/2
n = int(1/bin_width)*int(HistogramMaxHV - HistogramMinHV)


if len(inputs) != len(HV):
    print "Inputs and HV parsed argument of different sizes...\nExiting .."
    sys.exit(0)    

map_HV_file = {}
for k in range(len(HV)):
    map_HV_file[HV[k]] = inputs[k]


### ROOT style settings
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=6001;") 
ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0)
colorEta = [ROOT.TColor.GetColorPalette(254/8*(i)) for i in range(9)]

c2 = setUpCanvas("Comparison",1400,900)
c2.SetTickx()
c2.SetTicky()
c2.SetGrid()

## GMM Plot
if args.GMM:
    c3 = setUpCanvas("AAA",1400,900)
    c3.SetTickx()
    c3.SetTicky()
    # c3.SetGrid()
    GMM_Plot_High = ["GE11-M-26L2-L","GE11-P-07L2-S","GE11-M-23L1-S","GE11-M-32L2-L","GE11-M-30L2-L"]
    GMM_Plot_Low  = ["GE11-P-32L1-L","GE11-P-18L1-L","GE11-M-05L1-S"]
    GMM_MultiGraph = ROOT.TMultiGraph()
    GMM_MultiGraph.GetXaxis().SetLabelSize(0.03)
    GMM_MultiGraph.GetYaxis().SetLabelSize(0.03)
    GMMindex = 0
    totalNumber = len(GMM_Plot_High) + len(GMM_Plot_Low)
    colorGMM = [ROOT.TColor.GetColorPalette(int(i*254./totalNumber)) for i in range(totalNumber+1)]

Plot_Container = {}
for label in [ "Eta"+str(i) for i in range(1,9)] + ["All"]:
    Plot_Container[label] = {}

subprocess.call(["mkdir", "-p", "./Output/HVScan"])
OutF = ROOT.TFile("./Output/HVScan/Summary.root","RECREATE")

list_of_chambers = []
for r in [1,-1]:
    for l in [1,2]:
        for ch in range(1,37):
            list_of_chambers.append(ReChLa2chamberName(r,ch,l))

#list_of_chambers = ["GE11-M-31L2-S"]

## Loop on chambers
for chamber in list_of_chambers:

    ## ROOT Object Definition loop
    for label in [ "Eta"+str(i) for i in range(1,9)] + ["All"]:
        Plot_Container[label] ={"num":ROOT.TH1F('num_'+chamber+"_"+label, 'num_'+chamber+"_"+label,n,HistogramMinHV,HistogramMaxHV),"den":ROOT.TH1F('den_'+chamber+"_"+label, 'den_'+chamber+"_"+label,n,HistogramMinHV,HistogramMaxHV),'eff':ROOT.TGraphAsymmErrors()}
    print "Processing: "+chamber
    
    ## Loop on files
    for HV_Point in HV:
        selected_bin = int(1/bin_width)*(HV_Point - int(HistogramMinHV))

        try:
            file_path = map_HV_file[HV_Point]
            df = pd.read_csv(file_path, sep=',')
        except IOError:
            print "Couldn't open file : "+file_path+"\nExiting .."
            sys.exit(0)
        
        ## Overall den and num
        den = float(df[df['chamberID']==chamber]["propHit"].sum())
        num = float(df[df['chamberID']==chamber]["matchedRecHit"].sum())
        Plot_Container["All"]['den'].SetBinContent(selected_bin,den)
        Plot_Container["All"]['num'].SetBinContent(selected_bin,num)

        if den!=0: print "\t",HV_Point,"uA","\t",num,"/",den," = ",round(float(num)/den,3)
        
        ## Eta by Eta den and num
        for etaP in range(1,9):
            den = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==etaP) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==etaP) ]["propHit"]) !=0 else 0
            num = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==etaP) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==etaP) ]["matchedRecHit"]) !=0 else 0
            
            Plot_Container["Eta"+str(etaP)]['den'].SetBinContent(selected_bin,den)
            Plot_Container["Eta"+str(etaP)]['num'].SetBinContent(selected_bin,num)
        
    ## Loop on files ends    

    ## Divide and style
    for index,label in enumerate(["Eta"+str(i) for i in range(1,9)] + ["All"]):
        Plot_Container[label]['eff'].SetLineColor(colorEta[index])
        Plot_Container[label]['eff'].SetMarkerColor(colorEta[index])   
        Plot_Container[label]['eff'].Divide(Plot_Container[label]['num'],Plot_Container[label]['den'],"B")
        Plot_Container[label]['eff'].SetMaximum(1.1)
        Plot_Container[label]['eff'].SetMinimum(0.)
        Plot_Container[label]['eff'].GetXaxis().SetTitle("Equivalent Divider Current (#muA)")
        Plot_Container[label]['eff'].GetXaxis().SetLabelSize(0.03)
        Plot_Container[label]['eff'].GetYaxis().SetLabelSize(0.03)
        Plot_Container[label]['eff'].GetYaxis().SetTitle("MIP Efficiency")
        Plot_Container[label]['eff'].SetTitle(chamber+"_"+label)
        Plot_Container[label]['eff'].SetName(chamber+"_"+label)
        Plot_Container[label]['eff'].GetXaxis().SetLimits(HistogramMinHV,HistogramMaxHV)
        Plot_Container[label]['eff'].SetMarkerStyle(20)
        Plot_Container[label]['eff'].SetLineWidth(1)
        
    c2.cd()    
    ## Different style for the main plot
    Plot_Container["All"]['eff'].SetMarkerSize(.8)
    Plot_Container["All"]['eff'].SetMarkerColor(ROOT.kBlack)
    Plot_Container["All"]['eff'].SetLineColor(ROOT.kBlack)
    Plot_Container["All"]['eff'].SetMarkerStyle(20)
    Plot_Container["All"]['eff'].SetLineWidth(4)
    Plot_Container["All"]['eff'].Draw("APEL")
    
    
    ## Drawing the plots
    for label in ["Eta"+str(i) for i in range(1,9)]:
        Plot_Container[label]['eff'].Draw("SAMEPEL")
    c2.BuildLegend(0.8,0.2,1.,.6,"","PLE").SetBorderSize(1)

    ## Some more text
    CMS_Prel_Pad=ROOT.TPad("tempPad","tempPad",0.,0.,1.,1.)
    CMS_Prel_Pad.SetFillStyle(4000)
    CMS_Prel_Pad.Draw()
    CMS_Prel_Pad.cd()
    s = "CMS Preliminary"
    CMS_Prel = ROOT.TLatex(0.14,0.93,s)
    CMS_Prel.SetTextSize(0.035)
    CMS_Prel.Draw()
    CMS_Prel_Pad.Modified()
    CMS_Prel_Pad.Update()
    
    tempPad=ROOT.TPad("tempPad","tempPad",0.15, 0.8,0.43,0.89)
    tempPad.Draw()
    tempPad.cd()
    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)
    latex.SetTextSize(0.25)
    latex.DrawTextNDC(0., 0.8, "CRUZET: GE1/1 HV Scan")
    latex.DrawTextNDC(0., 0.3, "Gas Mixture Ar/CO2 (70/30%)")
    latex.Draw()
    tempPad.Modified()
    tempPad.Update()

    
    c2.Modified()
    c2.Update()


    for label in ["Eta"+str(i) for i in range(1,9)] + ["All"]:
        writeToTFile(OutF,Plot_Container[label]['eff'],chamber)
    writeToTFile(OutF,c2,chamber)
    c2.SaveAs("./Output/HVScan/"+chamber+".pdf")


    ## GMM Plot
    if args.GMM and (chamber in GMM_Plot_High or chamber in GMM_Plot_Low):
        Plot_Container["All"]['eff'].SetLineWidth(2)
        Plot_Container["All"]['eff'].SetLineColor(colorGMM[GMMindex])
        Plot_Container["All"]['eff'].SetMarkerColor(colorGMM[GMMindex])
        Plot_Container["All"]['eff'].SetTitle(chamber)
        Plot_Container["All"]['eff'].SetName(chamber)

        GMMindex +=1
        if chamber in GMM_Plot_Low:
            Plot_Container["All"]['eff'].SetLineStyle(2)
        GMM_MultiGraph.Add(Plot_Container["All"]['eff'])
## Loop on chambers ends

## GMM Plot
if args.GMM:
    c3.cd()
    GMM_MultiGraph.SetMinimum(0.)
    GMM_MultiGraph.SetTitle("GMM Plot")
    GMM_MultiGraph.SetName("GMM Plot")
    GMM_MultiGraph.SetMaximum(1.1)
    GMM_MultiGraph.GetXaxis().SetTitle("Equivalent Divider Current (#muA)")
    GMM_MultiGraph.GetXaxis().SetLimits(HistogramMinHV,HistogramMaxHV)
    GMM_MultiGraph.GetYaxis().SetTitle("MIP Efficiency")
    GMM_MultiGraph.Draw("0ALP")
    c3.BuildLegend(0.7,0.15,0.9,0.5,"","0PLE").SetBorderSize(0)

    CMS_Prel_Pad.Draw()
    CMS_Prel_Pad.Modified()
    CMS_Prel_Pad.Update()
    tempPad.Draw()
    tempPad.Modified()
    tempPad.Update()
    c3.Modified()
    c3.Update()
    c3.SaveAs("./Output/HVScan/GMM.pdf")
    OutF.cd()
    writeToTFile(OutF,GMM_MultiGraph,"GMM")
    writeToTFile(OutF,c3,"GMM")
    OutF.Close()

print "\nYour files are stored under ./Output/HVScan/\nBye now.."
