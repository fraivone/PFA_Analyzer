import ROOT
import argparse
import os.path
import sys
import pandas as pd
from argparse import RawTextHelpFormatter
### Let's add some more from different folder
lib_folder = os.path.expandvars('$myLIB')
sys.path.insert(1, lib_folder)
try:
    from ROOT_Utils import *
    

except:
  print("ERROR:\n\tCan't find the package CMS_lumi and tdrstlye\n\tPlease verify that this file are placed in the path $myLIB/ROOT_Utils/ \n\tAdditionally keep in mind to export the environmental variable $myLIB\nEXITING...\n") 
  sys.exit(0)
try:
    from PFA_Analyzer_Utils import *
except:
    print ("ERROR:\n\tCan't find the package PFA_Analyzer_Utils\nEXITING...\n")
    sys.exit(0)

parser = argparse.ArgumentParser(
        description='''Scripts that produces the HV Scan plot for the selected list_of_chambers. Takes as input the output csv from PFA_Analyzer.py''',
        epilog="""Typical exectuion\n\t  python HVScan.py --inputs  ./Plot/CRUZET_690_HVScan --HV 690""",
        formatter_class=RawTextHelpFormatter
        )


parser.add_argument('--inputs', type=str , help="file_path to the files to be compared",required=True,nargs='*')
parser.add_argument('--days', type=int , help="data taking day",required=True,nargs='*')
parser.add_argument('--ShowEta', default=False, action='store_true',help="Shows efficiency for each EtaPartition",required=False)

args = parser.parse_args()
inputs = args.inputs
days = args.days
ShowEta = args.ShowEta

HistogramMinX = min(days) - min(days)%10  - 0.5 
HistogramMaxX = max(days) + 10 - 0.5 #10 more days
n = int(HistogramMaxX - HistogramMinX)


if len(inputs) != len(days):
    print "Inputs and days parsed argument of different sizes...\nExiting .."
    sys.exit(0)    

map_day_file = {}
for k in range(len(days)):
    map_day_file[days[k]] = inputs[k]


### ROOT style settings
ROOT.gROOT.SetBatch(True)


ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=6001;") ##Ignore errors
ROOT.gStyle.SetLineScalePS(1) ##Nice PDF print
ROOT.gStyle.SetPalette(ROOT.kRainBow)
colorEta = [ROOT.TColor.GetColorPalette(255/8*(i)) for i in range(8)]
colorHV = [ROOT.TColor.GetColorPalette(255/8*(i)) for i in range(5)]
c2 = setUpCanvas("Comparison",1400,900)
c2.SetGrid()
c2.SetTickx()
c2.SetTicky()
c2.cd()
Plot_Container = {}
Plot_ContainerEta1 = {}
Plot_ContainerEta2 = {}
Plot_ContainerEta3 = {}
Plot_ContainerEta4 = {}
Plot_ContainerEta5 = {}
Plot_ContainerEta6 = {}
Plot_ContainerEta7 = {}
Plot_ContainerEta8 = {}

StyleHisto = []
Histonames = ["660uA - 100 Hz THR","680uA - 100 Hz THR","690uA - 100 Hz THR","690uA - 100 Hz THR + Trimming","690uA - 10 kHz THR + Trimming","700uA - 100 Hz THR","700uA - 100 Hz THR + Trimming"]
DayAtFixedHV = [[15,17,18,19,20,21],        #660 - 100 Hz THR
                [11,12,13,14],              #680 - 100 Hz THR
                [1,2,3,4,5,6,7,8,9,10,16],  #690 - 100 Hz THR
                [26,27,28,29],              #690 - 100 Hz THR + Trimming
                [51,52,53,54,55,56,57,58],  #690 - 10 kHz THR + Trimming
                [22,23,24,25],              #700 - 100 Hz THR
                [30,31,32,33]               #700 - 100 Hz THR + Trimming
                ]             
ColorForFixedHV = [ROOT.TColor.GetColorPalette(255/6*(i)) for i in range(7)]
print [(255/6*(i)) for i in range(7)]
for index in range(7):
    h = ROOT.TH1F(Histonames[index],Histonames[index],65,0.5,65.5)
    
    h.SetFillColorAlpha(ColorForFixedHV[index], 0.3)
    h.SetLineWidth(0)
    
    for bin_value in DayAtFixedHV[index]:
        h.SetBinContent(bin_value,4) ## Exceeds Eff = 1.1 so it is ok
    
    StyleHisto.append(h)
    



## chambers up to 690 uA
#list_of_chambers = ["GE11-M-08L1-L","GE11-M-10L1-L","GE11-M-25L1-S","GE11-M-07L2-S","GE11-M-10L2-L","GE11-M-25L2-S","GE11-M-34L2-L","GE11-P-10L1-L","GE11-P-24L1-L","GE11-P-07L2-S","GE11-P-10L2-L","GE11-P-25L2-S"]

## chambers up to 700 uA
# list_of_chambers = ["GE11-M-25L1-S","GE11-M-34L2-L","GE11-P-24L1-L","GE11-P-07L2-S","GE11-P-25L2-S"]

# All chambers
list_of_chambers = []
for r in [-1,1]:
    for l in [1,2]:
        for ch in range(1,37):
            list_of_chambers.append(ReChLa2chamberName(r,ch,l))

for chamber in list_of_chambers:
    Plot_Container[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta1[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta2[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta3[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta4[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta5[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta6[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta7[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}
    Plot_ContainerEta8[chamber]={"num":ROOT.TH1F('num_'+chamber, 'num_'+chamber,n,HistogramMinX,HistogramMaxX),"den":ROOT.TH1F('den_'+chamber, 'den_'+chamber,n,HistogramMinX,HistogramMaxX),'eff':ROOT.TGraphAsymmErrors()}

    print chamber
    for index,day in enumerate(days):
        try:
            file_path = map_day_file[day]+"/MatchingSummary_glb_rdphi.csv"
            df = pd.read_csv(file_path, sep=',')
        except IOError:
            print "Couldn't open file : "+file_path+"\nExiting .."
            sys.exit(0)

        den = float(df[df['chamberID']==chamber]["propHit"].sum())
        num = float(df[df['chamberID']==chamber]["matchedRecHit"].sum())
        
        denEta1 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==1) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==1) ]["propHit"]) !=0 else 0
        denEta2 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==2) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==2) ]["propHit"]) !=0 else 0
        denEta3 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==3) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==3) ]["propHit"]) !=0 else 0
        denEta4 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==4) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==4) ]["propHit"]) !=0 else 0
        denEta5 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==5) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==5) ]["propHit"]) !=0 else 0
        denEta6 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==6) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==6) ]["propHit"]) !=0 else 0
        denEta7 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==7) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==7) ]["propHit"]) !=0 else 0
        denEta8 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==8) ]["propHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==8) ]["propHit"]) !=0 else 0
        numEta1 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==1) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==1) ]["matchedRecHit"]) !=0 else 0
        numEta2 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==2) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==2) ]["matchedRecHit"]) !=0 else 0
        numEta3 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==3) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==3) ]["matchedRecHit"]) !=0 else 0
        numEta4 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==4) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==4) ]["matchedRecHit"]) !=0 else 0
        numEta5 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==5) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==5) ]["matchedRecHit"]) !=0 else 0
        numEta6 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==6) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==6) ]["matchedRecHit"]) !=0 else 0
        numEta7 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==7) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==7) ]["matchedRecHit"]) !=0 else 0
        numEta8 = float(df [ (df['chamberID']==chamber) & (df['etaPartition']==8) ]["matchedRecHit"]) if len(df [ (df['chamberID']==chamber) & (df['etaPartition']==8) ]["matchedRecHit"]) !=0 else 0


        selected_bin = day - int(HistogramMinX+0.5) + 1 

    
        Plot_Container[chamber]['den'].SetBinContent(selected_bin,den)
        Plot_Container[chamber]['num'].SetBinContent(selected_bin,num)


        Plot_ContainerEta1[chamber]['den'].SetBinContent(selected_bin,denEta1)
        Plot_ContainerEta1[chamber]['num'].SetBinContent(selected_bin,numEta1)
        Plot_ContainerEta2[chamber]['den'].SetBinContent(selected_bin,denEta2)
        Plot_ContainerEta2[chamber]['num'].SetBinContent(selected_bin,numEta2)
        Plot_ContainerEta3[chamber]['den'].SetBinContent(selected_bin,denEta3)
        Plot_ContainerEta3[chamber]['num'].SetBinContent(selected_bin,numEta3)
        Plot_ContainerEta4[chamber]['den'].SetBinContent(selected_bin,denEta4)
        Plot_ContainerEta4[chamber]['num'].SetBinContent(selected_bin,numEta4)
        Plot_ContainerEta5[chamber]['den'].SetBinContent(selected_bin,denEta5)
        Plot_ContainerEta5[chamber]['num'].SetBinContent(selected_bin,numEta5)
        Plot_ContainerEta6[chamber]['den'].SetBinContent(selected_bin,denEta6)
        Plot_ContainerEta6[chamber]['num'].SetBinContent(selected_bin,numEta6)
        Plot_ContainerEta7[chamber]['den'].SetBinContent(selected_bin,denEta7)
        Plot_ContainerEta7[chamber]['num'].SetBinContent(selected_bin,numEta7)
        Plot_ContainerEta8[chamber]['den'].SetBinContent(selected_bin,denEta8)
        Plot_ContainerEta8[chamber]['num'].SetBinContent(selected_bin,numEta8)

    
    Plot_Container[chamber]['den'].Fill(HistogramMaxX-0.0001) ## Adding fake point to preserve axis range
    Plot_Container[chamber]['den'].Fill(HistogramMinX+0.5) ## Adding fake point to preserve axis range
    Plot_Container[chamber]['eff'].Divide(Plot_Container[chamber]['num'],Plot_Container[chamber]['den'],"B")
    Plot_Container[chamber]['eff'].GetXaxis().SetRangeUser(HistogramMinX,HistogramMaxX)
    LastFakePointIndex = Plot_Container[chamber]['eff'].GetN()-1 ## Getting Fake point index
    Plot_Container[chamber]['eff'].RemovePoint(LastFakePointIndex) ## Removing Fake point
    Plot_Container[chamber]['eff'].RemovePoint(0) ## Removing Fake point
    
    
    Plot_ContainerEta1[chamber]['eff'].Divide(Plot_ContainerEta1[chamber]['num'],Plot_ContainerEta1[chamber]['den'],"B")
    Plot_ContainerEta2[chamber]['eff'].Divide(Plot_ContainerEta2[chamber]['num'],Plot_ContainerEta2[chamber]['den'],"B")
    Plot_ContainerEta3[chamber]['eff'].Divide(Plot_ContainerEta3[chamber]['num'],Plot_ContainerEta3[chamber]['den'],"B")
    Plot_ContainerEta4[chamber]['eff'].Divide(Plot_ContainerEta4[chamber]['num'],Plot_ContainerEta4[chamber]['den'],"B")
    Plot_ContainerEta5[chamber]['eff'].Divide(Plot_ContainerEta5[chamber]['num'],Plot_ContainerEta5[chamber]['den'],"B")
    Plot_ContainerEta6[chamber]['eff'].Divide(Plot_ContainerEta6[chamber]['num'],Plot_ContainerEta6[chamber]['den'],"B")
    Plot_ContainerEta7[chamber]['eff'].Divide(Plot_ContainerEta7[chamber]['num'],Plot_ContainerEta7[chamber]['den'],"B")
    Plot_ContainerEta8[chamber]['eff'].Divide(Plot_ContainerEta8[chamber]['num'],Plot_ContainerEta8[chamber]['den'],"B")
    ## Plots style
    Plot_Container[chamber]['eff'].SetMaximum(1.1)
    Plot_Container[chamber]['eff'].SetMinimum(0.)
    Plot_Container[chamber]['eff'].GetXaxis().SetTitle("CRUZET Day") #("Equivalent Divider Current (uA)")
    Plot_Container[chamber]['eff'].GetYaxis().SetTitle("MIP Efficiency")
    Plot_Container[chamber]['eff'].SetTitle(chamber)
    Plot_Container[chamber]['eff'].SetName(chamber)
    Plot_Container[chamber]['eff'].SetMarkerStyle(20)
    Plot_Container[chamber]['eff'].SetMarkerSize(.8)
    Plot_Container[chamber]['eff'].SetLineWidth(4)


    Plot_ContainerEta1[chamber]['eff'].SetTitle("Eta1")
    Plot_ContainerEta2[chamber]['eff'].SetTitle("Eta2")
    Plot_ContainerEta3[chamber]['eff'].SetTitle("Eta3")
    Plot_ContainerEta4[chamber]['eff'].SetTitle("Eta4")
    Plot_ContainerEta5[chamber]['eff'].SetTitle("Eta5")
    Plot_ContainerEta6[chamber]['eff'].SetTitle("Eta6")
    Plot_ContainerEta7[chamber]['eff'].SetTitle("Eta7")
    Plot_ContainerEta8[chamber]['eff'].SetTitle("Eta8")


    Plot_ContainerEta1[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta2[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta3[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta4[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta5[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta6[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta7[chamber]['eff'].SetMarkerStyle(4)
    Plot_ContainerEta8[chamber]['eff'].SetMarkerStyle(4)

    
    Plot_ContainerEta1[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta2[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta3[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta4[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta5[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta6[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta7[chamber]['eff'].SetLineWidth(1)
    Plot_ContainerEta8[chamber]['eff'].SetLineWidth(1)

    Plot_ContainerEta1[chamber]['eff'].SetLineColor(colorEta[1-1])
    Plot_ContainerEta2[chamber]['eff'].SetLineColor(colorEta[2-1])
    Plot_ContainerEta3[chamber]['eff'].SetLineColor(colorEta[3-1])
    Plot_ContainerEta4[chamber]['eff'].SetLineColor(colorEta[4-1])
    Plot_ContainerEta5[chamber]['eff'].SetLineColor(colorEta[5-1])
    Plot_ContainerEta6[chamber]['eff'].SetLineColor(colorEta[6-1])
    Plot_ContainerEta7[chamber]['eff'].SetLineColor(colorEta[7-1])
    Plot_ContainerEta8[chamber]['eff'].SetLineColor(colorEta[8-1])
    
    Plot_ContainerEta1[chamber]['eff'].SetMarkerColor(colorEta[1-1])
    Plot_ContainerEta2[chamber]['eff'].SetMarkerColor(colorEta[2-1])
    Plot_ContainerEta3[chamber]['eff'].SetMarkerColor(colorEta[3-1])
    Plot_ContainerEta4[chamber]['eff'].SetMarkerColor(colorEta[4-1])
    Plot_ContainerEta5[chamber]['eff'].SetMarkerColor(colorEta[5-1])
    Plot_ContainerEta6[chamber]['eff'].SetMarkerColor(colorEta[6-1])
    Plot_ContainerEta7[chamber]['eff'].SetMarkerColor(colorEta[7-1])
    Plot_ContainerEta8[chamber]['eff'].SetMarkerColor(colorEta[8-1])

    Plot_Container[chamber]['eff'].Draw("APE")
    legend = ROOT.TLegend (0.5 ,0.15 ,0.8 ,0.45)
    
    for hst in StyleHisto:
        hst.Draw("SAME")
        legend.AddEntry(hst,hst.GetTitle(),"F")
    
    if ShowEta:
        Plot_ContainerEta1[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta2[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta3[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta4[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta5[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta6[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta7[chamber]['eff'].Draw("SAMEPE")
        Plot_ContainerEta8[chamber]['eff'].Draw("SAMEPE")
        c2.BuildLegend(0.7,0.21,0.9,0.5,"","PLE").SetBorderSize(1)
    else:
        legend.Draw("SAME")
        
    Plot_Container[chamber]['eff'].Draw("SAME PE")


    c2.Modified()
    c2.Update()
    
    c2.SaveAs("./TimeSeries/"+chamber+".pdf")
    # c2.SaveAs("./AllRunsAt690uA/"+chamber+".png")
