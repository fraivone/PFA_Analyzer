import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import ROOT
import math
import numpy as np
import os.path
import sys
from ..lib.DAQ_Chamber_Mapping import *

def ReChLa2chamberName(re,ch,la):
    endcap = "M" if re == -1 else "P"
    size = "S" if ch%2 == 1 else "L"
    chID = 'GE11-'+endcap+'-%02d' % ch +"L"+str(la)+"-"+size 
    return chID


### Script to colate TRIM DAC data to DataFrame

ROOT.gROOT.SetBatch(True)
folder_TRIM = "TRIM_ARM_DAC"
folder_labels = ["Trimming"]
TRIMDAC2fC = 2.1/63  #63 TRIM DAC = 2.1 fC

color = [ROOT.TColor.GetColorPalette(255/3*(i)) for i in range(3)]

c2 = ROOT.TCanvas("TRIM DAC Comparison")
c2.Divide(2,2)
ROOT.gStyle.SetOptTitle(0) #Don't print titles
ROOT.gStyle.SetLineScalePS(1) 

master_df = pd.DataFrame(columns=["Chamber_ID","VFATN","Channel","TRIM_DAC","TRIM_fC"])


myPlot = {}
## tupl = crate, amc, OH
chamber_list = chamber_mapping.keys()
counter = 1
for re in [-1,1]:
    for layer in [1,2]:

        x_data = []
        y_data = []
        ey_data = []
        ex_data = []
        for ch in range(1,37):
            chamber = ReChLa2chamberName(re,ch,layer)
            tupl = chamber_mapping[chamber]
            crate = int(tupl[0])
            amc = int(tupl[1])
            OHLink = int(tupl[2])
            print chamber,tupl
            ## TRIM ARM DAC
            meanChannels = []
            stdChannels = []
            for VFATN in range(24):
                             
        
                file_path = folder_TRIM+"/shelf%02d" % crate + "/amc%02d" % amc +"/config_OH%d" % OHLink +"_VFAT%d" %VFATN +"_append.cfg"
                try:
                    data = pd.read_csv(file_path, sep=' ',header=None)
                except:
                    print "File "+file_path+"\t not found... skipping\n"
                    continue
        
                data[0] = data[0].apply(lambda x: (x.replace("VFAT_CHANNELS.","")).replace(".ARM_TRIM_"," ").replace("CHANNEL",""))
                data.columns = ["Channel","Value","Quantity"]
        
                data["Quantity"] = data["Channel"].apply(lambda x: x.split(" ")[1].replace(" ",""))
                data["Channel"] = data["Channel"].apply(lambda x: int(x.split(" ")[0]))
                data["Value"] = data["Value"].astype(int)
        
                ##adding polarity as singed int
                for channel in range(128):
                    sub_data = data[data["Channel"]==channel]
                    polarity = sub_data[sub_data["Quantity"]=="POLARITY"]["Value"].values[0]
                    polarity = 1 if polarity==1 else -1
                    index_of_interest = data.index[((data["Channel"]==channel) & (data["Quantity"]=="AMPLITUDE"))][0]
                    unsigned_amplitude = data["Value"].loc[index_of_interest]
                    signed_amplitude = unsigned_amplitude * polarity
                    data.loc[index_of_interest,'Value']  = signed_amplitude

                data = data[data["Quantity"]=="AMPLITUDE"]         
                data["Chamber_ID"]= [ chamber for i in range( len(data["Value"]) )]
                data["VFATN"]= [ VFATN for i in range( len(data["Value"]) )]
                data.rename(columns={'Value': 'TRIM_DAC'}, inplace=True)
                data["TRIM_fC"] = data["TRIM_DAC"]*TRIMDAC2fC
                data.drop(['Quantity'], axis=1,inplace=True)
                data.reset_index(inplace=True, drop=True) 
                data = data[["Chamber_ID","VFATN","Channel","TRIM_DAC","TRIM_fC"]]
                
                meanChannels.append(data["TRIM_fC"].mean())
                stdChannels.append(data["TRIM_fC"].std())
                master_df = master_df.append(data,ignore_index=True)

            
            x_data.append( ch )
            y_data.append( np.mean(np.asarray(meanChannels,dtype=float)) )        
            ey_data.append(0)
            ex_data.append(0)

        print x_data

        n = len(x_data)
        x = np.asarray(x_data,dtype=float)
        y = np.asarray(y_data,dtype=float)
        ex = np.asarray(ex_data,dtype=float)
        ey = np.asarray(ey_data,dtype=float)

        myPlot[counter] = ROOT.TGraphErrors(n,x,y,ex,ey)
        myPlot[counter].SetMarkerColor(color[0])
        myPlot[counter].SetLineColor(color[0])
        myPlot[counter].GetXaxis().SetTitle("Chamber")
        myPlot[counter].GetYaxis().SetTitle("AVG TRIM fC")
        myPlot[counter].SetMinimum(-2.1)
        myPlot[counter].SetMaximum(2.1)
        myPlot[counter].SetMarkerStyle(21)
        myPlot[counter].SetTitle("Trimming Data")
        myPlot[counter].SetName("Trimming Data")
        c2.cd(counter)
        
        myPlot[counter].Draw("APE")
        counter += 1

c2.cd(1).BuildLegend(0.3,0.21,0.3,0.21,"","PE").SetBorderSize(1)
c2.cd(2).BuildLegend(0.3,0.21,0.3,0.21,"","PE").SetBorderSize(1)
c2.cd(3).BuildLegend(0.3,0.21,0.3,0.21,"","PE").SetBorderSize(1)
c2.cd(4).BuildLegend(0.3,0.21,0.3,0.21,"","PE").SetBorderSize(1)


## Text
c2.cd(1)
t1 = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 ML1", "brNDC")
t1.SetBorderSize(0)
t1.SetFillColor(ROOT.gStyle.GetTitleFillColor())
t1.SetFillStyle(0)
t1.Draw()
c2.cd(2)
t2 = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 ML2", "brNDC")
t2.SetBorderSize(0)
t2.SetFillColor(ROOT.gStyle.GetTitleFillColor())
t2.SetFillStyle(0)
t2.Draw()
c2.cd(3)
t3 = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 PL1", "brNDC")
t3.SetBorderSize(0)
t3.SetFillColor(ROOT.gStyle.GetTitleFillColor())
t3.SetFillStyle(0)
t3.Draw()
c2.cd(4)
t4 = ROOT.TPaveLabel(0.3, 0.9, 0.7, 1.0, "GE11 PL2", "brNDC")
t4.SetBorderSize(0)
t4.SetFillColor(ROOT.gStyle.GetTitleFillColor())
t4.SetFillStyle(0)
t4.Draw()


c2.Modified()
c2.Update()

master_df.to_csv('TRIM_ARM_DAC/TrimmingData.csv', index=False)