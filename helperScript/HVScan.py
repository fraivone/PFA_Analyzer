import ROOT
import argparse
import sys
from argparse import RawTextHelpFormatter
import pandas as pd
from ROOT_Utils import *
from PFA_Analyzer_Utils import *
import array
import numpy as np
import TDR_Approval_Style

ver=""

if ver=="":
    hv_scan_output = "/eos/user/f/fivone/www/P5_Operations/Run3/RunMerge/HVScan/"
    fa1 = ROOT.TF1("erf(x)","[2] + [2]*ROOT::Math::erf((x-[0])/(sqrt(2)*[1]))",0,720)
    
if ver == "rev1":
    hv_scan_output = "/eos/user/f/fivone/www/P5_Operations/Run3/RunMerge/HVScan_rev1/"
    fa1 = ROOT.TF1("erf(x)","2*[2]/(1 + exp( -(x-[0])/[1]  ) )",0,720)

fa1.SetParNames("Mean","Sigma","Norm")
fa1.SetParameters(650,30,50)
fa1.SetLineColor(ROOT.kBlue)
fa1.SetLineStyle(1)
fa1.SetParLimits(2,0,50)
fa1.SetParLimits(1,5,50000)
fa1.SetParLimits(0,500,2000)


parser = argparse.ArgumentParser(
        description='''Scripts that prouces HV Scan plotsstarting from csv files output by PFA_Analyzer.py''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.HVScan --inputs 360019_700 360019_600 --HV 700 600 """,
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged",required=True,nargs='*')
parser.add_argument('--HV', type=int , help="HV values",required=True,nargs="*")
parser.add_argument('--batch', default=False, action='store_true',help="ROOT in batch mode",required=False)


args = parser.parse_args()
inputs = [f"{OUTPUT_PATH}//PFA_Analyzer_Output/CSV/"+i+"/" for i in args.inputs]
hv_list = [ k for k in args.HV ]

if len(hv_list) != len(inputs):
    print("Parsed inputs and labels are different in number...\nExiting ..")
    print(f"{hv_list}")
    print(f"{inputs}")
    sys.exit(0)    


QC8_THR = []
PRIMARY_XAXIS_MIN = 580
PRIMARY_XAXIS_MAX = 720
PRIMARY_YAXIS_MIN = 0
PRIMARY_YAXIS_MAX = 105


### ROOT style settings
ROOT.gROOT.SetBatch(args.batch)

ROOT.gStyle.SetLineScalePS(1)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetOptTitle(0) #Don't print titles

## ROOTStyle
TDR_Approval_Style.setTDRStyle()
color = [ROOT.TColor.GetColorPalette( int(255/(len(inputs))*(i))) for i in range(len(inputs))]

H_ref = 800
W_ref = 800
W = W_ref
H = H_ref

T = 0.08*H
B = 0.16*H
L = 0.16*W
R = 0.08*W

c2 = ROOT.TCanvas("c1", "c1", W, H)
c2.SetFillColor(0)
c2.SetBorderMode(0)
c2.SetFrameFillStyle(0)
c2.SetFrameBorderMode(0)
c2.SetLeftMargin( L/W )
c2.SetRightMargin( R/W )
c2.SetTopMargin( T/H )
c2.SetBottomMargin( B/H )
c2.SetTickx(0)
c2.SetTicky(0)
c2.SetGrid()

## ROOT Objects
TH1_OnsetDistr = ROOT.TH1F("mean Erd","mean Erf",100,550,850)
TH1_chi2Distr = ROOT.TH1F("chi2 Fit","chi2 Fit",40,0,0)
TH1_sigmaDistr = ROOT.TH1F("sigma Erf","sigma Erf",40,0,0)
TH1_normDistr = ROOT.TH1F("norm Erf","norm Erf",40,0,0)

TGraphError_Dict = {}
Text_Dict = {}
TGraph_THR_Dict = {}
fit_results = {}

for re,la in [(-1,1),(1,1),(-1,2),(1,2)]: 
    label = EndcapLayer2label(re,la)
    TGraphError_Dict[label] = {}

## Operative loop

for re in [-1,1]:
    for chamber in range(1,37):
        for la in [1,2]:
            label = EndcapLayer2label(re,la)           
            x,y,exl,exh,eyl,eyh = [],[],[],[],[],[]
            chamberID = getChamberName(re,chamber,la)
            muons_used,muons_matched = [],[]
            for index,file_path in enumerate(inputs):
                HV = hv_list[index]
                matched,propagated,efficiency_CL68 = ChamberEfficiencyFromCSV(file_path,chamberID)
                muons_used.append(propagated)
                muons_matched.append(matched)
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
            TGraphError_Dict[label][chamberID].Fit(fa1,"RMBQ")
            chi2 = fa1.GetChisquare()
            mean = fa1.GetParameter(0)
            sigma = fa1.GetParameter(1)
            norm = fa1.GetParameter(2)
            TH1_OnsetDistr.Fill(mean)
            TH1_chi2Distr.Fill(chi2)
            TH1_sigmaDistr.Fill(sigma)
            TH1_normDistr.Fill(norm*2)
            fit_results[chamberID] = {
                                    "mean":mean,
                                    "sigma":sigma,
                                    "norm":norm,
                                    "chi2":chi2}
            

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
            TGraphError_Dict[label][chamberID].GetXaxis().SetTitleSize(0.04)
            TGraphError_Dict[label][chamberID].SetMaximum(PRIMARY_YAXIS_MAX)
            TGraphError_Dict[label][chamberID].SetMinimum(PRIMARY_YAXIS_MIN)
            TGraphError_Dict[label][chamberID].GetYaxis().SetLabelSize(0.03)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTickLength(0.015)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTitleOffset(1.2)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTitleSize(0.04)
            TGraphError_Dict[label][chamberID].GetYaxis().SetTitle("Efficiency (%)")

            # c2.cd()
            TGraphError_Dict[label][chamberID].Draw("APE")

            ## Add labels on point
            
            for x_val,y_val in zip(x,y):
                norm_ypos =  ((7+y_val)//5)*5 if y_val < 93 else 100
                prop = muons_used[ hv_list.index(x_val) ]
                matc = muons_matched[ hv_list.index(x_val) ]
                pointAnnotation = ROOT.TLatex(x_val,norm_ypos,f"{int(prop)}")
                pointAnnotation.SetTextSize(0.02)
                TGraphError_Dict[label][chamberID].GetListOfFunctions().Add(pointAnnotation) 
            TGraphError_Dict[label][chamberID].Draw("APE")           

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextAngle(0)
            latex.SetTextColor(ROOT.kBlack)
            # Left aligned
            latex.SetTextAlign(12)
            latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())
            latex.SetTextFont(61)
            latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.65*ROOT.gPad.GetTopMargin(), "CMS")
            latex.SetTextFont(52)
            latex.SetTextSize(0.45*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.6*ROOT.gPad.GetLeftMargin(),1 - 0.74*ROOT.gPad.GetTopMargin(), "Preliminary")

            latex.SetTextFont(42)
            latex.SetTextSize(0.3*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 1.2*c2.GetTopMargin(), "#chi^{2}: "+f"{chi2:1.3f}")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 1.6*c2.GetTopMargin(), f"Mean: {mean:3.2f} uA")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.4*c2.GetTopMargin(), f"Norm: {2*norm:2.2f}")
            if ver=="":
                latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.8*c2.GetTopMargin(), f"Slope: {np.sqrt(2/np.pi) * norm /(sigma):2.2f}")
                latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2*c2.GetTopMargin(), f"Sigma: {sigma:2.2f}")
            if ver=="rev1":
                latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.8*c2.GetTopMargin(), f"Slope: {norm /(2*sigma):2.2f}")
                latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2*c2.GetTopMargin(), f"Lambda: {1/sigma:2.4f}")
            
        
            ## Right aligned
            latex.SetTextAlign(32)
            
            latex.SetTextSize(0.44*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1-1.1*c2.GetRightMargin(), 1 - 0.7*ROOT.gPad.GetTopMargin(), f"{chamberID}")

            # c2.BuildLegend()
            c2.Modified()
            c2.Update()

            ## Output files
            folder_name = hv_scan_output + label
            CreatEOSFolder(folder_name)
            outputPath = folder_name+"/"+chamberID +".pdf"
            c2.SaveAs(outputPath)
            Convert2png(outputPath)

        ## PLOT SUPERCHAMBER
        ch1,l1 = getChamberName(re,chamber,1),EndcapLayer2label(re,1)
        ch2,l2 = getChamberName(re,chamber,2),EndcapLayer2label(re,2)
        if TGraphError_Dict[l1].get(ch1) == None or TGraphError_Dict[l2].get(ch2) == None:
            continue
        else:
            TGraphError_Dict[l1][ch1].SetMarkerColor(ROOT.kBlue)
            TGraphError_Dict[l1][ch1].SetLineColor(ROOT.kBlue)
            TGraphError_Dict[l1][ch1].GetListOfFunctions().FindObject("erf(x)").SetLineColor(ROOT.kBlue)

            TGraphError_Dict[l2][ch2].SetLineColor(ROOT.kRed)
            TGraphError_Dict[l2][ch2].SetMarkerColor(ROOT.kRed)
            TGraphError_Dict[l2][ch2].GetListOfFunctions().FindObject("erf(x)").SetLineColor(ROOT.kRed)

            TGraphError_Dict[l1][ch1].Draw("APE")
            TGraphError_Dict[l2][ch2].Draw("SAME PE")

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextAngle(0)
            latex.SetTextColor(ROOT.kBlack)

            latex.SetTextAlign(12)
            latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())
            latex.SetTextFont(61)
            latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.65*ROOT.gPad.GetTopMargin(), "CMS")
            latex.SetTextFont(52)
            latex.SetTextSize(0.45*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.6*ROOT.gPad.GetLeftMargin(),1 - 0.74*ROOT.gPad.GetTopMargin(), "Preliminary")
            latex.SetTextFont(42)

            latex.SetTextSize(0.35*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 1.2*c2.GetTopMargin(), "#color[4]{"+f"{ch1}"+"}")
            latex.SetTextSize(0.3*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 1.6*c2.GetTopMargin(), "#chi^{2}: "+f"{fit_results[ch1]['chi2']:1.3f}")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2*c2.GetTopMargin(), f"Mean: {fit_results[ch1]['mean']:3.2f} uA")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.4*c2.GetTopMargin(), f"Std Dev: {fit_results[ch1]['sigma']:2.2f}")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.8*c2.GetTopMargin(), f"Norm: {2*fit_results[ch1]['norm']:2.2f}")

            latex.SetTextSize(0.35*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 3.3*c2.GetTopMargin(), "#color[2]{"+f"{ch2}"+"}")
            latex.SetTextSize(0.3*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 3.7*c2.GetTopMargin(), "#chi^{2}: "+f"{fit_results[ch2]['chi2']:1.3f}")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.1*c2.GetTopMargin(), f"Mean: {fit_results[ch2]['mean']:3.2f} uA")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.5*c2.GetTopMargin(), f"Std Dev: {fit_results[ch2]['sigma']:2.2f}")
            latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.9*c2.GetTopMargin(), f"Norm: {2*fit_results[ch2]['norm']:2.2f}")


            latex.SetTextAlign(32)
            latex.SetTextSize(0.45*ROOT.gPad.GetTopMargin())
            latex.DrawLatex(1-1.1*ROOT.gPad.GetRightMargin(),1 - 0.74*ROOT.gPad.GetTopMargin(), "#sqrt{s} = 13.6 TeV (2022)")


            legend = ROOT.TLegend (0.6 ,0.2 ,0.9 ,0.4)
            legend.AddEntry ( TGraphError_Dict[l1][ch1] , ch1, "pl" )
            legend.AddEntry ( TGraphError_Dict[l2][ch2] , ch2, "pl" )
            legend.SetLineWidth (0)
            legend.Draw ("same") 
            c2.Modified()
            c2.Update()


            file_name = f"{hv_scan_output}/GE{re}1_{chamber:02d}.pdf"
            c2.SaveAs(file_name)
            Convert2png(file_name)
        ## PLOT SUPERCHAMBER
            


file_name = hv_scan_output+"HVOnsetDistribution.pdf"
TH1_OnsetDistr.Draw("HIST")
c2.SaveAs(file_name)
Convert2png(file_name)

file_name = hv_scan_output+"chi2Distribution.pdf"
TH1_chi2Distr.Draw("HIST")
c2.SaveAs(file_name)
Convert2png(file_name)

file_name = hv_scan_output+"sigmaDistribution.pdf"
TH1_sigmaDistr.Draw("HIST")
c2.SaveAs(file_name)
Convert2png(file_name)

file_name = hv_scan_output+"normDistribution.pdf"
TH1_normDistr.Draw("HIST")
c2.SaveAs(file_name)
Convert2png(file_name)


#### FOR PIET 
# ch1,l1 = "GE11-M-02L1-L","ML1"
# ch2,l2 = "GE11-M-02L2-L","ML2"
# fit_results = {ch1:{},ch2:{}}

# TGraphError_Dict[l1][ch1].Fit(fa1,"RMBQ")
# fit_results[ch1]["chi2"] = fa1.GetChisquare()
# fit_results[ch1]["mean"] = fa1.GetParameter(0)
# fit_results[ch1]["sigma"] = fa1.GetParameter(1)
# fit_results[ch1]["norm"] = fa1.GetParameter(2)

# fa1.SetLineColor(ROOT.kRed)
# TGraphError_Dict[l2][ch2].Fit(fa1,"RMBQ")
# fit_results[ch2]["chi2"] = fa1.GetChisquare()
# fit_results[ch2]["mean"] = fa1.GetParameter(0)
# fit_results[ch2]["sigma"] = fa1.GetParameter(1)
# fit_results[ch2]["norm"] = fa1.GetParameter(2)


# TGraphError_Dict[l1][ch1].SetMarkerColor(ROOT.kBlue)
# TGraphError_Dict[l2][ch2].SetMarkerColor(ROOT.kRed)
# TGraphError_Dict[l1][ch1].SetLineColor(ROOT.kBlue)
# TGraphError_Dict[l2][ch2].SetLineColor(ROOT.kRed)
# TGraphError_Dict[l1][ch1].Draw("APE")
# TGraphError_Dict[l2][ch2].Draw("SAME PE")
            

# latex = ROOT.TLatex()
# latex.SetNDC()
# latex.SetTextAngle(0)
# latex.SetTextColor(ROOT.kBlack)

# latex.SetTextAlign(12)
# latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())
# latex.SetTextFont(61)
# latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.65*ROOT.gPad.GetTopMargin(), "CMS")
# latex.SetTextFont(52)
# latex.SetTextSize(0.45*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1.6*ROOT.gPad.GetLeftMargin(),1 - 0.74*ROOT.gPad.GetTopMargin(), "Preliminary")
# latex.SetTextFont(42)

# latex.SetTextSize(0.35*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 1.2*c2.GetTopMargin(), "#color[4]{GE11-M-02L1-L}")
# latex.SetTextSize(0.3*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 1.6*c2.GetTopMargin(), "#chi^{2}: "+f"{fit_results[ch1]['chi2']:1.3f}")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2*c2.GetTopMargin(), f"Mean: {fit_results[ch1]['mean']:3.2f} uA")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.4*c2.GetTopMargin(), f"Std Dev: {fit_results[ch1]['sigma']:2.2f}")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.8*c2.GetTopMargin(), f"Norm: {2*fit_results[ch1]['norm']:2.2f}")

# latex.SetTextSize(0.35*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 3.3*c2.GetTopMargin(), "#color[2]{GE11-M-02L2-L}")
# latex.SetTextSize(0.3*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(),1 - 3.7*c2.GetTopMargin(), "#chi^{2}: "+f"{fit_results[ch2]['chi2']:1.3f}")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.1*c2.GetTopMargin(), f"Mean: {fit_results[ch2]['mean']:3.2f} uA")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.5*c2.GetTopMargin(), f"Std Dev: {fit_results[ch2]['sigma']:2.2f}")
# latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4.9*c2.GetTopMargin(), f"Norm: {2*fit_results[ch2]['norm']:2.2f}")


# latex.SetTextAlign(32)
# latex.SetTextSize(0.45*ROOT.gPad.GetTopMargin())
# latex.DrawLatex(1-1.1*ROOT.gPad.GetRightMargin(),1 - 0.74*ROOT.gPad.GetTopMargin(), "#sqrt{s} = 13.6 TeV (2022)")
            

# legend = ROOT.TLegend (0.6 ,0.2 ,0.9 ,0.4)
# legend.AddEntry ( TGraphError_Dict[l1][ch1] , ch1, "pl" )
# legend.AddEntry ( TGraphError_Dict[l2][ch2] , ch2, "pl" )
# legend.SetLineWidth (0)
# legend.Draw ("same") 
# c2.Modified()
# c2.Update()


# file_name = hv_scan_output + "ForPiet/EffCurve.pdf"
# c2.SaveAs(file_name)
# Convert2png(file_name)
