import ROOT
import csv
import os.path
from os import mkdir
import subprocess
import numpy as np
import math
import sys
import time
import argparse
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
        description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Plot Sanity Checks\n\t-Plot Residuals (takes the cut as parameter)\n\t-Plot efficiency\nCurrently allows the track matching on glb_phi and glb_rdphi''',
        epilog="""Typical exectuion\n\t python Analysis_lxplus.py  --phi_cut 0.001 --rdphi_cut 0.15""",
        formatter_class=RawTextHelpFormatter
)

parser.add_argument('-pc','--phi_cut', type=float,help="Maximum allowed dphi between RecoHit and PropHit to be counted as matched hit",required=False)
parser.add_argument('-rdpc','--rdphi_cut', type=float,help="Maximum allowed rdphi between RecoHit and PropHit to be counted as matched hit",required=False)
parser.add_argument('--chi2cut', type=float,help="Maximum normalized chi2 for which accept propagated tracks",required=False)
parser.add_argument('--chamberOFF', type=str , help="file_path to the file containing a list of the chambers you want to exclude for the run (i.e. GE11-M-29L2)",required=False)
parser.add_argument('--outputname', type=str, help="output file name",required=False)
parser.add_argument('--fiducialR','-fR', type=float , help="fiducial cut along R axis",required=False)
parser.add_argument('--fiducialPhi','-fP', type=float , help="fiducial cut along phi axis",required=False)
parser.add_argument('--DLE', default=False, action='store_true',help="Swtiches on the Double Layer Efficiency (DLE) analisys. False by default",required=False)

parser.add_argument('--minME1', type=int, help="Min number of ME1 hits",required=False)
parser.add_argument('--minME2', type=int, help="Min number of ME2 hits",required=False)
parser.add_argument('--minME3', type=int, help="Min number of ME3 hits",required=False)
parser.add_argument('--minME4', type=int, help="Min number of ME4 hits",required=False)
parser.set_defaults(phi_cut=0.001)
parser.set_defaults(rdphi_cut=0.15)
parser.set_defaults(chi2cut=9999999)
parser.set_defaults(minME1=0)
parser.set_defaults(minME2=0)
parser.set_defaults(minME3=0)
parser.set_defaults(minME4=0)
parser.set_defaults(fiducialR=1)
parser.set_defaults(fiducialPhi=0.005)
args = parser.parse_args()



ROOT.gROOT.SetBatch(True)


start_time = time.time()

files = files_in_folder("/eos/user/f/fivone/GEMNTuples/MC/Output/MC_Cosmics2021/wMEXHits/")
files = [f for f in files if ".root" in f]

matching_variables = ['glb_phi','glb_rdphi']
matching_variable_units = {'glb_phi':'rad','glb_rdphi':'cm'}
ResidualCutOff= {'glb_phi':args.phi_cut,'glb_rdphi':args.rdphi_cut}

DLE = args.DLE
fiducialCut = True
maxErrOnPropR = 100
maxErrOnPropPhi = 100
fiducialR = args.fiducialR
fiducialPhi = args.fiducialPhi
CutminPt = 0.
maxSTA_NormChi2 = args.chi2cut
minME1Hit = args.minME1
minME2Hit = args.minME2
minME3Hit = args.minME3
minME4Hit = args.minME4

noisyEtaPID = []
chamberNumberOFF = [] if args.chamberOFF is None else importOFFChamber(args.chamberOFF)

chamberOFFCanvas = setUpCanvas("ExcludedChambers",1200,1200)
chamberOFFCanvas.Divide(1,4)

for counter,plot in enumerate(ChambersOFFHisto(chamberNumberOFF)):
    chamberOFFCanvas.cd(counter+1)
    plot.Draw("")

chamberOFFCanvas.Modified()
chamberOFFCanvas.Update()


TH1nbins = 120
TH2nbins = 200
TH2min = -80

EfficiencyDictGlobal = dict((m,{}) for m in matching_variables)
EfficiencyDictDirection = dict((m,{}) for m in matching_variables)
EfficiencyDictLayer = dict((m,{}) for m in matching_variables)

## ROOT Object declaration
test = ROOT.TH2F("DirX_CLS","DirX_CLS",40,-1,1,128,0,128)
test.GetXaxis().SetTitle("dirX")
test.GetYaxis().SetTitle("CLS")
DLE_ErrPhi = ROOT.TH1F("DLE_ErrPhi","DLE_ErrPhi",100,0,0.0025)
DLE_ErrR = ROOT.TH1F("DLE_ErrR","DLE_ErrR",100,0,5)

TH1Fresidual_collector = {}
TH2Fresidual_collector = generate2DResidualContainer(matching_variables,TH2nbins,TH2min)  
THSanityChecks = {'Occupancy':{}, 
                  'NHits':{},
                  'PropagationError':{},
                  'etaP_vs_pt':[],
                  'Residual_Correlation':{},
                  'PropHit_DirLoc_xOnGE11':{'BeforeMatching':{'Long':{},'Short':{}},
                                            'AfterMatching':{'Long':{},'Short':{}}
                                            },
                  'RecHitperStrip':{}
                }

TH1MetaData = { 'isFiducialCut':[],
                'PropRErr':[],
                'PropPhiErr':[],
                'fiducialR':[],
                'fiucialPhi':[],
                'glb_phi':[],
                'glb_rdphi':[],
                'pt':[],
                'maxSTA_NormChi2':[],
                'minME1Hit':[],
                'minME2Hit':[],
                'minME3Hit':[],
                'minME4Hit':[]}

## Initialize Collectors

for key in TH1MetaData.keys():
    TH1MetaData[key] = ROOT.TH1F("CUT_on_"+key,"CUT_on_"+key,10,1,1)

TH1MetaData['isFiducialCut'].Fill(fiducialCut)
TH1MetaData['PropRErr'].Fill(maxErrOnPropR)
TH1MetaData['PropPhiErr'].Fill(maxErrOnPropPhi)
TH1MetaData['fiducialR'].Fill(fiducialR)
TH1MetaData['fiucialPhi'].Fill(fiducialPhi)
TH1MetaData['glb_phi'].Fill(ResidualCutOff['glb_phi'])
TH1MetaData['glb_rdphi'].Fill(ResidualCutOff['glb_rdphi'])
TH1MetaData['pt'].Fill(CutminPt)
TH1MetaData['maxSTA_NormChi2'].Fill(maxSTA_NormChi2)
TH1MetaData['minME1Hit'].Fill(minME1Hit)
TH1MetaData['minME2Hit'].Fill(minME2Hit)
TH1MetaData['minME3Hit'].Fill(minME3Hit)
TH1MetaData['minME4Hit'].Fill(minME4Hit)
TH1MetaData['chamberOFFCanvas']=chamberOFFCanvas




THSanityChecks['NHits']['BeforeMatching'] = {'PerEVT':{'Reco':ROOT.TH1F("NRecoHitsPerEVT","NRecoHitsPerEVT",200,0,200),'Prop':ROOT.TH1F("NPropHitsPerEVT","NPropHitsPerEVT",20,0,20)},
                                             'ML1':ROOT.TH1F("ML1_NRecoHitsPerEVT","ML1_NRecoHitsPerEVT",200,0,200),
                                             'ML2':ROOT.TH1F("ML2_NRecoHitsPerEVT","ML1_NRecoHitsPerEVT",200,0,200),
                                             'PL1':ROOT.TH1F("PL1_NRecoHitsPerEVT","PL1_NRecoHitsPerEVT",200,0,200),
                                             'PL2':ROOT.TH1F("PL2_NRecoHitsPerEVT","PL2_NRecoHitsPerEVT",200,0,200)
                                         }

THSanityChecks['NHits']['AfterMatching'] = {'All':ROOT.TH1F("N_MatchedRecoHitsPerEVT","N_MatchedRecoHitsPerEVT",20,0,20),
                                             'ML1':ROOT.TH1F("ML1_N_MatchedRecoHitsPerEVT","ML1_N_MatchedRecoHitsPerEVT",20,0,20),
                                             'ML2':ROOT.TH1F("ML2_N_MatchedRecoHitsPerEVT","ML2_N_MatchedRecoHitsPerEVT",20,0,20),
                                             'PL1':ROOT.TH1F("PL1_N_MatchedRecoHitsPerEVT","PL1_N_MatchedRecoHitsPerEVT",20,0,20),
                                             'PL2':ROOT.TH1F("PL2_N_MatchedRecoHitsPerEVT","PL2_N_MatchedRecoHitsPerEVT",20,0,20)
                                         }

THSanityChecks['NHits']['PerEVT_PerEtaPartitionID'] = {'Reco':ROOT.TH1F("NRecoHitsPerEVTPetEtaPartitionID","NRecoHitsPerEVTPetEtaPartitionID",10,0,10),
                                                       'Prop':ROOT.TH1F("NPropHitsPerEVTPetEtaPartitionID","NPropHitsPerEVTPetEtaPartitionID",10,0,10)}

THSanityChecks['Occupancy'].setdefault('BeforeMatching',{'Reco':ROOT.TH2F("RecoHitOccupancyBeforeMatching","RecoHitOccupancyBeforeMatching",200,-300,300,200,-300,300),
                                                         'Prop':ROOT.TH2F("PropHitOccupancyBeforeMatching","PropHitOccupancyBeforeMatching",200,-300,300,200,-300,300),
                                                         'PropLocalLong':ROOT.TH2F("PropLocLongHitBeforeMatching","PropLocLongHitBeforeMatching",TH2nbins,TH2min,-TH2min,TH2nbins,TH2min,-TH2min),
                                                         'PropLocalShort':ROOT.TH2F("PropLocShortHitBeforeMatching","PropLocShortHitBeforeMatching",TH2nbins,TH2min,-TH2min,TH2nbins,TH2min,-TH2min),
                                                         'ML1':{},
                                                         'ML2':{},
                                                         'PL1':{},
                                                         'PL2':{}})

for key in ['ML1','ML2','PL1','PL2']:
    THSanityChecks['Occupancy']['BeforeMatching'][key]['RecHits'] = ROOT.TH2F(key+"_RecHitsOccupancy",key+"_RecHitsOccupancy",38,-0.5,37.5,10,-0.5,9.5)
    THSanityChecks['Occupancy']['BeforeMatching'][key]['PropHits'] = ROOT.TH2F(key+"_PropHitsOccupancy",key+"_PropHitsOccupancy",38,-0.5,37.5,10,-0.5,9.5)
    THSanityChecks['Occupancy']['BeforeMatching'][key]['RecHits'].GetXaxis().SetTitle("Chamber Number")
    THSanityChecks['Occupancy']['BeforeMatching'][key]['RecHits'].GetYaxis().SetTitle("EtaPartition")
    THSanityChecks['Occupancy']['BeforeMatching'][key]['PropHits'].GetXaxis().SetTitle("Chamber Number")
    THSanityChecks['Occupancy']['BeforeMatching'][key]['PropHits'].GetYaxis().SetTitle("EtaPartition")

    THSanityChecks['RecHitperStrip'][key] = {}
    for ch in range(1,37):
        size = "S" if ch%2 == 1 else "L"
        chID = 'GE11-'+key[0]+'-%02d' % ch + key[1:]+"-"+size 
        THSanityChecks['RecHitperStrip'][key][ch] = ROOT.TH2F(chID,chID,384,-0.5,383.5,10,-0.5,9.5)
        THSanityChecks['RecHitperStrip'][key][ch].SetStats(0)
        THSanityChecks['RecHitperStrip'][key][ch].GetXaxis().SetTitle("StripNumber")
        THSanityChecks['RecHitperStrip'][key][ch].GetYaxis().SetTitle("EtaPartition")


THSanityChecks['PropagationError']['glb_phi_error'] = {'all':ROOT.TH1F("All_errProp_glb_phi","All_errProp_glb_phi",100,0,0.0025),
                                                        'long':ROOT.TH1F("Long_errProp_glb_phi","Long_errProp_glb_phi",100,0,0.0025),
                                                        'long_isGEM':ROOT.TH1F("Long_errProp_glb_phi && isGEM==1","Long_errProp_glb_phi && isGEM==1",100,0,0.0025),
                                                        'long_noGEM':ROOT.TH1F("Long_errProp_glb_phi && isGEM==0","Long_errProp_glb_phi && isGEM==0",100,0,0.0025),
                                                        'short':ROOT.TH1F("Short_errProp_glb_phi","Short_errProp_glb_phi",100,0,0.0025),
                                                        'short_isGEM':ROOT.TH1F("Short_errProp_glb_phi && isGEM==1","Short_errProp_glb_phi && isGEM==1",100,0,0.0025),
                                                        'short_noGEM':ROOT.TH1F("Short_errProp_glb_phi && isGEM==0","Short_errProp_glb_phi && isGEM==0",100,0,0.0025),
                                                        'eta1':ROOT.TH1F("eta1_errProp_glb_phi","eta1_errProp_glb_phi",100,0,0.0025),
                                                        'eta2':ROOT.TH1F("eta2_errProp_glb_phi","eta2_errProp_glb_phi",100,0,0.0025),
                                                        'eta3':ROOT.TH1F("eta3_errProp_glb_phi","eta3_errProp_glb_phi",100,0,0.0025),
                                                        'eta4':ROOT.TH1F("eta4_errProp_glb_phi","eta4_errProp_glb_phi",100,0,0.0025),
                                                        'eta5':ROOT.TH1F("eta5_errProp_glb_phi","eta5_errProp_glb_phi",100,0,0.0025),
                                                        'eta6':ROOT.TH1F("eta6_errProp_glb_phi","eta6_errProp_glb_phi",100,0,0.0025),
                                                        'eta7':ROOT.TH1F("eta7_errProp_glb_phi","eta7_errProp_glb_phi",100,0,0.0025),
                                                        'eta8':ROOT.TH1F("eta8_errProp_glb_phi","eta8_errProp_glb_phi",100,0,0.0025)}
THSanityChecks['PropagationError']['glb_r_error'] = {'all':ROOT.TH1F("All_errProp_glb_r","All_errProp_glb_r",100,0,5),
                                                        'long':ROOT.TH1F("Long_errProp_glb_r","Long_errProp_glb_r",100,0,5),
                                                        'long_isGEM':ROOT.TH1F("Long_errProp_glb_r && isGEM==1","Long_errProp_glb_r && isGEM==1",100,0,5),
                                                        'long_noGEM':ROOT.TH1F("Long_errProp_glb_r && isGEM==0","Long_errProp_glb_r && isGEM==0",100,0,5),
                                                        'short':ROOT.TH1F("Short_errProp_glb_r","Short_errProp_glb_r",100,0,5),
                                                        'short_isGEM':ROOT.TH1F("Short_errProp_glb_r && isGEM==1","Short_errProp_glb_r && isGEM==1",100,0,5),
                                                        'short_noGEM':ROOT.TH1F("Short_errProp_glb_r && isGEM==0","Short_errProp_glb_r && isGEM==0",100,0,5),
                                                        'eta1':ROOT.TH1F("eta1_errProp_glb_r","eta1_errProp_glb_r",100,0,5),
                                                        'eta2':ROOT.TH1F("eta2_errProp_glb_r","eta2_errProp_glb_r",100,0,5),
                                                        'eta3':ROOT.TH1F("eta3_errProp_glb_r","eta3_errProp_glb_r",100,0,5),
                                                        'eta4':ROOT.TH1F("eta4_errProp_glb_r","eta4_errProp_glb_r",100,0,5),
                                                        'eta5':ROOT.TH1F("eta5_errProp_glb_r","eta5_errProp_glb_r",100,0,5),
                                                        'eta6':ROOT.TH1F("eta6_errProp_glb_r","eta6_errProp_glb_r",100,0,5),
                                                        'eta7':ROOT.TH1F("eta7_errProp_glb_r","eta7_errProp_glb_r",100,0,5),
                                                        'eta8':ROOT.TH1F("eta8_errProp_glb_r","eta8_errProp_glb_r",100,0,5)}
THSanityChecks['etaP_vs_pt'] = ROOT.TH2F("allChmbrs_etaP_pt","allChmbrs_etaP_pt",8,0,8,11,0,110)
THSanityChecks['etaP_vs_pt'].GetXaxis().SetTitle("i#eta")
THSanityChecks['etaP_vs_pt'].GetYaxis().SetTitle("pt (GeV)")
THSanityChecks['etaP_vs_pt'].SetStats(0)

THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'] = ROOT.TH2F("Residual_Correlation #Delta#phi vs R#Delta#phi","Residual_Correlation #Delta#phi vs R#Delta#phi",100,-3*ResidualCutOff['glb_phi'],3*ResidualCutOff['glb_phi'],100,-3*ResidualCutOff['glb_rdphi'],3*ResidualCutOff['glb_rdphi'])
THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'] = ROOT.TH2F("Residual_Correlation R#Delta#phi vs Dir_x","Residual_Correlation R#Delta#phi vs Dir_x",100,-3*ResidualCutOff['glb_rdphi'],3*ResidualCutOff['glb_rdphi'],100,0,3.1415)
THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'].GetXaxis().SetTitle("#Delta#phi (rad)")
THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'].GetYaxis().SetTitle("R#Delta#phi (cm)")
THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'].SetStats(0)
THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'].GetXaxis().SetTitle("R#Delta#phi (cm)")
THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'].GetYaxis().SetTitle("Dir_x (as Cos(#alpha) )")

THSanityChecks['STA_Normchi2'] = ROOT.TH1F("STA_NormChi2","STA_NormChi2",200,0,20)
THSanityChecks['nME1Hits'] = ROOT.TH1F("nME1Hits in STA","nME1Hits in STA",20,0,20)
THSanityChecks['nME2Hits'] = ROOT.TH1F("nME2Hits in STA","nME2Hits in STA",20,0,20)
THSanityChecks['nME3Hits'] = ROOT.TH1F("nME3Hits in STA","nME3Hits in STA",20,0,20)
THSanityChecks['nME4Hits'] = ROOT.TH1F("nME4Hits in STA","nME4Hits in STA",20,0,20)
THSanityChecks['nCSCHits'] = ROOT.TH1F("nCSCHits in STA","nCSCHits in STA",40,0,40)

for key_1 in matching_variables:
    THSanityChecks['Occupancy'].setdefault(key_1,{'AfterMatching':{'Reco':ROOT.TH2F("RecoHitAfterMatching_"+key_1,"RecoHitAfterMatching_"+key_1,200,-300,300,200,-300,300),
                                                                   'Prop':ROOT.TH2F("PropHitAfterMatching_"+key_1,"PropHitAfterMatching_"+key_1,200,-300,300,200,-300,300),
                                                                   'PropLocalLong':ROOT.TH2F("PropLocLongHitAfterMatching_"+key_1,"PropLocLongHitAfterMatching_"+key_1,TH2nbins,TH2min,-TH2min,TH2nbins,TH2min,-TH2min),
                                                                   'PropLocalShort':ROOT.TH2F("PropLocShortHitAfterMatching_"+key_1,"PropLocShortHitAfterMatching_"+key_1,TH2nbins,TH2min,-TH2min,TH2nbins,TH2min,-TH2min)}})

    TH1Fresidual_collector.setdefault(key_1,{'all':{},'short':{},'long':{},'ML1':{},'ML2':{},'PL1':{},'PL2':{}})
    for key_2 in ['all','short','long']+["eta"+str(k) for k in range(1,9)]:
        TH1Fresidual_collector[key_1][key_2] = {}
        for key_3 in ['glb_phi','glb_rdphi']:
            titleTH1 = key_2 + "_MatchedBy_"+key_1+"_"+key_3+"_residuals" 
            min_x = -ResidualCutOff[key_3]*2
            TH1Fresidual_collector[key_1][key_2][key_3]= ROOT.TH1F(titleTH1,titleTH1,TH1nbins,min_x,-min_x)
            TH1Fresidual_collector[key_1][key_2][key_3].GetXaxis().SetTitle(key_3+" Residuals ("+matching_variable_units[key_3]+")")
    if key_1 == 'glb_rdphi':
        for key_4 in ['ML1','ML2','PL1','PL2']:
            for ch in range(1,37):
                size = "S" if ch%2 == 1 else "L"
                chID = 'GE11-'+key_4[0]+'-%02d' % ch + key_4[1:]+"-"+size

                titleTH1 = "ResidualsOn_"+chID+"_MatchedBy_"+key_1
                min_x = -ResidualCutOff['glb_rdphi']*2
                TH1Fresidual_collector[key_1][key_4][ch]= ROOT.TH1F(titleTH1,titleTH1,TH1nbins,min_x,-min_x)

for key_1 in ['Long','Short']:
    for key_2 in ["eta"+str(k) for k in range(1,9)]:
        THSanityChecks['PropHit_DirLoc_xOnGE11']['BeforeMatching'][key_1][key_2] = ROOT.TH1F('BeforeMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,'BeforeMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,200,-1,1)
        THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching'][key_1][key_2] = ROOT.TH1F('AfterMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,'AfterMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,200,-1,1)

## Chain files
chain = ROOT.TChain("muNtupleProducer/MuDPGTree")

for fl in files:
    chain.Add(fl)

chainEntries = chain.GetEntries()


for chain_index,evt in enumerate(chain):
    
    if chain_index % 40000 ==0:
        print round(float(chain_index)/float(chainEntries),3)*100,"%"


    n_gemprop = len(evt.mu_propagated_chamber)
    n_gemrec = len(evt.gemRecHit_chamber)
    
    
    THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Prop'].Fill(n_gemprop)
    THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Reco'].Fill(n_gemrec)
    

    ## Discard null entries
    if n_gemprop==0:
        continue


    ML1_NGEMRecoHits = sum( [1 for RecHit_index in range(0,n_gemrec) if evt.gemRecHit_region[RecHit_index] == -1 and evt.gemRecHit_layer[RecHit_index] == 1 and evt.gemRecHit_region[RecHit_index]*(100*evt.gemRecHit_chamber[RecHit_index]+10*evt.gemRecHit_layer[RecHit_index]+evt.gemRecHit_etaPartition[RecHit_index]) not in noisyEtaPID ] )
    ML2_NGEMRecoHits = sum( [1 for RecHit_index in range(0,n_gemrec) if evt.gemRecHit_region[RecHit_index] == -1 and evt.gemRecHit_layer[RecHit_index] == 2 and evt.gemRecHit_region[RecHit_index]*(100*evt.gemRecHit_chamber[RecHit_index]+10*evt.gemRecHit_layer[RecHit_index]+evt.gemRecHit_etaPartition[RecHit_index]) not in noisyEtaPID ] )
    PL1_NGEMRecoHits = sum( [1 for RecHit_index in range(0,n_gemrec) if evt.gemRecHit_region[RecHit_index] == 1 and evt.gemRecHit_layer[RecHit_index] == 1 and evt.gemRecHit_region[RecHit_index]*(100*evt.gemRecHit_chamber[RecHit_index]+10*evt.gemRecHit_layer[RecHit_index]+evt.gemRecHit_etaPartition[RecHit_index]) not in noisyEtaPID ] )
    PL2_NGEMRecoHits = sum( [1 for RecHit_index in range(0,n_gemrec) if evt.gemRecHit_region[RecHit_index] == 1 and evt.gemRecHit_layer[RecHit_index] == 2 and evt.gemRecHit_region[RecHit_index]*(100*evt.gemRecHit_chamber[RecHit_index]+10*evt.gemRecHit_layer[RecHit_index]+evt.gemRecHit_etaPartition[RecHit_index]) not in noisyEtaPID ] )

    THSanityChecks['NHits']['BeforeMatching']['ML1'].Fill(ML1_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['ML2'].Fill(ML2_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['PL1'].Fill(PL1_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['PL2'].Fill(PL2_NGEMRecoHits)

    
    RecHit_Dict = {}
    PropHit_Dict = {}

    for RecHit_index in range(0,n_gemrec):
        region = evt.gemRecHit_region[RecHit_index]
        chamber = evt.gemRecHit_chamber[RecHit_index]
        layer = evt.gemRecHit_layer[RecHit_index]
        etaP = evt.gemRecHit_etaPartition[RecHit_index]
        RecHitEtaPartitionID = region*(100*chamber+10*layer+etaP)
        endcapKey = "PL"+str(layer) if region > 0 else "ML"+str(layer)

        ## discard chambers that were kept OFF from the analysis
        if [region,chamber,layer] in chamberNumberOFF:
            continue


        RecHit_Dict.setdefault(RecHitEtaPartitionID, {'loc_x':[],'glb_x':[],'glb_y':[],'glb_z':[],'glb_r':[],'glb_phi':[],'firstStrip':[],'cluster_size':[]})
        RecHit_Dict[RecHitEtaPartitionID]['loc_x'].append(evt.gemRecHit_loc_x[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_x'].append(evt.gemRecHit_g_x[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_y'].append(evt.gemRecHit_g_y[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_z'].append(evt.gemRecHit_g_z[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_r'].append(evt.gemRecHit_g_r[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_phi'].append(evt.gemRecHit_g_phi[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['firstStrip'].append(evt.gemRecHit_firstClusterStrip[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['cluster_size'].append(evt.gemRecHit_cluster_size[RecHit_index])

        THSanityChecks['Occupancy']['BeforeMatching']['Reco'].Fill(evt.gemRecHit_g_x[RecHit_index],evt.gemRecHit_g_y[RecHit_index])
        THSanityChecks['Occupancy']['BeforeMatching'][endcapKey]['RecHits'].Fill(chamber,etaP)
        
        for j in range(0,RecHit_Dict[RecHitEtaPartitionID]['cluster_size'][-1]):
            strip = RecHit_Dict[RecHitEtaPartitionID]['firstStrip'][-1] + j
            THSanityChecks['RecHitperStrip'][endcapKey][chamber].Fill(strip,etaP)
                    

        
    
    for PropHit_index in range(0,n_gemprop):
        
        region = evt.mu_propagated_region[PropHit_index]
        chamber = evt.mu_propagated_chamber[PropHit_index]
        layer = evt.mu_propagated_layer[PropHit_index]
        etaP = evt.mu_propagated_etaP[PropHit_index]
        PropHitChamberID = region*(100*chamber+10*layer+etaP)
        endcapKey = "PL"+str(layer) if region > 0 else "ML"+str(layer)

        outermost_z = evt.mu_propagated_Outermost_z[PropHit_index]
        is_incoming = evt.mu_isincoming[PropHit_index]

        ## discard chambers that were kept OFF from the analysis
        if [region,chamber,layer] in chamberNumberOFF:
            continue

        # if region == 1 and outermost_z < 0:
        #     continue
        # if region == -1 and outermost_z > 0:
        #     continue

        propHitFromME11 = bool(evt.mu_propagated_isME11[PropHit_index])
        if propHitFromME11:
            PropHit_Dict.setdefault(PropHitChamberID,{'loc_x':[],'loc_y':[],'glb_x':[],'glb_y':[],'glb_z':[],'glb_r':[],'glb_phi':[],'pt':[],'etaP':[],'err_glb_r':[],'err_glb_phi':[],'Loc_dirX':[],'Loc_dirY':[],'Loc_dirZ':[],'mu_propagated_isME11':[],'mu_propagated_EtaPartition_rMax':[],'mu_propagated_EtaPartition_rMin':[],'mu_propagated_isGEM':[],'mu_propagated_EtaPartition_phiMin':[],'mu_propagated_EtaPartition_phiMax':[],'STA_Normchi2':[],'nME1Hits':[],'nME2Hits':[],'nME3Hits':[],'nME4Hits':[]})    
            PropHit_Dict[PropHitChamberID]['loc_x'].append(evt.mu_propagatedLoc_x[PropHit_index])
            PropHit_Dict[PropHitChamberID]['loc_y'].append(evt.mu_propagatedLoc_y[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_x'].append(evt.mu_propagatedGlb_x[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_y'].append(evt.mu_propagatedGlb_y[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_z'].append(evt.mu_propagatedGlb_z[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_r'].append(evt.mu_propagatedGlb_r[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_phi'].append(evt.mu_propagatedGlb_phi[PropHit_index])
            PropHit_Dict[PropHitChamberID]['err_glb_r'].append(evt.mu_propagatedGlb_errR[PropHit_index])
            PropHit_Dict[PropHitChamberID]['err_glb_phi'].append(evt.mu_propagatedGlb_errPhi[PropHit_index])
            PropHit_Dict[PropHitChamberID]['Loc_dirX'].append(evt.mu_propagatedLoc_dirX[PropHit_index])
            PropHit_Dict[PropHitChamberID]['Loc_dirY'].append(evt.mu_propagatedLoc_dirY[PropHit_index])
            PropHit_Dict[PropHitChamberID]['Loc_dirZ'].append(evt.mu_propagatedLoc_dirZ[PropHit_index])
            PropHit_Dict[PropHitChamberID]['pt'].append(evt.mu_propagated_pt[PropHit_index])
            PropHit_Dict[PropHitChamberID]['etaP'].append(etaP)
            PropHit_Dict[PropHitChamberID]['mu_propagated_isME11'].append(evt.mu_propagated_isME11[PropHit_index])
            PropHit_Dict[PropHitChamberID]['mu_propagated_isGEM'].append(evt.mu_propagated_isGEM[PropHit_index])
            PropHit_Dict[PropHitChamberID]['mu_propagated_EtaPartition_rMax'].append(evt.mu_propagated_EtaPartition_rMax[PropHit_index])
            PropHit_Dict[PropHitChamberID]['mu_propagated_EtaPartition_rMin'].append(evt.mu_propagated_EtaPartition_rMin[PropHit_index])
            PropHit_Dict[PropHitChamberID]['mu_propagated_EtaPartition_phiMax'].append(evt.mu_propagated_EtaPartition_phiMax[PropHit_index])
            PropHit_Dict[PropHitChamberID]['mu_propagated_EtaPartition_phiMin'].append(evt.mu_propagated_EtaPartition_phiMin[PropHit_index])
            PropHit_Dict[PropHitChamberID]['STA_Normchi2'].append(evt.mu_propagated_TrackNormChi2[PropHit_index])            
            PropHit_Dict[PropHitChamberID]['nME1Hits'].append(evt.mu_propagated_nME1hits[PropHit_index])
            PropHit_Dict[PropHitChamberID]['nME2Hits'].append(evt.mu_propagated_nME2hits[PropHit_index])
            PropHit_Dict[PropHitChamberID]['nME3Hits'].append(evt.mu_propagated_nME3hits[PropHit_index])
            PropHit_Dict[PropHitChamberID]['nME4Hits'].append(evt.mu_propagated_nME4hits[PropHit_index])

            THSanityChecks['Occupancy']['BeforeMatching']['Prop'].Fill(evt.mu_propagatedGlb_x[PropHit_index],evt.mu_propagatedGlb_y[PropHit_index])
            THSanityChecks['Occupancy']['BeforeMatching'][endcapKey]['PropHits'].Fill(chamber,etaP)
            THSanityChecks['etaP_vs_pt'].Fill(PropHit_Dict[PropHitChamberID]['etaP'][-1]-1,10*pt_index(evt.mu_propagated_pt[PropHit_index]))
            THSanityChecks['STA_Normchi2'].Fill(evt.mu_propagated_TrackNormChi2[PropHit_index])


            if chamber % 2 == 0:
                THSanityChecks['Occupancy']['BeforeMatching']['PropLocalLong'].Fill(evt.mu_propagatedLoc_x[PropHit_index],evt.mu_propagatedLoc_y[PropHit_index])
                THSanityChecks['PropHit_DirLoc_xOnGE11']['BeforeMatching']['Long']['eta'+str(etaP)].Fill(evt.mu_propagatedLoc_dirX[PropHit_index])
            if chamber % 2 == 1:
                THSanityChecks['Occupancy']['BeforeMatching']['PropLocalShort'].Fill(evt.mu_propagatedLoc_x[PropHit_index],evt.mu_propagatedLoc_y[PropHit_index])
                THSanityChecks['PropHit_DirLoc_xOnGE11']['BeforeMatching']['Short']['eta'+str(etaP)].Fill(evt.mu_propagatedLoc_dirX[PropHit_index])

    ML1_N_MatchedGEMRecoHits = 0
    ML2_N_MatchedGEMRecoHits = 0
    PL1_N_MatchedGEMRecoHits = 0
    PL2_N_MatchedGEMRecoHits = 0

    ### Matching criteria between propagated hits from ME11 and RecHits : 
    ##      1.SAME REGION,SC,LAYER,ETA -->SAME etaPartitionID
    ##      When using DLE, only evts w/ exactly 2 PropHit in a SC: 1 hit per Layer with Delta(etaP) < 4

    if DLE and (len(PropHit_Dict.keys()) != 2 or abs(PropHit_Dict.keys()[0] - PropHit_Dict.keys()[1] ) > 13) :
        continue

    layer1Match = False
    layer2Match = False
    layer1PassedCut = False
    layer2PassedCut = False
    layer1pt = 0
    layer1etapID = 0
    layer2pt = 0
    layer2etapID = 0
    for etaPartitionID in PropHit_Dict.keys():
        
        region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)
        PropHitonEta = PropHit_Dict[etaPartitionID]


        nPropHitsOnEtaID = len(PropHitonEta['glb_phi'])
        THSanityChecks['NHits']['PerEVT_PerEtaPartitionID']['Prop'].Fill(nPropHitsOnEtaID)
        cos_of_alpha_list = [np.sqrt(PropHitonEta['Loc_dirX'][i]**2 + PropHitonEta['Loc_dirY'][i]**2) for i in range(nPropHitsOnEtaID)]        

        ## Defining Efficiency dict global: [matchingVar][etaPartitionID][pt]
        for mv in matching_variables:
            EfficiencyDictGlobal[mv].setdefault(etaPartitionID,{})
            EfficiencyDictDirection[mv].setdefault(etaPartitionID,{})
            EfficiencyDictLayer[mv].setdefault(etaPartitionID,{})
            for pt in range(0,11):
                EfficiencyDictGlobal[mv][etaPartitionID].setdefault(pt,{'num':0,'den':0})
                EfficiencyDictLayer[mv][etaPartitionID].setdefault(pt,{'num':0,'den':0})
            for j in range(0,10):
                EfficiencyDictDirection[mv][etaPartitionID].setdefault(j,{'num':0,'den':0})

        isGoodTrack = []
        passedCutProp = {key:[] for key in PropHitonEta.keys()}
        ## Applying cuts on the propagated tracks to be used
        for index in range(nPropHitsOnEtaID):
            if fiducialCut and passCut(PropHitonEta,index,maxPropR_Err=maxErrOnPropR,maxPropPhi_Err=maxErrOnPropPhi,fiducialCutR=fiducialR,fiducialCutPhi=fiducialPhi,minPt=CutminPt,maxChi2=maxSTA_NormChi2,minME1Hit=minME1Hit,minME2Hit=minME2Hit,minME3Hit=minME3Hit,minME4Hit=minME4Hit) == False:
                isGoodTrack.append(False)
            else:
                EfficiencyDictGlobal['glb_phi'][etaPartitionID][pt_index(PropHitonEta['pt'][index])]['den'] += 1
                EfficiencyDictGlobal['glb_rdphi'][etaPartitionID][pt_index(PropHitonEta['pt'][index])]['den'] += 1

                angle_index = int ( (cos_of_alpha_list[index] * 10 ) )
                EfficiencyDictDirection['glb_phi'][etaPartitionID][angle_index]['den'] += 1 
                EfficiencyDictDirection['glb_rdphi'][etaPartitionID][angle_index]['den'] += 1 

                isGoodTrack.append(True)
                for key in PropHitonEta.keys():
                    passedCutProp[key].append(PropHitonEta[key][index])
        
        #any is the logical or across all elements of a list
        if any(isGoodTrack) == False:
            #print "No good STA propagation for etaPartitionID =  ",etaPartitionID
            continue
    
        PropHitonEta = passedCutProp
        nGoodPropagation = len(PropHitonEta['glb_phi'])        
        if DLE:    
            if layer == 1:
                layer1PassedCut = True
                layer1etapID = etaPartitionID
                layer1pt = 5 ## Fake pt value in GeV ... no B field
            if layer == 2:
                layer2PassedCut = True
                layer2etapID = etaPartitionID
                layer2pt = 5 ## Fake pt value in GeV ... no B field

        ## Filling STA properties in histos
        for k in range(nGoodPropagation):
            THSanityChecks['PropagationError']['glb_phi_error']['all'].Fill(PropHitonEta['err_glb_phi'][k])
            THSanityChecks['PropagationError']['glb_r_error']['all'].Fill(PropHitonEta['err_glb_r'][k])
            THSanityChecks['nME1Hits'].Fill(PropHitonEta['nME1Hits'][k])
            THSanityChecks['nME2Hits'].Fill(PropHitonEta['nME2Hits'][k])
            THSanityChecks['nME3Hits'].Fill(PropHitonEta['nME3Hits'][k])
            THSanityChecks['nME4Hits'].Fill(PropHitonEta['nME4Hits'][k])
            THSanityChecks['nCSCHits'].Fill( PropHitonEta['nME1Hits'][k] + PropHitonEta['nME2Hits'][k] + PropHitonEta['nME3Hits'][k] + PropHitonEta['nME4Hits'][k] )

            if chamber%2 == 0:
                THSanityChecks['PropagationError']['glb_phi_error']['long'].Fill(PropHitonEta['err_glb_phi'][k])
                THSanityChecks['PropagationError']['glb_r_error']['long'].Fill(PropHitonEta['err_glb_r'][k])
                if PropHitonEta['mu_propagated_isGEM'][k] == True:
                    THSanityChecks['PropagationError']['glb_phi_error']['long_isGEM'].Fill(PropHitonEta['err_glb_phi'][k])
                    THSanityChecks['PropagationError']['glb_r_error']['long_isGEM'].Fill(PropHitonEta['err_glb_r'][k])
                elif PropHitonEta['mu_propagated_isGEM'][k] == False:
                    THSanityChecks['PropagationError']['glb_phi_error']['long_noGEM'].Fill(PropHitonEta['err_glb_phi'][k])
                    THSanityChecks['PropagationError']['glb_r_error']['long_noGEM'].Fill(PropHitonEta['err_glb_r'][k])
            if chamber % 2 == 1:
                THSanityChecks['PropagationError']['glb_phi_error']['short'].Fill(PropHitonEta['err_glb_phi'][k])
                THSanityChecks['PropagationError']['glb_r_error']['short'].Fill(PropHitonEta['err_glb_r'][k])
                if PropHitonEta['mu_propagated_isGEM'][k] == True:
                    THSanityChecks['PropagationError']['glb_phi_error']['short_isGEM'].Fill(PropHitonEta['err_glb_phi'][k])
                    THSanityChecks['PropagationError']['glb_r_error']['short_isGEM'].Fill(PropHitonEta['err_glb_r'][k])
                elif PropHitonEta['mu_propagated_isGEM'][k] == False:
                    THSanityChecks['PropagationError']['glb_phi_error']['short_noGEM'].Fill(PropHitonEta['err_glb_phi'][k])
                    THSanityChecks['PropagationError']['glb_r_error']['short_noGEM'].Fill(PropHitonEta['err_glb_r'][k])
            
            THSanityChecks['PropagationError']['glb_phi_error']['eta'+str(eta)].Fill(PropHitonEta['err_glb_phi'][k])
            THSanityChecks['PropagationError']['glb_r_error']['eta'+str(eta)].Fill(PropHitonEta['err_glb_r'][k])


        if etaPartitionID not in RecHit_Dict:
            # print "Nothing to match on ", etaPartitionID
            # print PropHit_Dict[etaPartitionID]['glb_phi']
            # raw_input()
            #print "No rechit in etaPartitionID =  ",etaPartitionID
            continue
        else: 
            RecHitonEta = RecHit_Dict[etaPartitionID]
    
        THSanityChecks['NHits']['PerEVT_PerEtaPartitionID']['Reco'].Fill(len(RecHitonEta['glb_phi']))


        ## Seek for 1-best match between rec and prop based on Matching Var and Cutoff
        for matchingVar in matching_variables:
                               
            residuals = []
            RecoMatched = []
            PropMatched = []
            if matchingVar == 'glb_rdphi':
                for RecHit_g_phi in RecHitonEta['glb_phi']:
                    deltardphis = [(PropHit_g_phi - RecHit_g_phi)*PropHitonEta['glb_r'][index] for index,PropHit_g_phi in enumerate(PropHitonEta['glb_phi'])]
                    temp_min = min(deltardphis,key=abs)
                    PropMatched.append(PropHitonEta['glb_phi'][deltardphis.index(temp_min)])
                    RecoMatched.append(RecHit_g_phi)
                    residuals.append(temp_min)
                
                min_residual = min(residuals,key=abs)
                min_residual_index = residuals.index(min_residual)
                prop_hit_index = PropHitonEta['glb_phi'].index(PropMatched[min_residual_index])
                reco_hit_index = RecHitonEta['glb_phi'].index(RecoMatched[min_residual_index])




            else:
                for RecHit_var in RecHitonEta[matchingVar]:
                    deltas = [PropHit_var - RecHit_var for PropHit_var in PropHitonEta[matchingVar]]
                    temp_min = min(deltas,key=abs)
                    PropMatched.append(PropHitonEta[matchingVar][deltas.index(temp_min)])
                    RecoMatched.append(RecHit_var)
                    residuals.append(temp_min)
                
                min_residual = min(residuals,key=abs)
                min_residual_index = residuals.index(min_residual)
                prop_hit_index = PropHitonEta[matchingVar].index(PropMatched[min_residual_index])
                reco_hit_index = RecHitonEta[matchingVar].index(RecoMatched[min_residual_index])
            
            glb_phi_residual = PropHitonEta['glb_phi'][prop_hit_index] - RecHitonEta['glb_phi'][reco_hit_index]
            glb_rdphi_residual = (PropHitonEta['glb_phi'][prop_hit_index] - RecHitonEta['glb_phi'][reco_hit_index])*PropHitonEta['glb_r'][prop_hit_index]
            loc_x_residual = PropHitonEta['loc_x'][prop_hit_index] - RecHitonEta['loc_x'][reco_hit_index]

            if matchingVar == 'glb_phi':
                THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'].Fill(glb_phi_residual,glb_rdphi_residual)
                THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'].Fill(glb_rdphi_residual,np.arccos(PropHitonEta['Loc_dirX'][prop_hit_index]))


            if abs(min_residual) < ResidualCutOff[matchingVar]:

                
                EfficiencyDictGlobal[matchingVar][etaPartitionID][pt_index(PropHitonEta['pt'][prop_hit_index])]['num'] += 1
                angle_index = int( cos_of_alpha_list[prop_hit_index] * 10)
                EfficiencyDictDirection[matchingVar][etaPartitionID][angle_index]['num'] += 1

                TH1Fresidual_collector[matchingVar]['all']['glb_phi'].Fill(glb_phi_residual)
                TH1Fresidual_collector[matchingVar]['all']['glb_rdphi'].Fill(glb_rdphi_residual)

                binx = int(round((PropHitonEta['loc_x'][prop_hit_index]-TH2min)*(TH2nbins-1)/(-2*TH2min)))+1
                biny = int(round((PropHitonEta['loc_y'][prop_hit_index]-TH2min)*(TH2nbins-1)/(-2*TH2min)))+1
                TH2Fresidual_collector[matchingVar]['all']['glb_phi'][binx][biny][0] += 1
                TH2Fresidual_collector[matchingVar]['all']['glb_rdphi'][binx][biny][0] += 1
                TH2Fresidual_collector[matchingVar]['all']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                TH2Fresidual_collector[matchingVar]['all']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)

                THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Reco'].Fill(RecHitonEta['glb_x'][reco_hit_index],RecHitonEta['glb_y'][reco_hit_index])
                THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Prop'].Fill(PropHitonEta['glb_x'][prop_hit_index],PropHitonEta['glb_y'][prop_hit_index])

                if chamber%2 == 0:
                    if matchingVar == 'glb_phi': THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching']['Long']['eta'+str(eta)].Fill(PropHitonEta['Loc_dirX'][prop_hit_index])
                    TH1Fresidual_collector[matchingVar]['long']['glb_phi'].Fill(glb_phi_residual)
                    TH1Fresidual_collector[matchingVar]['long']['glb_rdphi'].Fill(glb_rdphi_residual)
                    TH2Fresidual_collector[matchingVar]['long']['glb_phi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['long']['glb_rdphi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['long']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                    TH2Fresidual_collector[matchingVar]['long']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)
                    THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalLong'].Fill(PropHitonEta['loc_x'][prop_hit_index],PropHitonEta['loc_y'][prop_hit_index])
                if chamber%2 == 1:
                    if matchingVar == 'glb_phi': THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching']['Short']['eta'+str(eta)].Fill(PropHitonEta['Loc_dirX'][prop_hit_index])
                    TH1Fresidual_collector[matchingVar]['short']['glb_phi'].Fill(glb_phi_residual)
                    TH1Fresidual_collector[matchingVar]['short']['glb_rdphi'].Fill(glb_rdphi_residual)  
                    TH2Fresidual_collector[matchingVar]['short']['glb_phi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['short']['glb_rdphi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['short']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                    TH2Fresidual_collector[matchingVar]['short']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)                    
                    THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalShort'].Fill(PropHitonEta['loc_x'][prop_hit_index],PropHitonEta['loc_y'][prop_hit_index])
                
                TH1Fresidual_collector[matchingVar]['eta'+str(eta)]['glb_phi'].Fill(glb_phi_residual)
                TH1Fresidual_collector[matchingVar]['eta'+str(eta)]['glb_rdphi'].Fill(glb_rdphi_residual)

                
                if matchingVar == 'glb_phi':                    
                    if region == -1 and layer == 1:
                        ML1_N_MatchedGEMRecoHits += 1
                    if region == -1 and layer == 2:
                        ML2_N_MatchedGEMRecoHits += 1
                    if region == 1 and layer == 1:
                        PL1_N_MatchedGEMRecoHits += 1
                    if region == 1 and layer == 2:
                        PL2_N_MatchedGEMRecoHits += 1
                if matchingVar == 'glb_rdphi':
                    endcap = "M" if region == -1 else "P"
                    TH1Fresidual_collector[matchingVar][endcap+"L"+str(layer)][chamber].Fill(glb_rdphi_residual)
                    test.Fill(PropHitonEta['Loc_dirX'][prop_hit_index],RecHitonEta['cluster_size'][reco_hit_index])
                    
                    if layer == 1:
                        layer1Match = True
                    if layer == 2:
                        layer2Match = True
                    
            else:
                # print "Matching failed for ", etaPartitionID
                # print PropHit_Dict[etaPartitionID]['glb_phi'], RecHit_Dict[etaPartitionID]['glb_phi']
                # raw_input()
                pass
                
        ## Loop over etaPID 

    ## Double Layer Efficiency (DLE): test layer1(2) with tracks that do have a match in layer1(2)

    if DLE and layer2Match and layer1PassedCut and layer2PassedCut:
        EfficiencyDictLayer['glb_rdphi'][layer1etapID][pt_index(layer1pt)]['den'] += 1
        DLE_ErrPhi.Fill(PropHitonEta['err_glb_phi'][prop_hit_index])
        DLE_ErrR.Fill(PropHitonEta['err_glb_r'][prop_hit_index])

        if layer1Match == True:
            EfficiencyDictLayer['glb_rdphi'][layer1etapID][pt_index(layer1pt)]['num'] += 1

    if DLE and layer1Match and layer2PassedCut and layer1PassedCut:
        EfficiencyDictLayer['glb_rdphi'][layer2etapID][pt_index(layer2pt)]['den'] += 1
        DLE_ErrPhi.Fill(PropHitonEta['err_glb_phi'][prop_hit_index])
        DLE_ErrR.Fill(PropHitonEta['err_glb_r'][prop_hit_index])

        if layer2Match:
            EfficiencyDictLayer['glb_rdphi'][layer2etapID][pt_index(layer1pt)]['num'] += 1

    ## Loop over evts

    THSanityChecks['NHits']['AfterMatching']['All'].Fill(ML1_N_MatchedGEMRecoHits + ML2_N_MatchedGEMRecoHits  +  PL1_N_MatchedGEMRecoHits + PL2_N_MatchedGEMRecoHits)
    THSanityChecks['NHits']['AfterMatching']['ML1'].Fill(ML1_N_MatchedGEMRecoHits)
    THSanityChecks['NHits']['AfterMatching']['ML2'].Fill(ML2_N_MatchedGEMRecoHits)
    THSanityChecks['NHits']['AfterMatching']['PL1'].Fill(PL1_N_MatchedGEMRecoHits)
    THSanityChecks['NHits']['AfterMatching']['PL2'].Fill(PL2_N_MatchedGEMRecoHits)
## End of the evts loop


TH2Fresidual_collector = fillPlot2DResidualContainer(TH2Fresidual_collector,matching_variables,TH2nbins)

print("--- %s seconds ---" % (time.time() - start_time))



## Storing the results

timestamp =  time.strftime("%-y%m%d_%H%M") if args.outputname is None else args.outputname
subprocess.call(["mkdir", "-p", "./Plot/"+timestamp+"/"+"glb_phi/"])
subprocess.call(["mkdir", "-p", "./Plot/"+timestamp+"/"+"glb_rdphi/"])



OutF = ROOT.TFile("./"+timestamp+".root","RECREATE")

writeToTFile(OutF,test,"SanityChecks/AngleDirCorrelation/")

writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Reco'],"SanityChecks/NHits/BeforeMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Prop'],"SanityChecks/NHits/BeforeMatching/")

writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['ML1'],"SanityChecks/NHits/BeforeMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['ML2'],"SanityChecks/NHits/BeforeMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['PL1'],"SanityChecks/NHits/BeforeMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['BeforeMatching']['PL2'],"SanityChecks/NHits/BeforeMatching/")

writeToTFile(OutF,THSanityChecks['NHits']['AfterMatching']['All'],"SanityChecks/NHits/AfterMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['AfterMatching']['ML1'],"SanityChecks/NHits/AfterMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['AfterMatching']['ML2'],"SanityChecks/NHits/AfterMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['AfterMatching']['PL1'],"SanityChecks/NHits/AfterMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['AfterMatching']['PL2'],"SanityChecks/NHits/AfterMatching/")

writeToTFile(OutF,THSanityChecks['NHits']['PerEVT_PerEtaPartitionID']['Reco'],"SanityChecks/NHits/BeforeMatching/")
writeToTFile(OutF,THSanityChecks['NHits']['PerEVT_PerEtaPartitionID']['Prop'],"SanityChecks/NHits/BeforeMatching/")

writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching']['Prop'],"SanityChecks/Occupancy/BeforeMatching")
writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching']['PropLocalLong'],"SanityChecks/Occupancy/BeforeMatching")
writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching']['PropLocalShort'],"SanityChecks/Occupancy/BeforeMatching")
writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching']['Reco'],"SanityChecks/Occupancy/BeforeMatching")

for key in ['ML1','ML2','PL1','PL2']:
    writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching'][key]['PropHits'],"SanityChecks/Occupancy/BeforeMatching/"+key)
    writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching'][key]['RecHits'],"SanityChecks/Occupancy/BeforeMatching/"+key)
    for ch in range(1,37):
        writeToTFile(OutF,THSanityChecks['RecHitperStrip'][key][ch],"SanityChecks/Occupancy/RecHitByStrip/"+key)
        writeToTFile(OutF,TH1Fresidual_collector['glb_rdphi'][key][ch],"Residuals/MatchingOn_glb_rdphi/Residual_glb_rdphi/"+key)

writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['all'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['long'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['long_isGEM'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['long_noGEM'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['short'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['short_isGEM'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error']['short_noGEM'],"SanityChecks/PropagationError/glb_phi")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['all'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['long'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['long_isGEM'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['long_noGEM'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['short'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['short_isGEM'],"SanityChecks/PropagationError/glb_r")
writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error']['short_noGEM'],"SanityChecks/PropagationError/glb_r")

writeToTFile(OutF,THSanityChecks['etaP_vs_pt'],"SanityChecks/etaP_vs_pt/")

writeToTFile(OutF,THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'],"SanityChecks/Residual_Correlation/")
writeToTFile(OutF,THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'],"SanityChecks/Residual_Correlation/")
writeToTFile(OutF,THSanityChecks['STA_Normchi2'],"SanityChecks/STA_NormChi2/")
writeToTFile(OutF,THSanityChecks['nME1Hits'],"SanityChecks/HitsCSC/")
writeToTFile(OutF,THSanityChecks['nME2Hits'],"SanityChecks/HitsCSC/")
writeToTFile(OutF,THSanityChecks['nME3Hits'],"SanityChecks/HitsCSC/")
writeToTFile(OutF,THSanityChecks['nME4Hits'],"SanityChecks/HitsCSC/")
writeToTFile(OutF,THSanityChecks['nCSCHits'],"SanityChecks/HitsCSC/")




for key_1 in ['eta'+str(j) for j in range(1,9)]:
    writeToTFile(OutF,THSanityChecks['PropagationError']['glb_phi_error'][key_1],"SanityChecks/PropagationError/glb_phi/byEta")
    writeToTFile(OutF,THSanityChecks['PropagationError']['glb_r_error'][key_1],"SanityChecks/PropagationError/glb_r/byEta")
    for key_2 in ['Long','Short']:
        writeToTFile(OutF,THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching'][key_2][key_1],"SanityChecks/PropHitDirection/AfterMatching/"+key_2)
        writeToTFile(OutF,THSanityChecks['PropHit_DirLoc_xOnGE11']['BeforeMatching'][key_2][key_1],"SanityChecks/PropHitDirection/BeforeMatching/"+key_2)

for matchingVar in matching_variables:
    for chambers in ['all','long','short']+['eta'+str(j) for j in range(1,9)]:
        writeToTFile(OutF,TH1Fresidual_collector[matchingVar][chambers]['glb_phi'],"Residuals/MatchingOn_"+matchingVar+"/Residual_glb_phi")
        writeToTFile(OutF,TH1Fresidual_collector[matchingVar][chambers]['glb_rdphi'],"Residuals/MatchingOn_"+matchingVar+"/Residual_glb_rdphi")
        
    for chambers in ['all','long','short']:
        writeToTFile(OutF,TH2Fresidual_collector[matchingVar][chambers]['glb_phi']['TH2F'],"Residuals/MatchingOn_"+matchingVar+"/Residual_glb_phi")
        writeToTFile(OutF,TH2Fresidual_collector[matchingVar][chambers]['glb_rdphi']['TH2F'],"Residuals/MatchingOn_"+matchingVar+"/Residual_glb_rdphi")

    efficiency2DPlotAll,Num2DAll,Den2DAll,SummaryAll = generateEfficiencyPlot2DGE11(EfficiencyDictGlobal[matchingVar],[-1,1],[1,2])
    EffiDistrAll = generateEfficiencyDistribution(EfficiencyDictGlobal[matchingVar])
    GE11efficiencyByEta_Short,GE11efficiencyByEta_Long,GE11efficiencyByEta_All = generateEfficiencyPlotbyEta(EfficiencyDictGlobal[matchingVar],[1,-1],[1,2])
    GE11efficiencyByPt_Short,GE11efficiencyByPt_Long,GE11efficiencyByPt_All = generateEfficiencyPlotbyPt(EfficiencyDictGlobal[matchingVar])
    num_angle,den_angle,angle_vs_eff = incidenceAngle_vs_Eff(EfficiencyDictDirection[matchingVar],[-1,1],[1,2])

    writeToTFile(OutF,den_angle,"Efficiency/"+matchingVar+"/Angle/")
    writeToTFile(OutF,num_angle,"Efficiency/"+matchingVar+"/Angle/")
    writeToTFile(OutF,angle_vs_eff,"Efficiency/"+matchingVar+"/Angle/")

    writeToTFile(OutF,efficiency2DPlotAll,"Efficiency/"+matchingVar+"/2DView/")
    writeToTFile(OutF,Num2DAll,"Efficiency/"+matchingVar+"/2DView/")
    writeToTFile(OutF,Den2DAll,"Efficiency/"+matchingVar+"/2DView/")
    writeToTFile(OutF,SummaryAll,"Efficiency/"+matchingVar+"/")
    writeToTFile(OutF,EffiDistrAll,"Efficiency/"+matchingVar+"/")
    writeToTFile(OutF,GE11efficiencyByEta_Short,"Efficiency/"+matchingVar+"/ByEta/")
    writeToTFile(OutF,GE11efficiencyByEta_Long,"Efficiency/"+matchingVar+"/ByEta/")
    writeToTFile(OutF,GE11efficiencyByEta_All,"Efficiency/"+matchingVar+"/ByEta/")
    writeToTFile(OutF,GE11efficiencyByPt_Short,"Efficiency/"+matchingVar+"/ByPt/")
    writeToTFile(OutF,GE11efficiencyByPt_Long,"Efficiency/"+matchingVar+"/ByPt/")
    writeToTFile(OutF,GE11efficiencyByPt_All,"Efficiency/"+matchingVar+"/ByPt/")

    if DLE and matchingVar=='glb_rdphi':
        DLE_2DPlotAll,DLE_Num2DAll,DLE_Den2DAll,DLE_SummaryAll = generateEfficiencyPlot2DGE11(EfficiencyDictLayer[matchingVar],[-1,1],[1,2])
        DLE_ByEta_Short,DLE_ByEta_Long,DLE_ByEta_All = generateEfficiencyPlotbyEta(EfficiencyDictLayer[matchingVar],[1,-1],[1,2])
        writeToTFile(OutF,DLE_2DPlotAll,"Efficiency/DLE/2DView/")
        writeToTFile(OutF,DLE_Num2DAll,"Efficiency/DLE/2DView/")
        writeToTFile(OutF,DLE_Den2DAll,"Efficiency/DLE/2DView/")
        writeToTFile(OutF,DLE_SummaryAll,"Efficiency/DLE/")
        writeToTFile(OutF,DLE_ByEta_Short,"Efficiency/DLE/ByEta/")
        writeToTFile(OutF,DLE_ByEta_Long,"Efficiency/DLE/ByEta/")
        writeToTFile(OutF,DLE_ByEta_All,"Efficiency/DLE/ByEta/")
        writeToTFile(OutF,DLE_ErrR,"Efficiency/DLE/PropError/")
        writeToTFile(OutF,DLE_ErrPhi,"Efficiency/DLE/PropError/")

    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Reco'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/")
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Prop'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalLong'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalShort'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 


    for r in [-1,1]:
        for l in [1,2]:
            reg_tag_string = "P" if r == 1 else "M"
            lay_tag_string = "L1" if l == 1 else "L2"

            
            efficiency2DPlot,Num2D,Den2D,Summary = generateEfficiencyPlot2DGE11(EfficiencyDictGlobal[matchingVar],r,l)
            efficiencyByEta_Short,efficiencyByEta_Long,efficiencyByEta_All =  generateEfficiencyPlotbyEta(EfficiencyDictGlobal[matchingVar],r,l)
                       
            c1 = setUpCanvas("c1",2400,900)
            c1.SetLeftMargin(0.07)
            c1.SetRightMargin(0.09)
            c2 = setUpCanvas("c2",2400,900)
            c2.SetLeftMargin(0.07)
            c2.SetRightMargin(0.09)
            c3 = setUpCanvas("c3",2400,900)
            c3.SetLeftMargin(0.07)
            c3.SetRightMargin(0.09)
            c4 = setUpCanvas("c4",2400,900)
            c4.SetLeftMargin(0.08)
            c4.SetRightMargin(0.04)

            writeToTFile(OutF,efficiency2DPlot,"Efficiency/"+matchingVar+"/2DView/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,Num2D,"Efficiency/"+matchingVar+"/2DView/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,Den2D,"Efficiency/"+matchingVar+"/2DView/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,Summary,"Efficiency/"+matchingVar+"/")


            if DLE and matchingVar=='glb_rdphi':
                DLE_2DPlot,DLE_Num2D,DLE_Den2D,DLE_Summary = generateEfficiencyPlot2DGE11(EfficiencyDictLayer[matchingVar],r,l)
                DLE_ByEta_Short,DLE_ByEta_Long,DLE_ByEta_All =  generateEfficiencyPlotbyEta(EfficiencyDictLayer[matchingVar],r,l)
                writeToTFile(OutF,DLE_2DPlot,"Efficiency/DLE/2DView/"+reg_tag_string+lay_tag_string+"/")
                writeToTFile(OutF,DLE_Num2D,"Efficiency/DLE/2DView/"+reg_tag_string+lay_tag_string+"/")
                writeToTFile(OutF,DLE_Den2D,"Efficiency/DLE/2DView/"+reg_tag_string+lay_tag_string+"/")
                writeToTFile(OutF,DLE_Summary,"Efficiency/DLE/"+reg_tag_string+lay_tag_string+"/")
            
            
            c1.cd()
            efficiency2DPlot.Draw("COLZ TEXT45")
            c1.Modified()
            c1.Update()
            c1.SaveAs("./Plot/"+timestamp+"/"+matchingVar+"/"+efficiency2DPlot.GetTitle()+".pdf")
            
            c2.cd()
            Num2D.Draw("COLZ TEXT45")
            c2.Modified()
            c2.Update()
            c2.SaveAs("./Plot/"+timestamp+"/"+matchingVar+"/"+Num2D.GetTitle()+".pdf")
            
            c3.cd()
            Den2D.Draw("COLZ TEXT45")
            c3.Modified()
            c3.Update()
            c3.SaveAs("./Plot/"+timestamp+"/"+matchingVar+"/"+Den2D.GetTitle()+".pdf")
            c4.cd()
            Summary.Draw("APE")
            c4.Modified()
            c4.Update()
            c4.SaveAs("./Plot/"+timestamp+"/"+matchingVar+"/"+Summary.GetTitle()+".pdf")

            writeToTFile(OutF,efficiencyByEta_Short,"Efficiency/"+matchingVar+"/ByEta/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,efficiencyByEta_Long,"Efficiency/"+matchingVar+"/ByEta/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,efficiencyByEta_All,"Efficiency/"+matchingVar+"/ByEta/"+reg_tag_string+lay_tag_string+"/")

            writeToTFile(OutF,DLE_ByEta_Short,"Efficiency/DLE/ByEta/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,DLE_ByEta_Long,"Efficiency/DLE/ByEta/"+reg_tag_string+lay_tag_string+"/")
            writeToTFile(OutF,DLE_ByEta_All,"Efficiency/DLE/ByEta/"+reg_tag_string+lay_tag_string+"/")



for key in TH1MetaData.keys():
    writeToTFile(OutF,TH1MetaData[key],"Metadata/")

printSummary(EfficiencyDictGlobal,matching_variables,ResidualCutOff,matching_variable_units)


for matchingVar in ['glb_phi','glb_rdphi']:
    tempList = []
    for etaPID,subDict in EfficiencyDictGlobal[matchingVar].items():
        region,chamber,layer,eta = getInfoFromEtaID(etaPID)
        
        matchedRecHit = sum([subDict[k]['num'] for k in subDict.keys()])
        propHit = sum([subDict[k]['den'] for k in subDict.keys()])
    
        size = "S" if chamber%2 == 1 else "L"
        miplus = "M" if region == -1 else "P"
        chID = 'GE11-'+miplus+'-%02d' % chamber +"L"+str(layer)+"-"+size 
    
        tempList.append([chID,region,chamber,layer,eta,matchedRecHit,propHit])

    data = pd.DataFrame(tempList,columns=['chamberID',"region","chamber","layer","etaPartition","matchedRecHit","propHit"])
    data.to_csv('./Plot/'+timestamp+'/MatchingSummary_'+matchingVar+'.csv', index=False)
