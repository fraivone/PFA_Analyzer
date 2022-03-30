import ROOT
import subprocess
import sys
import numpy as np
import time
import argparse
import pandas as pd
from argparse import RawTextHelpFormatter
from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *


parser = argparse.ArgumentParser(
        description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Plot Sanity Checks\n\t-Plot Residuals (takes the cut as parameter)\n\t-Plot efficiency\nCurrently allows the hits matching on glb_phi and glb_rdphi''',
        epilog="""Typical exectuion\n\t python PFA_Analyzer.py  --phi_cut 0.001 --rdphi_cut 0.15""",
        formatter_class=RawTextHelpFormatter
)

parser.add_argument('-pc','--phi_cut', type=float,help="Maximum allowed dphi between RecoHit and PropHit to be counted as matched hit",required=False)
parser.add_argument('-rdpc','--rdphi_cut', type=float,help="Maximum allowed rdphi between RecoHit and PropHit to be counted as matched hit",required=False)
parser.add_argument('--chi2cut', type=float,help="Maximum normalized chi2 for which accept propagated tracks",required=False)
parser.add_argument('--minPt', type=float,help="Minimum pt for which accept propagated tracks",required=False)
parser.add_argument('--chamberOFF', type=str , help="file_path to the file containing a list of the chambers you want to exclude for the run (i.e. GE11-M-29L2)",required=False,nargs='*')
parser.add_argument('--VFATOFF', type=str , help="file_path to the file containing a list of the VFAT you want to exclude from efficiency evaluation. The file must be tab separated with \tregion,layer,chamber,VFAT,reason_mask",required=False,nargs='*')
parser.add_argument('--outputname', type=str, help="output file name",required=False)
parser.add_argument('--fiducialR','-fR', type=float , help="fiducial cut along R axis",required=False)
parser.add_argument('--fiducialPhi','-fP', type=float , help="fiducial cut along phi axis",required=False)
parser.add_argument('--maxErrPropR', type=float , help="max error on propagated R in order to accept the muon (cm)",required=False)
parser.add_argument('--maxErrPropPhi', type=float , help="max error on propagated phi in order to accept the muon (rad)",required=False)
parser.add_argument('--DLE', default=False, action='store_true',help="Swtiches on the Double Layer Efficiency (DLE) analisys. False by default",required=False)
parser.add_argument('--FD', default=False, action='store_true',help="When enabled, allows the storage of all GEM RecHit Digis. False by default, which means that GEM RecHit Digis are stored only for EVTs in which STA propagation hits GEM",required=False)
parser.add_argument('--verbose', default=False, action='store_true',help="Verbose printing",required=False)
parser.add_argument('--dataset','-ds', type=str,help="TAG to the folder containing the NTuples to be analyzed",required=True,nargs='*')

parser.add_argument('--minME1', type=int, help="Min number of ME1 hits",required=False)
parser.add_argument('--minME2', type=int, help="Min number of ME2 hits",required=False)
parser.add_argument('--minME3', type=int, help="Min number of ME3 hits",required=False)
parser.add_argument('--minME4', type=int, help="Min number of ME4 hits",required=False)
parser.set_defaults(phi_cut=0.001)
parser.set_defaults(rdphi_cut=0.15)
parser.set_defaults(chi2cut=9999999)
parser.set_defaults(minPt=0)
parser.set_defaults(minME1=0)
parser.set_defaults(minME2=0)
parser.set_defaults(minME3=0)
parser.set_defaults(minME4=0)
parser.set_defaults(fiducialR=1)
parser.set_defaults(fiducialPhi=0.005)
parser.set_defaults(maxErrPropR=1)
parser.set_defaults(maxErrPropPhi=0.005)
parser.set_defaults(outputname=time.strftime("%-y%m%d_%H%M"))
args = parser.parse_args()

chamberForEventDisplay = ["GE11-P-28L2-L"]


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetLineScalePS(1)
if not args.verbose: ROOT.gROOT.ProcessLine("gErrorIgnoreLevel=2001;") #suppressed everything less-than-or-equal-to kWarning



start_time = time.time()

files = []


for folder in args.dataset:
    files = files + files_in_folder(folder)

matching_variables = ['glb_phi','glb_rdphi']
matching_variable_units = {'glb_phi':'rad','glb_rdphi':'cm'}
ResidualCutOff= {'glb_phi':args.phi_cut,'glb_rdphi':args.rdphi_cut}

DLE = args.DLE
FD = args.FD
fiducialCut = True
maxErrOnPropR = args.maxErrPropR
maxErrOnPropPhi = args.maxErrPropPhi
fiducialR = args.fiducialR
fiducialPhi = args.fiducialPhi
CutminPt = args.minPt
maxSTA_NormChi2 = args.chi2cut
minME1Hit = args.minME1
minME2Hit = args.minME2
minME3Hit = args.minME3
minME4Hit = args.minME4
outputname = args.outputname

noisyEtaPID = []
VFATOFFDict = {} if args.VFATOFF is None else importOFFVFAT(args.VFATOFF)
chamberOFFLS = {} if args.chamberOFF is None else ChamberOFF_byLS(args.chamberOFF)

TH1nbins = 120
TH2nbins = 200
TH2min = -80

EfficiencyDictGlobal = dict((m,{}) for m in matching_variables)
EfficiencyDictVFAT = generateVFATDict(matching_variables)
EfficiencyDictDirection = dict((m,{}) for m in matching_variables)
EfficiencyDictLayer = dict((m,{}) for m in matching_variables)

## ROOT Object declaration
STAdirX_vs_CLS = ROOT.TH2F("DirX_CLS","DirX_CLS",40,-1,1,128,0,128)
STAdirX_vs_CLS.GetXaxis().SetTitle("dirX")
STAdirX_vs_CLS.GetYaxis().SetTitle("CLS")

DLE_ErrPhi = ROOT.TH1F("DLE_ErrPhi","DLE_ErrPhi",100,0,0.0025)
DLE_ErrR = ROOT.TH1F("DLE_ErrR","DLE_ErrR",100,0,5)
DLE_pt = ROOT.TH1F("pT of STA muons used to probe DLE","pT of STA muons used to probe DLE",200,0,100)

TH1Fresidual_collector = generate1DResidualContainer(matching_variables,TH1nbins,ResidualCutOff)
TH1FpropError_collector = generatePropagationErrorContainer(maxErrOnPropR, maxErrOnPropPhi)
TH2Fresidual_collector = generate2DResidualContainer(matching_variables,TH2nbins,TH2min)  
THSanityChecks = {'Occupancy':{}, 
                  'NHits':{},
                  'PropagationError':{},
                  'etaP_vs_pt':[],
                  'Residual_Correlation':{},
                  'PropHit_DirLoc_xOnGE11':{'BeforeMatching':{'Long':{},'Short':{}},
                                            'AfterMatching':{'Long':{},'Short':{}}
                                            },
                  'RecHitperStrip':{},
                  'NEvts':ROOT.TH1F("NumberOfAnalyzedEVTs","NumberOfAnalyzedEVTs",100,1,1),
                  'CLS':{'glb_phi':{},'glb_rdphi':{}},
                  'pt':{'glb_phi':{},'glb_rdphi':{},'All':ROOT.TH1F("pT of STA muons used to probe efficiency","pT of STA muons used to probe efficiency",200,0,100)}
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
TH1MetaData['chamberOFFCanvas']=setUpCanvas("ChamberOFF_byLS",1200,1200)
TH1MetaData['VFATOFFCanvas']=setUpCanvas("ExcludedVFAT",1200,1200)
TH1MetaData['ExclusionSummaryCanvas']=setUpCanvas("ChamberMask_Summary",1200,1200)



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
for t_re in [-1,1]:
    for t_la in [1,2]:
        for t_ch in range(1,37):
            ch_id = ReChLa2chamberName(t_re,t_ch,t_la)
            THSanityChecks['CLS']['glb_rdphi'].setdefault(ch_id,{})
            THSanityChecks['CLS']['glb_phi'].setdefault(ch_id,{})
            THSanityChecks['pt'].setdefault(ch_id,{})
            for t_eta in range(1,9)+["All"]:
                THSanityChecks['CLS']['glb_rdphi'][ch_id].setdefault(t_eta,ROOT.TH1F(ch_id+"_etaP"+str(t_eta)+" ClusterSize rdphi",ch_id+"_etaP"+str(t_eta)+" ClusterSize rdphi",20,0,20))
                THSanityChecks['CLS']['glb_phi'][ch_id].setdefault(t_eta,ROOT.TH1F(ch_id+"_etaP"+str(t_eta)+" ClusterSize phi",ch_id+"_etaP"+str(t_eta)+" ClusterSize phi",20,0,20))
                THSanityChecks['pt'][ch_id].setdefault(t_eta,ROOT.TH1F(ch_id+"_etaP"+str(t_eta)+" pT",ch_id+"_etaP"+str(t_eta)+" pT",200,0,100))

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
        THSanityChecks['RecHitperStrip'][key][ch].SetMaximum(600)
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

for key_1 in ['Long','Short']:
    for key_2 in ["eta"+str(k) for k in range(1,9)]:
        THSanityChecks['PropHit_DirLoc_xOnGE11']['BeforeMatching'][key_1][key_2] = ROOT.TH1F('BeforeMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,'BeforeMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,200,-1,1)
        THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching'][key_1][key_2] = ROOT.TH1F('AfterMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,'AfterMatchPropHit_DirLoc_xOnGE11_'+key_1+key_2,200,-1,1)

# matchingFile = ROOT.TFile("./Output/PFA_Analyzer_Output/ROOT_File/MatchingTree_"+outputname+".root","RECREATE")
# tree_out,generalMatching,recHit_Matching,propHit_Matching = initializeMatchingTree()

## Chain files
chain = ROOT.TChain("muNtupleProducer/MuDPGTree")

print args.dataset, "TChaining ",len(files)," files..."
print'\n'.join(files)
print
for fl in files:
    chain.Add(fl)
## Enabling only the branches which are actually in use
# 1. Disabling them all
chain.SetBranchStatus("*",0);     


branchList=["event_eventNumber","event_lumiBlock","event_runNumber","gemRecHit_region", "gemRecHit_chamber", "gemRecHit_layer", "gemRecHit_etaPartition", "gemRecHit_g_r", "gemRecHit_loc_x", "gemRecHit_g_x", "gemRecHit_g_y", "gemRecHit_g_z", "gemRecHit_g_phi", "gemRecHit_firstClusterStrip", "gemRecHit_cluster_size", "mu_propagated_region", "mu_propagated_chamber", "mu_propagated_layer", "mu_propagated_etaP", "mu_propagated_Outermost_z",  "mu_propagated_isME11", "mu_propagatedGlb_r", "mu_propagatedLoc_x", "mu_propagatedLoc_y", "mu_propagatedGlb_x", "mu_propagatedGlb_y", "mu_propagatedGlb_z", "mu_propagatedGlb_phi", "mu_propagatedGlb_errR", "mu_propagatedGlb_errPhi", "mu_propagatedLoc_dirX", "mu_propagatedLoc_dirY", "mu_propagatedLoc_dirZ", "mu_propagated_pt", "mu_propagated_isGEM", "mu_propagated_TrackNormChi2", "mu_propagated_nME1hits", "mu_propagated_nME2hits", "mu_propagated_nME3hits", "mu_propagated_nME4hits"]

# 2. Enabling the useful ones
for b in branchList:
    chain.SetBranchStatus(b,1)


chainEntries = chain.GetEntries()
maxLS = 0
print "\n#############\nStarting\n#############"

try:
    print "Analysing run(s): \t", [int(GetRunNumber(i)) for i in args.dataset]
except:
    pass

print "Number of evts \t\t%.2fM\n" %(round(float(chainEntries)/10**6,2))

THSanityChecks['NEvts'].Fill(chainEntries)
for chain_index,evt in enumerate(chain):
    
    if chain_index % 40000 ==0:
        print "[",time.strftime("%B %d - %H:%M:%S"),"]\t",round(float(chain_index)/float(chainEntries),3)*100,"%"
        
    n_gemprop = len(evt.mu_propagated_chamber)
    n_gemrec = len(evt.gemRecHit_chamber)
    
    EventNumber = evt.event_eventNumber
    LumiSection = evt.event_lumiBlock
    maxLS = max(maxLS,LumiSection)
    RunNumber = evt.event_runNumber
    
    THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Prop'].Fill(n_gemprop)
    THSanityChecks['NHits']['BeforeMatching']['PerEVT']['Reco'].Fill(n_gemrec)
    

    # If FullDigis == True, never skip evts
    # If FullDigis == False, skip evts with 0 propagations
    if FD == False and n_gemprop==0:
        continue


    ML1_NGEMRecoHits = 0
    ML2_NGEMRecoHits = 0
    PL1_NGEMRecoHits = 0
    PL2_NGEMRecoHits = 0    
    
    RecHit_Dict = {}
    PropHit_Dict = {}

    for RecHit_index in range(0,n_gemrec):
        region = evt.gemRecHit_region[RecHit_index]
        chamber = evt.gemRecHit_chamber[RecHit_index]
        layer = evt.gemRecHit_layer[RecHit_index]
        etaP = evt.gemRecHit_etaPartition[RecHit_index]
        RecHitEtaPartitionID =  region*(100*chamber+10*layer+etaP)
        endcapKey = EndcapLayer2label(region,layer)
        chamberID = ReChLa2chamberName(region,chamber,layer)

        ## discard chambers that were kept OFF from the analysis
        if chamberID in chamberOFFLS.keys() and (LumiSection in chamberOFFLS[chamberID] or -1 in chamberOFFLS[chamberID] ):
            continue

        rec_glb_r = evt.gemRecHit_g_r[RecHit_index]
        rec_loc_x = evt.gemRecHit_loc_x[RecHit_index]

        if RecHitEtaPartitionID in VFATOFFDict:
            if propHit2VFAT(rec_glb_r,rec_loc_x,etaP,region,chamber) in VFATOFFDict[RecHitEtaPartitionID]:
                continue

        if region == 1:
            if layer == 1:
                PL1_NGEMRecoHits += 1 
            else:
                PL2_NGEMRecoHits += 1 
        if region == -1:
            if layer == 1:
                ML1_NGEMRecoHits += 1 
            else:
                ML2_NGEMRecoHits += 1 

        RecHit_Dict.setdefault(RecHitEtaPartitionID, {'loc_x':[],'glb_x':[],'glb_y':[],'glb_z':[],'glb_r':[],'glb_phi':[],'firstStrip':[],'cluster_size':[]})
        RecHit_Dict[RecHitEtaPartitionID]['loc_x'].append(rec_loc_x)
        RecHit_Dict[RecHitEtaPartitionID]['glb_x'].append(evt.gemRecHit_g_x[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_y'].append(evt.gemRecHit_g_y[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_z'].append(evt.gemRecHit_g_z[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_r'].append(rec_glb_r)
        RecHit_Dict[RecHitEtaPartitionID]['glb_phi'].append(evt.gemRecHit_g_phi[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['firstStrip'].append(evt.gemRecHit_firstClusterStrip[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['cluster_size'].append(evt.gemRecHit_cluster_size[RecHit_index])

        THSanityChecks['Occupancy']['BeforeMatching']['Reco'].Fill(evt.gemRecHit_g_x[RecHit_index],evt.gemRecHit_g_y[RecHit_index])
        THSanityChecks['Occupancy']['BeforeMatching'][endcapKey]['RecHits'].Fill(chamber,etaP)
        
        for j in range(0,RecHit_Dict[RecHitEtaPartitionID]['cluster_size'][-1]):
            strip = RecHit_Dict[RecHitEtaPartitionID]['firstStrip'][-1] + j
            THSanityChecks['RecHitperStrip'][endcapKey][chamber].Fill(strip,etaP)
                    
    THSanityChecks['NHits']['BeforeMatching']['ML1'].Fill(ML1_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['ML2'].Fill(ML2_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['PL1'].Fill(PL1_NGEMRecoHits)
    THSanityChecks['NHits']['BeforeMatching']['PL2'].Fill(PL2_NGEMRecoHits)    
    
    for PropHit_index in range(0,n_gemprop):
        
        region = evt.mu_propagated_region[PropHit_index]
        chamber = evt.mu_propagated_chamber[PropHit_index]
        layer = evt.mu_propagated_layer[PropHit_index]
        etaP = evt.mu_propagated_etaP[PropHit_index]
        if etaP >= 9: ## from CMSSW_12_2_1 GE21 demonstrator is also included in the propagated chambers
            continue

        PropHitChamberID = region*(100*chamber+10*layer+etaP)
        endcapKey = EndcapLayer2label(region,layer)

        

        outermost_z = evt.mu_propagated_Outermost_z[PropHit_index]
        # is_incoming = evt.mu_propagated_isincoming[PropHit_index]

        if region == 1 and outermost_z < 0:
            continue
        if region == -1 and outermost_z > 0:
            continue

        chamberID = ReChLa2chamberName(region,chamber,layer)
        ## discard chambers that were kept OFF from the analysis
        if chamberID in chamberOFFLS.keys() and (LumiSection in chamberOFFLS[chamberID] or -1 in chamberOFFLS[chamberID] ):
            continue

        propHitFromME11 = bool(evt.mu_propagated_isME11[PropHit_index])
        if propHitFromME11:            
            PropHit_Dict.setdefault(PropHitChamberID,{'loc_x':[],'loc_y':[],'glb_x':[],'glb_y':[],'glb_z':[],'glb_r':[],'glb_phi':[],'pt':[],'etaP':[],'err_glb_r':[],'err_glb_phi':[],'Loc_dirX':[],'Loc_dirY':[],'Loc_dirZ':[],'mu_propagated_isME11':[],'mu_propagated_isGEM':[],'STA_Normchi2':[],'nME1Hits':[],'nME2Hits':[],'nME3Hits':[],'nME4Hits':[]})
            prop_glb_r = evt.mu_propagatedGlb_r[PropHit_index]
            prop_loc_x = evt.mu_propagatedLoc_x[PropHit_index]

            if PropHitChamberID in VFATOFFDict:
                if propHit2VFAT(prop_glb_r,prop_loc_x,etaP,region,chamber) in VFATOFFDict[PropHitChamberID]:
                    continue

            
            PropHit_Dict[PropHitChamberID]['loc_x'].append(prop_loc_x)
            PropHit_Dict[PropHitChamberID]['loc_y'].append(evt.mu_propagatedLoc_y[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_x'].append(evt.mu_propagatedGlb_x[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_y'].append(evt.mu_propagatedGlb_y[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_z'].append(evt.mu_propagatedGlb_z[PropHit_index])
            PropHit_Dict[PropHitChamberID]['glb_r'].append(prop_glb_r)
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
        endcapTag = EndcapLayer2label(region,layer)
        current_chamber_ID = ReChLa2chamberName(region,chamber,layer)
        
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
            for j in range(0,20):
                EfficiencyDictDirection[mv][etaPartitionID].setdefault(j,{'num':0,'den':0})

        isGoodTrack = []
        passedCutProp = {key:[] for key in PropHitonEta.keys()}
        ## Applying cuts on the propagated tracks to be used
        for index in range(nPropHitsOnEtaID):
            if fiducialCut and passCut(PropHitonEta,etaPartitionID,index,maxPropR_Err=maxErrOnPropR,maxPropPhi_Err=maxErrOnPropPhi,fiducialCutR=fiducialR,fiducialCutPhi=fiducialPhi,minPt=CutminPt,maxChi2=maxSTA_NormChi2,minME1Hit=minME1Hit,minME2Hit=minME2Hit,minME3Hit=minME3Hit,minME4Hit=minME4Hit) == False:
                isGoodTrack.append(False)
            else:
                EfficiencyDictGlobal['glb_phi'][etaPartitionID][pt_index(PropHitonEta['pt'][index])]['den'] += 1
                EfficiencyDictGlobal['glb_rdphi'][etaPartitionID][pt_index(PropHitonEta['pt'][index])]['den'] += 1
                
                VFAT_propagated = propHit2VFAT(PropHitonEta['glb_r'][index],PropHitonEta['loc_x'][index],eta,region,chamber)
                EfficiencyDictVFAT['glb_phi'][endcapTag][current_chamber_ID][VFAT_propagated]['den'] += 1
                EfficiencyDictVFAT['glb_rdphi'][endcapTag][current_chamber_ID][VFAT_propagated]['den'] += 1

                angle_index = int ( (cos_of_alpha_list[index] * 20 ) )
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
           
            TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID]["All"]["ErrPhi"].Fill(PropHitonEta['err_glb_phi'][k])
            TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID][eta]["ErrPhi"].Fill(PropHitonEta['err_glb_phi'][k])
            TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID]["All"]["ErrR"].Fill(PropHitonEta['err_glb_r'][k])
            TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID][eta]["ErrR"].Fill(PropHitonEta['err_glb_r'][k])

            THSanityChecks['pt'][current_chamber_ID]["All"].Fill(PropHitonEta['pt'][k])
            THSanityChecks['pt'][current_chamber_ID][eta].Fill(PropHitonEta['pt'][k])
            THSanityChecks['pt']['All'].Fill(PropHitonEta['pt'][k])
            
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
            # if current_chamber_ID in chamberForEventDisplay: store4evtDspl(outputname,RunNumber,LumiSection,EventNumber)
            # if current_chamber_ID in chamberForEventDisplay: 
            #     print "event",EventNumber
            #     print "\t\tNothing to match on ", eta,"\t Match looked for phi = ",PropHitonEta['glb_phi']
            #     for s in range(1,9):
            #         sf =  region*(100*chamber+10*layer+s)
            #         if sf in RecHit_Dict.keys(): print sf,"\t",RecHit_Dict[sf]
                
            # print PropHit_Dict[etaPartitionID]['glb_phi']
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
            propagation_error_phi = PropHitonEta['err_glb_phi'][prop_hit_index]
            propagation_error_r = PropHitonEta['err_glb_r'][prop_hit_index]

            if matchingVar == 'glb_phi':
                THSanityChecks['Residual_Correlation']['glb_phi_vs_glb_rdphi'].Fill(glb_phi_residual,glb_rdphi_residual)
                THSanityChecks['Residual_Correlation']['glb_rdphi_dir_x'].Fill(glb_rdphi_residual,np.arccos(PropHitonEta['Loc_dirX'][prop_hit_index]))


            if abs(min_residual) < ResidualCutOff[matchingVar]:


                # generalMatching['event'][0] = EventNumber
                # generalMatching['lumiblock'][0] = LumiSection
                # generalMatching['chamber'][0] = chamber
                # generalMatching['region'][0] = region
                # generalMatching['layer'][0] = layer
                # generalMatching['etaPartition'][0] = etaP
                # recHit_Matching,propHit_Matching = fillMatchingTreeArray(PropHitonEta,prop_hit_index,RecHitonEta,reco_hit_index,propHit_Matching,recHit_Matching)
                # tree_out.Fill()

                EfficiencyDictGlobal[matchingVar][etaPartitionID][pt_index(PropHitonEta['pt'][prop_hit_index])]['num'] += 1
                
                VFAT_propagated = propHit2VFAT(PropHitonEta['glb_r'][prop_hit_index],PropHitonEta['loc_x'][prop_hit_index],eta,region,chamber)
                EfficiencyDictVFAT[matchingVar][endcapTag][current_chamber_ID][VFAT_propagated]['num'] += 1

                angle_index = int( np.sqrt(PropHitonEta['Loc_dirX'][prop_hit_index]**2 + PropHitonEta['Loc_dirY'][prop_hit_index]**2) * 20)
                EfficiencyDictDirection[matchingVar][etaPartitionID][angle_index]['num'] += 1


                if matchingVar == "glb_rdphi":
                    TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID]["All"]["ErrPhi"].Fill(propagation_error_phi)
                    TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID][eta]["ErrPhi"].Fill(propagation_error_phi)
                    TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID]["All"]["ErrR"].Fill(propagation_error_r)
                    TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID][eta]["ErrR"].Fill(propagation_error_r)
                    TH1Fresidual_collector[matchingVar][endcapTag][current_chamber_ID]["All"]["Residual"].Fill(glb_rdphi_residual)
                    TH1Fresidual_collector[matchingVar][endcapTag][current_chamber_ID][eta]["Residual"].Fill(glb_rdphi_residual)
                else:
                    TH1Fresidual_collector[matchingVar][endcapTag][current_chamber_ID]["All"]["Residual"].Fill(glb_phi_residual)
                    TH1Fresidual_collector[matchingVar][endcapTag][current_chamber_ID][eta]["Residual"].Fill(glb_phi_residual)

                binx = int(round((PropHitonEta['loc_x'][prop_hit_index]-TH2min)*(TH2nbins-1)/(-2*TH2min)))+1
                biny = int(round((PropHitonEta['loc_y'][prop_hit_index]-TH2min)*(TH2nbins-1)/(-2*TH2min)))+1
                TH2Fresidual_collector[matchingVar]['all']['glb_phi'][binx][biny][0] += 1
                TH2Fresidual_collector[matchingVar]['all']['glb_rdphi'][binx][biny][0] += 1
                TH2Fresidual_collector[matchingVar]['all']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                TH2Fresidual_collector[matchingVar]['all']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)

                THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Reco'].Fill(RecHitonEta['glb_x'][reco_hit_index],RecHitonEta['glb_y'][reco_hit_index])
                THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Prop'].Fill(PropHitonEta['glb_x'][prop_hit_index],PropHitonEta['glb_y'][prop_hit_index])

                THSanityChecks['CLS'][matchingVar][current_chamber_ID]["All"].Fill(RecHitonEta['cluster_size'][reco_hit_index])
                THSanityChecks['CLS'][matchingVar][current_chamber_ID][eta].Fill(RecHitonEta['cluster_size'][reco_hit_index])

                if chamber%2 == 0:
                    if matchingVar == 'glb_phi': THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching']['Long']['eta'+str(eta)].Fill(PropHitonEta['Loc_dirX'][prop_hit_index])
                    TH2Fresidual_collector[matchingVar]['long']['glb_phi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['long']['glb_rdphi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['long']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                    TH2Fresidual_collector[matchingVar]['long']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)
                    THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalLong'].Fill(PropHitonEta['loc_x'][prop_hit_index],PropHitonEta['loc_y'][prop_hit_index])
                if chamber%2 == 1:
                    if matchingVar == 'glb_phi': THSanityChecks['PropHit_DirLoc_xOnGE11']['AfterMatching']['Short']['eta'+str(eta)].Fill(PropHitonEta['Loc_dirX'][prop_hit_index])
                    TH2Fresidual_collector[matchingVar]['short']['glb_phi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['short']['glb_rdphi'][binx][biny][0] += 1
                    TH2Fresidual_collector[matchingVar]['short']['glb_phi'][binx][biny][1] += abs(glb_phi_residual)
                    TH2Fresidual_collector[matchingVar]['short']['glb_rdphi'][binx][biny][1] += abs(glb_rdphi_residual)                    
                    THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalShort'].Fill(PropHitonEta['loc_x'][prop_hit_index],PropHitonEta['loc_y'][prop_hit_index])

                
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
                    STAdirX_vs_CLS.Fill(PropHitonEta['Loc_dirX'][prop_hit_index],RecHitonEta['cluster_size'][reco_hit_index])
                    
                    if layer == 1:
                        layer1Match = True
                        layer1pt = PropHitonEta['pt'][prop_hit_index]
                    if layer == 2:
                        layer2Match = True
                        layer2pt = PropHitonEta['pt'][prop_hit_index]
                    
            else:
                # if current_chamber_ID in chamberForEventDisplay: store4evtDspl(outputname,RunNumber,LumiSection,EventNumber)
                # print "Matching failed for ", etaPartitionID
                # print PropHit_Dict[etaPartitionID]['glb_phi'], RecHit_Dict[etaPartitionID]['glb_phi']
                # raw_input()
                pass
                
        ## Loop over etaPID 

    ## Double Layer Efficiency (DLE): test layer1(2) with tracks that have matched in layer1(2)

    if DLE and layer2Match and layer1PassedCut and layer2PassedCut:
        EfficiencyDictLayer['glb_rdphi'][layer1etapID][pt_index(layer2pt)]['den'] += 1
        DLE_pt.Fill(layer2pt)
        # DLE_ErrPhi.Fill(PropHitonEta['err_glb_phi'][prop_hit_index])
        # DLE_ErrR.Fill(PropHitonEta['err_glb_r'][prop_hit_index])

        if layer1Match == True:
            EfficiencyDictLayer['glb_rdphi'][layer1etapID][pt_index(layer2pt)]['num'] += 1

    if DLE and layer1Match and layer2PassedCut and layer1PassedCut:
        EfficiencyDictLayer['glb_rdphi'][layer2etapID][pt_index(layer1pt)]['den'] += 1
        DLE_pt.Fill(layer1pt)
        # DLE_ErrPhi.Fill(PropHitonEta['err_glb_phi'][prop_hit_index])
        # DLE_ErrR.Fill(PropHitonEta['err_glb_r'][prop_hit_index])

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
# matchingFile.Write()
OutF = ROOT.TFile("./Output/PFA_Analyzer_Output/ROOT_File/"+outputname+".root","RECREATE")

subprocess.call(["mkdir", "-p", "./Output/PFA_Analyzer_Output/CSV/"+outputname])
subprocess.call(["mkdir", "-p", "./Output/PFA_Analyzer_Output/CSV/"+outputname])
subprocess.call(["mkdir", "-p", "./Output/PFA_Analyzer_Output/Plot/"+outputname+"/"+"glb_phi/"])
subprocess.call(["mkdir", "-p", "./Output/PFA_Analyzer_Output/Plot/"+outputname+"/"+"glb_rdphi/"])


### Masking
TH1MetaData['chamberOFFCanvas'].Divide(2,2)
TH1MetaData['VFATOFFCanvas'].Divide(2,2)
TH1MetaData['ExclusionSummaryCanvas'].Divide(2,2)

OFFChambers_plots = ChambersOFFHisto(chamberOFFLS,maxLS)
for counter,plot in enumerate(OFFChambers_plots):
    TH1MetaData['chamberOFFCanvas'].cd(counter+1)
    plot.Draw()

OFFVFATs_plots = VFATOFFHisto(VFATOFFDict)
for counter,plot in enumerate(OFFVFATs_plots):
    TH1MetaData['VFATOFFCanvas'].cd(counter+1)
    plot.Draw("COLZ")

GE11Discarded_plots = GE11DiscardedSummary(chamberOFFLS,maxLS,VFATOFFDict)
for counter,plot in enumerate(GE11Discarded_plots):
    TH1MetaData['ExclusionSummaryCanvas'].cd(counter+1)
    plot.Draw("COLZ")
    TH1MetaData['ExclusionSummaryCanvas'].cd(counter+1).Update()
    palette = plot.GetListOfFunctions().FindObject("palette")
    palette.SetX2NDC(0.93)
    palette.Draw()
setUpCanvas("GE11 Masked",1200,1200).cd()

TH1MetaData['chamberOFFCanvas'].Modified()
TH1MetaData['chamberOFFCanvas'].Update()
TH1MetaData['VFATOFFCanvas'].Modified()
TH1MetaData['VFATOFFCanvas'].Update()
TH1MetaData['ExclusionSummaryCanvas'].Modified()
TH1MetaData['ExclusionSummaryCanvas'].Update()
TH1MetaData['ExclusionSummaryCanvas'].SaveAs("./Output/PFA_Analyzer_Output/Plot/"+outputname+"/ExclusionCanvas.pdf")

writeToTFile(OutF,THSanityChecks['NEvts'],"SanityChecks/NumberOfEVTs/")

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

writeToTFile(OutF,STAdirX_vs_CLS,"SanityChecks/CLS_vs_Direction/")

writeToTFile(OutF,THSanityChecks['pt']['All'],"SanityChecks/pt/")
if DLE: writeToTFile(OutF,DLE_pt,"SanityChecks/pt/")

for key in ['ML1','ML2','PL1','PL2']:
    writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching'][key]['PropHits'],"SanityChecks/Occupancy/BeforeMatching/"+key)
    writeToTFile(OutF,THSanityChecks['Occupancy']['BeforeMatching'][key]['RecHits'],"SanityChecks/Occupancy/BeforeMatching/"+key)
    for ch in range(1,37):
        writeToTFile(OutF,THSanityChecks['RecHitperStrip'][key][ch],"SanityChecks/Occupancy/RecHitByStrip/"+key)

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
    for chambers in ['all','long','short']:
        writeToTFile(OutF,TH2Fresidual_collector[matchingVar][chambers]['glb_phi']['TH2F'],"Residuals/MatchingOn_"+matchingVar+"/2D_glb_phi")
        writeToTFile(OutF,TH2Fresidual_collector[matchingVar][chambers]['glb_rdphi']['TH2F'],"Residuals/MatchingOn_"+matchingVar+"/2D_glb_rdphi")

    efficiency2DPlotAll,Num2DAll,Den2DAll,SummaryAll = generateEfficiencyPlot2DGE11(EfficiencyDictGlobal[matchingVar],[-1,1],[1,2],debug=args.verbose)
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

    for r in [-1,1]:
        for l in [1,2]:
            endcapTag = EndcapLayer2label(r,l)

            efficiency2DPlot,Num2D,Den2D,Summary = generateEfficiencyPlot2DGE11(EfficiencyDictGlobal[matchingVar],r,l)
            efficiencyByEta_Short,efficiencyByEta_Long,efficiencyByEta_All =  generateEfficiencyPlotbyEta(EfficiencyDictGlobal[matchingVar],r,l)
                       
            c1,c2 = setUpCanvas(matchingVar+"_"+endcapTag,2400,900),setUpCanvas(matchingVar+"_VFAT2D_"+endcapTag,2000,2000)
            c1.SetLeftMargin(0.07)
            c1.SetRightMargin(0.09)                     
            c2.SetLeftMargin(0.1)
            c2.SetRightMargin(0.1)  
            c1.cd()                   

            writeToTFile(OutF,efficiency2DPlot,"Efficiency/"+matchingVar+"/2DView/"+endcapTag+"/")
            writeToTFile(OutF,Num2D,"Efficiency/"+matchingVar+"/2DView/"+endcapTag+"/")
            writeToTFile(OutF,Den2D,"Efficiency/"+matchingVar+"/2DView/"+endcapTag+"/")
            writeToTFile(OutF,generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT[matchingVar],endcapTag),"Efficiency/"+matchingVar+"/ByVFAT/"+endcapTag)
            writeToTFile(OutF,Summary,"Efficiency/"+matchingVar+"/")


            if DLE and matchingVar=='glb_rdphi':
                DLE_2DPlot,DLE_Num2D,DLE_Den2D,DLE_Summary = generateEfficiencyPlot2DGE11(EfficiencyDictLayer[matchingVar],r,l)
                DLE_ByEta_Short,DLE_ByEta_Long,DLE_ByEta_All =  generateEfficiencyPlotbyEta(EfficiencyDictLayer[matchingVar],r,l)
                writeToTFile(OutF,DLE_2DPlot,"Efficiency/DLE/2DView/"+endcapTag+"/")
                writeToTFile(OutF,DLE_Num2D,"Efficiency/DLE/2DView/"+endcapTag+"/")
                writeToTFile(OutF,DLE_Den2D,"Efficiency/DLE/2DView/"+endcapTag+"/")
                writeToTFile(OutF,DLE_Summary,"Efficiency/DLE/"+endcapTag+"/")
                writeToTFile(OutF,DLE_ByEta_Short,"Efficiency/DLE/ByEta/"+endcapTag+"/")
                writeToTFile(OutF,DLE_ByEta_Long,"Efficiency/DLE/ByEta/"+endcapTag+"/")
                writeToTFile(OutF,DLE_ByEta_All,"Efficiency/DLE/ByEta/"+endcapTag+"/")
            
            ## Save plots
            for plot_obj in [efficiency2DPlot,Num2D,Den2D]:
                plot_obj.Draw("COLZ TEXT45")
                c1.Modified()
                c1.Update()
                c1.SaveAs("./Output/PFA_Analyzer_Output/Plot/"+outputname+"/"+matchingVar+"/"+plot_obj.GetTitle()+".pdf")
            
            Summary.Draw("APE")
            c1.Modified()
            c1.Update()
            c1.SaveAs("./Output/PFA_Analyzer_Output/Plot/"+outputname+"/"+matchingVar+"/"+Summary.GetTitle()+".pdf")

            c2.cd()
            endcapVFAT2D = generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT[matchingVar],endcapTag)
            endcapVFAT2D.Draw("COLZ TEXT")
            c2.Modified()
            c2.Update()
            c2.SaveAs("./Output/PFA_Analyzer_Output/Plot/"+outputname+"/"+matchingVar+"/"+endcapVFAT2D.GetTitle()+".pdf")

            writeToTFile(OutF,efficiencyByEta_Short,"Efficiency/"+matchingVar+"/ByEta/"+endcapTag+"/")
            writeToTFile(OutF,efficiencyByEta_Long,"Efficiency/"+matchingVar+"/ByEta/"+endcapTag+"/")
            writeToTFile(OutF,efficiencyByEta_All,"Efficiency/"+matchingVar+"/ByEta/"+endcapTag+"/")

            for t_ch in range(1,37):
                current_chamber_ID = ReChLa2chamberName(r,t_ch,l)
                writeToTFile(OutF,generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT[matchingVar],current_chamber_ID),"Efficiency/"+matchingVar+"/ByVFAT/"+endcapTag+"/"+current_chamber_ID)
                writeToTFile(OutF,generate1DEfficiencyPlotbyVFAT(EfficiencyDictVFAT[matchingVar],current_chamber_ID),"Efficiency/"+matchingVar+"/ByVFAT/"+endcapTag+"/"+current_chamber_ID)
                for t_eta in range(1,9) + ["All"]:
                    writeToTFile(OutF,THSanityChecks['CLS'][matchingVar][current_chamber_ID][t_eta],"SanityChecks/CLS/"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID)
                    writeToTFile(OutF,THSanityChecks['pt'][current_chamber_ID][t_eta],"SanityChecks/pt/"+endcapTag+"/"+current_chamber_ID)
                    writeToTFile(OutF,TH1Fresidual_collector[matchingVar][endcapTag][current_chamber_ID][t_eta]["Residual"],"Residuals/MatchingOn_"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID+"/Residual")
                    writeToTFile(OutF,TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID][t_eta]["ErrPhi"],"Residuals/MatchingOn_"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID+"/PropagationError/AllHits")
                    writeToTFile(OutF,TH1FpropError_collector["AllHits"][endcapTag][current_chamber_ID][t_eta]["ErrR"],"Residuals/MatchingOn_"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID+"/PropagationError/AllHits")

                    ## feature unavailable for glb_phi
                    if matchingVar == "glb_rdphi":
                        writeToTFile(OutF,TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID][t_eta]["ErrR"],"Residuals/MatchingOn_"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID+"/PropagationError/MatchedHits")
                        writeToTFile(OutF,TH1FpropError_collector["MatchedHits"][endcapTag][current_chamber_ID][t_eta]["ErrPhi"],"Residuals/MatchingOn_"+matchingVar+"/"+endcapTag+"/"+current_chamber_ID+"/PropagationError/MatchedHits")

    writeToTFile(OutF,SummaryAll,"Efficiency/"+matchingVar+"/")
    writeToTFile(OutF,EffiDistrAll,"Efficiency/"+matchingVar+"/")


    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Reco'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/")
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['Prop'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalLong'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 
    writeToTFile(OutF,THSanityChecks['Occupancy'][matchingVar]['AfterMatching']['PropLocalShort'],"SanityChecks/Occupancy/AfterMatching_"+matchingVar+"/") 


for key in TH1MetaData.keys():
    writeToTFile(OutF,TH1MetaData[key],"Metadata/")

printSummary(EfficiencyDictGlobal,matching_variables,ResidualCutOff,matching_variable_units,debug=args.verbose)

## CSV
for matchingVar in matching_variables:
    tempList = []
    tempList_byVFAT = []
    tempList_byDLE = []
    for etaPID,subDict in EfficiencyDictGlobal[matchingVar].items():
        region,chamber,layer,eta = getInfoFromEtaID(etaPID)
        
        matchedRecHit = sum([subDict[k]['num'] for k in subDict.keys()])
        propHit = sum([subDict[k]['den'] for k in subDict.keys()])
        chID = ReChLa2chamberName(region,chamber,layer)

        AVG_CLS = THSanityChecks['CLS'][matchingVar][chID][eta].GetMean()
        AVG_pt = THSanityChecks['pt'][chID][eta].GetMean()
        tempList.append([chID,region,chamber,layer,eta,matchedRecHit,propHit,AVG_CLS,AVG_pt])
            
    for (re,la) in [(1,1),(1,2),(-1,1),(-1,2)]:
        endcap_key = EndcapLayer2label(re,la)
        for chID,subDict in EfficiencyDictVFAT[matchingVar][endcap_key].items():
            for VFATN in subDict.keys():
                tempList_byVFAT.append([chID,endcap_key,VFATN,subDict[VFATN]['num'],subDict[VFATN]['den']])

    data = pd.DataFrame(tempList,columns=['chamberID',"region","chamber","layer","etaPartition","matchedRecHit","propHit","AVG_CLS","AVG_pt"])
    data.to_csv('./Output/PFA_Analyzer_Output/CSV/'+outputname+'/MatchingSummary_'+matchingVar+'.csv', index=False)
    
    data_byVFAT = pd.DataFrame(tempList_byVFAT,columns=['chamberID',"EndcapTag","VFATN","matchedRecHit","propHit"])
    data_byVFAT.to_csv('./Output/PFA_Analyzer_Output/CSV/'+outputname+'/MatchingSummary_'+matchingVar+'_byVFAT.csv', index=False)
if DLE:
    for etaPID,subDict in EfficiencyDictLayer['glb_rdphi'].items():
        region,chamber,layer,eta = getInfoFromEtaID(etaPID)

        matchedRecHit = sum([subDict[k]['num'] for k in subDict.keys()])
        propHit = sum([subDict[k]['den'] for k in subDict.keys()])
        chID = ReChLa2chamberName(region,chamber,layer)
        tempList_byDLE.append([chID,region,chamber,layer,eta,matchedRecHit,propHit])
    data_byDLE = pd.DataFrame(tempList_byDLE,columns=['chamberID',"region","chamber","layer","etaPartition","matchedRecHit","propHit"])
    data_byDLE.to_csv('./Output/PFA_Analyzer_Output/CSV/'+outputname+'/MatchingSummary_glb_rdphi_byDLE.csv', index=False)

print "\n#############\nOUTPUT\n#############"
print "\tCSVs in \t"+'./Output/PFA_Analyzer_Output/CSV/'+outputname
print "\tROOT_File \t"+"./Output/PFA_Analyzer_Output/ROOT_File/"+outputname+".root"
# print "\tMatchingTTree File \t"+"./Output/PFA_Analyzer_Output/ROOT_File/MatchingTree_"+outputname+".root"
print "\tPlots in\t"+'./Output/PFA_Analyzer_Output/Plot/'+outputname


