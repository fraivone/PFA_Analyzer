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
sys.path.insert(1, lib_folder+'/ROOT_Utils/')
sys.path.insert(2, lib_folder+'/GE11Geometry/')
try:
    from ROOT_Utils import *
    from PFA_Analyzer_Utils import *

except:
  print("ERROR:\n\tCan't find the package CMS_lumi and tdrstlye\n\tPlease verify that this file are placed in the path $myLIB/ROOT_Utils/ \n\tAdditionally keep in mind to export the environmental variable $myLIB\nEXITING...\n") 
  sys.exit(0)



ROOT.gROOT.SetBatch(True)


start_time = time.time()

TH1STAangleAll=ROOT.TH1F("STAangleAll","STAangleAll",100,0,1)
TH1STAangleisME11=ROOT.TH1F("STAangleisME11","STAangleisME11",100,0,1)
TH1STAangle1xME1=ROOT.TH1F("STAangle1xME1","STAangle1xME1",100,0,1)
TH1STAangle4xME1=ROOT.TH1F("STAangle4xME1","STAangle4xME1",100,0,1)
TH1STAangle4xME12=ROOT.TH1F("STAangle4xME12","STAangle4xME12",100,0,1)
TH1STAangle4xME123=ROOT.TH1F("STAangle4xME123","STAangle4xME123",100,0,1)
TH1STAangle4xME1234=ROOT.TH1F("STAangle4xME1234","STAangle4xME1234",100,0,1)
TH2Chi2vsAngle = ROOT.TH2F("Chi2vsAngle","Chi2vsAngle",200,0,10,200,0,1)

## Input files
files = files_in_folder("/eos/user/f/fivone/GEMNTuples/MC/Output/MC_Cosmics2021/wMEXHits")
files = [f for f in files if ".root" in f]
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

    if n_gemprop == 0:
        continue

    RecHit_Dict = {}
    PropHit_Dict = {}


    for RecHit_index in range(0,n_gemrec):
        region = evt.gemRecHit_region[RecHit_index]
        chamber = evt.gemRecHit_chamber[RecHit_index]
        layer = evt.gemRecHit_layer[RecHit_index]
        etaP = evt.gemRecHit_etaPartition[RecHit_index]
        RecHitEtaPartitionID = region*(100*chamber+10*layer+etaP)
        endcapKey = "PL"+str(layer) if region > 0 else "ML"+str(layer)


        RecHit_Dict.setdefault(RecHitEtaPartitionID, {'loc_x':[],'glb_x':[],'glb_y':[],'glb_z':[],'glb_r':[],'glb_phi':[],'firstStrip':[],'cluster_size':[]})
        RecHit_Dict[RecHitEtaPartitionID]['loc_x'].append(evt.gemRecHit_loc_x[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_x'].append(evt.gemRecHit_g_x[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_y'].append(evt.gemRecHit_g_y[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_z'].append(evt.gemRecHit_g_z[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_r'].append(evt.gemRecHit_g_r[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['glb_phi'].append(evt.gemRecHit_g_phi[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['firstStrip'].append(evt.gemRecHit_firstClusterStrip[RecHit_index])
        RecHit_Dict[RecHitEtaPartitionID]['cluster_size'].append(evt.gemRecHit_cluster_size[RecHit_index])

    for PropHit_index in range(0,n_gemprop):
        
        region = evt.mu_propagated_region[PropHit_index]
        chamber = evt.mu_propagated_chamber[PropHit_index]
        layer = evt.mu_propagated_layer[PropHit_index]
        etaP = evt.mu_propagated_etaP[PropHit_index]
        PropHitChamberID = region*(100*chamber+10*layer+etaP)
        endcapKey = "PL"+str(layer) if region > 0 else "ML"+str(layer)

        outermost_z = evt.mu_propagated_Outermost_z[PropHit_index]
        is_incoming = evt.mu_isincoming[PropHit_index]

        
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
        
            

    ### Matching criteria between propagated hits from ME11 and RecHits : 
    ##      1.SAME REGION,SC,LAYER,ETA -->SAME etaPartitionID
    for etaPartitionID in PropHit_Dict.keys():
        
        region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)
        PropHitonEta = PropHit_Dict[etaPartitionID]

        nPropHitsOnEtaID = len(PropHitonEta['glb_phi'])
        i = 0
        # for i in range(nPropHitsOnEtaID):
        TrackChi2 = PropHitonEta['STA_Normchi2'][i]
        isME11 = PropHitonEta['mu_propagated_isME11'][i]
        nME1 = PropHitonEta['nME1Hits'][i]
        nME2 = PropHitonEta['nME2Hits'][i]
        nME3 = PropHitonEta['nME3Hits'][i]
        nME4 = PropHitonEta['nME4Hits'][i]
        
        cosx = PropHitonEta['Loc_dirX'][i]
        cosy = PropHitonEta['Loc_dirY'][i]
        dirChamber = ( cosx**2 + cosy**2 )**(0.5)

        
        TH1STAangleAll.Fill(dirChamber)
        TH2Chi2vsAngle.Fill(TrackChi2,dirChamber)
        if isME11:
            TH1STAangleisME11.Fill(dirChamber)

        if nME1 > 0:
            TH1STAangle1xME1.Fill(dirChamber)
        if nME1 > 3:
            TH1STAangle4xME1.Fill(dirChamber)
        if nME1>3 and nME2 > 3:
            TH1STAangle4xME12.Fill(dirChamber)
        if nME1 > 3 and nME2 > 3 and nME3 > 3:
            TH1STAangle4xME123.Fill(dirChamber)
        if nME1 > 3 and nME2 > 3 and nME3 > 3 and nME4 > 3:
            TH1STAangle4xME1234.Fill(dirChamber)
        

        
            
    ## End on the matching/filling loop  

print("--- %s seconds ---" % (time.time() - start_time))


OutF = ROOT.TFile("./Extractor.root","RECREATE")

writeToTFile(OutF,TH1STAangleAll)
writeToTFile(OutF,TH2Chi2vsAngle)

writeToTFile(OutF,TH1STAangleisME11)
writeToTFile(OutF,TH1STAangle1xME1)
writeToTFile(OutF,TH1STAangle4xME1)
writeToTFile(OutF,TH1STAangle4xME12)
writeToTFile(OutF,TH1STAangle4xME123)
writeToTFile(OutF,TH1STAangle4xME1234)        
