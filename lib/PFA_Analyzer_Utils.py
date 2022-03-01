from numpy.core.numeric import NaN
import ROOT
import numpy
import pandas as pd
from ROOT_Utils import *
import sys
import json
from array import array


def getInfoFromEtaID(id):
    etaPartition = abs(id)%10
    layer = ((abs(id)-etaPartition)%100)/10
    chamber = (abs(id)-etaPartition-10*layer)/100
    region = id/abs(id)
    return region,chamber,layer,etaPartition

## Expects to have a list of file paths each containing a list of "GE11-*-XXLY" to be excluded
def importOFFChamber(list_of_file_path):
    chamberNumberOFF = []
    for file_path in list_of_file_path:   
        try:
            with open(file_path, "r") as my_file:
                content = my_file.read()
                chamber_OFF_list = content.split("\n")
                try: ## Will remove 1 empty entry from the list, if any
                    chamber_OFF_list.remove('')
                except:
                    pass
                try:
                    findChambers = [ chamberName2ReChLa(k) for k in chamber_OFF_list]
                except:
                    print "Error in the formatting of Chamber names in the file " + file_path+"\nExiting .."
                    sys.exit(0)
                ## Avoiding duplicates
                for chamber in findChambers:
                    if chamber not in chamberNumberOFF:
                        chamberNumberOFF.append(chamber)
        except IOError:
            print "Couldn't open file : "+file_path+"\nExiting .."
            sys.exit(0)
    
    return chamberNumberOFF



## Expects to have a list of file paths each containing a list of "GE11-*-XXLY" to be excluded
def ChamberOFF_byLS(file_path):
    try:
        io = open(file_path[0],"r")
        dictionary = json.load(io)  
    except IOError:
            print "Couldn't open file : "+file_path+"\nExiting .."
            sys.exit(0)

    for key,item in dictionary.items():
        if item == []:
            dictionary.pop(key)
    return dictionary

def VFAT2iEta_iPhi(VFATN):
    try:
        vfatPosition = int(VFATN)
    except:
        print "VFAT Number provided is not a number.\nExiting..."
        sys.exit(0)

    if vfatPosition <0 or vfatPosition>23:
        print "Invalid VFAT position.\nExiting..."
        sys.exit(0)

    iEta = (8 - vfatPosition%8)
    iPhi = vfatPosition/8
    return iEta,iPhi

def iEta_iPhi2VFAT(iEta,iPhi):
    try:
        etaP = int(iEta)
        phi = int(iPhi)
    except:
        print "Provided iEta and/or iPhi are not numbers.\nExiting..."
        sys.exit(0)
    
    if iEta <1 or iEta>8 or iPhi <0 or iPhi>2:
        print "Invalid iEta and/or iPhi provided position.\nExiting..."
        sys.exit(0)
    
    VFAT = iPhi*8 + (8-iEta)
    return VFAT

## Associates a propagated hit to a VFAT
## It is based on default GE11 geometry
## Might need refinements after alignment 
def propHit2VFAT(glb_r,loc_x,etaP,region,chamber_number):
    ## Magic numbers taken from strip geometry
    #  They are loc_x/glb_r (basically cos phi) for strip 127 and 255 respectively
    LeftMost_iPhi_boundary = -0.029824804
    RightMost_iPhi_boundary = 0.029362374

    prophit_cosine = loc_x/glb_r
    iPhi =  0 if  prophit_cosine<=LeftMost_iPhi_boundary else (2 if prophit_cosine>RightMost_iPhi_boundary else 1)
        
    VFAT = iEta_iPhi2VFAT(etaP,iPhi)

    return VFAT

## returns a dict with a list of OFF VFAT for each  etapartitionID (key of the dict)
def importOFFVFAT(list_of_file_path):
    off_VFAT = {}
    for file_path in list_of_file_path:
        try:
            df = pd.read_csv(file_path, sep='\t')
        except IOError:
            print "Couldn't open file : "+file_path+"\nExiting .."
            sys.exit(0)

        off_VFAT = {}
        for index, row in df.iterrows():
            region =  -1 if ( row['region']=='N' or row['region']=='M') else (1 if row['region']=='P' else None) ## In the past N was used... currently M indicates the minus endcap
            layer = int(row['layer'])
            chamber = int(row['chamber'])
            VFAT = int(row['VFAT'])
            maskReason = int(row['reason_mask'])

            iEta , iPhi = VFAT2iEta_iPhi(VFAT)
            etaPartitionID = region*(100*chamber+10*layer+iEta)

            ## If the key doesn't exsist create it, then fill the list
            off_VFAT.setdefault(etaPartitionID,[])
            if maskReason == 1:
                off_VFAT[etaPartitionID].append(VFAT)

    return off_VFAT

def chamberName2ReChLa(chamberName):
    ## Accepts as input either 
        # GE11-M-03L1 or GE11-M-03L1-S
    if len(chamberName)==13:
        chamberName = chamberName[:11]
    re = -1 if "M" in chamberName else 1
    ch = int( chamberName.split("-")[-1][:2] )
    la = int( chamberName.split("-")[-1][-1] )

    return [re,ch,la]

def ReChLa2chamberName(re,ch,la):
    endcap = "M" if re == -1 else "P"
    size = "S" if ch%2 == 1 else "L"
    chID = 'GE11-'+endcap+'-%02d' % ch +"L"+str(la)+"-"+size 
    return chID

def EndcapLayer2label(re,layer):
    label = "ML" if re == -1 else "PL"
    return label+str(layer)
    

def ChambersOFFHisto(chamberNumberOFF_byLS,maxLS):
    GE11_OFF = {-1:{},1:{}}
    GE11_OFF[-1][1] = ROOT.TH1F("GE11-M-L1  Masked Chambers","GE11-M-L1  Masked Chambers",36,0.5,36.5)
    GE11_OFF[-1][2] = ROOT.TH1F("GE11-M-L2  Masked Chambers","GE11-M-L2  Masked Chambers",36,0.5,36.5)
    GE11_OFF[1][1] = ROOT.TH1F("GE11-P-L1  Masked Chambers","GE11-P-L1  Masked Chambers",36,0.5,36.5)
    GE11_OFF[1][2] = ROOT.TH1F("GE11-P-L2  Masked Chambers","GE11-P-L2  Masked Chambers",36,0.5,36.5)

    for region in [-1,1]:
        for layer in [1,2]:
            for chamber in range(1,37):
                ch_ID = ReChLa2chamberName(region,chamber,layer) 
                n_OFF = 0 if ch_ID not in chamberNumberOFF_byLS else (len(chamberNumberOFF_byLS[ch_ID]) if -1 not in chamberNumberOFF_byLS[ch_ID] else maxLS)
                percentageOFF = 100*float(n_OFF)/maxLS
                GE11_OFF[region][layer].SetBinContent(chamber,percentageOFF)
    for key_1 in GE11_OFF.keys():
        for key_2 in GE11_OFF[key_1].keys():
            GE11_OFF[key_1][key_2].GetYaxis().SetTickLength(0.005)
            GE11_OFF[key_1][key_2].GetYaxis().SetTitle("Masked Lumisection (%)")
            GE11_OFF[key_1][key_2].SetStats(False)
            GE11_OFF[key_1][key_2].SetMinimum(0)
            GE11_OFF[key_1][key_2].SetMaximum(110)
    return GE11_OFF[-1][1],GE11_OFF[1][1],GE11_OFF[-1][2],GE11_OFF[1][2]



def VFATOFFHisto(VFATOFF_dictionary):
    VFAT_OFFTH2D = {-1:{},1:{}}
    VFAT_OFFTH2D[-1][1] = ROOT.TH2F("GE11-M-L1  Masked VFAT","GE11-M-L1  Masked VFAT",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[-1][2] = ROOT.TH2F("GE11-M-L2  Masked VFAT","GE11-M-L2  Masked VFAT",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[1][1] = ROOT.TH2F("GE11-P-L1  Masked VFAT","GE11-P-L1  Masked VFAT",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[1][2] = ROOT.TH2F("GE11-P-L2  Masked VFAT","GE11-P-L2  Masked VFAT",36,0.5,36.5,24,-0.5,23.5)

    for key_1 in VFAT_OFFTH2D.keys():
        for key_2 in VFAT_OFFTH2D[key_1].keys():
            VFAT_OFFTH2D[key_1][key_2].GetYaxis().SetTitle("VFAT Number")
            VFAT_OFFTH2D[key_1][key_2].GetXaxis().SetTitle("Chamber")
            VFAT_OFFTH2D[key_1][key_2].SetStats(False)
            VFAT_OFFTH2D[key_1][key_2].SetMinimum(0)

    for key, value in VFATOFF_dictionary.items():
        region,chamber,layer,etaPartition = getInfoFromEtaID(key)
        for VFATN in value:
            VFAT_OFFTH2D[region][layer].Fill(chamber,VFATN)
    
    return VFAT_OFFTH2D[-1][1],VFAT_OFFTH2D[1][1],VFAT_OFFTH2D[-1][2],VFAT_OFFTH2D[1][2]


def GE11DiscardedSummary(chamberNumberOFF_byLS,maxLS,VFATOFF_dictionary):
    VFAT_OFFTH2D = {-1:{},1:{}}
    VFAT_OFFTH2D[-1][1] = ROOT.TH2F("GE11-M-L1  Masked","GE11-M-L1  Masked",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[-1][2] = ROOT.TH2F("GE11-M-L2  Masked","GE11-M-L2  Masked",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[1][1] = ROOT.TH2F("GE11-P-L1  Masked","GE11-P-L1  Masked",36,0.5,36.5,24,-0.5,23.5)
    VFAT_OFFTH2D[1][2] = ROOT.TH2F("GE11-P-L2  Masked","GE11-P-L2  Masked",36,0.5,36.5,24,-0.5,23.5)


    for key_1 in VFAT_OFFTH2D.keys():
        for key_2 in VFAT_OFFTH2D[key_1].keys():
            VFAT_OFFTH2D[key_1][key_2].GetYaxis().SetTitle("VFAT Number")
            VFAT_OFFTH2D[key_1][key_2].GetYaxis().SetTickLength(0.005)
            VFAT_OFFTH2D[key_1][key_2].GetXaxis().SetTickLength(0.005)
            VFAT_OFFTH2D[key_1][key_2].GetXaxis().SetTitle("Chamber")
            VFAT_OFFTH2D[key_1][key_2].GetZaxis().SetTitle("Masked Lumisection (%)")
            VFAT_OFFTH2D[key_1][key_2].GetZaxis().SetTitleSize(0.03)
            VFAT_OFFTH2D[key_1][key_2].GetZaxis().SetTitleOffset(1.1)
            VFAT_OFFTH2D[key_1][key_2].GetZaxis().SetLabelSize(0.025)
            VFAT_OFFTH2D[key_1][key_2].SetStats(False)
            VFAT_OFFTH2D[key_1][key_2].SetMinimum(0)
            VFAT_OFFTH2D[key_1][key_2].SetMaximum(100)
    
    for region in [-1,1]:
        for layer in [1,2]:
            for chamber in range(1,37):
                current_chamber_ID = ReChLa2chamberName(region,chamber,layer)
                n_OFF = 0 if current_chamber_ID not in chamberNumberOFF_byLS else (len(chamberNumberOFF_byLS[current_chamber_ID]) if -1 not in chamberNumberOFF_byLS[current_chamber_ID] else maxLS)
                percentageOFF = 100*float(n_OFF)/maxLS
                for VFATN in range(0,24):
                    eta,phi = VFAT2iEta_iPhi(VFATN)
                    etaPartitionID = region*(100*chamber+10*layer+eta)

                    if etaPartitionID in VFATOFF_dictionary.keys() and VFATN in VFATOFF_dictionary[etaPartitionID]:
                        VFAT_OFFTH2D[region][layer].SetBinContent(chamber,VFATN+1,100.)
                    else:
                        VFAT_OFFTH2D[region][layer].SetBinContent(chamber,VFATN+1,percentageOFF)

    return [VFAT_OFFTH2D[-1][1],VFAT_OFFTH2D[1][1],VFAT_OFFTH2D[-1][2],VFAT_OFFTH2D[1][2]]

def incidenceAngle_vs_Eff(sourceDict,input_region=1,input_layer=1):
    ## Transforming in list
    reg_tag_string = "All" if isinstance(input_region, list) else "P" if input_region == 1 else "M"
    lay_tag_string = "" if isinstance(input_layer, list) else "L1" if input_layer == 1 else "L2"
    input_region = [input_region] if not isinstance(input_region, list) else input_region
    input_layer = [input_layer] if not isinstance(input_layer, list) else input_layer

    title = "GE11"+reg_tag_string+lay_tag_string
    angle_nbins,angle_min,angle_max = 20,0,1.
    
    
    Eff_Plot = ROOT.TGraphAsymmErrors()
    NumTH1F = ROOT.TH1F(title+"_incidenceAngle_Num",title+"_incidenceAngle_Num",angle_nbins,angle_min,angle_max)
    DenTH1F = ROOT.TH1F(title+"_incidenceAngle_Den",title+"_incidenceAngle_Den",angle_nbins,angle_min,angle_max)
        

    for j in range(0,20):
        etaPartitionRecHits = 0
        etaPartitionPropHits = 0
        for etaPartitionID,value in sourceDict.items():
            region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)
            
            if layer not in input_layer or region not in input_region:
                continue

            etaPartitionRecHits  += value[j]['num']
            etaPartitionPropHits += value[j]['den']
        NumTH1F.SetBinContent((j+1),etaPartitionRecHits)
        DenTH1F.SetBinContent((j+1),etaPartitionPropHits)
    
    Eff_Plot.Divide(NumTH1F,DenTH1F,"B")
    Eff_Plot.SetTitle(title+"_incidenceAngle_Eff")
    Eff_Plot.SetName(title+"_incidenceAngle_Eff")
    Eff_Plot.GetXaxis().SetTitle("Cos(#alpha)")
    Eff_Plot.GetYaxis().SetTitle("Efficiency")
    
    return NumTH1F, DenTH1F, Eff_Plot



## Structure:: TH2Fresidual_collector[matchingvar][chambers][residual_of_what][Plot] == Contains the TH2F
## Structure:: TH2Fresidual_collector[matchingvar][chambers][residual_of_what][binx][biny] == list(n entries, sum(abs(residual)))
## Usage:: TH2Fresidual_collector[matchingvar][chambers][residual_of_what][binx][biny] == Fill at each iteration
## Usage:: TH2Fresidual_collector[matchingvar][chambers][residual_of_what][Plot] == Fill at the end of the loop
def generate2DResidualContainer(matching_variables,nbins,minB):
    TH2Fresidual_collector = {}

    for key_1 in matching_variables:
        TH2Fresidual_collector.setdefault(key_1,{'all':{},'short':{},'long':{}})
        for key_2 in ['all','short','long']:
            TH2Fresidual_collector[key_1][key_2] = {}
            for key_3 in ['glb_phi','glb_rdphi']:
                titleTH2 = key_2 + "Chmbrs_MatchedBy_"+key_1+"_2DMap_of_"+key_3+"_res" 
                TH2Fresidual_collector[key_1][key_2][key_3] = {}
                TH2Fresidual_collector[key_1][key_2][key_3]['TH2F'] = ROOT.TH2F(titleTH2,titleTH2,nbins,minB,-minB,nbins,minB,-minB)
                TH2Fresidual_collector[key_1][key_2][key_3]['TH2F'].GetXaxis().SetTitle("Loc_x (cm)")
                TH2Fresidual_collector[key_1][key_2][key_3]['TH2F'].GetYaxis().SetTitle("Loc_y (cm)")
    
                for x_bin in range(1,nbins+1):
                    TH2Fresidual_collector[key_1][key_2][key_3][x_bin] = {}
                    for y_bin in range(1,nbins+1):
                        TH2Fresidual_collector[key_1][key_2][key_3][x_bin].setdefault(y_bin,[0,0])

    return TH2Fresidual_collector


def generate1DResidualContainer(matching_variables,nbins,ResidualCutOff):
    output_dict = {}
    for mv in matching_variables:
        output_dict[mv] = {}
        minB = -ResidualCutOff["glb_rdphi"]*2

        for (re,la) in [(1,1),(1,2),(-1,1),(-1,2)]:
            endcap_key = EndcapLayer2label(re,la)
            output_dict[mv][endcap_key] = {}

            for chamber in range(1,37):
                chID = ReChLa2chamberName(re,chamber,la)
                output_dict[mv][endcap_key][chID] = {}
                for eta in range(1,9)+["All"]:
                    titleTH1 = chID+"_eta"+str(eta)+"_"
                    output_dict[mv][endcap_key][chID][eta]={"Residual":ROOT.TH1F(titleTH1+"Residual",titleTH1+"Residual",nbins,minB,-minB)}
                    output_dict[mv][endcap_key][chID][eta]["Residual"].GetXaxis().SetTitle("#Deltardphi (cm)")
    return output_dict

# dict[matchingVar][endcaptag][chamberID][VFAT]
def generateVFATDict(matching_variables):
    output_dict = {}
    for mv in matching_variables:
        output_dict[mv] = {}

        for (re,la) in [(1,1),(1,2),(-1,1),(-1,2)]:
            endcap_key = EndcapLayer2label(re,la)
            output_dict[mv][endcap_key] = {}
            for chamber in range(1,37):
                chID = ReChLa2chamberName(re,chamber,la)
                output_dict[mv][endcap_key][chID] = {}
                for VFAT in range(24):
                    output_dict[mv][endcap_key][chID][VFAT]={"num":0,"den":0}
    return output_dict



def generatePropagationErrorContainer(maxErrOnPropR, maxErrOnPropPhi):
    output_dict = {}

    for which_hit in ["MatchedHits","AllHits"]:
        output_dict[which_hit] = {}
        for (re,la) in [(1,1),(1,2),(-1,1),(-1,2)]:
            endcap_key = EndcapLayer2label(re,la)
            output_dict[which_hit][endcap_key] = {}

            for chamber in range(1,37):
                chID = ReChLa2chamberName(re,chamber,la)
                output_dict[which_hit][endcap_key][chID] = {}
                for eta in range(1,9)+["All"]:
                    titleTH1 = chID+"_eta"+str(eta)+"_"
                    output_dict[which_hit][endcap_key][chID][eta] = {"ErrPhi":ROOT.TH1F(titleTH1+"PropagationError on Phi",titleTH1+"PropagationError on Phi",200,0,2*maxErrOnPropPhi),
                                                          "ErrR":ROOT.TH1F(titleTH1+"PropagationError on R",titleTH1+"PropagationError on R",100,0,2*maxErrOnPropR)}
                    output_dict[which_hit][endcap_key][chID][eta]["ErrPhi"].GetXaxis().SetTitle("ErrPhi (rad)")
                    output_dict[which_hit][endcap_key][chID][eta]["ErrR"].GetXaxis().SetTitle("ErrR (cm)")

    return output_dict
        

def fillPlot2DResidualContainer(TH2Fresidual_collector,matching_variables,nbins):
    
    for key_1 in matching_variables:
        for key_2 in ['all','short','long']:
            for key_3 in ['glb_phi','glb_rdphi']:
                for x_bin in range(1,nbins+1):
                    for y_bin in range(1,nbins+1):
                        AVG_Residual = TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][1]/TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][0] if TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][0]!=0 else 0
                        TH2Fresidual_collector[key_1][key_2][key_3]['TH2F'].SetBinContent(x_bin,y_bin,AVG_Residual)
    return TH2Fresidual_collector

def passCut(PropHitonEta,prop_hit_index,maxPropR_Err=0.7,maxPropPhi_Err=0.001,fiducialCutR=0.5,fiducialCutPhi=0.002,minPt=0.,maxChi2=9999999,minME1Hit=0,minME2Hit=0,minME3Hit=0,minME4Hit=0):
    passedCut = True
    if PropHitonEta['err_glb_phi'][prop_hit_index] > maxPropPhi_Err:
        passedCut = False
    if PropHitonEta['err_glb_r'][prop_hit_index] > maxPropR_Err:
        passedCut = False

    if PropHitonEta['STA_Normchi2'][prop_hit_index] > maxChi2 or PropHitonEta['STA_Normchi2'][prop_hit_index] < 0.5:
        passedCut = False

    if PropHitonEta['nME1Hits'][prop_hit_index] < minME1Hit:
        passedCut = False
    if PropHitonEta['nME2Hits'][prop_hit_index] < minME2Hit:
        passedCut = False
    if PropHitonEta['nME3Hits'][prop_hit_index] < minME3Hit:
        passedCut = False
    if PropHitonEta['nME4Hits'][prop_hit_index] < minME4Hit:
        passedCut = False

    
    PhiMin = PropHitonEta['mu_propagated_EtaPartition_phiMin'][prop_hit_index]
    PhiMax = PropHitonEta['mu_propagated_EtaPartition_phiMax'][prop_hit_index]
    PropHitPhi = PropHitonEta['glb_phi'][prop_hit_index]
    PropHitPt = PropHitonEta['pt'][prop_hit_index]

    if PhiMin > PhiMax: # Happens for chamber 19 cause 181 degrees becomes -179 degrees. So chamber 19 has phiMin = 174 degrees phiMax = -174
        PhiMax = 2*numpy.pi + PhiMax 
        if PropHitPhi<0:
            PropHitPhi = 2*numpy.pi + PropHitPhi


    if PropHitPhi < (PhiMin+fiducialCutPhi) or PropHitPhi > (PhiMax-fiducialCutPhi):
        passedCut = False

    if PropHitPt < minPt:
        passedCut = False
    
    ## Fiducial cut on chamber perimeter
    # if PropHitonEta['etaP'][prop_hit_index] == 1 and PropHitonEta['glb_r'][prop_hit_index] > (PropHitonEta['mu_propagated_EtaPartition_rMax'][prop_hit_index]-fiducialCutR):
    #     passedCut = False
    # if PropHitonEta['etaP'][prop_hit_index] == 8 and PropHitonEta['glb_r'][prop_hit_index] < (PropHitonEta['mu_propagated_EtaPartition_rMin'][prop_hit_index]+fiducialCutR):
    #     passedCut = False    

    ## Fiducial cut on etaP perimeter
    if PropHitonEta['glb_r'][prop_hit_index] > (PropHitonEta['mu_propagated_EtaPartition_rMax'][prop_hit_index]-fiducialCutR):
        passedCut = False
    if PropHitonEta['glb_r'][prop_hit_index] < (PropHitonEta['mu_propagated_EtaPartition_rMin'][prop_hit_index]+fiducialCutR):
        passedCut = False
    return passedCut

## Generate confidence level limits for value obtained from ratio of Poissonian
## Knowing that MUST BE Num < Den
## Not needed...ROOT can do it automagically with TGraphAsymmErrors::DIVIDE
def generateClopperPeasrsonInterval(num,den):
    confidenceLevel = 0.68
    alpha = 1 - confidenceLevel
    
    lowerLimit = ROOT.Math.beta_quantile(alpha/2,num,den-num + 1)
    if num==den:
        upperLimit=1
    else:
        upperLimit = ROOT.Math.beta_quantile(1-alpha/2,num + 1,den-num)
    return lowerLimit,upperLimit

def generateEfficiencyPlotbyEta(sourceDict,input_region=1,input_layer=1):
    ## Transforming in list
    reg_tag_string = "All" if isinstance(input_region, list) else "P" if input_region == 1 else "M"
    lay_tag_string = "" if isinstance(input_layer, list) else "L1" if input_layer == 1 else "L2"
    input_region = [input_region] if not isinstance(input_region, list) else input_region
    input_layer = [input_layer] if not isinstance(input_layer, list) else input_layer
    
    TH1F_TempContainer = {}
    Plot_Container = {}
    ColorAssociation = {'All':ROOT.kBlack,'Long':ROOT.kGreen+1,'Short':ROOT.kRed}
    
    title = "GE11"+reg_tag_string+lay_tag_string+"_EffbyEta"
    for chambers in ['All','Long','Short']:   
        Plot_Container[chambers] = ROOT.TGraphAsymmErrors()
        Plot_Container[chambers].GetXaxis().SetTitle("GE11 Eta Partition")
        Plot_Container[chambers].GetYaxis().SetTitle("Efficiency")
        Plot_Container[chambers].SetMaximum(1.1)
        Plot_Container[chambers].SetTitle(chambers+"_"+title)
        Plot_Container[chambers].SetName(chambers+"_"+title)
        Plot_Container[chambers].SetLineColor(ColorAssociation[chambers])
        Plot_Container[chambers].SetMarkerColor(ColorAssociation[chambers])
        Plot_Container[chambers].SetFillColorAlpha(ColorAssociation[chambers],.4)
        Plot_Container[chambers].SetMarkerStyle(20)
        Plot_Container[chambers].SetMarkerSize(.8)

        TH1F_TempContainer.setdefault(chambers,{'num':ROOT.TH1F('num_'+chambers+title, title,8,0.5,8.5),'den':ROOT.TH1F('den_'+chambers+title, title,8,0.5,8.5)})

    for etaPartitionID in sourceDict.keys():
        region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)
    
        if layer not in input_layer or region not in  input_region:
            continue
        
        TH1F_TempContainer['All']['num'].SetBinContent( eta, TH1F_TempContainer['All']['num'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['num'] for pt in range(0,11)]) )
        TH1F_TempContainer['All']['den'].SetBinContent( eta, TH1F_TempContainer['All']['den'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['den'] for pt in range(0,11)]) )

        if chamber%2==0:
            TH1F_TempContainer['Long']['num'].SetBinContent( eta, TH1F_TempContainer['Long']['num'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['num'] for pt in range(0,11)]) )
            TH1F_TempContainer['Long']['den'].SetBinContent( eta, TH1F_TempContainer['Long']['den'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['den'] for pt in range(0,11)]) )
        else:
            TH1F_TempContainer['Short']['num'].SetBinContent( eta, TH1F_TempContainer['Short']['num'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['num'] for pt in range(0,11)]) )
            TH1F_TempContainer['Short']['den'].SetBinContent( eta, TH1F_TempContainer['Short']['den'].GetBinContent(eta) + sum([sourceDict[etaPartitionID][pt]['den'] for pt in range(0,11)]) )
    
    
    for chambers in ['All','Long','Short']: 
        Plot_Container[chambers].Divide(TH1F_TempContainer[chambers]['num'],TH1F_TempContainer[chambers]['den'],"B")

    return Plot_Container['Short'],Plot_Container['Long'],Plot_Container['All']

def generate1DEfficiencyPlotbyVFAT(EfficiencyDictVFAT,chamber_ID):
    
    re_ch_la_list = chamberName2ReChLa(chamber_ID)
    re = re_ch_la_list[0]
    ch = re_ch_la_list[1]
    la = re_ch_la_list[2]
    endcapTag = EndcapLayer2label(re,la)

    title = chamber_ID+"_1DEffVFAT"
    Summary = ROOT.TGraphAsymmErrors()
    Summary.GetXaxis().SetTitle("VFAT Number")
    Summary.GetYaxis().SetTitle("Efficiency")
    Summary.SetMaximum(1.1)
    Summary.SetMinimum(0)
    Summary.SetTitle(title)
    Summary.SetName(title)
    Summary.SetMarkerStyle(20)
    Summary.SetMarkerSize(.8)

    TH1Num = ROOT.TH1F(title, title,24,-0.5,23.5)
    TH1Den = ROOT.TH1F(title, title,24,-0.5,23.5)

    for VFATN in range(24):
        TH1Num.SetBinContent(VFATN,EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['num'])
        TH1Den.SetBinContent(VFATN,EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['den'])
    
    Summary.Divide(TH1Num,TH1Den,"B")
    Summary.GetXaxis().SetLimits(-0.5,23.5)

    return Summary


def generate2DEfficiencyPlotbyVFAT(EfficiencyDictVFAT,efficiency_target):
    ## GEOMETRY DEFINITION
    VFATGeometryDict = {"Long":dict((VFATN,{}) for VFATN in range(0,24)),"Short":dict((VFATN,{}) for VFATN in range(0,24))}

    VFATGeometryDict['Short'][0]['x'] = np.asarray([-3.83938562003834,-4.14080692494689,-12.4513908518566,-11.5450181214881],dtype=float)
    VFATGeometryDict['Short'][0]['y'] = np.asarray([0,10.206,10.206,0],dtype=float)
    VFATGeometryDict['Short'][1]['x'] = np.asarray([-4.14080692494689,-4.44222822985544,-13.3577635822252,-12.4513908518566],dtype=float)
    VFATGeometryDict['Short'][1]['y'] = np.asarray([10.206,20.412,20.412,10.206],dtype=float)
    VFATGeometryDict['Short'][2]['x'] = np.asarray([-4.44222822985544,-4.79536310569235,-14.4196388259069,-13.3577635822252],dtype=float)
    VFATGeometryDict['Short'][2]['y'] = np.asarray([20.412,32.369,32.369,20.412],dtype=float)
    VFATGeometryDict['Short'][3]['x'] = np.asarray([-4.79536310569235,-5.14849798152926,-15.4815140695887,-14.4196388259069],dtype=float)
    VFATGeometryDict['Short'][3]['y'] = np.asarray([32.369,44.326,44.326,32.369],dtype=float)
    VFATGeometryDict['Short'][4]['x'] = np.asarray([-5.14849798152926,-5.56365370199756,-16.7298857598484,-15.4815140695887],dtype=float)
    VFATGeometryDict['Short'][4]['y'] = np.asarray([44.326,58.383,58.383,44.326],dtype=float)
    VFATGeometryDict['Short'][5]['x'] = np.asarray([-5.56365370199756,-5.97880942246586,-17.9782574501081,-16.7298857598484],dtype=float)
    VFATGeometryDict['Short'][5]['y'] = np.asarray([58.383,72.44,72.44,58.383],dtype=float)
    VFATGeometryDict['Short'][6]['x'] = np.asarray([-5.97880942246586,-6.47069378786385,-19.4573518871341,-17.9782574501081],dtype=float)
    VFATGeometryDict['Short'][6]['y'] = np.asarray([72.4399999999999,89.0949999999999,89.0949999999999,72.4399999999999],dtype=float)
    VFATGeometryDict['Short'][7]['x'] = np.asarray([-6.47069378786385,-6.96257815326184,-20.9364463241602,-19.4573518871341],dtype=float)
    VFATGeometryDict['Short'][7]['y'] = np.asarray([89.0949999999999,105.75,105.75,89.0949999999999],dtype=float)
    VFATGeometryDict['Short'][8]['x'] = np.asarray([-3.83938562003834,-4.14080692494689,4.14080692494689,3.83938562003834],dtype=float)
    VFATGeometryDict['Short'][8]['y'] = np.asarray([0,10.206,10.206,0],dtype=float)
    VFATGeometryDict['Short'][9]['x'] = np.asarray([-4.14080692494689,-4.44222822985544,4.44222822985544,4.14080692494689],dtype=float)
    VFATGeometryDict['Short'][9]['y'] = np.asarray([10.206,20.412,20.412,10.206],dtype=float)
    VFATGeometryDict['Short'][10]['x'] = np.asarray([-4.44222822985544,-4.79536310569235,4.79536310569235,4.44222822985544],dtype=float)
    VFATGeometryDict['Short'][10]['y'] = np.asarray([20.412,32.369,32.369,20.412],dtype=float)
    VFATGeometryDict['Short'][11]['x'] = np.asarray([-4.79536310569235,-5.14849798152926,5.14849798152926,4.79536310569235],dtype=float)
    VFATGeometryDict['Short'][11]['y'] = np.asarray([32.369,44.326,44.326,32.369],dtype=float)
    VFATGeometryDict['Short'][12]['x'] = np.asarray([-5.14849798152926,-5.56365370199756,5.56365370199756,5.14849798152926],dtype=float)
    VFATGeometryDict['Short'][12]['y'] = np.asarray([44.326,58.383,58.383,44.326],dtype=float)
    VFATGeometryDict['Short'][13]['x'] = np.asarray([-5.56365370199756,-5.97880942246586,5.97880942246586,5.56365370199756],dtype=float)
    VFATGeometryDict['Short'][13]['y'] = np.asarray([58.383,72.44,72.44,58.383],dtype=float)
    VFATGeometryDict['Short'][14]['x'] = np.asarray([-5.97880942246586,-6.47069378786385,6.47069378786385,5.97880942246586],dtype=float)
    VFATGeometryDict['Short'][14]['y'] = np.asarray([72.4399999999999,89.0949999999999,89.0949999999999,72.4399999999999],dtype=float)
    VFATGeometryDict['Short'][15]['x'] = np.asarray([-6.47069378786385,-6.96257815326184,6.96257815326184,6.47069378786385],dtype=float)
    VFATGeometryDict['Short'][15]['y'] = np.asarray([89.0949999999999,105.75,105.75,89.0949999999999],dtype=float)
    VFATGeometryDict['Short'][16]['x'] = np.asarray([3.83938562003834,4.14080692494689,12.4513908518566,11.5450181214881],dtype=float)
    VFATGeometryDict['Short'][16]['y'] = np.asarray([0,10.206,10.206,0],dtype=float)
    VFATGeometryDict['Short'][17]['x'] = np.asarray([4.14080692494689,4.44222822985544,13.3577635822252,12.4513908518566],dtype=float)
    VFATGeometryDict['Short'][17]['y'] = np.asarray([10.206,20.412,20.412,10.206],dtype=float)
    VFATGeometryDict['Short'][18]['x'] = np.asarray([4.44222822985544,4.79536310569235,14.4196388259069,13.3577635822252],dtype=float)
    VFATGeometryDict['Short'][18]['y'] = np.asarray([20.412,32.369,32.369,20.412],dtype=float)
    VFATGeometryDict['Short'][19]['x'] = np.asarray([4.79536310569235,5.14849798152926,15.4815140695887,14.4196388259069],dtype=float)
    VFATGeometryDict['Short'][19]['y'] = np.asarray([32.369,44.326,44.326,32.369],dtype=float)
    VFATGeometryDict['Short'][20]['x'] = np.asarray([5.14849798152926,5.56365370199756,16.7298857598484,15.4815140695887],dtype=float)
    VFATGeometryDict['Short'][20]['y'] = np.asarray([44.326,58.383,58.383,44.326],dtype=float)
    VFATGeometryDict['Short'][21]['x'] = np.asarray([5.56365370199756,5.97880942246586,17.9782574501081,16.7298857598484],dtype=float)
    VFATGeometryDict['Short'][21]['y'] = np.asarray([58.383,72.44,72.44,58.383],dtype=float)
    VFATGeometryDict['Short'][22]['x'] = np.asarray([5.97880942246586,6.47069378786385,19.4573518871341,17.9782574501081],dtype=float)
    VFATGeometryDict['Short'][22]['y'] = np.asarray([72.4399999999999,89.0949999999999,89.0949999999999,72.4399999999999],dtype=float)
    VFATGeometryDict['Short'][23]['x'] = np.asarray([6.47069378786385,6.96257815326184,20.9364463241602,19.4573518871341],dtype=float)
    VFATGeometryDict['Short'][23]['y'] = np.asarray([89.0949999999999,105.75,105.75,89.0949999999999],dtype=float)
    VFATGeometryDict['Long'][0]['x'] = np.asarray([-3.83938562003834,-4.17332356777506,-12.5491682745625,-11.5450181214881],dtype=float)
    VFATGeometryDict['Long'][0]['y'] = np.asarray([0,11.307,11.307,0],dtype=float)
    VFATGeometryDict['Long'][1]['x'] = np.asarray([-4.17332356777506,-4.50726151551178,-13.5533184276368,-12.5491682745625],dtype=float)
    VFATGeometryDict['Long'][1]['y'] = np.asarray([11.307,22.614,22.614,11.307],dtype=float)
    VFATGeometryDict['Long'][2]['x'] = np.asarray([-4.50726151551178,-4.90466746092129,-14.7483166110425,-13.5533184276368],dtype=float)
    VFATGeometryDict['Long'][2]['y'] = np.asarray([22.614,36.07,36.07,22.614],dtype=float)
    VFATGeometryDict['Long'][3]['x'] = np.asarray([-4.90466746092129,-5.30207340633079,-15.9433147944483,-14.7483166110425],dtype=float)
    VFATGeometryDict['Long'][3]['y'] = np.asarray([36.07,49.526,49.526,36.07],dtype=float)
    VFATGeometryDict['Long'][4]['x'] = np.asarray([-5.30207340633079,-5.77777328465355,-17.3737425397006,-15.9433147944483],dtype=float)
    VFATGeometryDict['Long'][4]['y'] = np.asarray([49.526,65.633,65.633,49.526],dtype=float)
    VFATGeometryDict['Long'][5]['x'] = np.asarray([-5.77777328465355,-6.2534731629763,-18.804170284953,-17.3737425397006],dtype=float)
    VFATGeometryDict['Long'][5]['y'] = np.asarray([65.633,81.74,81.74,65.633],dtype=float)
    VFATGeometryDict['Long'][6]['x'] = np.asarray([-6.2534731629763,-6.82657530110586,-20.5274862591644,-18.804170284953],dtype=float)
    VFATGeometryDict['Long'][6]['y'] = np.asarray([81.74,101.145,101.145,81.74],dtype=float)
    VFATGeometryDict['Long'][7]['x'] = np.asarray([-6.82657530110586,-7.39967743923543,-22.2508022333757,-20.5274862591644],dtype=float)
    VFATGeometryDict['Long'][7]['y'] = np.asarray([101.145,120.55,120.55,101.145],dtype=float)
    VFATGeometryDict['Long'][8]['x'] = np.asarray([-3.83938562003834,-4.17332356777506,4.17332356777506,3.83938562003834],dtype=float)
    VFATGeometryDict['Long'][8]['y'] = np.asarray([0,11.307,11.307,0],dtype=float)
    VFATGeometryDict['Long'][9]['x'] = np.asarray([-4.17332356777506,-4.50726151551178,4.50726151551178,4.17332356777506],dtype=float)
    VFATGeometryDict['Long'][9]['y'] = np.asarray([11.307,22.614,22.614,11.307],dtype=float)
    VFATGeometryDict['Long'][10]['x'] = np.asarray([-4.50726151551178,-4.90466746092129,4.90466746092129,4.50726151551178],dtype=float)
    VFATGeometryDict['Long'][10]['y'] = np.asarray([22.614,36.07,36.07,22.614],dtype=float)
    VFATGeometryDict['Long'][11]['x'] = np.asarray([-4.90466746092129,-5.30207340633079,5.30207340633079,4.90466746092129],dtype=float)
    VFATGeometryDict['Long'][11]['y'] = np.asarray([36.07,49.526,49.526,36.07],dtype=float)
    VFATGeometryDict['Long'][12]['x'] = np.asarray([-5.30207340633079,-5.77777328465355,5.77777328465355,5.30207340633079],dtype=float)
    VFATGeometryDict['Long'][12]['y'] = np.asarray([49.526,65.633,65.633,49.526],dtype=float)
    VFATGeometryDict['Long'][13]['x'] = np.asarray([-5.77777328465355,-6.2534731629763,6.2534731629763,5.77777328465355],dtype=float)
    VFATGeometryDict['Long'][13]['y'] = np.asarray([65.633,81.74,81.74,65.633],dtype=float)
    VFATGeometryDict['Long'][14]['x'] = np.asarray([-6.2534731629763,-6.82657530110586,6.82657530110586,6.2534731629763],dtype=float)
    VFATGeometryDict['Long'][14]['y'] = np.asarray([81.74,101.145,101.145,81.74],dtype=float)
    VFATGeometryDict['Long'][15]['x'] = np.asarray([-6.82657530110586,-7.39967743923543,7.39967743923543,6.82657530110586],dtype=float)
    VFATGeometryDict['Long'][15]['y'] = np.asarray([101.145,120.55,120.55,101.145],dtype=float)
    VFATGeometryDict['Long'][16]['x'] = np.asarray([3.83938562003834,4.17332356777506,12.5491682745625,11.5450181214881],dtype=float)
    VFATGeometryDict['Long'][16]['y'] = np.asarray([0,11.307,11.307,0],dtype=float)
    VFATGeometryDict['Long'][17]['x'] = np.asarray([4.17332356777506,4.50726151551178,13.5533184276368,12.5491682745625],dtype=float)
    VFATGeometryDict['Long'][17]['y'] = np.asarray([11.307,22.614,22.614,11.307],dtype=float)
    VFATGeometryDict['Long'][18]['x'] = np.asarray([4.50726151551178,4.90466746092129,14.7483166110425,13.5533184276368],dtype=float)
    VFATGeometryDict['Long'][18]['y'] = np.asarray([22.614,36.07,36.07,22.614],dtype=float)
    VFATGeometryDict['Long'][19]['x'] = np.asarray([4.90466746092129,5.30207340633079,15.9433147944483,14.7483166110425],dtype=float)
    VFATGeometryDict['Long'][19]['y'] = np.asarray([36.07,49.526,49.526,36.07],dtype=float)
    VFATGeometryDict['Long'][20]['x'] = np.asarray([5.30207340633079,5.77777328465355,17.3737425397006,15.9433147944483],dtype=float)
    VFATGeometryDict['Long'][20]['y'] = np.asarray([49.526,65.633,65.633,49.526],dtype=float)
    VFATGeometryDict['Long'][21]['x'] = np.asarray([5.77777328465355,6.2534731629763,18.804170284953,17.3737425397006],dtype=float)
    VFATGeometryDict['Long'][21]['y'] = np.asarray([65.633,81.74,81.74,65.633],dtype=float)
    VFATGeometryDict['Long'][22]['x'] = np.asarray([6.2534731629763,6.82657530110586,20.5274862591644,18.804170284953],dtype=float)
    VFATGeometryDict['Long'][22]['y'] = np.asarray([81.74,101.145,101.145,81.74],dtype=float)
    VFATGeometryDict['Long'][23]['x'] = np.asarray([6.82657530110586,7.39967743923543,22.2508022333757,20.5274862591644],dtype=float)
    VFATGeometryDict['Long'][23]['y'] = np.asarray([101.145,120.55,120.55,101.145],dtype=float)

    ## Additional polygon to standardize Axis ranges
    shape_x = np.asarray([-120.55/2,-120.55/2,120.55/2,120.55/2],dtype=float)
    shape_y = np.asarray([125.,125.01,125.01,125.],dtype=float)
    ## END OF GEOMETRY DEFINITION
    h2p = ROOT.TH2Poly()
    if efficiency_target in ["ML1","ML2","PL1","PL2"]:
        endcapTag = efficiency_target
        for chamber_number in range(1,37):
            chamberSize = "Short" if chamber_number%2 == 1 else "Long"
            
            region = 1 if endcapTag[0]=="P" else -1
            chamber_ID = ReChLa2chamberName(region,chamber_number,endcapTag[-1])           
            
            angle_deg = -90 + (chamber_number - 1)*10## degrees
            angle_rad = np.radians(angle_deg)
            for VFATN in range(24):

                original_x = VFATGeometryDict[chamberSize][VFATN]['x']
                original_y = VFATGeometryDict[chamberSize][VFATN]['y']

                ## fix needed to respect GE11 installation, long SCs have Drift facing the Interaction Point. Short SCs have ReadOut Board facing the IP
                ## Applying reflection along y axis
                if (region == 1 and chamber_number % 2 == 0) or (region == -1 and chamber_number % 2 == 1):
                    original_x = -original_x

                translated_x = original_x
                translated_y = original_y + 130 ## GE11 ~ R = 130 cm                
                rotated_and_translated_x = np.asarray([ translated_x[i]*np.cos(angle_rad) - translated_y[i]*np.sin(angle_rad) for i in range(len(translated_x)) ],dtype=float)
                rotated_and_translated_y = np.asarray([ translated_x[i]*np.sin(angle_rad) + translated_y[i]*np.cos(angle_rad) for i in range(len(translated_x)) ],dtype=float)

                VFAT_den = EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['den'] 
                VFAT_num = EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['num'] 
                
                if VFAT_den == 0:
                    efficiency = 0 ## Bin has no text in case of no prop hit
                elif VFAT_num == 0:
                    efficiency = 0.01 ## Avoid empty label when printing it with text
                else:
                    efficiency = round(float(VFAT_num)/VFAT_den,2)


                selected_bin = h2p.AddBin(4,rotated_and_translated_x,rotated_and_translated_y)
                h2p.SetBinContent(selected_bin,efficiency)
                h2p.SetMarkerSize(.3) ## Reduces the text size when drawing the option "COLZ TEXT"
                
                h2p.GetXaxis().SetTitle("Glb x (cm)")
                h2p.GetYaxis().SetTitle("Glb y (cm)")
    else:
        chamber_ID = efficiency_target
        re_ch_la_list = chamberName2ReChLa(chamber_ID)
        re = re_ch_la_list[0]
        ch = re_ch_la_list[1]
        la = re_ch_la_list[2]
        chamberSize = "Short" if ch%2 == 1 else "Long"
        endcapTag = EndcapLayer2label(re,la)
        

        for VFATN in range(24):
            selected_bin = h2p.AddBin(4,VFATGeometryDict[chamberSize][VFATN]['x'],VFATGeometryDict[chamberSize][VFATN]['y'])
                
            VFAT_den = EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['den'] 
            VFAT_num = EfficiencyDictVFAT[endcapTag][chamber_ID][VFATN]['num'] 
                
            if VFAT_den == 0:
                efficiency = 0 ## Bin has no text in case of no prop hit
            elif VFAT_num == 0:
                    efficiency = 0.01 ## Avoid empty label when printing it with text
            else:
                efficiency = round(float(VFAT_num)/VFAT_den,2)

            h2p.SetBinContent(selected_bin,efficiency)

        h2p.AddBin(4,shape_x,shape_y)

    # h2p.SetMaximum(1.1)
    h2p.SetMinimum(0)
    h2p.SetTitle(efficiency_target+"_2DEffVFAT")
    h2p.SetName(efficiency_target+"_2DEffVFAT")
    h2p.GetYaxis().SetTickLength(0.005)
    h2p.GetYaxis().SetTitleOffset(1.3)
    h2p.GetXaxis().SetTickLength(0.005)
    h2p.GetYaxis().SetLabelSize(0.03)
    h2p.GetXaxis().SetLabelSize(0.03)
    h2p.SetStats(False)
    return h2p


def generateEfficiencyPlotbyPt(sourceDict,input_region=[-1,1],input_layer=[1,2]):
    ## Transforming in list
    reg_tag_string = "All" if isinstance(input_region, list) else "P" if input_region == 1 else "M"
    lay_tag_string = "" if isinstance(input_layer, list) else "L1" if input_layer == 1 else "L2"
    input_region = [input_region] if not isinstance(input_region, list) else input_region
    input_layer = [input_layer] if not isinstance(input_layer, list) else input_layer
    
    TH1F_TempContainer = {}
    Plot_Container = {}
    ColorAssociation = {'All':ROOT.kBlack,'Long':ROOT.kGreen+1,'Short':ROOT.kRed}
    
    title = "GE11"+reg_tag_string+lay_tag_string+"_EffbyPt"
    for chambers in ['All','Long','Short']:   
        Plot_Container[chambers] = ROOT.TGraphAsymmErrors()
        Plot_Container[chambers].GetXaxis().SetTitle("pt (GeV)")
        Plot_Container[chambers].GetYaxis().SetTitle("Efficiency")
        Plot_Container[chambers].SetMaximum(1.1)
        Plot_Container[chambers].SetMinimum(0)
        Plot_Container[chambers].SetTitle(chambers+"_"+title)
        Plot_Container[chambers].SetName(chambers+"_"+title)
        Plot_Container[chambers].SetLineColor(ColorAssociation[chambers])
        Plot_Container[chambers].SetMarkerColor(ColorAssociation[chambers])
        Plot_Container[chambers].SetFillColorAlpha(ColorAssociation[chambers],.4)
        Plot_Container[chambers].SetMarkerStyle(20)
        Plot_Container[chambers].SetMarkerSize(.8)

        TH1F_TempContainer.setdefault(chambers,{'num':ROOT.TH1F('num_'+chambers+title, title,11,0,110),'den':ROOT.TH1F('den_'+chambers+title, title,11,0,110)})

    for pt in range(0,11):
        
        TH1F_TempContainer['All']['num'].SetBinContent(pt+1, TH1F_TempContainer['All']['num'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['num'] for etaPartitionID in sourceDict.keys()]))
        TH1F_TempContainer['All']['den'].SetBinContent(pt+1, TH1F_TempContainer['All']['den'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['den'] for etaPartitionID in sourceDict.keys()]))
        
        long_chambers_etaPartitionID = [etaPartitionID for etaPartitionID in sourceDict.keys() if getInfoFromEtaID(etaPartitionID)[1] % 2 == 0]
        short_chambers_etaPartitionID = [etaPartitionID for etaPartitionID in sourceDict.keys() if getInfoFromEtaID(etaPartitionID)[1] % 2 == 1]
        
        TH1F_TempContainer['Long']['num'].SetBinContent(pt+1, TH1F_TempContainer['Long']['num'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['num'] for etaPartitionID in long_chambers_etaPartitionID]))
        TH1F_TempContainer['Long']['den'].SetBinContent(pt+1, TH1F_TempContainer['Long']['den'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['den'] for etaPartitionID in long_chambers_etaPartitionID]))

        TH1F_TempContainer['Short']['num'].SetBinContent(pt+1, TH1F_TempContainer['Short']['num'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['num'] for etaPartitionID in short_chambers_etaPartitionID]))
        TH1F_TempContainer['Short']['den'].SetBinContent(pt+1, TH1F_TempContainer['Short']['den'].GetBinContent(pt) + sum([sourceDict[etaPartitionID][pt]['den'] for etaPartitionID in short_chambers_etaPartitionID]))
    
    
    for chambers in ['All','Long','Short']:
        Plot_Container[chambers].Divide(TH1F_TempContainer[chambers]['num'],TH1F_TempContainer[chambers]['den'],"B")

    return Plot_Container['Short'],Plot_Container['Long'],Plot_Container['All']


def generateEfficiencyDistribution(sourceDict):
    EfficiencyDistribution = ROOT.TH1F("EfficiencyDistribution","EfficiencyDistribution",100,0,1.)

    Num = {}
    Den = {}
    Eff = {}
    
    for region in [-1,1]:
        for chamber in range(1,37):
            for layer in [1,2]:
                key = region*(100*chamber+10*layer)
                
                Num[key] = 0
                Den[key] = 0
                Eff[key] = -9

    for etaPartitionID,value in sourceDict.items():
        region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)

        key = region*(100*chamber+10*layer)

        Num[key] += sum([value[k]['num'] for k in value.keys()])
        Den[key] += sum([value[k]['den'] for k in value.keys()])

    for k in Num.keys():
        if (Den[k] != 0):
            Eff[k] = float(Num[k])/float(Den[k])
        EfficiencyDistribution.Fill(Eff[k])
    OverFlow = EfficiencyDistribution.GetBinContent(101)
    UnderFlow = EfficiencyDistribution.GetBinContent(0)

    return EfficiencyDistribution


            

def generateEfficiencyPlot2DGE11(sourceDict,input_region=1,input_layer=1,debug=False):
    ## Transforming in list
    reg_tag_string = "All" if isinstance(input_region, list) else "P" if input_region == 1 else "M"
    lay_tag_string = "" if isinstance(input_layer, list) else "L1" if input_layer == 1 else "L2"
    input_region = [input_region] if not isinstance(input_region, list) else input_region
    input_layer = [input_layer] if not isinstance(input_layer, list) else input_layer

    title = "GE11"+reg_tag_string+lay_tag_string
    phi_nbins,phi_min,phi_max = 36,0.5,36.5
    etaP_nbins, etaP_min,etaP_max = 10, -0.5, 9.5    

    EfficiencyTH2D = ROOT.TH2F(title+"_ChambersEfficiency", title+"_ChambersEfficiency", phi_nbins,phi_min,phi_max,etaP_nbins,etaP_min,etaP_max)
    NumTH2D = ROOT.TH2F(title+"_ChambersNum", title+"_ChambersNum", phi_nbins,phi_min,phi_max,etaP_nbins,etaP_min,etaP_max)
    DenTH2D = ROOT.TH2F(title+"_ChambersDen", title+"_ChambersDen", phi_nbins,phi_min,phi_max,etaP_nbins,etaP_min,etaP_max)


    Summary = ROOT.TGraphAsymmErrors()
    Summary.GetXaxis().SetTitle("Chamber Number")
    Summary.GetXaxis().SetTitleSize(0.05)
    Summary.GetYaxis().SetTitle("Efficiency")
    Summary.GetYaxis().SetTitleSize(0.05)
    Summary.SetMaximum(1.1)
    Summary.SetMinimum(0.)
    Summary.SetTitle(title+"_Efficiency")
    Summary.SetName(title+"_Efficiency")
    Summary.SetLineColor(ROOT.kBlack)
    Summary.SetMarkerColor(ROOT.kBlack)
    Summary.SetFillColorAlpha(ROOT.kBlack,.4)
    Summary.SetMarkerStyle(20)
    Summary.SetMarkerSize(.8)

    SummaryNum = ROOT.TH1F("n","n",phi_nbins,phi_min,phi_max)
    SummaryDen = ROOT.TH1F("d","d",phi_nbins,phi_min,phi_max)
    
    N = 1                                                                        
    for etaPartitionID,value in sourceDict.items():
        region,chamber,layer,eta = getInfoFromEtaID(etaPartitionID)
        
        if layer not in input_layer or region not in input_region:
            continue

        etaPartitionRecHits = sum([value[k]['num'] for k in value.keys()])
        etaPartitionPropHits = sum([value[k]['den'] for k in value.keys()])
        
        try:
            eta_efficiency = round(float(etaPartitionRecHits)/float(etaPartitionPropHits),2)
        except:
            if debug:print "Warning on Re,Ch,La,etaP = ", region,chamber,layer,eta, "\tDenominator is 0"
            eta_efficiency = 0
            N+=1

        binx = chamber
        biny = eta  + 1
        
        
        existingNumSummary = SummaryNum.GetBinContent(binx)
        existingDenSummary = SummaryDen.GetBinContent(binx)
        existingNumTH2D = NumTH2D.GetBinContent(binx,biny)
        existingDenTH2D = DenTH2D.GetBinContent(binx,biny)
        

        SummaryNum.SetBinContent(binx,existingNumSummary+etaPartitionRecHits)
        SummaryDen.SetBinContent(binx,existingDenSummary+etaPartitionPropHits)

        NumTH2D.SetBinContent(binx,biny,existingNumTH2D+etaPartitionRecHits)
        DenTH2D.SetBinContent(binx,biny,existingDenTH2D+etaPartitionPropHits)
        #EfficiencyTH2D.SetBinContent(binx,biny,eta_efficiency)

    EfficiencyTH2D = NumTH2D.Clone()
    EfficiencyTH2D.SetTitle(title+"_ChambersEfficiency")
    EfficiencyTH2D.SetName(title+"_ChambersEfficiency")
    EfficiencyTH2D.Divide(DenTH2D)
    for x in range(1,phi_nbins+1):
        for y in range(1,etaP_nbins+1):
            EfficiencyTH2D.SetBinContent(x,y,  round(EfficiencyTH2D.GetBinContent(x,y),2) )


    EfficiencyTH2D.SetStats(False)
    EfficiencyTH2D.GetYaxis().SetTickLength(0.005)
    EfficiencyTH2D.GetXaxis().SetTitle("Chamber number")
    EfficiencyTH2D.GetXaxis().SetTitleSize(0.05)
    EfficiencyTH2D.GetYaxis().SetTitle("GEM EtaPartition")
    EfficiencyTH2D.GetYaxis().SetTitleSize(0.05)

    NumTH2D.SetStats(False)
    NumTH2D.GetYaxis().SetTickLength(0.005)
    NumTH2D.GetXaxis().SetTitle("Chamber number")
    NumTH2D.GetXaxis().SetTitleSize(0.05)
    NumTH2D.GetYaxis().SetTitle("GEM EtaPartition")
    NumTH2D.GetYaxis().SetTitleSize(0.05)

    DenTH2D.SetStats(False)
    DenTH2D.GetYaxis().SetTickLength(0.005)
    DenTH2D.GetXaxis().SetTitle("Chamber number")
    DenTH2D.GetXaxis().SetTitleSize(0.05)
    DenTH2D.GetYaxis().SetTitle("GEM EtaPartition")
    DenTH2D.GetYaxis().SetTitleSize(0.05)

    Summary.Divide(SummaryNum,SummaryDen,"B")
    Summary.GetXaxis().SetRangeUser(0.,36.5)
    Summary.GetYaxis().SetTickLength(0.005)
    Summary.GetXaxis().SetTickLength(0.005)
    Summary.SetMarkerColor(ROOT.kBlue)
    Summary.SetLineColor(ROOT.kBlue)
    return EfficiencyTH2D,NumTH2D,DenTH2D,Summary

def pt_index(num):
    if num <10:
        index = 0
    elif num <20:
        index = 1
    elif num <30:
        index = 2
    elif num <40:
        index = 3
    elif num <50:
        index = 4
    elif num <60:
        index = 5
    elif num <70:
        index = 6
    elif num <80:
        index = 7
    elif num <90:
        index = 8
    elif num <100:
        index = 9
    else:
        index = 10
    
    return index

## returns a dict containing the strip pitch for a given ch,etaP
def GetStripGeometry():
    stripGeometryDict = {}
    df = pd.read_csv("/afs/cern.ch/user/f/fivone/Documents/myLIB/GE11Geometry/GE11StripSpecs.csv")
        
    for ch in range(1,37):
        stripGeometryDict[ch] = {}
        for etaP in range(1,9):
            if ch%2 == 0:
                chID = "Long"
            if ch%2 == 1:
                chID = "Short"
            df_temp = df.loc[df['Chamber'] == chID]
            df_temp = df_temp.loc[df_temp['EtaP'] == etaP]
            firstStripPosition = df_temp["Loc_x"].min()
            lastStripPosition = df_temp["Loc_x"].max()
            stripPitch = (lastStripPosition - firstStripPosition)/float(len(df_temp["Loc_x"]))
            stripGeometryDict[ch].setdefault(etaP,{'stripPitch':stripPitch,'firstStrip':firstStripPosition,"lastStrip":lastStripPosition})

    return stripGeometryDict


def printSummary(sourceDict,matching_variables,ResidualCutOff,matching_variable_units,debug=False):
    for matching_variable in matching_variables:
        print "\n\n#############\nSUMMARY\n#############\nMatchingVariable = "+matching_variable+"\nCutoff = ",ResidualCutOff[matching_variable],matching_variable_units[matching_variable],"\n"
        for eta in range(1,9):
            num = sum([sourceDict[matching_variable][key][k]['num'] for key in sourceDict[matching_variable].keys() for k in range(0,11) if abs(key)%10 == eta])
            den = sum([sourceDict[matching_variable][key][k]['den'] for key in sourceDict[matching_variable].keys() for k in range(0,11) if abs(key)%10 == eta])
            try:
                print "Efficiency GE11 ETA"+str(eta)+"  ==> ", num , "/",den, "\t = ", round(float(num)/float(den),3)
            except:
                print "EtaP = ",str(eta)," has no propagated hits..."

        num = sum([sourceDict[matching_variable][key][k]['num'] for key in sourceDict[matching_variable].keys() for k in range(0,11)])
        den = sum([sourceDict[matching_variable][key][k]['den'] for key in sourceDict[matching_variable].keys() for k in range(0,11)])
        try:
            print "Efficiency GE11       ==> ", num , "/",den, "\t = ", round(float(num)/float(den),3)
        except:
            if debug: print "WARNING"
        print "#############"

def initializeMatchingTree():
    generalMatching = {}
    recHit_Matching = {}
    propHit_Matching = {}
    
    generalMatching['event'] = array('i', [0])
    generalMatching['lumiblock'] = array('i', [0])
    generalMatching['chamber'] = array('i', [0])
    generalMatching['region'] = array('i', [0])
    generalMatching['layer'] = array('i', [0])
    generalMatching['etaPartition'] = array('i', [0])

    recHit_Matching['loc_x'] = array('f', [0.])
    recHit_Matching['glb_x'] = array('f', [0.])
    recHit_Matching['glb_y'] = array('f', [0.])
    recHit_Matching['glb_z'] = array('f', [0.])
    recHit_Matching['glb_r'] = array('f', [0.])
    recHit_Matching['glb_phi'] = array('f', [0.])
    recHit_Matching['firstStrip'] = array('i', [0])
    recHit_Matching['cluster_size'] = array('i', [0])

    propHit_Matching['loc_x'] = array('f', [0.])
    propHit_Matching['loc_y'] = array('f', [0.])
    propHit_Matching['glb_x'] = array('f', [0.])
    propHit_Matching['glb_y'] = array('f', [0.])
    propHit_Matching['glb_z'] = array('f', [0.])
    propHit_Matching['glb_r'] = array('f', [0.])
    propHit_Matching['glb_phi'] = array('f', [0.])
    propHit_Matching['err_glb_r'] = array('f', [0.])
    propHit_Matching['err_glb_phi'] = array('f', [0.])
    propHit_Matching['pt'] = array('f', [0.])
    propHit_Matching['mu_propagated_isGEM'] = array('i', [0])
    propHit_Matching['STA_Normchi2'] = array('f', [0.])
    propHit_Matching['nME1Hits'] = array('i', [0])
    propHit_Matching['nME2Hits'] = array('i', [0])
    propHit_Matching['nME3Hits'] = array('i', [0])
    propHit_Matching['nME4Hits'] = array('i', [0])

    t_out = ROOT.TTree("MatchingTree", "MatchingTree")
    for key in generalMatching.keys():
        if type(generalMatching[key][0]) == float: label="F"
        elif type(generalMatching[key][0]) == int: label="I"
        t_out.Branch(key, generalMatching[key], key+"/"+label)
    for key in recHit_Matching.keys():
        if type(recHit_Matching[key][0]) == float: label="F"
        elif type(recHit_Matching[key][0]) == int: label="I"
        t_out.Branch("recHit_"+key, recHit_Matching[key], "recHit_"+key+"/"+label)
    for key in propHit_Matching.keys():
        if type(propHit_Matching[key][0]) == float: label="F"
        elif type(propHit_Matching[key][0]) == int: label="I"
        t_out.Branch("propHit_"+key, propHit_Matching[key], "propHit_"+key+"/"+label)
    
    return t_out,generalMatching,recHit_Matching,propHit_Matching

def fillMatchingTreeArray(PropHitonEta,prop_hit_index,RecHitonEta,reco_hit_index,propHit_Matching,recHit_Matching):
    for key in recHit_Matching.keys():
        recHit_Matching[key][0] =RecHitonEta[key][reco_hit_index]
    for key in propHit_Matching.keys():
        propHit_Matching[key][0] =PropHitonEta[key][prop_hit_index]
    
    return recHit_Matching,propHit_Matching

def store4evtDspl(name,run,lumi,evt):
    evtDspl_dir = "/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/EvtDspl/"
    with open(evtDspl_dir+name+".txt", 'a+') as f:
        f.write(str(run)+":"+str(lumi)+":"+str(evt)+"\n")

def printEfficiencyFromCSV(path):
    file_path = path + "/MatchingSummary_glb_rdphi.csv"
    print file_path,"\n"
    print "Label","\t","Analyzed Chambers","\t","AVG pt","\t","Matched","\t","Propagated","\t","Eff. 68%CL"

    df = pd.read_csv(file_path, sep=',')
    
    Matched = df["matchedRecHit"].sum()
    Propagated = df["propHit"].sum()
    if "AVG_pt" in df.columns: 
        AVGpt  =  round(df["AVG_pt"].mean(),4)
    else:
        AVGpt = 0
    analyzed_Chambers = float(len(df))/8
    
    for eta in range(1,9):
        df_temp = df[ df["etaPartition"] == eta ]
        etaMatched = df_temp["matchedRecHit"].sum()
        etaProp = df_temp["propHit"].sum()
        if "AVG_pt" in df.columns: 
            AVGpteta  = round(df_temp["AVG_pt"].mean(),4)
        else:
            AVGpteta  = 0


        # print "Eta",eta,"\t",AVGpteta,"\t",etaMatched,"\t",etaProp,"\t",generateClopperPeasrsonInterval(etaMatched,etaProp)
    print "GE11","\t",analyzed_Chambers,"\t", AVGpt,"\t",Matched,"\t",Propagated,"\t",generateClopperPeasrsonInterval(Matched,Propagated)
    print "##############\n"  

if __name__ == '__main__':
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_20/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_15/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_10/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_8/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_6/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_4/")
    printEfficiencyFromCSV("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/Output/PFA_Analyzer_Output/CSV/MergedRuns_Prompt_700uA_100HZ_Trim_chi2_2/")

    pass
