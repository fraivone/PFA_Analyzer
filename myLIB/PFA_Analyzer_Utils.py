import ROOT
import numpy
import pandas as pd
from ROOT_Utils import *



def getInfoFromEtaID(id):
    etaPartition = abs(id)%10
    layer = ((abs(id)-etaPartition)%100)/10
    chamber = (abs(id)-etaPartition-10*layer)/100
    region = id/abs(id)
    return region,chamber,layer,etaPartition

def chamberName2ReChLa(chamberName):
    re = -1 if "M" in chamberName else 1
    ch = int( chamberName.split("-")[-1][:2] )
    la = int( chamberName.split("-")[-1][-1] )

    return [re,ch,la]

def ReChLa2chamberName(re,ch,la):
    endcap = "M" if re == -1 else "P"
    size = "S" if ch%2 == 1 else "L"
    chID = 'GE11-'+endcap+'-%02d' % ch +"L"+str(la)+"-"+size 

def ChambersOFFHisto(chamberNumberOFF):
    GE11_OFF = {-1:{},1:{}}
    GE11_OFF[-1][1] = ROOT.TH1F("GE11-M-L1_isAnalyzed","GE11-M-L1_isAnalyzed",36,0.5,36.5)
    GE11_OFF[-1][2] = ROOT.TH1F("GE11-M-L2_isAnalyzed","GE11-M-L2_isAnalyzed",36,0.5,36.5)
    GE11_OFF[1][1] = ROOT.TH1F("GE11-P-L1_isAnalyzed","GE11-P-L_isAnalyzed",36,0.5,36.5)
    GE11_OFF[1][2] = ROOT.TH1F("GE11-P-L2_isAnalyzed","GE11-P-L2_isAnalyzed",36,0.5,36.5)

    for region in [-1,1]:
        for layer in [1,2]:
            for chamber in range(1,37):
                if [region,chamber,layer] in chamberNumberOFF:
                    GE11_OFF[region][layer].SetBinContent(chamber,0)
                else:
                    GE11_OFF[region][layer].SetBinContent(chamber,1)
    for key_1 in GE11_OFF.keys():
        for key_2 in GE11_OFF[key_1].keys():
            GE11_OFF[key_1][key_2].GetYaxis().SetTickLength(0.005)
            GE11_OFF[key_1][key_2].SetStats(False)

    return GE11_OFF[-1][1],GE11_OFF[-1][2],GE11_OFF[1][1],GE11_OFF[1][2]

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

def fillPlot2DResidualContainer(TH2Fresidual_collector,matching_variables,nbins):
    
    for key_1 in matching_variables:
        for key_2 in ['all','short','long']:
            for key_3 in ['glb_phi','glb_rdphi']:
                for x_bin in range(1,nbins+1):
                    for y_bin in range(1,nbins+1):
                        AVG_Residual = TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][1]/TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][0] if TH2Fresidual_collector[key_1][key_2][key_3][x_bin][y_bin][0]!=0 else 0
                        TH2Fresidual_collector[key_1][key_2][key_3]['TH2F'].SetBinContent(x_bin,y_bin,AVG_Residual)
    return TH2Fresidual_collector

def passCut(PropHitonEta,prop_hit_index,maxPropR_Err=0.7,maxPropPhi_Err=0.001,fiducialCutR=0.5,fiducialCutPhi=0.002,minPt=0.,maxChi2=9999999):
    passedCut = True
    if PropHitonEta['err_glb_phi'][prop_hit_index] > maxPropPhi_Err:
        passedCut = False
    if PropHitonEta['err_glb_r'][prop_hit_index] > maxPropR_Err:
        passedCut = False

    if PropHitonEta['STA_Normchi2'][prop_hit_index] > maxChi2 or PropHitonEta['STA_Normchi2'][prop_hit_index] < 0.7:
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

        TH1F_TempContainer.setdefault(chambers,{'num':ROOT.TH1F('num_'+chambers+title, title,8,0,8),'den':ROOT.TH1F('den_'+chambers+title, title,8,0,8)})

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
        Plot_Container[chambers].Divide(TH1F_TempContainer[chambers]['num'],TH1F_TempContainer[chambers]['den'])

    return Plot_Container['Short'],Plot_Container['Long'],Plot_Container['All']

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
        Plot_Container[chambers].Divide(TH1F_TempContainer[chambers]['num'],TH1F_TempContainer[chambers]['den'])

    return Plot_Container['Short'],Plot_Container['Long'],Plot_Container['All']


def generateEfficiencyDistribution(sourceDict):
    EfficiencyDistribution = ROOT.TH1F("EfficiencyDistribution","EfficiencyDistribution",100,0.,1.)

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
    
    return EfficiencyDistribution


            

def generateEfficiencyPlot2DGE11(sourceDict,input_region=1,input_layer=1):
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
            print "Warning on Re,Ch,La,etaP = ", region,chamber,layer,eta, "\tDenominator is 0"
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

    Summary.Divide(SummaryNum,SummaryDen)
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


def printSummary(sourceDict,matching_variables,ResidualCutOff,matching_variable_units):
    for matching_variable in matching_variables:
        print "\n\n#############\nSUMMARY\n#############\nMatchingVariable = "+matching_variable+"\nCutoff = ",ResidualCutOff[matching_variable],matching_variable_units[matching_variable],"\n"
        for eta in range(1,9):
            num = sum([sourceDict[matching_variable][key][k]['num'] for key in sourceDict[matching_variable].keys() for k in range(0,11) if abs(key)%10 == eta])
            den = sum([sourceDict[matching_variable][key][k]['den'] for key in sourceDict[matching_variable].keys() for k in range(0,11) if abs(key)%10 == eta])
            print "Efficiency GE11 ETA"+str(eta)+"  ==> ", num , "/",den, " = ", float(num)/float(den)

        num = sum([sourceDict[matching_variable][key][k]['num'] for key in sourceDict[matching_variable].keys() for k in range(0,11)])
        den = sum([sourceDict[matching_variable][key][k]['den'] for key in sourceDict[matching_variable].keys() for k in range(0,11)])
        try:
            print "Efficiency GE11  ==> ", num , "/",den, " = ", float(num)/float(den)
        except:
            print "WARNING"
        print "#############"

