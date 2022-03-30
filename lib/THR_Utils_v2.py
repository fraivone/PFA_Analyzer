'''
Version for 2022 data
'''
import pandas as pd
import numpy as np
import sys
import os
import json
 
sys.path.append(os.path.abspath("/eos/project-c/cmsgemonline/public/"))
from chamber_mapping import chamber_mapping as mapping

TRIMDAC2fC = 2.1/63     # VFAT Design --> 63 TRIM DAC = 2.1 fC
THRDAC2fC = 12. / 100   # VFAT Design --> 100 THR DAC = 12 fC

## inverting the map so that the chamber name is the key
mapping = {v: k for k, v in mapping.items()}

## 
def GetCHDict(chamber_ID,run=348832):
        tupl = mapping[chamber_ID]
        crate = int(tupl[0])
        amc = int(tupl[1])
        OHLink = int(tupl[2])
        fed = 1466 + crate
        
        file_path = "/eos/project-c/cmsgemonline/public/runs/run"+str(run)+"/fed%d" % fed + "-amc%02d" % amc +"_ConfigInfo.json"

        try:
                with open(file_path) as json_file:
                        data_dict = json.load(json_file)

        except:
                print "Can't open ",file_path
                print "Exiting..."
                sys.exit(0)
        ## check exsistance off all keys
        try:
                data_dict = data_dict['fed'][str(fed)]['slot'][str(amc)]['link'][str(OHLink)]
        except:
                data_dict = {}
        return data_dict

## Returns the THR_DAC set for a given VFAT in a Chamber
## run can be provided
def GetVFAT_THRDAC(chamber_ID,VFATN,run=348832):

        data = GetCHDict(chamber_ID,run)

        try:
                threshold = data['vfat'][str(VFATN)]['THRESHOLD_DAC']
                latency = data['vfat'][str(VFATN)]['LATENCY']
        except:
                threshold = -1
                latency = -1

        return threshold

def GetVFAT_LATENCY(chamber_ID,VFATN,run=348832):

        data = GetCHDict(chamber_ID,run)

        try:
                threshold = data['vfat'][str(VFATN)]['THRESHOLD_DAC']
                latency = data['vfat'][str(VFATN)]['LATENCY']
        except:
                threshold = -1
                latency = -1

        return latency

## Returns the AVG THR of a chamber       
def GetOverallChamberThreshold(chamberID,run=348832):

        if chamberID not in mapping.keys():
                print "Invalid chamberID",chamberID
                print "Exiting..."
                sys.exit(0)
        

        data = GetCHDict(chamberID,run)

        temp_array = np.empty(0,dtype=float)
    
        for VFAT_N in range(24):
                # THR data
                try:
                        threshold = data['vfat'][str(VFAT_N)]['THRESHOLD_DAC']
                        temp_array = np.append(temp_array,threshold)
                except:
                        print "Can't find THR for\t",chamberID,"\tVFAT:",VFAT_N,"\tSkipping...\n"
                
        # if overall_vfat_trheshold == 10**6:
        #         continue
        # elif abs(overall_vfat_trheshold) > 255:
        #     if verbose: print "No-sense DAC value for VFAT",VFAT_N," on chamber ",chamberID,"\nPlease check, skipping VFAT value..."
        #     continue
        # else:


        if len(temp_array) != 0:
                return np.mean(temp_array)
        else:
                return None
        

if __name__ == "__main__":
        counter = 0
        chamber = "GE11-P-04L2-L"
        for VFATN in range(24):
                th = GetVFAT_THRDAC(chamber,VFATN,348773)
                lt = GetVFAT_LATENCY(chamber,VFATN)
                print chamber,VFATN, th, lt

