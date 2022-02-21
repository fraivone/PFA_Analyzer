import pandas as pd
import numpy as np
import sys

TRIMDAC2fC = 2.1/63     # VFAT Design --> 63 TRIM DAC = 2.1 fC
THRDAC2fC = 12. / 100   # VFAT Design --> 100 THR DAC = 12 fC


chamber_mapping = {
            "GE11-M-01L1-S":(1,1,0), 
            "GE11-M-01L2-S":(1,1,1), 
            "GE11-M-02L1-L":(1,1,2), 
            "GE11-M-02L2-L":(1,1,3), 
            "GE11-M-03L1-S":(1,1,4), 
            "GE11-M-03L2-S":(1,1,5), 
            "GE11-M-04L1-L":(1,1,6), 
            "GE11-M-04L2-L":(1,1,7), 
            "GE11-M-05L1-S":(1,1,8), 
            "GE11-M-05L2-S":(1,1,9), 
            "GE11-M-06L1-L":(1,1,10),
            "GE11-M-06L2-L":(1,1,11),
            "GE11-M-07L1-S":(1,3,0), 
            "GE11-M-07L2-S":(1,3,1), 
            "GE11-M-08L1-L":(1,3,2), 
            "GE11-M-08L2-L":(1,3,3), 
            "GE11-M-09L1-S":(1,3,4), 
            "GE11-M-09L2-S":(1,3,5), 
            "GE11-M-10L1-L":(1,3,6), 
            "GE11-M-10L2-L":(1,3,7), 
            "GE11-M-11L1-S":(1,3,8), 
            "GE11-M-11L2-S":(1,3,9), 
            "GE11-M-12L1-L":(1,3,10),
            "GE11-M-12L2-L":(1,3,11),
            "GE11-M-13L1-S":(1,5,0), 
            "GE11-M-13L2-S":(1,5,1), 
            "GE11-M-14L1-L":(1,5,2), 
            "GE11-M-14L2-L":(1,5,3), 
            "GE11-M-15L1-S":(1,5,4), 
            "GE11-M-15L2-S":(1,5,5), 
            "GE11-M-16L1-L":(1,5,6), 
            "GE11-M-16L2-L":(1,5,7), 
            "GE11-M-17L1-S":(1,5,8), 
            "GE11-M-17L2-S":(1,5,9), 
            "GE11-M-18L1-L":(1,5,10),
            "GE11-M-18L2-L":(1,5,11),
            "GE11-M-19L1-S":(1,7,0), 
            "GE11-M-19L2-S":(1,7,1), 
            "GE11-M-20L1-L":(1,7,2), 
            "GE11-M-20L2-L":(1,7,3), 
            "GE11-M-21L1-S":(1,7,4), 
            "GE11-M-21L2-S":(1,7,5), 
            "GE11-M-22L1-L":(1,7,6), 
            "GE11-M-22L2-L":(1,7,7), 
            "GE11-M-23L1-S":(1,7,8), 
            "GE11-M-23L2-S":(1,7,9), 
            "GE11-M-24L1-L":(1,7,10),
            "GE11-M-24L2-L":(1,7,11),
            "GE11-M-25L1-S":(1,9,0), 
            "GE11-M-25L2-S":(1,9,1), 
            "GE11-M-26L1-L":(1,9,2), 
            "GE11-M-26L2-L":(1,9,3), 
            "GE11-M-27L1-S":(1,9,4), 
            "GE11-M-27L2-S":(1,9,5), 
            "GE11-M-28L1-L":(1,9,6), 
            "GE11-M-28L2-L":(1,9,7), 
            "GE11-M-29L1-S":(1,9,8), 
            "GE11-M-29L2-S":(1,9,9), 
            "GE11-M-30L1-L":(1,9,10),
            "GE11-M-30L2-L":(1,9,11),
            "GE11-M-31L1-S":(1,11,0),
            "GE11-M-31L2-S":(1,11,1),
            "GE11-M-32L1-L":(1,11,2),
            "GE11-M-32L2-L":(1,11,3),
            "GE11-M-33L1-S":(1,11,4),
            "GE11-M-33L2-S":(1,11,5),
            "GE11-M-34L1-L":(1,11,6),
            "GE11-M-34L2-L":(1,11,7),
            "GE11-M-35L1-S":(1,11,8),
            "GE11-M-35L2-S":(1,11,9),
            "GE11-M-36L1-L":(1,11,10),
            "GE11-M-36L2-L":(1,11,11),
            "GE11-P-01L1-S":(2,1,0), 
            "GE11-P-01L2-S":(2,1,1), 
            "GE11-P-02L1-L":(2,1,2), 
            "GE11-P-02L2-L":(2,1,3), 
            "GE11-P-03L1-S":(2,1,4), 
            "GE11-P-03L2-S":(2,1,5), 
            "GE11-P-04L1-L":(2,1,6), 
            "GE11-P-04L2-L":(2,1,7), 
            "GE11-P-05L1-S":(2,1,8), 
            "GE11-P-05L2-S":(2,1,9), 
            "GE11-P-06L1-L":(2,1,10),
            "GE11-P-06L2-L":(2,1,11),
            "GE11-P-07L1-S":(2,3,0), 
            "GE11-P-07L2-S":(2,3,1), 
            "GE11-P-08L1-L":(2,3,2), 
            "GE11-P-08L2-L":(2,3,3), 
            "GE11-P-09L1-S":(2,3,4), 
            "GE11-P-09L2-S":(2,3,5), 
            "GE11-P-10L1-L":(2,3,6), 
            "GE11-P-10L2-L":(2,3,7), 
            "GE11-P-11L1-S":(2,3,8), 
            "GE11-P-11L2-S":(2,3,9), 
            "GE11-P-12L1-L":(2,3,10),
            "GE11-P-12L2-L":(2,3,11),
            "GE11-P-13L1-S":(2,5,0), 
            "GE11-P-13L2-S":(2,5,1), 
            "GE11-P-14L1-L":(2,5,2), 
            "GE11-P-14L2-L":(2,5,3), 
            "GE11-P-15L1-S":(2,5,4), 
            "GE11-P-15L2-S":(2,5,5), 
            "GE11-P-16L1-L":(2,5,6), 
            "GE11-P-16L2-L":(2,5,7), 
            "GE11-P-17L1-S":(2,5,8), 
            "GE11-P-17L2-S":(2,5,9), 
            "GE11-P-18L1-L":(2,5,10),
            "GE11-P-18L2-L":(2,5,11),
            "GE11-P-19L1-S":(2,7,0), 
            "GE11-P-19L2-S":(2,7,1), 
            "GE11-P-20L1-L":(2,7,2), 
            "GE11-P-20L2-L":(2,7,3), 
            "GE11-P-21L1-S":(2,7,4), 
            "GE11-P-21L2-S":(2,7,5), 
            "GE11-P-22L1-L":(2,7,6), 
            "GE11-P-22L2-L":(2,7,7), 
            "GE11-P-23L1-S":(2,7,8), 
            "GE11-P-23L2-S":(2,7,9), 
            "GE11-P-24L1-L":(2,7,10),
            "GE11-P-24L2-L":(2,7,11),
            "GE11-P-25L1-S":(2,9,0), 
            "GE11-P-25L2-S":(2,9,1), 
            "GE11-P-26L1-L":(2,9,2), 
            "GE11-P-26L2-L":(2,9,3), 
            "GE11-P-27L1-S":(2,9,4), 
            "GE11-P-27L2-S":(2,9,5), 
            "GE11-P-28L1-L":(2,9,6), 
            "GE11-P-28L2-L":(2,9,7), 
            "GE11-P-29L1-S":(2,9,8), 
            "GE11-P-29L2-S":(2,9,9), 
            "GE11-P-30L1-L":(2,9,10),
            "GE11-P-30L2-L":(2,9,11),
            "GE11-P-31L1-S":(2,11,0),
            "GE11-P-31L2-S":(2,11,1),
            "GE11-P-32L1-L":(2,11,2),
            "GE11-P-32L2-L":(2,11,3),
            "GE11-P-33L1-S":(2,11,4),
            "GE11-P-33L2-S":(2,11,5),
            "GE11-P-34L1-L":(2,11,6),
            "GE11-P-34L2-L":(2,11,7),
            "GE11-P-35L1-S":(2,11,8),
            "GE11-P-35L2-S":(2,11,9),
            "GE11-P-36L1-L":(2,11,10),
            "GE11-P-36L2-L":(2,11,11)
}

## Returns the THR_DAC set for a given VFAT in a Chamber
## THR_ARM_Folder better to be absolute
def GetVFAT_THRDAC(THR_ARM_folder,chamber_ID,VFATN):
        tupl = chamber_mapping[chamber_ID]
        crate = int(tupl[0])
        amc = int(tupl[1])
        OHLink = int(tupl[2])
        file_path = THR_ARM_folder+"/crate%02d" % crate + "-amc%02d" % amc +".txt"
        THR_Data = pd.read_csv(file_path, sep=':')
        Threshold_ARM_DAC = THR_Data[( (THR_Data["OH/I"] == OHLink) &  (THR_Data["vfatN/I"] == VFATN) )]["threshold/I"].values
        
        Threshold_ARM_DAC = Threshold_ARM_DAC[0] if len(Threshold_ARM_DAC) != 0 else None
        return Threshold_ARM_DAC

## Returns the overall THR for a given VFAT in a chamber, combining THR_DAC and TRIM_DAC (if the trimming was applied)        
def GetOverallVFATThreshold(source_folder,chamberID,VFAT_N,verbose=False):
    if chamberID not in chamber_mapping.keys():
        print "Invalid chamberID",chamberID
        print "Exiting..."
        sys.exit(0)
    if VFAT_N not in range(0,24):
        print "Invalid VFAT number",VFAT_N
        print "Exiting..."
        sys.exit(0)

    if "Trimming" in source_folder:
        Trimming = True
        data = pd.read_csv("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/THR_Data/TRIM_ARM_DAC/TrimmingData.csv")
    else:
        Trimming = False

    # THR data
    Threshold_ARM_DAC = GetVFAT_THRDAC(source_folder,chamberID,VFAT_N)
    if Threshold_ARM_DAC == None:
            if verbose: print "No THR data for",chamberID," VFAT",VFAT_N
            return 10**6

    if Trimming:
        pandas_column_with_trim_fc = data[( (data["Chamber_ID"] == chamberID) & (data["VFATN"] == VFAT_N) )]["TRIM_fC"]
        overall_thresholds_dac =  pandas_column_with_trim_fc.mean() / THRDAC2fC + Threshold_ARM_DAC  ## Converting TRIM fC into DAC THR (conversion factors in PFA_Analyzer_Utils)
        if len(pandas_column_with_trim_fc) == 0:
            if verbose: print "No TRIM data for",chamberID," VFAT",VFAT_N
            return 10**6
    else:
        overall_thresholds_dac = Threshold_ARM_DAC

    return np.mean(overall_thresholds_dac)

## Returns the overall THR of a chamber, combining THR_DAC and TRIM_DAC (if the trimming was applied)        
def GetOverallChamberThreshold(source_folder,chamberID,verbose=False):
    if chamberID not in chamber_mapping.keys():
        print "Invalid chamberID",chamberID
        print "Exiting..."
        sys.exit(0)
        
    temp_array = np.empty(0,dtype=float)
    
    for VFAT_N in range(24):
        # THR data
        overall_vfat_trheshold = GetOverallVFATThreshold(source_folder,chamberID,VFAT_N,verbose=verbose)
        if overall_vfat_trheshold == 10**6:
                continue
        elif abs(overall_vfat_trheshold) > 255:
            if verbose: print "No-sense DAC value for VFAT",VFAT_N," on chamber ",chamberID,"\nPlease check, skipping VFAT value..."
            continue
        else:
            temp_array = np.append(temp_array,overall_vfat_trheshold)

    if len(temp_array) != 0:
        return np.mean(temp_array)
    else:
        return None
        

if __name__ == "__main__":
        for VFAT_N in range(24):
                print VFAT_N,GetOverallVFATThreshold("/afs/cern.ch/user/f/fivone/Test/PFA_Analyzer/THR_Data/THR_ARM_DAC/SBit10000_Trimming","GE11-M-16L1-L",VFAT_N,verbose=True)
        pass

