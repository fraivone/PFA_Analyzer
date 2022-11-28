"""
Generate a unique dataframe containing VFAT efficiency info for all the runs.
Adding a function to dump the dataframe to csv file.
"""
import pandas as pd
from PFA_Analyzer_Utils import *
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt
from scipy.stats import t as tstudent
import time
import subprocess

parser = argparse.ArgumentParser(
        description='''Scripts that generates a dataframe containing efficiency per VFATs over input runs''',
        epilog="""Typical exectuion  (from top folder)\n\t  python -m helperScript.CreateDataset --input 357735_ZMu_pT20 357700_ZMu_pT20""",
        formatter_class=RawTextHelpFormatter
        )

parser.add_argument('--input', type=str , nargs="*", help="Tag of the run to be plotted",required=True)





def calculateTstudent(a,b): #a=(var1,errvar1,nvar1) b=(var2,errvar2,nvar2)
    var1,errvar1,nvar1 = a
    var2,errvar2,nvar2 = b
    DOF = nvar1 + nvar2 - 2
    A = (nvar1 + nvar2) / (nvar1*nvar2)
    B = ( (nvar1 - 1)*errvar1**2 + (nvar2 - 1)*errvar2**2 ) / ( nvar1 + nvar2 - 2 )

    t = np.abs(var1 - var2) / np.sqrt(A*B)
    return t,DOF


## as a standard, I keep all VFATs there. When nprop == 0, efficiency will be nan
def num_den(row):
    num,den = row['matchedRecHit'],row['propHit']
    return num,den
def calcEfficiency(num,den):
    if den == 0:
        return np.nan,np.nan
    else:
        lowL,upL = generateClopperPeasrsonInterval(num,den)
        eff,effErr = float(num)/den,(upL - lowL)/2
        return eff,effErr

def generate(inputTags,save=True):
    mergedf = pd.DataFrame()
    outputFolder = f"{OUTPUT_PATH}/VFAT_Efficiency/{ time.strftime('%-y%m%d_%H%M')}"
    subprocess.call(["mkdir", "-p", outputFolder])

    for tag in inputTags:
        run_number = GetRunNumber(tag) 
        inputFile = f"{OUTPUT_PATH}/PFA_Analyzer_Output/CSV/{tag}/MatchingSummary_glb_rdphi_byVFAT.csv"
        df = pd.read_csv(inputFile, sep=',')
        df['efficiency'] = df.apply (lambda row: calcEfficiency(num_den(row)[0],num_den(row)[1])[0], axis=1)
        df['efficiencyError'] = df.apply (lambda row: calcEfficiency(num_den(row)[0],num_den(row)[1])[1], axis=1)
        df['run'] = run_number

        mergedf = mergedf.append(df)
    
    if save:
        for j in ["chamber","layer","region"]:
            mergedf.drop(labels=j,inplace=True,axis=1)
        mergedf.to_csv(f"{outputFolder}/EfficiencyPerVFAT.csv", index=False)
        whichrunsfile = f"{outputFolder}/RunUsed.txt"
        with open(whichrunsfile, "w") as outfile:
            outfile.write("\n".join(inputTags))

    return mergedf



def checkVFATOutliers(df):
    reference_eff_chamber = {}
    ChamberIDs = ["GE11-M-02L1-L"]#set(df["chamberID"].tolist())
    for ch in ChamberIDs:
        ## Step 1 reference efficiency value
        df_sel = df[df["chamberID"] == ch ]
        ## Runs in which efficiency for chamber ch was calculated
        runs = set(df_sel.dropna()["run"].tolist())        
        n_runs = len(runs)

        matched, prop = df_sel["matchedRecHit"].sum(), df_sel["propHit"].sum()
        eff, effErr = calcEfficiency(matched,prop)
        reference_eff_chamber[ch] = (eff,effErr,n_runs)


        for run in runs:
            current_run_df = df_sel[ df_sel["run"] == run]
            num, den = sum(current_run_df["matchedRecHit"].tolist()) , sum(current_run_df["propHit"].tolist())

            eff,effErr = calcEfficiency(num,den)
            t ,DOF = calculateTstudent( (eff,effErr,1) , reference_eff_chamber[ch] )
            print(run,t,tstudent.cdf(t, DOF),eff,effErr,reference_eff_chamber[ch])
            #print(run,eff,effErr,reference_eff_chamber[ch], np.abs(eff - reference_eff_chamber[ch][0])/reference_eff_chamber[ch][1])
            
    
if __name__ == '__main__':
    args = parser.parse_args()
    inputTags = args.input
    df = generate(inputTags)
    # checkVFATOutliers(df)