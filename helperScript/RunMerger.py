import ROOT
import argparse
import os.path
import sys
import pandas as pd
import subprocess
from argparse import RawTextHelpFormatter
import json
import re as regularExpression

from lib.ROOT_Utils import *
from lib.PFA_Analyzer_Utils import *


parser = argparse.ArgumentParser(
        description='''Scripts that merges the csv output from PFA_Analyzer.py for different runs, so that many runs can be merged without re-analyzing the NTuples.\nSupports chamber exclusion per run.\nREQUIRES inputs to be in the standard format day<N>_<runNumber>_<EqDivdCurr>uA''',
        epilog="""Typical exectuion (from top folder)\n\t  python -m helperScript.RunMerger --inputs day51_344817_690uA day52_344824_690uA day53_344859_690uA day5_342966_690uA --exclusion '{"GE11-P-12L2-L":[344824,344817],"GE11-P-03L1-S":[344859]}' --output myOutput""",
        formatter_class=RawTextHelpFormatter
        )


parser.add_argument('--inputs', type=str , help="Tag of the runs to be merged.\nMUST be in the standard format day<N>_<runNumber>_<EqDivdCurr>uA",required=True,nargs='*')
parser.add_argument('--output', type=str , help="output name",required=True)
parser.add_argument('--exclusion', type=json.loads,default={}, help="dictionary[chamber_ID] = list_of_runs_to_be_excluded from the merging. [-1] means 'do not merge at all' ",required=False)
parser.add_argument('--DLE', type=bool,default=False, help="Switch to true to merge Double Layer Efficiency csv data",required=False)

args = parser.parse_args()
inputs =  ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi.csv" for i in args.inputs]
inputs_byVFAT =  ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi_byVFAT.csv" for i in args.inputs]
inputs_byDLE =  ["./Output/PFA_Analyzer_Output/CSV/"+i+"/MatchingSummary_glb_rdphi_byDLE.csv" for i in args.inputs]
exclusion_dict = args.exclusion
output = args.output
DLE = args.DLE
run_list = [GetRunNumber(i) for i in args.inputs]


print "\n########\nRun(s) to be merged:\t\t\t",str(run_list)[1:-1]

## safety check the exclusion list
if type(exclusion_dict) != dict:
    print ("ERROR:\n\tParsed exclusion dict is not a python dictionary\nEXITING...\n")
    sys.exit(0)

for key,item in exclusion_dict.items():
    re,chamber,la = chamberName2ReChLa(key)
    if re not in [-1,1] or la not in [1,2] or chamber not in range(1,37) or len(key)!=13: ## A standard Chamber_ID has 13 chars
        print "ERROR:\n\tParsed exclusion dict containing bad chamber_ID:",key,"\nEXITING...\n"
        sys.exit(0)

    if (chamber%2==0 and key[-1]!="L") or (chamber%2==1 and key[-1]!="S"): ## If even must be Long, odd must be Short
        print "ERROR:\n\tBad chamber_ID: chamber number and chamber size do not match for",key,"\nEXITING...\n"
        sys.exit(0)
    print "Excluding ",key," from run(s):\t",str(exclusion_dict[key])[1:-1]
## safety check the exclusion list
print "########\n"


try:
    df_list = [ pd.read_csv(i, sep=',') for i in inputs]
except IOError:
            print "Pandas couldn't open one or more input files in ",inputs,"\nExiting .."
            sys.exit(0)
try:
    df_list_byVFAT = [ pd.read_csv(i, sep=',') for i in inputs_byVFAT]
except IOError:
            print "Pandas couldn't open one or more input files in ",inputs_byVFAT,"\nExiting .."
            sys.exit(0)
if DLE:
    try:
        df_list_byDLE = [ pd.read_csv(i, sep=',') for i in inputs_byDLE]
    except IOError:
        print "Pandas couldn't open one or more input files in ",inputs_byDLE,"\nExiting .."
        sys.exit(0)

tempList = []
tempList_byVFAT = []
tempList_byDLE = []


print "Merging ",len(df_list)," files for standard procedure"
print "Merging ",len(df_list_byVFAT)," files for efficiency by VFAT"
if DLE: print "Merging ",len(df_list_byDLE)," files for efficiency by DLE"

## Loop over all eta partitions
for region in [-1,1]:
    for layer in [1,2]:
        for ch in range(1,37):
            chamberID = ReChLa2chamberName(region,ch,layer)
            if chamberID in exclusion_dict.keys() and -1 in exclusion_dict[chamberID]:
                print chamberID+": completely excluded from the merge"
                continue

            for eta in range(1,9):
                matchedRecHit = 0
                propHit = 0
                matchedRecHit_DLE = 0
                propHit_DLE = 0

                ## Loop over all runs/files
                for file_index,df in enumerate(df_list):
                    if DLE: df_byDLE = df_list_byDLE[file_index]
                    if chamberID in exclusion_dict.keys() and int(run_list[file_index]) in exclusion_dict[chamberID]:
                        print "ChamberID ",chamberID,"\t",args.inputs[file_index],"\tExclusion DICT ", exclusion_dict[chamberID]
                        continue

                    temp_df = df[ (df['chamberID']==chamberID) & (df['etaPartition']==eta)]
                    if DLE: temp_df_byDLE = df_byDLE[ (df_byDLE['chamberID']==chamberID) & (df_byDLE['etaPartition']==eta)]
                    
                    if len(temp_df) != 0:
                        matchedRecHit += float(temp_df["matchedRecHit"])
                        propHit += float(temp_df["propHit"])
                    if DLE and len(temp_df_byDLE) != 0:
                        matchedRecHit_DLE += float(temp_df_byDLE["matchedRecHit"])
                        propHit_DLE += float(temp_df_byDLE["propHit"])
                
                tempList.append([chamberID,region,ch,layer,eta,matchedRecHit,propHit])
                if DLE: tempList_byDLE.append([chamberID,region,ch,layer,eta,matchedRecHit_DLE,propHit_DLE])


            for VFATN in range(24):
                matchedRecHit_VFAT = 0
                propHit_VFAT = 0

                for file_index,df_byVFAT in enumerate(df_list_byVFAT):

                    if chamberID in exclusion_dict.keys() and run_list[file_index] in exclusion_dict[chamberID]:
                        continue
                    temp_df_byVFAT = df_byVFAT[ (df_byVFAT['chamberID']==chamberID) & (df_byVFAT['VFATN']==VFATN)]

                    if len(temp_df_byVFAT) != 0:
                        matchedRecHit_VFAT += float(temp_df_byVFAT["matchedRecHit"])
                        propHit_VFAT += float(temp_df_byVFAT["propHit"])
                    
                ## End loop over runs/files
                
                tempList_byVFAT.append([chamberID,region,ch,layer,VFATN,matchedRecHit_VFAT,propHit_VFAT])
## End of the loop over all eta partitions


subprocess.call(["mkdir", "-p", "./Output/PFA_Analyzer_Output/CSV/"+output])
data = pd.DataFrame(tempList,columns=['chamberID',"region","chamber","layer","etaPartition","matchedRecHit","propHit"])
data.to_csv("./Output/PFA_Analyzer_Output/CSV/"+output+'/MatchingSummary_glb_rdphi.csv', index=False)

data_byVFAT = pd.DataFrame(tempList_byVFAT,columns=['chamberID',"region","chamber","layer","VFATN","matchedRecHit","propHit"])
data_byVFAT.to_csv("./Output/PFA_Analyzer_Output/CSV/"+output+'/MatchingSummary_glb_rdphi_byVFAT.csv', index=False)

if DLE:
    data_byDLE = pd.DataFrame(tempList_byDLE,columns=['chamberID',"region","chamber","layer","etaPartition","matchedRecHit","propHit"])
    data_byDLE.to_csv("./Output/PFA_Analyzer_Output/CSV/"+output+'/MatchingSummary_glb_rdphi_byDLE.csv', index=False)


print 
jsonFile = open("./Output/PFA_Analyzer_Output/CSV/"+output+"/Exclusion_Dict.json", "w")
json_data = json.dumps(exclusion_dict)
jsonFile.write(json_data)

print "Merge stored in\t\t\t"+"./Output/PFA_Analyzer_Output/CSV/"+output
print "Exclusion Dict stored in"+"\t./Output/PFA_Analyzer_Output/CSV/"+output+'/Exclusion_Dict.json'
print
