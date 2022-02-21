import ROOT
import time
from PFA_Analyzer_Utils import *
from ROOT_Utils import *
import argparse
from argparse import RawTextHelpFormatter
import subprocess
ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(
        description='''Scripts that: \n\t-Reads the GEMMuonNtuple\n\t-Skims them to take only the useful branches''',
        epilog="""Typical exectuion\n\t python GEMNTuplesSkimmer.py  --dataset 2Jul2021""",
        formatter_class=RawTextHelpFormatter
)
parser.add_argument('--dataset','-ds', type=str,help="TAG to the folder containing the NTuples to be analyzed",required=True,nargs='*')

args = parser.parse_args()

start_time = time.time()


## Input data
files = []

for folder in args.dataset:
    
    temp_files = []
    all_folders = subprocess.check_output(["find", "/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021", "-type", "d", "-name", "*"+str(folder)+"*"]).split("\n")
    all_folders.remove('')
    all_folders = [folder for folder in all_folders if "L1A" not in folder and "new" not in folder] 
    for s in all_folders: 
        temp_files += files_in_folder(s)
    files += [f for f in temp_files if ".root" in f]


n = len(files)
print (files)
chain = ROOT.TChain("muNtupleProducer/MuDPGTree")
print "Chaining ",n," files"
for index,fl in enumerate(files):
    chain.Add(fl)
    print "Chained ",index+1,"/",n, " files"


## Prepping output file
OutF = ROOT.TFile("/eos/user/f/fivone/GEMNTuples/"+args.dataset[0]+".root","RECREATE")

## Selecting useful branches, skim the others
chain.SetBranchStatus("*", 0)
chain.SetBranchStatus("event_eventNumber", 1)
chain.SetBranchStatus("mu_propagatedGlb_x", 1)
chain.SetBranchStatus("mu_propagatedGlb_y", 1)
chain.SetBranchStatus("mu_propagated_isME11", 1)
chain.SetBranchStatus("mu_propagated_layer", 1)
chain.SetBranchStatus("mu_propagated_chamber", 1)
chain.SetBranchStatus("mu_propagated_region", 1)






OutF.cd()
tree = ROOT.TTree()
print "Created output TTree"
tree = chain.CloneTree()
OutF.Write()
OutF.Close()

## Summary Plots
print("--- %s seconds ---" % (time.time() - start_time))
print("\n"+"/eos/user/f/fivone/GEMNTuples/"+args.dataset[0]+".root")
