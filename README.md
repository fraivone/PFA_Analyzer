# PFA_Analyzer

Repository aimed to host the code to analyze the [GEM Common Muon NTuples](https://github.com/gmilella12/MuonDPGNTuples)

# PFA_Analyzer

### Installation
```
git clone https://github.com/fraivone/PFA_Analyzer.git
cd PFA_Analyzer
source setup.sh
```
### Usage
Let s suppose you have a bunch of  [compatible](#Compatibility) NTuples stored in the folder */A/B/C/*

1. Change the line containing `files = files_in_folder("<>") ` in `files = files_in_folder("/A/B/C/")`
1. Run  ``` python PFA_Analyzer.py -pc 0.01 -rdpc 4```

where
* pc = "φ cut" max angular distance (in rad) between RecHit and PropHit  to allow matching
* rdpc = "RΔφ cut" max distance (in cm) between RecHit and PropHit  to allow matching

Many other options can be provided as input (e.g STA chi2 cut, output name, fiducial cuts and etc...). You are encouraged to have a look at them `python PFA_Analyzer.py --help`

### Nuts and bolts
1. Runs through all the events included in the source NTuples
1. Fetches all the propagated hits on GEM that do pass the cut selection aka *Matchable PropHits*
1. For each *Matchable PropHits* checks if, in the same eta partition of the same GEM chamber, there is a GEM RecHits closer than **matching variable** cut

The analysis performs these operations for 2 different **matching variables**: 
* **φ** -->  glb_phi
* **RΔφ** --> glb_rdphi

The output data are always divided in 2 independent branches, named accordingly.

More info on the Analysis workflow --> [MWGR4 PFA Report](https://indico.cern.ch/event/1048923/contributions/4406801/attachments/2264472/3844543/PFA_FIvone_MWGR4_v1.pdf#page=33)

### Special Feature: Double Layer Efficiency (DLE)
DLE stands for double layer efficiency. In short, this method adds a tighter selection criteria on STA tracks to be used for efficiency evaluation.
* Consider only STA tracks with 2 PropHits, 1 for each layer of the same SC
* For efficiency on L1(2) consider only STA tracks with matched RecHitin L2(1)

When the boolean option DLE is provided, the analysis will produce an additional set of plots under "Efficiency/DLE".

**Why do I care?**

When the option DLE is selected, the efficiency is still evaluated in the classical way. However only STA tracks with 2 PropHits, 1 for each layer of the same SC are considered. So you do care because it lowers down the statistics.

### Output
The output file is named based on the outputname provided as input. If no outputname was provided the analysis date will be used:
```
outputname.root // yyMMdd_hhmm.root
```
where outputname 
Additionally, further details will be saved in the sub-folder **./Plot/outputname/** // **./Plot/yyMMdd_hhmm/**:

* Two csv files, containing the number of RecHit and PropHit for each unique etaP analyzed; one for glb_phi and one for RΔφ
* Two subfolders containing the most relevant plots as pdf, one for glb_phi and one for glb_rdphi

### Compatibility 
Compatible with GEM Common Ntuples produced with the release 
```
2021_MWGR4_v2
```