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
1. Change the line  `chambersOFF = ["GE11-P-14L1", "GE11-M-20L1"] ` with the list of  chambers you want to exclude from the analysis, mantaining the same naming convention
1. Run  ``` python PFA_Analyzer.py -pc 0.01 -rdpc 4```

where
* pc = "φ cut" max angular distance (in rad) between RecHit and PropHit  to allow matching
* rdpc = "RΔφ cut" max distance (in cm) between RecHit and PropHit  to allow matching

Other 2nd order parameter as the fiducial cuts, propagation error cuts, muon pt cut  can be adjusted by direct editing of the PFA_Analyzer.py file.

### Nuts and bolts
1. Runs through all the events included in the source NTuples
1. Fetches all the propagated hits on GEM that do pass the cut selection aka *Matchable PropHits*
1. For each *Matchable PropHits* checks if, in the same eta partition of the same GEM chamber, there is a GEM RecHits closer than **matching variable** cut

The analysis performs these operations for 2 different **matching variables**: 
* **φ** -->  glb_phi
* **RΔφ** --> glb_rdphi

The output data are always divided in 2 independent branches, named accordingly.

More info on the Analysis workflow --> [MWGR4 PFA Report](https://indico.cern.ch/event/1048923/contributions/4406801/attachments/2264472/3844543/PFA_FIvone_MWGR4_v1.pdf#page=33)


### Output
The output will be stored in the root file named:
```
PFA_Output_yyMMdd_hhmm.root
```
Additionally, further details will be saved in the sub-folder **./Plot/yyMMdd_hhmm/**:

* Two csv files, containing the number of RecHit and PropHit for each unique etaP analyzed; one for glb_phi and one for RΔφ
* Two subfolders containing the most relevant plots as pdf, one for glb_phi and one for glb_rdphi

### Compatibility 
Compatible with:
 1. [MuonDPGNtuples 2021_MWGR4_v1](https://github.com/gmilella12/MuonDPGNTuples/releases) 