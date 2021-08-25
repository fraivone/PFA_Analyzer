- [PFA_Analyzer](#pfa-analyzer)
    + [Installation](#installation)
    + [Usage](#usage)
      - [Single run](#single-run)
        * [Required Inputs](#required-inputs)
        * [Optional Inputs (some of them)](#optional-inputs--some-of-them-)
      - [Merging runs](#merging-runs)
        * [Required Inputs](#required-inputs-1)
        * [Optional Inputs (some of them)](#optional-inputs--some-of-them--1)
    + [Nuts and bolts](#nuts-and-bolts)
    + [Special Features](#special-features)
      - [Double Layer Efficiency (DLE)](#double-layer-efficiency--dle-)
      - [Full Digis (FD)](#full-digis--fd-)
    + [Output](#output)
    + [Compatibility](#compatibility)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# PFA_Analyzer
**NOTE : Not Final**
Repository aimed to host the code to analyze the [GEM Common Muon NTuples](https://github.com/gmilella12/MuonDPGNTuples)



### Installation
```
git clone https://github.com/gem-sw/pfa.git
cd PFA_Analyzer
source setup.sh
```
### Usage
#### Single run
##### Required Inputs
* Set of [compatible](#Compatibility) GEM NTuples related to your run, stored under */A/B/C/*

##### Optional Inputs (some of them)
* List of chambers to be masked ([./ChamberOFF_Example.txt](./ExcludeMe/ChamberOFF_CRUZET_342728.txt)) cause were OFF/Tripped/inError during the run
* Tab separated list of VFAT to be masked ([./ListOfDeadVFAT_Example.txt](/ExcludeMe/ListOfDeadVFAT_run343034.txt)) cause were inError/Noisy/NoData during the run
* *myOutput* as label for the output data 

By executing
```
python PFA_Analyzer.py --dataset /A/B/C/ --chamberOFF ./ChamberOFF_Example.txt  --VFATOFF ./ListOfDeadVFAT_Example.txt -pc 0.02 -rdpc 4  --outputname myOutput
```

the GEM NTuples stored in */A/B/C/* will be anlyzed using 
* pc = 0.02 rad  &emsp;&emsp; is the "φ cut" value, the max angular distance between RecHit and PropHit to allow matching
* rdpc = 4 cm   &emsp;&emsp;is the "RΔφ cut" value, the max distance (in cm) between RecHit and PropHit  to allow matching

The following surfaces will be ignored:
* Entire surface of the chambers listed in ./ChamberOFF_Example.txt
* Entire surface of the VFATs listed in ./ListOfDeadVFAT_Example.txt

Additionaly (set by default but can be parsed as option):
* Fiducial Cut on R = 1 cm  (suggested for cosmics data)
* Fiducial Cut on Phi = 5 mrad (suggested for cosmics data)
* Max Error On Propagated Hit in R = 1 cm (suggested for cosmics data)
* Max Error On Propagated Hit in φ = 10 mrad (suggested for cosmics data)

Many other options can be provided as input (e.g STA chi2 cut, fiducial cuts values, number of MEX hits etc...). You are encouraged to have a look at them `python PFA_Analyzer.py --help`


#### Merging runs

To achieve better statistics you want to analyze many sets of GEM NTuples coming from different runs. Let's suppose they are all [compatible](#Compatibility) with this analyzer release and stored in:
##### Required Inputs
* /A/B/C/Run1/
* /A/B/C/Run2/
* /A/B/C/Run3/

##### Optional Inputs (some of them)
Let's suppose you have a list of chamber OFF for each of this run:
* ./ChamberOFF_Run_1.txt, * ./ChamberOFF_Run_2.txt, * ./ChamberOFF_Run_3.txt
Let's suppose you have a list of VFAT OFF for each of this run:
* ./ListOfDeadVFAT_Run_1.txt, * ./ListOfDeadVFAT_Run_2.txt, * ./ListOfDeadVFAT_Run_3.txt
* *myOutput* as label for the output data 

By executing
```
python PFA_Analyzer.py --dataset /A/B/C/Run1/ /A/B/C/Run2/ /A/B/C/Run3/ --chamberOFF ./ChamberOFF_Run_1.txt ./ChamberOFF_Run_2.txt ./ChamberOFF_Run_3.txt --VFATOFF ./ListOfDeadVFAT_Run_1.txt ./ListOfDeadVFAT_Run_2.txt ./ListOfDeadVFAT_Run_3.txt -pc 0.02 -rdpc 4  --outputname myOutput
```

the GEM NTuples stored in /A/B/C/Run1/, /A/B/C/Run2/, /A/B/C/Run3/ will be merged together and anlyzed using :
* pc = 0.02 rad  &emsp;&emsp; is the "φ cut" value, the max angular distance between RecHit and PropHit to allow matching
* rdpc = 4 cm   &emsp;&emsp;is the "RΔφ cut" value, the max distance (in cm) between RecHit and PropHit  to allow matching

The following surfaces will be ignored:
* The surface that comes from the union of the chambers listed in ./ChamberOFF_Run_1.txt, ./ChamberOFF_Run_2.txt, ./ChamberOFF_Run_3.txt
* The surface that comes from the union of the VFATs listed in ./ListOfDeadVFAT_Run_1.txt, ./ListOfDeadVFAT_Run_2.txt, ./ListOfDeadVFAT_Run_3.txt

Additionaly (set by default but can be parsed as option):
* Fiducial Cut on R = 1 cm  (suggested for cosmics data)
* Fiducial Cut on Phi = 5 mrad (suggested for cosmics data)
* Max Error On Propagated Hit in R = 1 cm (suggested for cosmics data)
* Max Error On Propagated Hit in φ = 10 mrad (suggested for cosmics data)


### Nuts and bolts
1. Runs through all the events included in the source NTuples
1. Fetches all the propagated hits on GEM that do pass the cut selection and are not ignored due to masking aka *Matchable PropHits*
1. For each *Matchable PropHits* checks if, in the same eta partition of the same GEM chamber, there is a GEM RecHits closer than **matching variable** cut

The analysis performs these operations for 2 different **matching variables**: 
* **φ** -->  glb_phi
* **RΔφ** --> glb_rdphi

The output data are always divided in 2 independent branches, named accordingly.

More info on the Analysis workflow --> [MWGR4 PFA Report](https://indico.cern.ch/event/1048923/contributions/4406801/attachments/2264472/3844543/PFA_FIvone_MWGR4_v1.pdf#page=33)

### Special Features

#### Double Layer Efficiency (DLE)
DLE stands for double layer efficiency. In short, this method adds a tighter selection criteria on STA tracks to be used for efficiency evaluation.
* Consider only STA tracks with 2 PropHits, 1 for each layer of the same SC
* For efficiency on L1(2) consider only STA tracks with matched RecHitin L2(1)

When the boolean option DLE is provided, the analysis will produce an additional set of plots under "Efficiency/DLE".

**Why do I care?**

When the option DLE is selected, the efficiency is still evaluated in the classical way. However only the events containing STA tracks with exactly 2 PropHits ( 1 for each layer of the same SC) are considered. So you do care because it lowers down the statistics.

#### Full Digis (FD)
This analysis uses the propagated hits coming from reconstructed STA muons to probe GEMs performance. 
Therefore events without propagated hits on GEM are not processed. This allows to speed up the analysis of the data, but ignores the RecHit data associated with the event.

When the boolean option --FD is provided, no events are ignored and GEM digis are collected and stored in the SanityChekc plots even for events without propagated hits. 
The processing time increases.

### Output
A typical output file is [./day19_344063_660uA.root](./day19_344063_660uA.root)

Three main groups of folders can be found in the output file:
* Performance Related Folder
  + **Residuals**
  + **Efficiency**
* Overview of other quantities
  + **SanityChecks**
* Description of the input parameters for the analysis
  + **Metadata**

The output file is named based on the outputname provided as input. If no outputname was provided, the analysis date will be used:
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