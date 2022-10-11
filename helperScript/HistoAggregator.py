import ROOT 
from ROOT_Utils import *
from PFA_Analyzer_Utils import *
from tqdm import tqdm

"""
Plots that loops through the ROOT files and merges together the plots for a certain quantity.
"""

#tags = ["357482_ZMu_cut.root", "357613_ZMu_cut.root", "357333_ZMu_cut.root", "357610_ZMu_cut.root", "357612_ZMu_cut.root", "357611_ZMu_cut.root", "357438_ZMu_cut.root", "357542_ZMu_cut.root", "357330_ZMu_cut.root", "357440_ZMu_cut.root", "357401_ZMu_cut.root", "357442_ZMu_cut.root", "357479_ZMu_cut.root", "357332_ZMu_cut.root", "357329_ZMu_cut.root"]
## run 11 17 august
#tags = ["357333_ZMu_pt_chi2_cut.root", "357613_ZMu_pt_chi2_cut.root", "357610_ZMu_pt_chi2_cut.root", "357482_ZMu_pt_chi2_cut.root", "357542_ZMu_pt_chi2_cut.root", "357330_ZMu_pt_chi2_cut.root", "357611_ZMu_pt_chi2_cut.root", "357440_ZMu_pt_chi2_cut.root", "357438_ZMu_pt_chi2_cut.root", "357612_ZMu_pt_chi2_cut.root", "357401_ZMu_pt_chi2_cut.root", "357442_ZMu_pt_chi2_cut.root", "357332_ZMu_pt_chi2_cut.root", "357479_ZMu_pt_chi2_cut.root", "357329_ZMu_pt_chi2_cut.root"]
tags = ["357479_Prompt_pt_chi2_cut.root"]

fileNames = [f"/afs/cern.ch/user/f/fivone/Documents/test/Output/PFA_Analyzer_Output/ROOT_File/{tag}" for tag in tags]

ROOT.gROOT.SetBatch(True)
c1 = setUpCanvas("temp",1200,1200)
OutF = ROOT.TFile("../Output/Plot/Aggregate.root","RECREATE")

etas = [k for k in range(1,9)] + ["All"]
endcaps = ["PL1_","PL2_","ML1_","ML2_",""]
sizes = ["long","short"]
CLSs = [1,2,3,4,5]
### Fill this one to add another histo
## {etaP} {endcaptag} and {chamber} will be replaced with value
## for each main entry, the search path must be unique
histo_directives = {"Residual":
                    {"Title":"Residual eta{etaP}",
                     "Xaxis Title":"R#Delta#phi (cm)",
                     "bins":120,
                     "range":[-0.5,0.5],
                     "search":"{chamber}*eta*{etaP}*Residual",
                     "path":"Residuals/MatchingOn_glb_rdphi",
                     "color":ROOT.kOrange+2,
                     "linewidth":2
                 },
                    "ResidualY":
                    {"Title":"Residual eta{etaP}",
                     "Xaxis Title":"#Delta#LocalY (cm)",
                     "bins":120,
                     "range":[-10,10],
                     "search":"{chamber}*eta*{etaP}*ResidualY",
                     "path":"Residuals/MatchingOn_glb_rdphi",
                     "color":ROOT.kOrange+2,
                     "linewidth":2
                 },
                    "CLSBeforeMatching":
                    {"Title":"CLS Before Matching eta{etaP}",
                     "Xaxis Title":"ClusterSize",
                     "bins":20,
                     "range":[0,20],
                     "search":"{chamber}*eta*{etaP}*ClusterSize*",
                     "path":"SanityChecks/CLS/beforematching",
                     "color":ROOT.kGreen+2,
                     "linewidth":2
                 },
                    "CLSAfterMatching":
                    {"Title":"CLS After Matching eta{etaP}",
                     "Xaxis Title":"Clustersize",
                     "bins":20,
                     "range":[0,20],
                     "search":"{chamber}*eta*{etaP}*ClusterSize*",
                     "path":"SanityChecks/CLS/aftermatching",
                     "color":ROOT.kBlue+2,
                     "linewidth":2
                 },
                    "NumberRecHitperEvt":
                    {"Title":"Number of GEMRecoHits per Evt {endcaptag}",
                     "Xaxis Title":"NGEM Rechits",
                     "bins":200,
                     "range":[0,200],
                     "search":"{endcaptag}NRecoHitsPerEVT",
                     "path":"SanityChecks/NHits/BeforeMatching",
                     "color":ROOT.kRed+1,
                     "linewidth":2
                 },
                    "LocalResidual":
                    {"Title":"RecHit residual by local position: {size} chamber",
                     "Xaxis Title":"Loc x (cm)",
                     "Yaxis Title":"Loc y (cm)",
                     "binx":200,
                     "biny":200,
                     "rangex":[-80,80],
                     "rangey":[-80,80],
                     "search":"{size}Chmbrs*res",
                     "path":"Residuals/MatchingOn_glb_rdphi/2D_glb_rdphi"},

                    "ResidualCLS":
                    {"Title":"Residuals ClusterSize{CLS}",
                     "Xaxis Title":"R#Delta#phi (cm)",
                     "bins":120,
                     "range":[-0.5,0.5],
                     "search":"Residual CLS {CLS}",
                     "path":"Residuals/MatchingOn_glb_rdphi",
                     "color":ROOT.kMagenta+1,
                     "linewidth":2
                 }
                }




def generateHisto(histo_directives):
    TH1_Plots = {}
    ## Loop over all quantities
    for k in histo_directives:
        TH1_Plots[k] = {}

        if "{etaP}" in histo_directives[k]["Title"]:
            for eta in etas:
                title = histo_directives[k]["Title"].format(etaP=eta)
                bins = histo_directives[k]["bins"]
                x_range = histo_directives[k]["range"]
                x_axis_title = histo_directives[k]["Xaxis Title"]
                color = histo_directives[k]["color"]
                linewidth = histo_directives[k]["linewidth"]

                TH1_Plots[k][eta] = ROOT.TH1F(title,title,bins,x_range[0],x_range[1])
                TH1_Plots[k][eta].GetXaxis().SetTitle(x_axis_title)
                TH1_Plots[k][eta].SetLineColor(color)
                TH1_Plots[k][eta].SetLineWidth(linewidth)
        elif "{CLS}" in histo_directives[k]["Title"]:
            for CLS in CLSs:
                title = histo_directives[k]["Title"].format(CLS=CLS)
                bins = histo_directives[k]["bins"]
                x_range = histo_directives[k]["range"]
                x_axis_title = histo_directives[k]["Xaxis Title"]
                color = histo_directives[k]["color"]
                linewidth = histo_directives[k]["linewidth"]

                TH1_Plots[k][CLS] = ROOT.TH1F(title,title,bins,x_range[0],x_range[1])
                TH1_Plots[k][CLS].GetXaxis().SetTitle(x_axis_title)
                TH1_Plots[k][CLS].SetLineColor(color)
                TH1_Plots[k][CLS].SetLineWidth(linewidth)


        elif "{endcaptag}" in histo_directives[k]["Title"]:
            for endcaptag in endcaps: 
                title = histo_directives[k]["Title"].format(endcaptag=endcaptag)
                bins = histo_directives[k]["bins"]
                x_range = histo_directives[k]["range"]
                x_axis_title = histo_directives[k]["Xaxis Title"]
                color = histo_directives[k]["color"]
                linewidth = histo_directives[k]["linewidth"]               

                TH1_Plots[k][endcaptag] = ROOT.TH1F(title,title,bins,x_range[0],x_range[1])
                TH1_Plots[k][endcaptag].GetXaxis().SetTitle(x_axis_title)
                TH1_Plots[k][endcaptag].SetLineColor(color)
                TH1_Plots[k][endcaptag].SetLineWidth(linewidth)
        elif "{size}" in histo_directives[k]["Title"]:
            for size in sizes:
                title = histo_directives[k]["Title"].format(size=size)
                binx = histo_directives[k]["binx"]
                biny = histo_directives[k]["biny"]
                x_range = histo_directives[k]["rangex"]
                y_range = histo_directives[k]["rangey"]
                x_axis_title = histo_directives[k]["Xaxis Title"]
                y_axis_title = histo_directives[k]["Yaxis Title"]
                
                TH1_Plots[k][size] = ROOT.TH2F(title,title,binx,x_range[0],x_range[1],biny,y_range[0],y_range[1])
                TH1_Plots[k][size].GetXaxis().SetTitle(x_axis_title)
                TH1_Plots[k][size].GetYaxis().SetTitle(y_axis_title)
        else:
            print("Generate histo: TO BE IMPLEMENTED")
            input("press a key")
            
    return TH1_Plots

def generateSkimmedDicts(histo_directives,mp):
    skimmed_dicts = {}
    for k in histo_directives:
        splitted_path = histo_directives[k]["path"].split("/")
        temp = mp
        for step in splitted_path:
            print(k,step)
            temp = temp[step][1]
        skimmed_dicts[k] = temp

    return skimmed_dicts

TH1Plots = generateHisto(histo_directives)

print(TH1Plots["ResidualCLS"].keys())
for fn in tqdm(fileNames):
    print(fn)
    mp=Map_TFile(fn)

    skimmed_dicts = generateSkimmedDicts(histo_directives,mp[fn][1])

    ## Fill histograms
    for k in histo_directives:
        search_tag = histo_directives[k]["search"]
        if "{chamber}" in search_tag and "{etaP}" in search_tag:
            for ch in GE11_ChamberIDs():
                for eta in etas:
                    ### Only one hist found during the search
                    query_result_dict = getObjects(skimmed_dicts[k], search_tag.format(etaP=eta,chamber=ch) )
                    if len( query_result_dict ) == 1: 
                            TH1Plots[k][eta].Add( list(query_result_dict.values())[0] )
        elif "{CLS}" in histo_directives[k]["Title"]:
            for CLS in CLSs: 
                ### Only one hist found during the search
                query_result_dict = getObjects(skimmed_dicts[k], search_tag.format(CLS=CLS) )
                if len( query_result_dict ) == 1: 
                    TH1Plots[k][CLS].Add( list(query_result_dict.values())[0] )
        elif "{endcaptag}" in histo_directives[k]["Title"]:
            for endcaptag in endcaps: 
                ### Only one hist found during the search
                query_result_dict = getObjects(skimmed_dicts[k], search_tag.format(endcaptag=endcaptag) )
                if len( query_result_dict ) == 1: 
                    TH1Plots[k][endcaptag].Add( list(query_result_dict.values())[0] )
        elif "{size}" in histo_directives[k]["Title"]:
            for size in sizes: 
                ### Only one hist found during the search
                query_result_dict = getObjects(skimmed_dicts[k], search_tag.format(size=size) )
                if len( query_result_dict ) == 1: 
                    TH1Plots[k][size].Add( list(query_result_dict.values())[0] )
        else:
            print("Filling Hist: TO BE IMPLEMENTED")
            input("press a key")


for k in TH1Plots:
    for subk,item in TH1Plots[k].items():
        c1.cd()
        item.Draw("HIST")
        if type(item) == ROOT.TH2F: item.Draw("COLZ")
        c1.Modified()
        c1.Update()
        c1.SaveAs(f"../Output/Plot/{item.GetTitle()}.png")
        writeToTFile(OutF,item,k+"/")

OutF.Close()
