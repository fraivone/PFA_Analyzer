#!/usr/bin/env python
import argparse
import os
import yaml
import pathlib


parser = argparse.ArgumentParser(description='generates config file for PFA_Analyzer ')
parser.add_argument('run_number', help='run number')
parser.add_argument('hv', help='hv')
parser.add_argument('lumisection', help='lumisection of start',type=int)

args = parser.parse_args()
config_filename = str(args.run_number)+"_"+str(args.hv)+"uA.yml"
config_file = pathlib.Path(os.path.abspath(config_filename))

## Run chamber masker
os.system("cd /afs/cern.ch/user/f/fivone/Test/Chamber_MaskMaker; /usr/bin/python PFA_MaskGenerator.py -rl "+ str(args.run_number)+" -iexpl "+str(args.hv)+" -o "+str(args.run_number)+"_"+str(args.hv)+"uA;")
## Copy config file template
os.system("cd /afs/cern.ch/user/f/fivone/Documents/test/config; cp test_config.yml "+config_filename )

chamber_off = f"/afs/cern.ch/user/f/fivone/Test/Analyzer/ExcludeMe/{args.run_number}_{args.hv}uA.json"
latest_folder = os.popen("ls -tr /eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2022/GEMCommonNtuples/Muon/"+str(args.run_number)+"_ZMu/ | tail  -1").read()
latest_folder = latest_folder.split("\n")[0]

with config_file.open() as ymlfile:
        cfg = yaml.full_load(ymlfile)
for k in cfg["data"]:
    cfg["data"][int(args.run_number)] = cfg["data"].pop(k)
    cfg["data"][int(args.run_number)]["tag"] = str(args.run_number)+"_ZMu"
    cfg["data"][int(args.run_number)]["path_tag"] = latest_folder
    cfg["data"][int(args.run_number)]["lumisection"] = (args.lumisection,10000)
    cfg["data"][int(args.run_number)]["chamberOFF"] = chamber_off
    cfg["parameters"]["outputname"] = f"{args.run_number}_{args.hv}uA"

with open(config_file, 'w') as outfile:
    yaml.safe_dump(cfg, outfile, default_flow_style=False,sort_keys=False)
