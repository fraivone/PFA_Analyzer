import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import re as regularExpression

base_folder = os.path.expandvars('$PFA')
sys.path.insert(1, base_folder+"/Analyzer/lib/")
try:
    from PFA_Analyzer_Utils import *
except:
   print ("ERROR:\n\tCan't find the package PFA_Analyzer_Utils in ",lib_folder,"\nEXITING...\n")
   sys.exit(0)


def generateSubFile(outputName,shell_name):
    job_flavour = "tomorrow"
    SubfileName = base_folder+"/Analyzer/CondorSub/JobFiles/SubmitFile_"+str(outputName)+".sub"
    with open(SubfileName, 'w') as sout:
        sout.write("executable              = "+shell_name+"\n")
        sout.write("getenv                  = true"+"\n")
        sout.write("arguments               = $(ClusterId) $(ProcId)"+"\n")
        sout.write("output                  = ./Logs/out.$(ClusterId)_Run"+str(outputName)+".dat"+"\n")
        sout.write("error                   = ./Logs/error.$(ClusterId)_Run"+str(outputName)+".err"+"\n")
        sout.write("log                     = ./Logs/log.$(ClusterId)_Run"+str(outputName)+".log"+"\n")
        sout.write("+JobFlavour             = \""+job_flavour+"\" "+"\n")
        sout.write("notify_user             = francesco.ivone@cern.ch\n")
        sout.write("notification            = always\n")
        sout.write("queue"+"\n")
    return SubfileName

def generateJobShell(run_number,outputName,pc,rdpc,minPt,max_NormChi2,minME1Hit,minME2Hit,minME3Hit,minME4Hit,maxErrOnPropR,maxErrOnPropPhi,maskChVFAT,doubleLayerEfficiency):
    run_number_int = GetRunNumber(run_number)
    main_command = "python PFA_Analyzer.py --dataset "+str(run_number)+" -pc " +str(pc) + " -rdpc "+str(rdpc)+" --outputname "+outputName+" --minPt "+str(minPt) +   " --chi2cut "+str(max_NormChi2) +" --minME1 "+str(minME1Hit) + " --minME2 "+str(minME2Hit) +  " --minME3 "+str(minME3Hit) + " --minME4 "+str(minME4Hit) + " --maxErrPropR "+str(maxErrOnPropR)+" --maxErrPropPhi "+str(maxErrOnPropPhi)
        
    #if maskChVFAT == True: 
    main_command = main_command + " --chamberOFF ./ExcludeMe/ChamberOFF_Run_"+str(run_number_int)+".json --VFATOFF ./ExcludeMe/ListOfDeadVFAT_run"+str(run_number_int)+".txt"
    if doubleLayerEfficiency == True:
        main_command = main_command + " --DLE"
        
    main_command = main_command + " \n"

    shell_script_name = base_folder+"/Analyzer/CondorSub/JobFiles/job_Run_"+str(outputName)+".sh"
    with open(shell_script_name, 'w') as fout:
        ####### Write the instruction for each job
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo %s_Run"+str(run_number)+"\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("cd "+base_folder+"/Analyzer/"+"\n")
        ## sourceing the right gcc version to compile the source code
        fout.write("source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9.1.0/x86_64-centos7/setup.sh"+"\n")
        fout.write("source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh"+"\n")
        ## Run the same job interval number of time
        fout.write("ClusterId=$1\n")
        fout.write("ProcId=$2\n")
        fout.write(main_command)
    return shell_script_name
if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='''Scripts runs PFA_Analyzer.py for many different input runs''',
        epilog="""Typical exectuion (from this folder)\n\t  python runJobs.py --runList 349527_Express 349528_Express --maskChVFAT""",
        formatter_class=RawTextHelpFormatter
        )


    parser.add_argument('--runList', type=str , help="List of the runs to be analyzed (make sure they're  ntuplized)",required=True,nargs='*')
    parser.add_argument('-pc','--phi_cut', type=float,help="Maximum allowed dphi between RecoHit and PropHit to be counted as matched hit",required=False)
    parser.add_argument('-rdpc','--rdphi_cut', type=float,help="Maximum allowed rdphi between RecoHit and PropHit to be counted as matched hit",required=False)
    parser.add_argument('--outputNames', type=str , help="List of output names",required=False,nargs='*')
    parser.add_argument('--maskChVFAT' ,dest='maskChVFAT', help="Activate chamber/VFAT masking",action='store_true')
    parser.add_argument('--chi2cut', type=float,help="Maximum normalized chi2 for which accept propagated tracks",required=False)
    parser.add_argument('--minPt', type=float,help="Minimum pt for which accept propagated tracks, in GeV",required=False)
    parser.add_argument('--minME1', type=int, help="Min number of ME1 hits, 0 by default",required=False)
    parser.add_argument('--minME2', type=int, help="Min number of ME2 hits, 0 by default",required=False)
    parser.add_argument('--minME3', type=int, help="Min number of ME3 hits, 0 by default",required=False)
    parser.add_argument('--minME4', type=int, help="Min number of ME4 hits, 0 by default",required=False)
    parser.add_argument('--maxErrPropR', type=float , help="max error on propagated R in order to accept the muon",required=False)
    parser.add_argument('--maxErrPropPhi', type=float , help="max error on propagated phi in order to accept the muon",required=False)
    parser.add_argument('--DLE', default=False, action='store_true',help="Swtiches on the Double Layer Efficiency (DLE) analisys. False by default",required=False)

    parser.set_defaults(phi_cut=0.02)
    parser.set_defaults(rdphi_cut=4)
    parser.set_defaults(minPt=0)
    parser.set_defaults(minME1=0)
    parser.set_defaults(minME2=0)
    parser.set_defaults(minME3=0)
    parser.set_defaults(minME4=0)
    parser.set_defaults(maxErrPropR=1)
    parser.set_defaults(maxErrPropPhi=0.005)
    parser.set_defaults(maskChVFAT=False)
    parser.set_defaults(chi2cut=9999999)

    args = parser.parse_args()
    inputs = args.runList
    pc = args.phi_cut
    rdpc = args.rdphi_cut
    if args.outputNames != None: outputs = args.outputNames
    else: outputs = inputs
    maskChVFAT = args.maskChVFAT
    DLE = args.DLE
    minPt = args.minPt
    maxSTA_NormChi2 = args.chi2cut
    minME1Hit = args.minME1
    minME2Hit = args.minME2
    minME3Hit = args.minME3
    minME4Hit = args.minME4
    maxErrOnPropR = args.maxErrPropR
    maxErrOnPropPhi = args.maxErrPropPhi


    if len(inputs) != len(outputs):
        print "Parsed runList and outputNames are different in number...\nExiting .."
        sys.exit(0)    

    if maskChVFAT:
        ## Mask chamber empty, in error or with eq. divider current not at 700
        run_string = " ".join( [str(GetRunNumber(i)) for i in inputs])
        eq_divid_current = ''.join(["690 " for i in range(len(inputs))])
        os.system("python "+base_folder+"/Chamber_MaskMaker/PFA_MaskGenerator.py -rl "+run_string+ " -iexpl "+eq_divid_current)


    for index in range(len(inputs)):
        run = inputs[index]
        name = outputs[index]

        shell_name = generateJobShell(run,name,pc,rdpc,minPt,maxSTA_NormChi2,minME1Hit,minME2Hit,minME3Hit,minME4Hit,maxErrOnPropR,maxErrOnPropPhi,maskChVFAT,DLE)
        SubfileName = generateSubFile(name,shell_name)
        condorDAG_file = "./condor_DAG_"+run+".dag"

        if maskChVFAT:
            ## Crate job files for vfat_masking
            os.system("python "+base_folder+"/VFAT_MaskMaker/run_step.py -r "+run)

            condorsubmit1_file = base_folder+"/VFAT_MaskMaker/CondorFiles/condor_step1_"+run+".submit"
            condorsubmit2_file = base_folder+"/VFAT_MaskMaker/CondorFiles/condor_step2_"+run+".submit"
        
            # prepare CONDOR DAG file:  will run step1 submission and then, at the step1 termination, will run step2 and finally the analysis        
            with open(condorDAG_file, "w") as DAG_file:
                DAG_file.write(
                    """
                    JOB A {step1VFAT_submit}
                    JOB B {step2VFAT_submit}
                    JOB C {analysis_submit}
                    PARENT A CHILD B
                    PARENT B CHILD C
                    """.format(step1VFAT_submit=condorsubmit1_file,step2VFAT_submit=condorsubmit2_file,analysis_submit=SubfileName))
            
            os.system("condor_submit_dag -dont_suppress_notification "+condorDAG_file)
        else:
            with open(condorDAG_file, "w") as DAG_file:
                DAG_file.write(
                    """
                    JOB A {analysis_submit}
                    """.format(analysis_submit=SubfileName))
            
            os.system("condor_submit_dag -dont_suppress_notification "+condorDAG_file)

