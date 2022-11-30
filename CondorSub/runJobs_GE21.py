import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import re as regularExpression
from PFA_Analyzer_Utils import *

thisDir = os.path.dirname(os.path.abspath(__file__))                                                                                                                                                        
base_folder = os.path.dirname(thisDir)
#base_folder = os.path.expandvars('$PY3PFA')


def generateSubFile(outputName,shell_name):
    job_flavour = "tomorrow"
    SubfileName = base_folder+"/CondorSub/JobFiles/SubmitFile_"+str(outputName)+".sub"
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

def generateJobShell(name,absPath_config):
    main_command = f"python PFA_Analyzer_GE21.py {absPath_config}"    
    main_command = main_command + " \n"

    shell_script_name = base_folder+"/CondorSub/JobFiles/job_Run_"+str(name)+".sh"
    with open(shell_script_name, 'w') as fout:
        ####### Write the instruction for each job
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo %s_Run"+str(name)+"\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write(f"cd {base_folder} \n")
        # ## sourceing the right gcc version to compile the source code
        # fout.write("source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9.1.0/x86_64-centos7/setup.sh"+"\n")
        # fout.write("source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh"+"\n")
        # activate poetry venv
        fout.write("source /afs/cern.ch/user/f/fivone/.cache/pypoetry/virtualenvs/pfa-analyzer-vSUoOE05-py3.6/bin/activate"+"\n")
        ## Run the same job interval number of time
        fout.write("ClusterId=$1\n")
        fout.write("ProcId=$2\n")
        fout.write(main_command)
    return shell_script_name
if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='''Scripts runs PFA_Analyzer_GE21.py for many different input runs''',
        epilog="""Typical exectuion (from this folder)\n\t  python runJobs.py config/file1.yml config/file2.yml""",
        formatter_class=RawTextHelpFormatter
        )
    parser.add_argument('config', help='Analysis description file', nargs="*")
    args = parser.parse_args()
    config_list = args.config
    
    for filename in config_list:
        path = os.path.abspath(filename)
        name = path.replace(".yml","")
        name = name.split("/")[-1]

        shell_name = generateJobShell(name,path)
        SubfileName = generateSubFile(name,shell_name)
        condorDAG_file = "./condor_DAG_"+name+".dag"

        # if maskChVFAT:
        #     ## Crate job files for vfat_masking
        #     os.system("python "+base_folder+"/VFAT_MaskMaker/run_step.py -r "+run)

        #     condorsubmit1_file = base_folder+"/VFAT_MaskMaker/CondorFiles/condor_step1_"+run+".submit"
        #     condorsubmit2_file = base_folder+"/VFAT_MaskMaker/CondorFiles/condor_step2_"+run+".submit"
        
        #     # prepare CONDOR DAG file:  will run step1 submission and then, at the step1 termination, will run step2 and finally the analysis        
        #     with open(condorDAG_file, "w") as DAG_file:
        #         DAG_file.write(
        #             """
        #             JOB A {step1VFAT_submit}
        #             JOB B {step2VFAT_submit}
        #             JOB C {analysis_submit}
        #             PARENT A CHILD B
        #             PARENT B CHILD C
        #             """.format(step1VFAT_submit=condorsubmit1_file,step2VFAT_submit=condorsubmit2_file,analysis_submit=SubfileName))
            
        #     os.system("condor_submit_dag -dont_suppress_notification "+condorDAG_file)
        with open(condorDAG_file, "w") as DAG_file:
            DAG_file.write(
                """
                JOB A {analysis_submit}
                """.format(analysis_submit=SubfileName))
        
        os.system("condor_submit_dag -dont_suppress_notification "+condorDAG_file)

