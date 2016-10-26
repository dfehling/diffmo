# -*- sh -*- # for font lock mode
# variable definitions
- env = cd /uscms_data/d3/dfehling/sl6/FINAL/CMSSW_5_3_28_patch1; eval `scramv1 runtime -sh`; cd -
- tag = 
- output = outputFile=
- tagmode = none
- tarfile = /uscms_data/d3/dfehling/sl6/FINAL/CMSSW_5_3_28_patch1/src/Analysis/Treemaker/test/CONDOR/tarball.tgz
- untardir = tardir
- copycommand = cp

# Sections listed
output_$(JID)        python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/ --outfile=test --sec=1 --totalSec=20 --isMC=1 --includePDF=1 --includeTrigger=1 --triggerFile=./tardir/trigger_eff_test --includePileup=1 --pileupFile=./tardir/hadPU100nom_test --doUnfold=1 --unfoldWeight=0.22339 -N 50000

