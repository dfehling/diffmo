#!/bin/bash

num_jobs=25

rm listofjobs.txt
tar czvfh tarball.tgz ../CONDOR/* ../run_top_xs_TreeMaker.py ../top_xs_TreeMaker*.py

for i in `seq $num_jobs`;
do

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/"               "--outfile=ttjets_powhegNC_muon_noTrig_noPU_sec_${i} --config=muon --isMC=1 --includeTrigger=0 --includePileup=0 --totalSec=$num_jobs --sec=${i} >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/  --outfile=ttjets_mcatnlo_muon_noTrig_noPU_sec_${i}"  "--config=muon --isMC=1 --includeTrigger=0 --includePileup=0 --totalSec=$num_jobs --sec=${i} --isMCatNLO=1 >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/"               "--outfile=ttjets_powhegNC_muon_Trig_PU_sec_${i}"     "--config=muon --isMC=1 --includeTrigger=1 --includePileup=1 --totalSec=$num_jobs --sec=${i} --triggerFile=./tardir/temp --pileupFile=./tardir/muon_pu_100 >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/  --outfile=ttjets_mcatnlo_muon_Trig_PU_sec_${i}"      "--config=muon --isMC=1 --includeTrigger=1 --includePileup=1 --totalSec=$num_jobs --sec=${i} --triggerFile=./tardir/temp --pileupFile=./tardir/muon_pu_100 --isMCatNLO=1 >> listofjobs.txt

done
exit
