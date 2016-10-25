#!/bin/bash

num_jobs=20

for i in `seq $num_jobs`;
do

python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets/' --outfile=output/ttjets_powhegNC_muon_noTrig_noPU_sec_${i} --config='muon' --sec=${i} --totalSec=$num_jobs --isMC --includeTrigger=0 --includePileup=0 &

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/'  --outfile=output/ttjets_mcatnlo_muon_noTrig_noPU_sec_${i} --config='muon' --sec=${i} --totalSec=$num_jobs --isMC=1 --includeTrigger=0 --includePileup=0 --isMCatNLO=1 &

python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets/' --outfile=output/ttjets_powhegNC_muon_Trig_PU_sec_${i} --config='muon' --sec=${i} --totalSec=$num_jobs --isMC=1 --includeTrigger=1 --triggerFile="temp" --includePileup=1 --pileupFile="muon_pu_100" &

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/'  --outfile=output/ttjets_mcatnlo_muon_Trig_PU_sec_${i} --config='muon' --sec=${i} --totalSec=$num_jobs --isMC=1 --includeTrigger=1 --triggerFile="temp" --includePileup=1 --pileupFile="muon_pu_100" --isMCatNLO=1 &

done
exit
