#!/bin/bash

#252.89×19700×282÷284÷32345466
weight=0.15294
# weight7=0.11623
# weight10=0.05427

num_jobs=20

for i in `seq $num_jobs`;

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/'  --outfile=output/ttjets_mcatnlo_252_pdfup_sec_${i} --sec=${i} --totalSec=$num_jobs --isMC=1 --includePDF=1 --unfoldWeight=$weight   --doUnfold=1 --isMCatNLO=1
python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/'  --outfile=output/ttjets_mcatnlo_252_pdfdn_sec_${i} --sec=${i} --totalSec=$num_jobs --isMC=1 --includePDF=-1 --unfoldWeight=$weight   --doUnfold=1 --isMCatNLO=1
do
for s in "" "JERup" "JERdn" "JESup" "JESdn"
do

python run_top_xs_TreeMaker.py --dirs='/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/'  --outfile=output/ttjets_mcatnlo_252_${s}_sec_${i} --sec=${i} --totalSec=$num_jobs --isMC=1 --includePDF=0 --unfoldWeight=$weight   --doUnfold=1 --useSyst=${s} --isMCatNLO=1

done
# python run_top_xs_TreeMaker.py --dirs='/uscms_data/d3/dfehling/NTUPLES/ttjets/' --outfile=output/ttjets_252_powhegLT7_${s}_sec_${i}    --sec=${i} --totalSec=$num_jobs --isMC=1 --includePDF=0 --includeTrigger=1 --triggerFile="trigger_eff_test" --includePileup=1 --pileupFile="hadPU100nom_test" --doUnfold=1 --unfoldWeight=$weight --invMassCut=700 --useSyst=${s} &
done
exit
