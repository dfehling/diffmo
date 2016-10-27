#!/bin/bash

#252.89×19700×282÷284÷32345466
weight=0.15294

num_jobs=25

rm listofjobs.txt
rm commands.cmd
tar czvfh tarball.tgz ../CONDOR/* ../run_top_xs_TreeMaker.py ../top_xs_TreeMaker*.py

for i in `seq $num_jobs`;
do

for s in "     " "JERup" "JERdn" "JESup" "JESdn"
do

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/  --outfile=ttjets_mcatnlo_252_noTrig_noPU_withPDF_unfold_sec_${i}_${s} --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=$weight --totalSec=$num_jobs --isMCatNLO=1 --useCondor=1 --sec=${i}  --useSyst=${s} >> listofjobs.txt

done

done

runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz listofjobs.txt commands.cmd

exit

#python run_top_xs_TreeMaker.py --dirs=/uscmst1b_scratch/lpc1/3DayLifetime/dfehling/mcatnlo/  --outfile=test.root --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=0.2 --totalSec=25 --isMCatNLO=1 --sec=1  --useSyst=JERup -N 50000