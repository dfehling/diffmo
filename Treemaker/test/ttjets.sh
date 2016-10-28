#!/bin/bash

#Weight for MC 
weight=0.22984      #Powheg <700     : 1.00/21675970*252.89*19700
weight7=0.11959     #Powheg 700-1000 : 0.074/3082812*252.89*19700    : 0.074 * 1.571
weight10=0.05584    #Powheg 1000+    : 0.014/1249111*252.89*19700    : 0.014 * 3.877

num_jobs=25

rm listofjobs.txt
rm commands.cmd
tar czvfh tarball.tgz ../CONDOR/* ../run_top_xs_TreeMaker.py ../top_xs_TreeMaker*.py

for i in `seq $num_jobs`;
do
# for s in "     " "JERup" "JERdn" "JESup" "JESdn"
for s in ""
do

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/  --outfile=ttjets_powhegLT7_pdf_${s}_sec_${i} --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=$weight --totalSec=$num_jobs --isMCatNLO=0 --useCondor=1 --invMassCut=700 --sec=${i}  --useSyst=${s} >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/  --outfile=ttjets_powhegNC_pdf_${s}_sec_${i} --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=$weight --totalSec=$num_jobs --isMCatNLO=0 --useCondor=1 --sec=${i}  --useSyst=${s} >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets7/  --outfile=ttjets_powheg7_pdf_${s}_sec_${i} --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=$weight7 --totalSec=$num_jobs --isMCatNLO=0 --useCondor=1 --sec=${i}  --useSyst=${s} >> listofjobs.txt

echo python ./tardir/run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets10/  --outfile=ttjets_powheg10_pdf_${s}_sec_${i} --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=$weight10 --totalSec=$num_jobs --isMCatNLO=0 --useCondor=1 --sec=${i}  --useSyst=${s} >> listofjobs.txt

done

done

runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz listofjobs.txt commands.cmd

exit

#Test job
#python run_top_xs_TreeMaker.py --dirs=/uscms_data/d3/dfehling/NTUPLES/ttjets/ --outfile=test.root --isMC=1 --includeTrigger=0 --includePileup=0 --includePDF=1 --doUnfold=1 --unfoldWeight=0.2 --totalSec=25 --isMCatNLO=0 --sec=1  --useSyst=JERup -N 50000