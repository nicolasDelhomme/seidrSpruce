#!/bin/bash -l

# paramss
proj=facility
mail=aaron.ayllon.benitez@umu.se
time="8:00:00"
in=$(realpath ../data/seidr/)
inFile="results/aggregated.sf"
backbone="backbone/aggregate-backbone-p"
outFile="aggregated"
runSeidrRoc=$(realpath ../UPSCb-common/pipeline/runSeidrRoc.sh)
positiveGold="/mnt/picea/storage/reference/goldStandard/Picea-abies_KEGG-based-positive-gold-standard.tsv"
negativeGold="/mnt/picea/storage/reference/goldStandard/Picea-abies_KEGG-based-negative-gold-standard.tsv"

module load bioinfo-tools seidr-devel


rocFolder=$(realpath $in/roc)

if [ ! -d $rocFolder ]; then
  mkdir $rocFolder
#else
 # rm -r $rocFolder
  #mkdir $rocFolder
fi


echo sbatch -A $proj --mail-user=$mail -t $time --mem=64G\
        -e $rocFolder/$outFile.err -o $rocFolder/$outFile.out \
        $runSeidrRoc $in/$inFile $positiveGold $negativeGold $rocFolder/$outFile.roc
        
        
for ((i=25;i<=50;i+5)); 
do 
  echo $i
  backboneFile=$backbone$i.sf
  sbatch -A $proj --mail-user=$mail -t $time\
        -e $rocFolder/$outFile-backbone-p$i.err -o $rocFolder/$outFile-backbone-p$i.out \
        $runSeidrRoc $in/$backboneFile $positiveGold $negativeGold $rocFolder/$outFile-backbone-p$i.roc
done
