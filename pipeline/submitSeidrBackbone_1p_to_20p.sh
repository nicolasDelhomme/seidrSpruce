#!/bin/bash -l
# paramss
proj=facility
mail=aaron.ayllon.benitez@umu.se
in=$(realpath ../data/seidr/)
inFile="results/aggregated.sf"
outFile="aggregate-backbone"
runSeidrBackbone=$(realpath ../UPSCb-common/pipeline/runSeidrBackbone.sh)

module load bioinfo-tools seidr-devel

backboneFolder=$in/backbone
if [ ! -d $backboneFolder ]; then
  mkdir $backboneFolder
else
  rm -r $backboneFolder
  mkdir $backboneFolder
fi


sdThresholds="2.33 2.05 1.88 1.75 1.64 1.55 1.48 1.41 1.34 1.28 1.23 1.17 1.13 1.08 1.04 0.99 0.95 0.92 0.88 0.84"
#sdThresholds="0.6744898 0.5244005 0.3853205 0.2533471 0.1256613 0.0000000"
count=1
for i in $sdThresholds; do

out=$outFile-p$count.sf

  sbatch -A $proj --mail-user=$mail -t 8:00:00\
        -e $backboneFolder/$out.err -o $backboneFolder/$out.out \
        $runSeidrBackbone $in/$inFile $i $backboneFolder/$out
count=$(( $count+1 ));
done


