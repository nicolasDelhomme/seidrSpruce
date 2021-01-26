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
fi


sdThresholds="2.33 2.05 1.88 1.75 1.64 1.55 1.48 1.41 1.34 1.28"
count=1
for i in $sdThresholds; do

out=$outFile-p$count.sf

sbatch -A $proj --mail-user=$mail -t 8:00:00\
        -e $backboneFolder/$out.err -o $backboneFolder/$out.out \
        $runSeidrBackbone $in/$inFile $i $backboneFolder/$out
count=$(( $count+1 ));
done


