#!/bin/bash -l

# params
proj=facility
mail=aaron.ayllon.benitez@umu.se 
in=$(realpath ../data/seidr/results)
outFile=$(realpath ../data/seidr/results/aggregated)
runSeidrAggregate=$(realpath ../UPSCb-common/pipeline/runSeidrAggregate.sh)
MEM=64G
files=""

module load bioinfo-tools seidr-devel


for f in $(find $outFile.*); do
  rm $f;
done


for folder in $(find $in -maxdepth 1 -mindepth 1 -type d); do
  if [ "$(basename $folder)" != "el-ensemble" ]; then
    files=$files" "$folder/$(basename $folder.sf);
  fi
done

## create the out dir

sbatch -A $proj --mail-user=$mail\
        --mem=$MEM -e $outFile.err -o $outFile.out \
        -J aaron-SeidrAggregate \
        $runSeidrAggregate $in $files