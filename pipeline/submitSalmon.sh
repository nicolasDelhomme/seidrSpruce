#!/bin/bash -l

## be verbose and print
set -eux

proj=u2019016
mail=mist0069@student.umu.se

## process the argument
ref=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/salmon/Pabies1.0-all-phase.gff3.CDSandLTR-TE
bind=/mnt:/mnt
img=/mnt/picea/projects/singularity/salmon-0.14.1.simg 

## Zygotic
in=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/RNA-Seq/trimmomatic
out=/mnt/picea/projects/spruce/uegertsdotter/ZE-developmental-series/RNA-Seq/salmon

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file

for f in $(find $in -name "*_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
sbatch -A $proj --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  ../UPSCb-common/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done