#!/bin/bash -l

## be verbose and print
set -eux

# params
proj=facility
mail=nicolas.delhomme@umu.se
ref=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/salmon/Pabies1.0-all-phase.gff3.CDSandLTR-TE
bind=/mnt:/mnt
img=/mnt/picea/projects/singularity/salmon-1.1.0.simg 
in=/mnt/picea/projects/spruce/facility/trimmomatic
out=$(realpath ../data/salmon)

declare -A DATASET=(
  [atlas]="atlas"
  [diurnal]="diurnal"
  [diurnal2]="diurnal"
  [flakaliden]="flakaliden"
  [seasonal]="seasonal")
  
declare -A RUN=(
  [atlas]=1
  [diurnal]=1
  [diurnal2]=1
  [flakaliden]=1
  [seasonal]=1
  )

default="*_trimmomatic_1.fq.gz"

declare -A PATTERN=(
  [atlas]=$default
  [diurnal]=$default
  [diurnal2]="*_trimmomatic.fq.gz"
  [flakaliden]=$default
  [seasonal]=$default
)

default="../UPSCb-common/pipeline/runSalmon.sh"

declare -A TOOL=(
  [atlas]=$default
  [diurnal]=$default
  [diurnal2]="../UPSCb-common/pipeline/runSalmonSE.sh ADD-EXTRA-ARGS-THERE fragment length mean and fragment length sd"
  [flakaliden]=$default
  [seasonal]=$default
)

default=1
declare -A PE=(
  [atlas]=$default
  [diurnal]=$default
  [diurnal2]=0
  [flakaliden]=$default
  [seasonal]=$default
)

# helper
source ../UPSCb-common/src/bash/functions.sh

# variables
FORCE=0

# usage
USAGETXT=\
"
$0 [options]
  
  Options -f force the rebuild
"

while getopts f option
do
    case "$option" in
    f) FORCE=1;;
    \?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

# Dataset
for DSET in "${!DATASET[@]}"; do
  if [ ${RUN[$DSET]} -eq 1 ]; then
  
    outDir=$out/${DATASET[$DSET]} 
    
    ## create the out dir
    if [ ! -d $outDir ]; then
      mkdir -p $outDir
    fi
  
    ## for every file
    REV=
    for f in $(find $in/${DATASET[$DSET]} -name ${PATTERN[$DSET]}); do
      if [ ${PE[$DSET]} -eq 1 ]; then
        fnam=$(basename ${f/_1.fq.gz/})
        REV=$in/$DSET/${fnam}_2.fq.gz
      else
        fnam=$(basename ${f/.fq.gz/})
      fi
      
      ## execute
      echo sbatch -A $proj --mail-user=$mail \
        -e $outDir/$fnam.err -o $outDir/$fnam.out -J salmon.$fnam \
        ${TOOL[$DSET]} -b $bind \
        -i $img $ref $f $REV $outDir
    done
  fi
done
  



