#!/bin/bash -l

# params
proj=facility
mail=aaron.ayllon.benitez@umu.se 
in=/mnt/picea/storage/reference/Picea-abies/v1.1/fasta/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa
out=$(realpath ../index)
outFile=$out/Pabies1.0-all.phase.gff3.CDSandLTR-TE_INDEX
img=/mnt/picea/projects/singularity/salmon-1.1.0.simg
runSalmonIndex=$(realpath ../UPSCb-common/pipeline/runSalmonIndex.sh)
## create the out dir
    if [ ! -d $out ]; then
      mkdir -p $out
    fi

## execute
      sbatch -A $proj --mail-user=$mail \
        -e $outFile.err -o $outFile.out -J aaron-generatingIndex \
        $runSalmonIndex -i $img $in $outFile
       
