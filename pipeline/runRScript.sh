#!/bin/bash
#SBATCH -A facility
#SBATCH -J ICcomputing
#SBATCH -t 12:00:00
#SBATCH -p core -n 1
#SBATCH --mem=64GB

set -ex

# usage
USAGETXT=\
"
  Usage: $0 [options] <script file> <var1> <var2>
  
"
if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

module load R/4.0.2

Rscript $1 $2 $3
