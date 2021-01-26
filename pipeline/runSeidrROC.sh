#!/bin/bash
#SBATCH -A facility
#SBATCH -J aaron-SeidrROC
#SBATCH -t 12:00:00
#SBATCH -p core -n 1
#SBATCH --mem=64GB

set -ex

# Check for options
NGS=

# usage
USAGETXT=\
"
  Usage: $0 [options] <seidr file> <gold-standard> <output filename>
  
  Options: 
                -x define your negative golden standard
"

source ../UPSCb-common/src/bash/functions.sh

isExec seidr

# Get the options
while getopts x: option
do
        case "$option" in
      x) NGS="-x $OPTARG";;
    \?) ## unknown flag
		abort;;
        esac
done
shift `expr $OPTIND - 1`


if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be an existing file"
fi

if [ ! -d $(dirname $3) ]; then
  abort "The third argument directory needs to exist"
fi

# run
seidr roc -p 1000 -a -E 1 -g $2 -n $1 $NGS -o $3