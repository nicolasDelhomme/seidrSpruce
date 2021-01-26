#!/bin/bash
set -ex


# usage
USAGETXT=\
"
  Usage: $0 [options] <seidr file> <infomap input folder> <output filename>
"

source ../UPSCb-common/src/bash/functions.sh

isExec seidr


if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi


if [ ! -d $2 ]; then
 mkdir $2
fi

if [ ! -d $(dirname $3) ]; then
  abort "The third argument directory needs to exist"
fi

# run
seidr view $1 -c -d $"\t" | cut -f 1,2,26 > $2/edgelist.txt
seidr view $1 -N -d $"\t" | cut -f 1,2,26 > $2/edgelistIndex.txt

Infomap $2/edgelistIndex.txt $2

seidr resolve $2/edgelistIndex.tree -s $1 > $3
