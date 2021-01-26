#!/bin/bash
set -ex

bb=$(realpath "../data/seidr/backbone/aggregate-backbone-p1.sf")
imFolder=$(realpath "../data/seidr/infomap-bb-p1")
out=$imFolder/clusterResolve.txt

module load bioinfo-tools seidr-devel
module load bioinfo-tools InfoMap

./runInfoMap.sh $bb $imFolder $out