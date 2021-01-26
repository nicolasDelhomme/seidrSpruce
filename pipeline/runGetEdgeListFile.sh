#!/bin/bash -l

backboneFolder=$(realpath ../data/seidr/backbone/)
backfile=aggregate-backbone-p1.sf
module load bioinfo-tools seidr-devel

# To get the specific columns by using awk

## columns: source, target, type-edge

 seidr view -c $backboneFolder/$backfile | \
awk '{print $1,$2,$3}' > $backboneFolder/backbone-edgelist.txt
 
## columns: source, target, type-edge, irp_score

 seidr view -c $backboneFolder/$backfile | \
awk '{      
            str1=$1
            str2=$2
            str3=$3
            str=$16
            split(str, arr,";")
            print str1,str2,str3,arr[1]
            }' > $backboneFolder/backbone-edgelist-weighted.txt
