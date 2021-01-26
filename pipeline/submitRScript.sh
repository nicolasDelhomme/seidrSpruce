# paramss
proj=facility
mail=aaron.ayllon.benitez@umu.se
time="2-00:00:00"
in=$(realpath "../R/computeICRes.R")
#annotation="tair goa_uniprot_all"
annotation="goa_uniprot_all"
outFolder=$(realpath "../../data")
runR=$(realpath "./runRScript.sh")

for i in $annotation; do

out=$outFolder/$i-IC.csv

sbatch -A $proj -w picea --mail-user=$mail -t $time\
        -e $out.err -o $out.out \
        $runR $in $i $out

done
