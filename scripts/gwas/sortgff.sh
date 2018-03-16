my_gff=/home/elias/analysis/data/genome/unsplit_ncbi.gff
my_out=/home/elias/analysis/data/genome/unsplit_sorted.gff

sort -k1,1V -k4,4V $my_gff | grep -Ev "^N[WC]_*"  > $my_out
