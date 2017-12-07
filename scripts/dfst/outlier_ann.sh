#!bin/bash

my_bed=/home/elias/program/bedtools2/bin/bedtools
my_gff=/home/elias/analysis/data/genome/unsplit_sorted.gff
my_snps=/home/elias/analysis/data/dfst/outliers/zshared_out_all
my_out=/home/elias/analysis/data/dfst/outliers/zshared_all_ann

$my_bed intersect \
-a $my_gff \
-b $my_snps  \
-wa | uniq > $my_out
