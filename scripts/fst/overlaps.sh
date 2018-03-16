#!bin/bash

my_bed=/home/elias/program/bedtools2/bin/bedtools
my_bb=/home/elias/analysis/data/fst/individual_pbs/bb_out.bed
my_vb=/home/elias/analysis/data/fst/individual_pbs/vb_out.bed
my_pb=/home/elias/analysis/data/fst/individual_pbs/pb_out.bed
my_sj=/home/elias/analysis/data/fst/individual_pbs/sj_out.bed
my_bnp=/home/elias/analysis/data/fst/individual_pbs/bnp_out.bed
my_out=/home/elias/analysis/data/fst/individual_pbs/overlaps

$my_bed multiinter \
-i $my_bb $my_vb $my_pb $my_sj $my_bnp \
-names BB VB PB SJ BNP > $my_out
