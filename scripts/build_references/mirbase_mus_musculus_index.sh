#!/bin/bash

mkdir -p  mmu
# bash species_index.sh ../v22/mature.fa "Mus musculus" mmu mature > mmu/bt2_mature.log

bash species_index.sh ../v22/hairpin.fa "Mus musculus" mmu hairpin > mmu/bt2_hairpin.log
