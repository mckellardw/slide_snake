# build_references
Scripts to build species-specific bowtie2 indices from miRbase

Example of how to do this for mouse:
```
mkdir -p  mmu
bash species_index.sh ../v22/mature.fa "Mus musculus" mmu mature > mmu/bt2_mature.log

bash species_hairpin_index.sh ../v22/hairpin.fa "Mus musculus" mmu hairpin > mmu/bt2_hairpin.log
``` 