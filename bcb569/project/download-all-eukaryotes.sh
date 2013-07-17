#!/bin/bash

#----- Animals

# human
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/protein/protein.fa.gz
mv protein.fa.gz H_sampiens.protein.fa.gz
gunzip H_sampiens.protein.fa.gz

# mouse
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/protein/protein.fa.gz
mv protein.fa.gz M_musculus.protein.fa.gz
gunzip M_musculus.protein.fa.gz

# zebrafish
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/D_rerio/protein/protein.fa.gz
mv protein.fa.gz D_rerio.protein.fa.gz
gunzip D_rerio.protein.fa.gz

# fruit fly
for dir in CHR_2 CHR_3 CHR_4 CHR_Un CHR_X; do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Drosophila_melanogaster/RELEASE_5_30/$dir/*.faa; done
cat N*.faa > D_melanogaster.protein.fa
rm N*.faa



#----- Fungi

# bread yeast
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Fungi/Saccharomyces_cerevisiae_uid128/*.faa
cat NC_001*.faa > S_cerevisiae.protein.fa
rm NC_001*.faa

# E. gossypii
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Fungi/Eremothecium_gossypii_uid10623/*.faa
cat *.faa > E_gossypii.protein.fa
rm *.faa

# A. nidulans
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Fungi/Aspergillus_nidulans_FGSC_A4_uid13961/*.faa
cat *.faa > A_nidulans.protein.fa
rm *.faa



#----- Plants

# maize
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/PLANTS/Zea_mays/Gnomon_v2/*_prots.fa
cat *prots.fa > Z_mays.protein.fa
rm *prots.fa

# arabidopsis
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Arabidopsis_thaliana/GNOMON/Arabidopsis_gnomon_prots.fa.gz
gunzip Arabidopsis_gnomon_prots.fa.gz
mv Arabidopsis_gnomon_prots.fa A_thaliana.protein.fa

# moss
wget ftp://ftp.plantgdb.org/download/Genomes/PpGDB/Ppatens_152_peptide.fa.gz
gunzip Ppatens_152_peptide.fa.gz
mv Ppatens_152_peptide.fa P_patens.protein.fa



#----- Protists

# P. falciparum
wget ftp://ftp.ncbi.nih.gov/genomes/Protozoa/Plasmodium_falciparum/*.faa
cat *.faa > P_falciparum.protein.fa
rm *.faa

# C. parvum
wget ftp://ftp.ncbi.nih.gov/genomes/Protozoa/Cryptosporidium_parvum/*.faa
cat *.faa > C_parvum.protein.fa
rm *.faa

