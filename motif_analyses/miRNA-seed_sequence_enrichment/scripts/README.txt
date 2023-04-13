# Here we want to get a coverage of k-mers and see to which miRNA the 6-mer (seed) belongs to. We first mapped 6-mers to miRBase and used the output to link miRNA names. 
# !!! We need to accept only reverse complements when mapping with bowtie !!!

# get seed motif, strand and miRNA name
cat bowtie-miRBase-ams-a2-hela_trimmed_single-sum-flanked30.6-mers_1.5_enriched.0_mismatche-nofw.out | awk '{print $1 "\t" $2 "\t" $3}' > miRNA-seed-targets-6-mers_1.5_enriched.0_mismatch-nofw.out.tab

# get coverage per seed
python miRNA-seed_coverage-norm-as_track.py ams-a2-hela_trimmed_single-sum-flanked30.fasta miRNA-seed-targets-6-mers_1.5_enriched.0_mismatch-nofw.out.tab ams-a2-hela_trimmed_single-sum-flanked30.miRNA-seed_enrichment.tab

# merge
cat *.coverage.txt > ams-a2-hela_trimmed_single-sum-flanked30.6-mers_1.5_enriched.ALL.seed.coverage.txt

# get RNA-maps and heatMaps
Motif-seed-enrichment.R
