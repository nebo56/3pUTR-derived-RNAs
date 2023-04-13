# get miRNA seed sequence enrichment target peaks
genome=./GRCh38.p10.genome.fa

target=$1	# in BED format
seeds=$2	# reverste complement seed sequence \t name 

# flank target and control positions for 30 nt
python ./scripts/flankBEDpositionsCustom.py ${target} ${target}.flank30.bed 30 30

# get fasta sequence
bedtools getfasta -s -fi $genome -bed ${target}.flank30.bed -fo ${target}.flank30.fasta

# get coverage per seed - this script will clip first and last nt of the 7-mer and keep the most enriched one
python ./scripts/miRNA-7mer-clipping_seed_coverage-norm-as_track.py ${target}.flank30.fasta ${seeds} ${target}.flank30-miRNA-seed_enrichment.tab

# merge
cat *.coverage.txt > ${target}.seed.coverage.tab
rm *.coverage.txt

# get RNA-maps and heatMaps
Rscript ./scripts/Motif-seed-enrichment.R ${target}.seed.coverage.tab

