VASP
====
VASP is a standalone perl software package designed to analyse sequenced pedigrees

To see available command line options type
./vasp.pl -h or ./vasp.pl -man

1) Three required input files:
i) ped files: ped file describing pedigree files (see http://www.sph.umich.edu/csg/abecasis/Pedstats/tour/input.html)

ii) vcf file: multisample vcf variant file
-if not inputting bam files (which is the preferred method) you must include the genotypes for each sample with the sample ordering the same as the ped file

iii) vep file: vep annotation file (see http://www.ensembl.org/info/docs/tools/vep/index.html)
-vep must have been run with the '--canonical' and '--gmaf' flags as canonical transcripts are utilised in all cases and allele frequency is used as a filter
-when multiple canonical transcripts overlap any causing protein changes are selected
-if only possess a vcf file you can generate vep file using the vcf as input (see vep documentation)
-if filtering on polyphen or sift you must run vep with '--poly b' and '--sift b'

Strongly encouraged to input these three arguments as leads to more accurate inheritance calls by using raw sequence calls; otherwise uses vcf genotype calls
iv) bam_list: File containing the bam file locations with format
Sample1_name  full_file_path1
Sample2_name  full_file_path2

v) ref: single reference fasta file used for alignment

vi) samtools: path to samtools binary (default=/usr/bin/samtools)

2) Filtering options for reducing search space. Most filters can be combined to together to further reduce search space 

i) vcf_cutoff (default=20): Variant quality cutoff to employ. This filter is ignored if no quality, i.e. quality is '.'. 

ii) min_num_aff (default=0): The minimum number of affected individuals containing the variant  

iii) polyphen (default=OFF): Filter on polyphen category (i.e. probably_damaging, possibly_damaging, or benign). Won't include variants without polyphen scores

iv) sift (default=OFF): Filter on polyphen category (i.e. deleterious or tolerated). Won't include variants without SIFT scores

v) comhet (default=OFF): Filter for variants found in genes determined to be either possible or definite compound hets (see vasp_output.pdf for full definition) 

vi) max_allele_freq (default=1): Filter for variants below the gmaf frequency.  Variants where the reference allele is the minor allele have their frequencies corrected to reflect this
 
vii) phase_var_num (default=2): Minimum number of consecutive variants to constitute genome block

viii) min_phase_block (default=10000bp): Minimum size in bp to constitute genome block

ix) denovo (default=OFF): Filter for candidate denovo mutations (defined as at least one case with ref parents and het child)

x) inheritance (default=OFF): Filter on inheritance type (options=ad,ar,xd,xr) for autosomal-dominant, autosomal-recessive, x-linked-dominant, and x-linked-recessive 

xi) coord (default=OFF):  Filter by genomics coordinate range (format=chr:start-end)

xii) gene_list (default=OFF): file containing list of genes, one on each line.  Genes can ensembl transcripts, ensembl genes, or hgnc names as long as vep was run with '--hgnc'

xiii) chrom (default=OFF): Filter by chromosome


3) Output options:
Columns reported depends on the constitution of the pedigree.
See vasp_output.pdf for details
Default output file = ./vasp.tsv but can be overidden with '-out' flag 
Output files are designed for loading into excel or libreoffice using a tab delimiter 

4) Add custom annotations and filter for matches (COMING...)

5) Sample usage:
To run with default parameters:
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep

Filter by chromosome:
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -chrom 22

Filter by genomic region:
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -coord 22:10000000-20000000

Filter by vcf cutoff:
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -vcf_cutoff 40

Filter for compound heterozygous genes:
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -comhet

Filter for polyphen 'probably_damaging'
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -polyphen probably_damaging

Filter for SIFT 'deleterious'
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -sift deleterious

Filter for max_allele_frequencies <=10%
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -max_allele_freq 0.1

Filter for genes in gene list file
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -gene_list sample/gene_list

Filter for inheritance type autosomal-dominant
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -inheritance ad

Change definition of genome block to need 5 variants and be 50000bp
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -phase_var_num 5 -min_phase_block 50000

To use bam file for generating zygosity, inheritance, parental alleles, etc (preferred method but beware TAKES LONG TIME TO GENERATE PILEUPS....)
./vasp.pl -ped sample/sample.ped -vcf sample/sample.vcf -vep sample/sample.vep -bam_list sample/bam_list -ref sample/GRCh37.fa -samtools /usr/bin/samtools







