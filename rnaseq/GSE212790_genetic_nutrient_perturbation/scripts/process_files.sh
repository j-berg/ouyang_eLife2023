#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/rutter-group2/j-berg/yeyun_2022_04_13_novogene/slurmjob-%j
#SBATCH --partition=notchpeak

# Rutter lab GPU node specs 
# - 40 cores
# - 192 GB of memory
# - 4x RTX2080TI GPUs

# Note: Change --gres command to up to `:4` to gain access to all 4 GPUs on node
# Note: --mem=0 reserves all node memory for the current job 
# See https://www.chpc.utah.edu/documentation/guides/gpus-accelerators.php for more information on options
echo "+ Checking CUDA devices..."
echo "Devices: $CUDA_VISIBLE_DEVICES"
nvidia-smi -L

NCPU=$(grep -c ^processor /proc/cpuinfo)
echo "No. CPUs: $NCPU"


# Activate conda environment
source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh

# conda create --name ahmad_seq
# conda install -c bioconda fastp multiqc kallisto
source activate xpresspipe2

echo "Testing installs..."
conda config --show channels
conda list

# Set directory variables
SCRUSER=/scratch/general/lustre/$USER
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
PROJDIR=/uufs/chpc.utah.edu/common/home/rutter-group2/j-berg/yeyun_2022_04_13_novogene/usftp21.novogene.com/raw_data
REF=/uufs/chpc.utah.edu/common/home/rutter-group2/j-berg/yeyun_2022_04_13_novogene/index_files

# Set experimental metadata variables
OVERHANG=149
ADAPTER5=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
ADAPTER3=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG
EXPNAME=yeyun_rnaseq

# release-105 (2021-06-20)
echo "Building genome index"
mkdir -p $REF
mkdir -p $REF/genome

cd $REF

# DNA FASTAs
mkdir -p $REF/genome_fastas
cd $REF/genome_fastas
for J in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito; do curl -OL http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.${J}.fa.gz; done
gzip -d *.gz

# Transcriptome GTF
cd $REF
curl -OL http://ftp.ensembl.org/pub/release-106/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.106.gtf.gz
gzip -d *gtf.gz
mv Saccharomyces_cerevisiae.R64-1-1.106.gtf transcripts.gtf
cd $SCRDIR

# Build STAR index files for XPRESSpipe
echo "Building STAR index..."
xpresspipe makeReference -o $REF -f $REF/genome_fastas -g $REF/transcripts.gtf --sjdbOverhang $OVERHANG --genome_size 12000000

# Process files
echo "Processing sequence files..."
mkdir -p $SCRDIR/input_files
mkdir -p $SCRDIR/output_files
cp $PROJDIR/**/*.fq.gz $SCRDIR/input_files

xpresspipe peRNAseq -i $SCRDIR/input_files -o $SCRDIR/output_files -r $REF --gtf $REF/transcripts.gtf -e $EXPNAME -a $ADAPTER5 $ADAPTER3 --sjdbOverhang $OVERHANG --quantification_method htseq --remove_rrna
