#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0690617/slurmjob-%j
#SBATCH --account=rutter-gpu-np
#SBATCH --partition=rutter-gpu-np
#SBATCH --gres=gpu:2080ti:4
#SBATCH --mem=0

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


# Activate conda environment
# $ conda create -n pe_seq_2021_02_03
# $ conda activate pe_seq_2021_02_03
# $ conda config --add channels bioconda
# $ conda config --add channels conda-forge
# $ conda install fastp star samtools bedtools htseq multiqc
echo "Loading conda environment..."
source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
source activate pe_seq_2021_02_03

# Test installs
echo "Testing installs..."
conda config --show channels
conda list

fastp -v
STAR --version
samtools --version
bedtools --version
htseq-count --version
multiqc -v

# Set instance variables
echo "Initializing environment..."
SCRUSER=/scratch/general/lustre/u0690617
SCRDIR=/scratch/general/lustre/u0690617/$SLURM_JOBID

# File pattern: {ID}_210128_A00421_0278_BHT2YWDSXY_{SXX}_L001_R{1/2}_001.fastq.gz
FILES=(\
18607X1_210128_A00421_0278_BHT2YWDSXY_S16_L001 \
18607X2_210128_A00421_0278_BHT2YWDSXY_S17_L001 \
18607X3_210128_A00421_0278_BHT2YWDSXY_S18_L001 \
18607X4_210128_A00421_0278_BHT2YWDSXY_S19_L001 \
18607X5_210128_A00421_0278_BHT2YWDSXY_S20_L001 \
18607X6_210128_A00421_0278_BHT2YWDSXY_S21_L001 \
18607X7_210128_A00421_0278_BHT2YWDSXY_S22_L001 \
18607X8_210128_A00421_0278_BHT2YWDSXY_S23_L001 \
18607X9_210128_A00421_0278_BHT2YWDSXY_S24_L001 \
18607X10_210128_A00421_0278_BHT2YWDSXY_S25_L001 \
18607X11_210128_A00421_0278_BHT2YWDSXY_S26_L001 \
18607X12_210128_A00421_0278_BHT2YWDSXY_S27_L001 \
18607X13_210128_A00421_0278_BHT2YWDSXY_S28_L001 \
18607X14_210128_A00421_0278_BHT2YWDSXY_S29_L001 \
18607X15_210128_A00421_0278_BHT2YWDSXY_S30_L001 \
18607X16_210128_A00421_0278_BHT2YWDSXY_S31_L001 \
18607X17_210128_A00421_0278_BHT2YWDSXY_S32_L001 \
18607X18_210128_A00421_0278_BHT2YWDSXY_S33_L001 \
18607X19_210128_A00421_0278_BHT2YWDSXY_S34_L001 \
18607X20_210128_A00421_0278_BHT2YWDSXY_S35_L001 \
18607X21_210128_A00421_0278_BHT2YWDSXY_S36_L001 \
18607X22_210128_A00421_0278_BHT2YWDSXY_S37_L001 \
18607X23_210128_A00421_0278_BHT2YWDSXY_S38_L001 \
18607X24_210128_A00421_0278_BHT2YWDSXY_S39_L001 \
18607X25_210128_A00421_0278_BHT2YWDSXY_S40_L001 \
18607X26_210128_A00421_0278_BHT2YWDSXY_S41_L001 \
18607X27_210128_A00421_0278_BHT2YWDSXY_S42_L001 \
18607X28_210128_A00421_0278_BHT2YWDSXY_S43_L001 \
18607X29_210128_A00421_0278_BHT2YWDSXY_S44_L001 \
18607X30_210128_A00421_0278_BHT2YWDSXY_S45_L001 \
)


REF=/scratch/general/lustre/u0690617/references/sce_2021_02_03
FASTA=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
GTF=Saccharomyces_cerevisiae.R64-1-1.100.gtf


mkdir -p $REF
cd $REF


# Get reference files
mkdir -p $REF/fasta
cd $REF/fasta
echo "Downloading source FASTA file(s)..."
curl -OL ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gzip -d *gz
cd $REF

echo "Downloading source GTF file..."
curl -OL ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz
gzip -d *gz


# Generate genome reference
echo "Building genome index..."
mkdir -p $REF/genome
STAR --runMode genomeGenerate \
  --genomeDir $REF/genome \
  --genomeFastaFiles $REF/fasta/$FASTA \
  --sjdbGTFfile $REF/$GTF \
  --runThreadN 40 \
  --sjdbOverhang 150 \
  --genomeSAindexNbases 10


# Prep output locations
echo "Prepping output locations..."
mkdir -p $SCRDIR/input
mkdir -p $SCRDIR/output

mkdir -p $SCRDIR/output/preprocess
mkdir -p $SCRDIR/output/alignment
mkdir -p $SCRDIR/output/postprocess
mkdir -p $SCRDIR/output/counts
mkdir -p $SCRDIR/output/counts_dedup


echo "Copying input files..."
cp $FILEDIR/*fastq.gz $SCRDIR/input


cd $SCRDIR
# Process files
for FILE in ${FILES[@]}; do

  echo "Processing $FILE"

  # Clean FASTQ files
  echo "Preprocessing..."
  fastp --thread 40 -l 30 -q 28 \
    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -i $SCRDIR/input/${FILE}_R1_001.fastq.gz \
    -I $SCRDIR/input/${FILE}_R2_001.fastq.gz \
    -o $SCRDIR/output/preprocess/${FILE}_R1_001_out.fastq.gz \
    -O $SCRDIR/output/preprocess/${FILE}_R2_001_out.fastq.gz \
    -j $SCRDIR/output/preprocess/${FILE}.json \
    -h $SCRDIR/output/preprocess/${FILE}.html

    # Align and quantify FASTQ files
    echo "Aligning & quantifying..."
    STAR --runThreadN 40 --sjdbOverhang 150 \
      --readFilesCommand zcat \
      --outSAMunmapped Within \
      --outSAMtype BAM Unsorted \
      --quantMode TranscriptomeSAM \
      --genomeDir $REF/genome \
      --sjdbGTFfile $REF/$GTF \
      --readFilesIn $SCRDIR/output/preprocess/${FILE}_R1_001_out.fastq.gz $SCRDIR/output/preprocess/${FILE}_R2_001_out.fastq.gz \
      --outFileNamePrefix $SCRDIR/output/alignment/${FILE}_

    # Postprocess
    echo "Postprocessing..."
    samtools sort --threads 40 -n \
            -o $SCRDIR/output/postprocess/${FILE}_nameSorted.bam \
            $SCRDIR/output/alignment/${FILE}_Aligned.out.bam
    samtools fixmate --threads 40 -m \
            $SCRDIR/output/postprocess/${FILE}_nameSorted.bam \
            $SCRDIR/output/postprocess/${FILE}_nameSorted_fixed.bam
    bedtools intersect -abam \
            $SCRDIR/output/postprocess/${FILE}_nameSorted_fixed.bam \
            -b $REF/rrna_ref.bed \
            -v > $SCRDIR/output/postprocess/${FILE}_nameSorted_fixed_rrnaDepl.bam
    samtools sort --threads 40 \
            -o $SCRDIR/output/postprocess/${FILE}_output.bam \
            $SCRDIR/output/postprocess/${FILE}_nameSorted_fixed_rrnaDepl.bam
    samtools index -@ 40 \
            $SCRDIR/output/postprocess/${FILE}_output.bam

    samtools markdup --threads 40 -s \
            $SCRDIR/output/postprocess/${FILE}_output.bam \
            $SCRDIR/output/postprocess/${FILE}_sorted_dedupMarked.bam
    samtools markdup --threads 40 -s -r \
            $SCRDIR/output/postprocess/${FILE}_output.bam \
            $SCRDIR/output/postprocess/${FILE}_sorted_dedupRemoved.bam
    samtools index -@ 40 \
            $SCRDIR/output/postprocess/${FILE}_sorted_dedupRemoved.bam

    # Quantification
    echo "Counting..."
    htseq-count -q -f bam -m intersection-nonempty -t exon \
      -i gene_id -r pos -s no \
      $SCRDIR/output/postprocess/${FILE}_output.bam \
      $REF/$GTF > $SCRDIR/output/counts/${FILE}.tsv

    # Quantify deduplicated alignments
    echo "Counting deduplicated reads..."
    htseq-count -q -f bam -m intersection-nonempty -t exon \
      -i gene_id -r pos -s no \
      $SCRDIR/output/postprocess/${FILE}_sorted_dedupRemoved.bam \
      $REF/$GTF > $SCRDIR/output/counts_dedup/${FILE}_deduplicated.tsv;

done


# Perform QC
echo "Running final QC..."
multiqc $SCRDIR/output -i job-$SLURM_JOBID -o $SCRDIR/output


# Clean-up
echo "Cleaning up files..."
OUTPUT=/uufs/chpc.utah.edu/common/home/rutter-group1/jordan/yeyun_18607R/processed
mkdir -p $OUTPUT
cd $OUTPUT
cp $SCRDIR/output/postprocess/*_output.bam $OUTPUT
cp $SCRDIR/output/postprocess/*_output.bam.bai $OUTPUT
cp $SCRDIR/output/postprocess/*_sorted_dedupMarked.bam $OUTPUT
cp $SCRDIR/output/postprocess/*_sorted_dedupRemoved.bam $OUTPUT
cp $SCRDIR/output/postprocess/*_sorted_dedupRemoved.bam.bai $OUTPUT
cp $SCRDIR/output/counts/*.tsv $OUTPUT
cp $SCRDIR/output/counts_dedup/*.tsv $OUTPUT
cp $SCRDIR/output/*.html $OUTPUT

# rm -rf $SCRDIR
cp /uufs/chpc.utah.edu/common/home/u0690617/slurmjob-$SLURM_JOBID $OUTPUT
