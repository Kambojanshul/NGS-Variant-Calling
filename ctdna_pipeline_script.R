# ctDNA Variant Analysis Pipeline in R

# Define directories (replace with actual paths)
input.dir <- "/home/anshul/software/ctdna/ctdna_data"
output.dir <- "/home/anshul/software/ctdna/ctdna_data/ctdnaoutput"
reference.genome <- "/home/anshul/software/ctdna/ctdna_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
snpeff.database <- "GRCh38.99"  # Replace with your actual snpEff database
gvcf.dir <- file.path(output.dir, "gvcfs")
joint.vcf <- file.path(output.dir, "joint_genotyped.vcf.gz")

# Create necessary directories
dir.create(file.path(output.dir, "fastqc_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(gvcf.dir, recursive = TRUE, showWarnings = FALSE)
system(paste("sudo chown anshul", shQuote(output.dir)))
system(paste("sudo chmod 755", shQuote(output.dir)))

# List all FASTQ files in the input directory
fastq.files <- list.files(input.dir, pattern = "\\.fastq$", full.names = TRUE)

# Function to execute commands with error handling
execute <- function(cmd, quitOnError = TRUE) {
  cat("Running command:", cmd, "\n")
  res <- system(cmd, intern = TRUE)
  if (res[1] >= 1 && quitOnError) {
    stop("Command failed: ", cmd)
  }
}

# Process each FASTQ file
for (fastq.file in fastq.files) {
  sample.name <- gsub("\\.fastq$", "", basename(fastq.file))

  # Define paths for outputs
  trimmed.fastq <- file.path(output.dir, paste0("trimmed_", sample.name, ".fastq"))
  sam.file <- file.path(output.dir, paste0(sample.name, "_aligned.sam"))
  sorted.bam <- file.path(output.dir, paste0(sample.name, "_sorted.bam"))
  dedup.bam <- file.path(output.dir, paste0(sample.name, "_deduplicated.bam"))
  gvcf.file <- file.path(gvcf.dir, paste0(sample.name, ".g.vcf.gz"))

  # Step 1: FastQC
  execute(paste("fastqc", fastq.file, "-o", file.path(output.dir, "fastqc_results")))

  # Step 2: Trim adapters using cutadapt
  execute(paste("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o", shQuote(trimmed.fastq), shQuote(fastq.file)))

  # Step 3: Align reads with BWA
  execute(paste("bwa mem -M -t 4", shQuote(reference.genome), shQuote(trimmed.fastq), ">", shQuote(sam.file)))

  # Step 4: Convert SAM to BAM and sort
  execute(paste("samtools view -Sb", shQuote(sam.file), "| samtools sort -o", shQuote(sorted.bam)))

  # Step 5: Index the sorted BAM
  execute(paste("samtools index", shQuote(sorted.bam)))

  # Step 6: Deduplicate BAM using GATK MarkDuplicates
  dedup.metrics <- file.path(output.dir, paste0("dedup_metrics_", sample.name, ".txt"))
  execute(paste(
    "gatk MarkDuplicates",
    "-I", shQuote(sorted.bam),
    "-O", shQuote(dedup.bam),
    "-M", shQuote(dedup.metrics),
    "--REMOVE_DUPLICATES true"
  ))

  # Step 7: Add read groups
  rg.bam <- file.path(output.dir, paste0(sample.name, "_with_RG.bam"))
  execute(paste(
    "samtools addreplacerg -r \"@RG\\tID:", sample.name, "\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1\\tSM:", sample.name, "\" -o",
    shQuote(rg.bam), shQuote(dedup.bam)
  ))

  # Step 8: Index BAM with read groups
  execute(paste("samtools index", shQuote(rg.bam)))

  # Step 9: Variant calling (GVCF mode) using GATK HaplotypeCaller
  execute(paste("gatk HaplotypeCaller -R", shQuote(reference.genome), "-I", shQuote(rg.bam), "-O", shQuote(gvcf.file), "-ERC GVCF"))

  # Index the GVCF file
  execute(paste("gatk IndexFeatureFile -I", shQuote(gvcf.file)))
}

# Combine GVCFs and perform joint genotyping
combined.gvcf <- file.path(output.dir, "combined.g.vcf.gz")
execute(paste(
  "gatk CombineGVCFs -R", shQuote(reference.genome),
  paste(sapply(list.files(gvcf.dir, pattern = "\\.g.vcf.gz$", full.names = TRUE), function(f) paste("-V", shQuote(f))), collapse = " "),
  "-O", shQuote(combined.gvcf)
))
execute(paste("gatk GenotypeGVCFs -R", shQuote(reference.genome), "-V", shQuote(combined.gvcf), "-O", shQuote(joint.vcf)))

# Filter variants
filtered.vcf <- file.path(output.dir, "joint_genotyped_filtered.vcf.gz")
execute(paste(
  "gatk VariantFiltration -R", shQuote(reference.genome), "-V", shQuote(joint.vcf),
  "-O", shQuote(filtered.vcf), "--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0' --filter-name BasicFilter"
))

# Check if any variants passed the filters
filtered.count <- as.numeric(system(paste("zgrep -v '^#'", shQuote(filtered.vcf), "| wc -l"), intern = TRUE))
if (filtered.count == 0) {
  cat("No variants passed the filters in the joint genotyped VCF. Skipping annotation.\n")
} else {
  # Annotate filtered variants
  annotated.vcf <- file.path(output.dir, "joint_genotyped_annotated.vcf")
  execute(paste("snpEff ann -v -o vcf", shQuote(snpeff.database), shQuote(filtered.vcf), ">", shQuote(annotated.vcf)))
}

cat("Pipeline Completed.\n")
