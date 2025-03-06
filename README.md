# NGS Variant Calling Pipeline

## Overview

This pipeline performs **NGS variant calling** using Next-Generation Sequencing (NGS) data. It follows a **GATK-based workflow** for quality control, read alignment, deduplication, variant calling, filtering, and annotation.

### Features:

- **Quality control** using FastQC
- **Adapter trimming** using Cutadapt
- **Read alignment** with BWA MEM
- **Sorting & Deduplication** using Samtools & GATK MarkDuplicates
- **Variant Calling** with GATK HaplotypeCaller
- **Joint Genotyping** with GATK GenotypeGVCFs
- **Variant Filtering** using GATK VariantFiltration
- **Variant Annotation** with snpEff

---

## Dependencies

Ensure the following software/tools are installed before running the pipeline:

| Tool         | Version | Installation Guide                                    |
| ------------ | ------- | ----------------------------------------------------- |
| **R**        | 4.3+    | [R Project](https://cran.r-project.org/)              |
| **FastQC**   | 0.11+   | `sudo apt install fastqc`                             |
| **Cutadapt** | 4.0+    | `pip install cutadapt`                                |
| **BWA**      | 0.7+    | `sudo apt install bwa`                                |
| **Samtools** | 1.15+   | `sudo apt install samtools`                           |
| **GATK**     | 4.2+    | [GATK Installation](https://gatk.broadinstitute.org/) |
| **snpEff**   | Latest  | [snpEff Guide](https://pcingola.github.io/SnpEff/)    |

---

## Usage

### 1.Clone the repository:

```sh
 git clone https://github.com/yourusername/NGS-Variant-Calling.git
 cd NGS-Variant-Calling
```

### 2.Modify input parameters:

- Update the **input directory, reference genome, and output path** inside `pipeline.R`.

### 3.Run the pipeline in R:

```r
 source("pipeline.R")
```

---

## Pipeline Workflow

###  **Step 1: Quality Control**

- Uses **FastQC** to check raw sequencing reads.

###  **Step 2: Adapter Trimming**

- Uses **Cutadapt** to remove sequencing adapters.

### **Step 3: Read Alignment**

- Reads are aligned to the **human genome (GRCh38)** using **BWA MEM**.

###  **Step 4: Sorting & Deduplication**

- Sorted using **Samtools**.
- PCR duplicates are removed using **GATK MarkDuplicates**.

### **Step 5: Variant Calling**

- **GATK HaplotypeCaller** is used in GVCF mode.

###  **Step 6: Joint Genotyping**

- Uses **GATK CombineGVCFs** and **GenotypeGVCFs** to call variants across multiple samples.

### **Step 7: Variant Filtering & Annotation**

- Filters variants with **GATK VariantFiltration**.
- Annotates variants using **snpEff**.

---

## Expected Output

- **`fastqc_results/`** → Quality check reports
- **`*.bam`** → Aligned and deduplicated BAM files
- **`*.vcf.gz`** → Variant Call Format files
- **`joint_genotyped_filtered.vcf.gz`** → Final filtered variants

---

## Citation

If you use this pipeline, please cite:

- **GATK Best Practices** [(Link)](https://gatk.broadinstitute.org/)
- **snpEff Documentation** [(Link)](https://pcingola.github.io/SnpEff/)

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Contact

**Anshul Kamboj**\
📧 [kambojanshul51@gmail.com](mailto\:kambojanshul51@gmail.com)

---


