# Bulk_RNA_Seq_data_analysis with multiple contrasts

# STEP 1: INDEXING THE REFERENCE GENOME;  --sjdbOverhang readLength-1
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/reference/genome.fasta --sjdbGTFfile /path/to/reference/genome.gtf --sjdbOverhang 149
```

# STEP 2: QUALITY ASSESSMENT, QUALITY CONTROL, ALIGNMENT AGAINST THE REFERENCE GENOME, SORTING THE BAM FILES
# Create a bash script for running the pipeline
```
vim rnaseq_pipeline.sh
```

#Inside rnaseq_pipeline.sh:

```
#!/bin/bash
echo "Analysis started"
for R1 in *_1.fq.gz; 
do
    # Get the corresponding R2 file
    R2="${R1/_1.fq.gz/_2.fq.gz}"
    samplename=$(basename "$R1" _1.fq.gz)  
    echo "Sample in analysis: $samplename"
    echo "Files: $R1, $R2"

    # Quality control with FastQC
    echo "Running FastQC for $R1"
    fastqc "$R1" --outdir ./
    echo "Running FastQC for $R2"
    fastqc "$R2" --outdir ./

    # Read trimming with Trimmomatic
    echo "Trimming reads with Trimmomatic for $samplename..."
    trimmomatic PE -phred33 "$R1" "$R2" \
        "${samplename}_1_paired.fq.gz" "${samplename}_1_unpaired.fq.gz" \
        "${samplename}_2_paired.fq.gz" "${samplename}_2_unpaired.fq.gz" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo "    "

    # Update variables to use trimmed files
    trimmed_R1="${samplename}_1_paired.fq.gz"
    trimmed_R2="${samplename}_2_paired.fq.gz"

    echo "Alignment with STAR"
    STAR --runMode alignReads --runThreadN 10 --genomeDir /path/to/STAR/genome/folder \
         --sjdbGTFfile /path/to/genomeAnnotation.gtf --sjdbOverhang 149 \
         --readFilesIn "$trimmed_R1" "$trimmed_R2" --outFileNamePrefix "${samplename}_" \
         --outSAMtype BAM Unsorted --readFilesCommand gunzip -c

    # Sorting the BAM file for further analysis or indexing
    echo "Sorting BAM file"
    samtools sort -o "${samplename}_sorted.bam" "${samplename}_Aligned.out.bam"

    echo "Samtools indexing"
    samtools index "${samplename}_sorted.bam"

    echo "Generating bedgraph with STAR"
    STAR --runMode inputAlignmentsFromBAM --runThreadN 10 --inputBAMfile "${samplename}_sorted.bam" \
         --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix "${samplename}_"
    echo "    "
done
echo "Analysis finished"

# Grant execution permissions to the script
chmod +x rnaseq.sh
```

#save the file
```
esc
:wq
```
#Check this script on small datasets

Once the alignemnt is complete get the count matrix using featureCounts. -s 0,1,2 non-stranded data, stranded, reverse stranded, respectively
```
featureCounts -t exon -g gene_id -p -s 0 -a /path/to/genome.gtf -T 12 -o /path/to/output_featureCounts.txt /path/to/bam_files/*_Aligned.sortedByCoord.out.bam
```

Once the count matrix is generated we can proceed for differential gene expression analysis, GSEA, etc and the corresponding R script can be found iun differential_gene_expression.R 
