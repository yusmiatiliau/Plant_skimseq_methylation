This is a compilation of the scripts used in this publication:

**Low pass Nanopore sequencing for measurement of global methylation levels in plants**

Some of these scripts were adapte based on this paper: Faulk, C. (2023). Genome skimming with nanopore sequencing precisely determines global and transposon DNA methylation in vertebrates. Genome Research, 33(6), 948–956. and its supplementary document at: https://faulk-lab.github.io/skimming/.



**Data pre-process (pod5)**

If output from ONT sequencer are in fast5, run pod5 convert to convert fast5 into pod5 to enable re-basecalling using dorado
    
    pod5 convert path_to_fast5_dir/ --output file.pod5


**Re-basecall (dorado)**

Modbasecalling on ONT sequencer does not have modbase model for non-CG contexts, therefore, raw outputs from instrument were re-basecalled using dorado (v.0.3.2) to detect both all context 5mC and 6mA methylation.
    
    dorado basecaller --device cuda:all --emit-moves \
    path/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
    path_pod5_dir/ \
    --modified-bases-models path/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mC@v2, \
    path/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2 > \
    file.bam


**QC (nanoplot)**

Run nanoplot using the sequencing summary file to get some QC metrics, such as N50, estimated number of bases, etc.
    
    NanoPlot --outdir path/dir/ \
    --threads 8 \ --tsv_stats \
    --info_in_report --N50 \
    --summary  sequencing_summary.txt \
    --barcoded \
    --no_static
    
    
**Filter and map (chopper, samtools, minimap2)**

Use samtools 
    
    #Convert bam into fastq, carrying all tags including MM,ML tags for methylation signal,
    #And filter based on reads quality >=10
    samtools bam2fq -T '*' file.bam | chopper -q 10 | gzip > \
    file.fastq.gz
    
    #Map fastq to reference, sort and index
    minimap2 -ax map-ont -y --secondary=no -t 16 \ 
    ref.fasta \ 
    file.fastq.gz | samtools sort -@8 \
    -o file.mapped.sorted.bam
    samtools index file.mapped.sorted.bam
    
   
**Calculate coverage (mosdepth**
    
    cd outdir
    mosdepth -t 8 --fast-mode file_name file.mapped.sorted.bam


**Downsample (samtools)**

Based on coverage as calculated by mosdepth, determine downsampling proportion needed to achieve 10x, 1x, 0.1x, 0.01x and 0.001x coverage, each replicated with 10x bootstraps.
For example, for the Vitis vinifera sample, the coverage was calculated at 168x, so to achieve 10x coverage, downsample with proportion of 0.059x.
    
    for i in {01..10}; do; \
    samtools view -@ 16 -b -s$RANDOM.059 file.mapped.sorted.bam > file.10x.${i}.bam;
    samtools sort file.10x.${i}.bam -o file.10x.${i}.sorted.bam;
    samtools index -@ 16 $DIR/bam/10x/1031.10x.${i}.sorted.bam; done
    

**Call methylation and split into each context (modkit + bedtools)**

Call the methylation signal written in bam file, as a methylation level per site in a bedmethyl file using modkit pileup. Modkit motif used to create a motif reference bed file containing the chromosomal position of each CG, CHG, CHH, and 6mA context. Bedtools intersect was then use to split the bedmethyl files into each context.

    #Run modkit motif
    #example for CG context
    modkit motif-bed ref.fasta CG 0 > ref.CG.bed
    
    #Run modkit pileup for each replicate of each downsample
    for i in {01..10}; do; \
    for j in 10x 1x 0.1x 0.01x 0.001x; do; \
    modkit pileup file.${j}.${i}.sorted.bam \ 
    file.${j}.${i}.bed \
    --log-filepath file.${j}.${i}.log \
    --filter-threshold C:0.85 --mod-thresholds m:0.85 --only-tabs; \
    done; done
    
    #Run bedtools intersect to get context-specific methylation data
    #Example for CG context
    #There can be "a" methylation recorded due to SNPs,
    #include only "m" methylation signal using awk
    for j in 1x 0.1x 0.01x; do; \
    for i in {01..10}; do; \
    bedtools intersect -a file.${j}.${i}.bed -b ref.CG.bed -sorted | \
    awk '$4=="m"' > file.${j}.${i}.CG.bed; \
    done; done
    
Alternatively, later versions of modkit pileup enable adding information of methylation context into bedmethyl file using the flag --motif, for example adding --motif CG 0 --motif CHG 0 --motif CHH 0 --motif A 0, to modkit pileup command, to get all CG, CHG, CHH, and 6mA methylation information.


**Global methylation calculation (awk)**
     
     for k in CG CHG CHH 6mA; do; \
     for j in 0.001x 0.01x 0.1x 1x; do; \
     for i in path/*.${j}.${k}.bed; do; \
     awk '{can+=$13; mod+=$12} END{print 100*(mod/(can+mod))}' $i >> \
     globalmethylation.${j}.${k}.txt; \
     done; done; done


**Grouping of different read length (chopper)**

For analysis of different read length, reads from Vitis vinifera data were grouped into different read length, <10kb, 10-50 kb, and >50kb, using chopper
     
     #example for reads 10-50kb
     gunzip -c file.fastq.gz| chopper --minlength 10000 --maxlength 50000 | \
     gzip > file_10_50kb.fastq.gz
     minimap2 -ax map-ont -y --secondary=no -t 16 file_10_50kb.fastq.gz | \
     samtools sort -@8 -o file_10_50kb_sorted.bam
     samtools index file_10_50kb_sorted.bam

     


**trimming reads (shred.sh from bbmap) followed by modkit repair**

    #split reads into 5kb sizes
    shred.sh in=original.fastq.gz \
    out=5kb.fastq.gz length=5000 minlength=1 overlap=0 overwrite=true usejni=t
    
    #map and sort the new 5kb reads
    minimap2 -ax map-ont -y --secondary=no -t 8 ref.fasta \
     5kb.fastq.gz | samtools sort -@8 -o 5kb.sorted.bam 
    samtools index 5kb.sorted.bam

    #sort the new 5kb reads and the orignal coverage bam file by read names
    samtools sort -n original.bam -O BAM > original_sortedname.bam -@8
    samtools sort -n 5kb.sorted.bam -O BAM > 5kb_sortedname.bam -@8
    
    #run modkit repair to repair the MM,ML tags based on the original bam file
    modkit repair --donor-bam donor_sortedname.bam \
    --acceptor-bam 5kb_sortedname.bam \
    --log-filepath /repair.log \
    --output-bam 5kb_repaired.bam
    
  
**phasing (whatshap)**

    #run Clair3 to get phased variants information
    run_clair3.sh \
    --bam_fn=file.sorted.bam \
    --ref_fn=ref.fasta \
    --threads=16 \
    --platform="ont" \
    --model_path="path/r1041_e82_400bps_sup_v420" \
    --output=/path/dir \
    --include_all_ctgs \
    --use_longphase_for_intermediate_phasing \
    --use_whatshap_for_final_output_phasing \
    --use_whatshap_for_final_output_haplotagging

    
    #create haplotype information in tsv
    whatshap haplotag \
    --output file_haplotagged.bam \
    --reference ref.fasta \
    --output-haplotag-list /file_haplotype.tsv \
    --output-threads 32 \ 
    phased.vcf.gz \
    file.sorted.bam \
    --ignore-read-groups
    
    #split bam file into each haplotype
    whatshap split \
    --output-h1 file_h1.bam \
    --output-h2 file_h2.bam \ 
    file.sorted.bam \ 
    file_haplotype.tsv


**methylation entrophy (DMEAS)**

For analysis of DNA methylation heterogeneity in CG context across the genome  among different plant species assessed in this study (*Vitis vinifera, Arabidopsis thaliana, and Actinia melanandra*) and human data

    #Example showed here for 
    #Extract methylation data per read using modkit extract

plotting (ggplot)






