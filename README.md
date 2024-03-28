This is a compilation of the scripts used in this publication:

**Low pass Nanopore sequencing for measurement of global methylation levels in plants**

Some of these scripts were adapted from code provided in this paper: Faulk, C. (2023). Genome skimming with nanopore sequencing precisely determines global and transposon DNA methylation in vertebrates. Genome Research, 33(6), 948–956. and its supplementary document at: https://faulk-lab.github.io/skimming/.



**Data pre-process (pod5)**

If output from ONT sequencer are in fast5, run pod5 convert to convert fast5 into pod5 to enable re-basecalling using dorado
    
    pod5 convert path_to_fast5_dir/ --output file.pod5


**Re-basecall (dorado)**

The MinKNOW interface on the ONT sequencer did not have modbase model for non-CG contexts available, therefore raw outputs from instrument were re-basecalled after sequencing using dorado (v.0.3.2) to detect both all context 5mC and 6mA methylation.
    
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
    samtools fastq -T '*' file.bam | chopper -q 10 | gzip > \
    file.fastq.gz
    
    #Map fastq to reference, sort and index
    minimap2 -ax map-ont -y --secondary=no -t 16 \ 
    ref.fasta \ 
    file.fastq.gz | samtools sort -@8 \
    -o file.mapped.sorted.bam
    samtools index file.mapped.sorted.bam
    
   
**Calculate coverage (mosdepth)**
    
    cd outdir
    mosdepth -t 8 --fast-mode file_name file.mapped.sorted.bam


**Downsample (samtools)**

Based on coverage as calculated by mosdepth, determine downsampling proportion needed to achieve 10x, 1x, 0.1x, 0.01x and 0.001x coverage, each replicated with 10x bootstraps.
For example, for the Vitis vinifera sample, the coverage was calculated at 168x, so to achieve 10x coverage, downsample with proportion of 0.059x.
    
    for i in {01..10}; do \
    samtools view -@ 16 -b -s$RANDOM.059 file.mapped.sorted.bam > file.10x.${i}.bam;
    samtools sort file.10x.${i}.bam -o file.10x.${i}.sorted.bam;
    samtools index -@ 16 $DIR/bam/10x/Vv.10x.${i}.sorted.bam; done
    

**Call methylation and split into each context (modkit + bedtools)**

Call the methylation signals in the bam files, as a methylation level per site in a bedmethyl file using modkit pileup. Modkit motif was used to create a motif reference bed file containing the chromosomal position of each CG, CHG, CHH, and 6mA context. Bedtools intersect was then use to split the bedmethyl files into each context.

    #Run modkit motif
    #example for CG context
    modkit motif-bed ref.fasta CG 0 > ref.CG.bed
    
    #Run modkit pileup for each replicate of each downsample
    for i in {01..10}; do \
    for j in 10x 1x 0.1x 0.01x 0.001x; do \
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


**Calculate global methylation (awk)**
     
     for k in CG CHG CHH 6mA; do \
     for j in 0.001x 0.01x 0.1x 1x; do \
     for i in path/*.${j}.${k}.bed; do \
     awk '{can+=$13; mod+=$12} END{print 100*(mod/(can+mod))}' $i >> \
     globalmethylation.${j}.${k}.txt; \
     done; done; done


**Group reads of different length (chopper or bbmap)**

For analysis of different read length, reads from Vitis vinifera data were grouped into different read length, <10kb, 10-50 kb, and >50kb, using chopper
     
     #example for reads 10-50kb
     gunzip -c file.fastq.gz| chopper --minlength 10000 --maxlength 50000 | \
     gzip > file_10_50kb.fastq.gz
     minimap2 -ax map-ont -y --secondary=no -t 16 file_10_50kb.fastq.gz | \
     samtools sort -@8 -o file_10_50kb_sorted.bam
     samtools index file_10_50kb_sorted.bam

As this Vitis database contains mainly long reads, generation of reads <= 5kb length was done by trimming the reads to 5kb  using reformat.sh from bbmap.
```
reformat.sh in=file.mapped.sorted.bam \ 
out=file.5kb.mapped.sorted.bam \
allowidenticalnames=t \
forcetrimright=5000 \ 
ref=ref.fasta \
usejni=t 
```


 
 **Phase reads (whatshap)**

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


**Calculate methylation entrophy (DMEAS)**

For analysis of DNA methylation heterogeneity in CG context across the genome  among different plant species assessed in this study (*Vitis vinifera, Arabidopsis thaliana*, and *Actinia melanandra*) and human data

    #Example showed here for 
    #Extract methylation data per read using modkit extract
    
    modkit extract \ 
    arabidopsis_DB00183.aligned.bam \
    arabidopsis.extract.tsv \
    --read-calls arabidopsis.readcalls.tsv \
    --log-filepath arabidopsis_extract.log \
    --mapped-only \
    --reference GCA_000001735.2.fasta \
    --cpg
    
    #Convert methylation per read into bis format according to DMEAS requirement (https://sourceforge.net/projects/dmeas/files/). 
    #The format is a file with 4 columns:  reads IDs, methylation state, chromosome ID, genome start position and methylation call.
    #In methylation state: “+” means methylated, “-” means unmethylated, 
    #in methylation call: “Z” means methylated in CpG dinucleotide,“z” means unmethylated in CpG dinucleotide
    
    awk '{print $1,$12,$4,$3,$12}' \ 
    arabidopsis.readcalls.tsv | \
    awk '{ if ($2 == "m") $2="+"; print $0 }' | \
    awk '{ if ($5 == "m") $5="Z"; print $0 }' | \
    awk '{ if ($5 == "-") $5="z"; print $0 }' > \
    arabidopsis.readcalls.bis
    
    #DMEAS expect bis file per chromosome in a directory, as well as  gp file per chromosome which is position of the CG in the genome.
    
    #Split bis file based on chromosome (column 3 in the bis file)
    awk -F\  '{print>"/path/dir_tocontain_bis_files/"$3}' arabidopsis.readcalls.bis    
    
    #use CG_motif.bed file (generated using modkit motif-bed) to create gp file of each chromosome
    awk -F\  '{print $2>"/path/dir_tocontain_gp/"$1}' arabidopsis_CG.bed
    
    #Rename the bis and gp file to contain the .bis and .gp extension:
    cd /path/dir_tocontain_bis_files/
    for file in *; do \
    mv -- "$file" "${file%}.bis"; done
    
    cd /path/dir_tocontain_gp/
    for file in *; do \
    mv -- "$file" "${file%}.gp"; done
    
    #run DMEAS
    #DMEAS.pl was downloaded from https://sourceforge.net/projects/dmeas/files/
    
    perl DMEAS.pl -d gw -s ms \
    /path/dir_tocontain_bis_files/ \
    /path/dir_tocontain_gp/ \
    -o /path_to_outdir/
    
    #Calculate the average entrophy across each chromosome
    for i in {1..5}; do \
    cd /path_to_outdir/; \
    for j in /path_to_outdir/*4CG.txt; do; \
    awk '{ sum += $4 } END { if (NR > 0) print sum \
    / (NR-2) }' $j >> entropy.arabidopsis.txt'; \
    done; done
    
    
**Annotate bedmethyl files with TE or gene information**






