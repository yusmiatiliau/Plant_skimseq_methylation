This is a compilation of the scripts used in this publication:

**Low pass Nanopore sequencing for measurement of global methylation levels in plants**


data pre-processing (pod5)

basecalling (dorado)

QC (nanoplot)

methylation entrophy (DMEAS)

mapping (samtools + minimap2)

coverage (mosdepth)

downsampling (samtools)

methylation calling (modkit + bedtools)

calculation (awk)

grouping of different read length (chopper)

trimming reads (shred.sh from bbmap followed by modkit repair)

    shred.sh in=original.fastq.gz \
    out=5kb.fastq.gz length=5000 minlength=1 overlap=0 overwrite=true usejni=t
    
    minimap2 -ax map-ont -y --secondary=no -t 8 ref.fasta \
     5kb.fastq.gz | samtools sort -@8 -o 5kb.sorted.bam 
    samtools index 5kb.sorted.bam

    samtools sort -n original.bam -O BAM > original_sortedname.bam -@8
    samtools sort -n 5kb.sorted.bam -O BAM > 5kb_sortedname.bam -@8
    
    modkit repair --donor-bam donor_sortedname.bam \
    --acceptor-bam 5kb_sortedname.bam \
    --log-filepath /repair.log \
    --output-bam 5kb_repaired.bam
    

    
    
    

phasing (whatshap)

plotting (ggplot)






