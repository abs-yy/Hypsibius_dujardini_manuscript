# Keio side commands
```
GENOME=nHd3.0.fsa
```
1. Genome validation
```{genome_validation.sh}
# DNA-Seq mapping
DNASEQ=*** /ex DRR055040
R1=$DNASEQ.R1.fq
R2=$DNASEQ.R2.fq

bwa index $GENOME
bwa mem $GENOME $R1 $R2 -t 30  > $DNASEQ.mem.sam
samtools view -@ 30 -bS $DNASEQ.mem.sam > $DNASEQ.mem.bam
samtools sort -@ 30 $DNASEQ.mem.bam $DNASEQ.mem.sorted
samtools index $DNASEQ.mem.sorted.bam
bin/qualimap_v2.2/qualimap bamqc -bam $DNASEQ.mem.sorted.bam -outformat pdf --java-mem-size=16G

# CEGMA
cegma -g $GENOME -T 30

# Pilon polishing
bwa index $GENOME
java -Xmx200g -jar ~/bin/pilon-1.17.jar --genome $GENOME --frags $DNASEQ.mem.sorted --output $GENOME --changes --vcf

```

2. Genome annotation
```{genome_annot.sh}
# ribosomal RNA prediction
perl /home/yuki.yoshida/bin/rnammer-1.2/rnammer -multi -S euk -m lsu,ssu,tsu -gff $GENOME.rnammer.gff < $GENOME

# Repeat detection
build_lmer_table -sequence $GENOME -freq lmer_table
RepeatScout -sequence $GENOME -output rs_output.fa -freq lmer_table
~/bin/maker/exe/RepeatMasker/RepeatMasker -parallel 30 -dir . -lib rs_output.fa $GENOME

# transfer RNA prediction
## Thank you sujai.  (https://github.com/sujaikumar/assemblage/blob/master/README-annotation.md)
~/bin/tRNAscan-SE-1.3.1/tRNAscan-SE $GENOME -o $GENOME.trna.raw
perl < $GENOME.trna.raw -lne '
    ($contig, $trnanum, $st, $en, $type, $codon, $intron_st, $intron_en, $score)=split /\s+/;
    next unless $trnanum =~ /^\d+$/; # get rid of trnascan headers
    print "$contig\ttRNAscan-SE-1.3\ttRNA\t" . ($st < $en ? "$st\t$en\t$score\t+\t" : "$en\t$st\t$score\t-\t") . ".\tID=$contig\_tRNA$trnanum;Name=tRNA_$type\_$codon"
' > hypsibius-georgios-extended-filtered.fasta.trna.gff
```


3. BRAKER
```{braker.sh}
SAMPLE=***
R1=${SAMPLE}.R1.fq
R2=${SAMPLE}.R2.fq
OUT=$SAMPLE

bowtie2-build ${GENOME} $GENOME
# Paired end
tophat -p 30 -o ${OUT} ${GENOME} ${R1} ${R2}
# Single read
tophat -p 30 -o ${OUT} ${GENOME} ${R1}

perl ~/bin/BRAKER_v1.9/braker.pl \
  --genome=$GENOME \
    --cores 30 --species Hd_braker_all --overwrite \
    --bam A1_L1/accepted_hits.bam \
    --bam A1_L2/accepted_hits.bam \
    --bam A2_L1/accepted_hits.bam \
    --bam A2_L2/accepted_hits.bam \
    --bam A3_L1/accepted_hits.bam \
    --bam A3_L2/accepted_hits.bam \
    --bam Tun_1_L1/accepted_hits.bam \
    --bam Tun_1_L2/accepted_hits.bam \
    --bam Tun_2_L1/accepted_hits.bam \
    --bam Tun_2_L2/accepted_hits.bam \
    --bam Tun_3_L1/accepted_hits.bam \
    --bam Tun_3_L2/accepted_hits.bam  \
    --bam H-B1-1/accepted_hits.bam \
    --bam H-B1-2/accepted_hits.bam \
    --bam H-B1-3/accepted_hits.bam \
    --bam H-B2-1/accepted_hits.bam \
    --bam H-B2-2/accepted_hits.bam \
    --bam H-B2-3/accepted_hits.bam \
    --bam H-B3-1/accepted_hits.bam \
    --bam H-B3-2/accepted_hits.bam \
    --bam H-B3-3/accepted_hits.bam \
    --bam H-B4-1/accepted_hits.bam \
    --bam H-B4-2/accepted_hits.bam \
    --bam H-B4-3/accepted_hits.bam \
    --bam H-B5-1/accepted_hits.bam \
    --bam H-B5-2/accepted_hits.bam \
    --bam H-B5-3/accepted_hits.bam \
    --bam H-B6-1/accepted_hits.bam \
    --bam H-B6-2/accepted_hits.bam \
    --bam H-B6-3/accepted_hits.bam \
    --bam H-B7-1/accepted_hits.bam \
    --bam H-B7-2/accepted_hits.bam \
    --bam H-B7-3/accepted_hits.bam \
    --bam H-E1-1/accepted_hits.bam \
    --bam H-E1-2/accepted_hits.bam \
    --bam H-E1-3/accepted_hits.bam \
    --bam H-E2-1/accepted_hits.bam \
    --bam H-E2-2/accepted_hits.bam \
    --bam H-E2-3/accepted_hits.bam \
    --bam H-E3-1/accepted_hits.bam \
    --bam H-E3-2/accepted_hits.bam \
    --bam H-E3-3/accepted_hits.bam \
    --bam H-E4-1/accepted_hits.bam \
    --bam H-E4-2/accepted_hits.bam \
    --bam H-E4-3/accepted_hits.bam \
    --bam H-E5-1/accepted_hits.bam \
    --bam H-E5-2/accepted_hits.bam \
    --bam H-E5-3/accepted_hits.bam

cd Hd_braker_all
getAnnoFasta.pl --seqfile $GENOME augustus.gff

```

3. Predicted CDS sequence and amino acid sequence annotation and validation
```{cds_annot.sh}
CDSSEQ=augustus.codingseq
AASEQ=augustus.aa

# RNA-Seq mapping
SAMPLE=***
R1=$SAMPLE.R1.fq
R2=$SAMPLE.R2.fq

## Raw mapping
bwa mem $GENOME $R1 $R2 -t 30 > $SAMPLE.mem.sam
samtools flagstat $SAMPLE.mem.sam

## Trinity assembly

bin/trinityrnaseq-Trinity-v2.3.2/Trinity --seqType fq --max_memory 100G --left $SAMPLE.R1.fq --right $SAMPLE.R2.fq --CPU 15 --output trinity_$SAMPLE
formatdb -i trinity_$SAMPLE/Trinity.fasta -p F
blastall -p blastn -i $CDSSEQ  -d trinity_$SAMPLE/Trinity.fasta -m 8 -a 30 -e 1e-50 -o trinity_$SAMPLE/Trinity.fasta.blastn.cdsseq.1e-50
blastall -p blastn -i $CDSSEQ  -d $GENOME -a -m 8 -a 30 -e 1e-50 -o trinity_$SAMPLE/Trinity.fasta.blastn.genome.1e-50
cut -f 1 trinity_$SAMPLE/Trinity.fasta.blastn.cdsseq.1e-50 | sort | uniq | wc -l
cut -f 1 trinity_$SAMPLE/Trinity.fasta.blastn.genome.1e-50 | sort | uniq | wc -l


# Database search
blastall -p blastp -i $AASEQ -d /path/to/db/uniprot_sprot.fa -m 8 -a 30 -e 1e-15 -o $AASEQ.blastp.swissprot.1e-15
blastall -p blastp -i $AASEQ -d /path/to/db/uniprot_trembl.fa -m 8 -a 30 -e 1e-15 -o $AASEQ.blastp.trembl.1e-15
hmmsearch --cpu 30 --domE 1e-3 --domtblout $AASEQ.hmmsearch.Pfam.1e-3 /path/to/Pfam/Pfam-A.hmm $AASEQ

hmmsearch --cpu 30 --domE 1e-3 --domtblout $CDSSEQ.hmmsearch.Dfam.1e-3 /path/to/Dfam/Dfam.hmm $CDSSEQ


```


4. miRNA
```{miRNA_analysis.sh}

miRDeep2.pl reads_collapsed.fa ~/final/hypsibius-georgios-extended-filtered.fasta reads_collapsed_vs_genome.arf none /path/to/miRdeep/mature.fa none
```

5. HGT index based HGT candidate identification
```{HGT.sh}
wget -R ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/




```





