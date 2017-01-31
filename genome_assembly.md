1. FALCON assembly of PacBio reads
* Error correction:
daligner cutoff of 4,000 bp and -k16 -e0.70 -s1000 -t16 -l1000 -h64 -w7
* Second daligner stage:
-k20 -e.96 -s1000 -t32 -l1500 -h256,
* Final assembly:
min coverage of 2, max coverage of 80, max diff coverage of 40

1. I assembled your MiSeq data with SPAdes (v 3.8.1) after merging them with usearch. (I believe that SPAdes is better than MaSuR-CA for the given data)

```{spades.sh}
spades.py \
--only-assembler \
--careful \
--s1 notmerged.fq \
--s2 merged.fq \
-k 21,33,55,77,99,127 \
-m 800 \
-t 16 \
-o spades_output
``` 
2. Mapped the UNC reads with bowtie2 (v 2.2.4) to the SPAdes assembly having filtered contigs below 250. I retained reads that at least one read mapped to the assembly.

3. Assembled the UNC reads with Platanus (v 1.2.4). I have used -u 0.2 because I didn't want to overcollapse bubbles. Similarly, I think that redundans over collapses based on prior experiences.

```{platanus.sh}
platanus assemble \
-o platanus \
-f \
Lib300_unc_1.fq \
Lib300_unc_2.fq \
Lib500_unc_1.fq \
Lib500_unc_2.fq \
Lib800_unc_1.fq \
Lib800_unc_2.fq \
-t 64 \
-u 0.2 \
-m 400

platanus scaffold \
-o platanus \
-c platanus_contig.fa \
-b platanus_contigBubble.fa \
-u 0.2 \
-IP1 Lib300_unc_1.fq Lib300_unc_2.fq \
-IP2 Lib500_unc_1.fq Lib500_unc_2.fq \
-IP3 Lib800_unc_1.fq Lib800_unc_2.fq \
-t 64

platanus gap_close \
-o platanus \
-c platanus_scaffold.fa \
-IP1 Lib300_unc_1.fq Lib300_unc_2.fq \
-IP2 Lib500_unc_1.fq Lib500_unc_2.fq \
-IP3 Lib800_unc_1.fq Lib800_unc_2.fq \
-t 64
```

4. Run SSPACE-LongRead and PBJelly on the gapclosed assembly.
I have only changed the -m (Minimum MapQ) parameter because I believe that 200 is too high.

```{sspace_pbjelly.sh}
FILEDIR=./pbjelly
Jelly.py setup $FILEDIR/Protocol_hd.xml
Jelly.py mapping $FILEDIR/Protocol_hd.xml
Jelly.py support $FILEDIR/Protocol_hd.xml -x "-m 50"
Jelly.py extraction $FILEDIR/Protocol_hd.xml
Jelly.py assembly $FILEDIR/Protocol_hd.xml -x "-n 64"
Jelly.py output $FILEDIR/Protocol_hd.xml

# Protocol_hd.xml
<jellyProtocol>
    <reference>platanus_gapClosed.fasta</reference>  
    <outputDir>pbjelly/</outputDir>
    <blasr>-minMatch 8 -minPctIdentity 70 -bestn 1 -nCandidates 20 -maxScore -500 -nproc 64 -noSplitSubreads</blasr>
    <input baseDir="arakawa_pacbio_raw">
        <job>PacBio-raw.fastq</job>
    </input>
</jellyProtocol>
```


5. Filter low quality contigs
```
DNASEQ=DRR055040
R1=$DNASEQ.R1.fq
R2=$DNASEQ.R2.fq

bwa index $GENOME
bwa mem $GENOME $R1 $R2 -t 30  > $DNASEQ.mem.sam
samtools view -@ 30 -bS $DNASEQ.mem.sam > $DNASEQ.mem.bam
samtools sort -@ 30 $DNASEQ.mem.bam $DNASEQ.mem.sorted
samtools index $DNASEQ.mem.sorted.bam
bin/qualimap_v2.2/qualimap bamqc -bam $DNASEQ.mem.sorted.bam -outformat pdf --java-mem-size=16G

# remove low quality and short contigs
perl remove_short_lowquality_contigs.pl $DNASEQ.mem.sorted_stats/genome_stats.txt $GENOME > $GENOME_filtered
cegma -g 



```





