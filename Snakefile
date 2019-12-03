
IP = ['SRR495368']
input = ['SRR495378']
SRR = IP + input


## has to run on master node
localrules: fastqdump, get_index_files

rule all:
    input:
        expand("{SRR}_out/{SRR}_deeptools.bw", SRR=SRR),
        "showcase.html"

rule knitR:
    input:
        expand("{SRR}_out/{SRR}.bam", SRR=SRR),
        "peak_calls/SRR495368.bed"
    output:
        "showcase.html"
    script:
        "showcase.Rmd"

rule deeptoolsCoverage:
    input:
        "{SRR}_out/{SRR}.bam"
    output:
        "{SRR}_out/{SRR}_deeptools.bw"
    threads:
        8
    shell:
        """
        bamCoverage --bam {input} -o {output} \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 142573017 \
        --numberOfProcessors {threads} \
        --extendReads 150
        """

rule homer:
    input:
        IP=expand("{IP}_out/{IP}.bam", IP=IP),
        input=expand("{input}_out/{input}.bam", input=input)
    output:
        "peak_calls/{IP}.bed"
    shell:
        """
        makeTagDirectory tags/input.{wildcards.IP}/ {input.input}
        makeTagDirectory tags/ip.{wildcards.IP}/ {input.IP}
        findPeaks tags/ip.{wildcards.IP}/ -i tags/input.{wildcards.IP}/ -fragLength 150 -inputFragLength 150 -C 0 > peak_calls/{wildcards.IP}.txt
        pos2bed.pl peak_calls/{wildcards.IP}.txt > peak_calls/{wildcards.IP}.bed
        rm -rf tags/input.{wildcards.IP}
        rm -rf tags/ip.{wildcards.IP}
        """

rule bamSortIndex:
    input:
        "{SRR}_out/{SRR}.sam"
    output:
        "{SRR}_out/{SRR}.bam"
    threads:
        8
    shell:
        """
        samtools view -bS -o {wildcards.SRR}_out/aligned.bam {input}
        samtools sort -m 2G -@ {threads} -o {output} {wildcards.SRR}_out/aligned.bam
        samtools index {output}
        rm {wildcards.SRR}_out/aligned.bam
        """

rule bowtie:
    input:
        files="raw/{SRR}.fastq",
        index="bowtie_index/genome.1.ebwt"
    output:
        temp("{SRR}_out/{SRR}.sam")
    threads:
        16
    shell:
        """
        bowtie bowtie_index/genome -p 8 -m 1 -S {input.files} > {wildcards.SRR}_out/{wildcards.SRR}.sam
        """

# take only 5 Mio reads, there are 10x more
rule fastqdump:
    output:
        temp("raw/{SRR}.fastq")
    shell:
        """
        fastq-dump --outdir raw -X 5000000  {wildcards.SRR}
        """

rule genomefile:
    input:
        "genome.fa"
    output:
        "genomeFile.txt"
    shell:
        """
        samtools faidx {input}
        awk -v OFS='\t' {{'print $1,$2'}} genome.fa.fai > {output}
        """

rule index_bowtie:
    input:
        "genome.fa"
    output:
        "bowtie_index/genome.1.ebwt"
    threads: 8
    shell:
        """
        bowtie-build --threads {threads} {input} bowtie_index/genome
        """

rule get_index_files:
    output:
        "genome.fa"
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
        gunzip Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
        ln -s Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa genome.fa
        """

onsuccess:
        shell("rm *.out; rm -rf tags; rm -rf raw")
