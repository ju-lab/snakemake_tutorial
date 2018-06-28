# Snakefile to run FASTQ -> Processed BAM from a given forward/reverse fastq files. 
# 2018.06.06 Jongsoo Yoon

configfile: 'pathConfig.yaml'
configfile: 'sampleConfig.yaml'

rule all:
    input:
        'done'


rule bwa_align:
    input:
        bwa = config['bwa'],
        samtools = config['samtools'],
        ref = config['reference'],
        fq1 = lambda wildcards: config['samples'][wildcards.sample]['fq1'],
        fq2 = lambda wildcards: config['samples'][wildcards.sample]['fq2'], 
    output:
        bam = temp("temp_dna/{sample}.temp.bam"),
        sortedbam = temp('temp_dna/{sample}.temp.sorted.bam')
    params:
        rg = "@RG\tID:{sample}\tSM:{sample}\tPL:Illumina", 
    threads: 4
    log:
        "logs/{sample}.bwa.log"
    shell:
        "({input.bwa} mem -t {threads} -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} |"
        "{input.samtools} view -Sb - > {output.bam}; {input.samtools} sort -@ {threads} -o {output.sortedbam} {output.bam}) "
        " &> {log} " # redirects all stdin/out


rule markdup:
    input:
        java = config['java8'], 
        picard = config['picard'],
        samtools = config['samtools'],
        sortedbam = "temp_dna/{sample}.temp.sorted.bam", 
    output:
        mdbam = temp("temp_dna/{sample}.temp.sorted.md.bam"),
        mdbai = temp("temp_dna/{sample}.temp.sorted.md.bam.bai")

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.md.log"
    shell:
        "({input.java} -XX:ParallelGCThreads={threads} -Xmx{resources.mem_mb}m -jar {input.picard} MarkDuplicates "
        "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={input.sortedbam} O={output.mdbam} "
        "M={output.mdbam}.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR=temp_dna/md_temp QUIET=true; "
        "{input.samtools} index {output.mdbam}) "
        "&> {log} " # redirects all stdin/out

rule realign:
    input:
        java = config['java8'],
        gatk = config['gatk'],
        samtools = config['samtools'],
        ref = config['reference'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}.temp.sorted.md.bam", 
        mdbai = "temp_dna/{sample}.temp.sorted.md.bam.bai"
    output:
        realignedbam = temp("temp_dna/{sample}.temp.sorted.md.ir.bam"),
        realignedbai = temp("temp_dna/{sample}.temp.sorted.md.ir.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.realign.log"
    shell:
        "({input.java} -Xmx{resources.mem_mb}m -jar {input.gatk} -T RealignerTargetCreator -R {input.ref} "
        " -I {input.bam} --known {input.knownindel} -o {output.realignedbam}.intervals; "
        "{input.java} -Xmx{resources.mem_mb}m -jar {input.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} "
        "-targetIntervals {output.realignedbam}.intervals -o {output.realignedbam}; "
        "{input.samtools} index {output.realignedbam}) "
        "&> {log}"


rule baserecal:
    input:
        java = config['java8'], 
        gatk = config['gatk'], 
        samtools = config['samtools'],
        ref = config['reference'], 
        dbsnp = config['dbsnp'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}.temp.sorted.md.ir.bam"
    output:
        recalTable = 'temp_dna/{sample}.recaltable', 
        recalbam = protected("dna_bam/{sample}.bam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.baserecal.log"
    shell:
        "({input.java} -Xmx{resources.mem_mb}m -jar {input.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -knownSites {input.dbsnp} --knownSites {input.knownindel} -o {output.recalTable}; "
        "{input.java} -Xmx{resources.mem_mb}m -jar {input.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {output.recalTable} -o {output.recalbam} -nct {threads}; "
        "{input.samtools} index {output.recalbam}) "
        " &> {log}"

rule combine:
    input:
        outfile = expand('dna_bam/{sample}.bam', sample=config['samples'])
    output:
        'done'
    log:
        "logs/done.log"
    shell:
        "(touch {output}) &> {log}"
