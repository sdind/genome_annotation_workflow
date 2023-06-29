
rule fastqc:
    """RNA-seq quality check using FastQC"""

    input:
        RNA_dir = config['RNA_dir']
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_fastQC"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/fastqc.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 20
    resources:
        mem_mb = 30000
    params:
        runtime = '30:00:00'
    shell:
        """
        mkdir {output.outdir}
        cd {input.RNA_dir}
        for file in $(ls {input.RNA_dir})
        do
            SAMPLE=$(basename $file)
            fastqc -t {threads} $SAMPLE -o {output.outdir}
        done    
        """


rule trimm:
    """Trim reads using trimmomatic"""

    input:
        RNA_dir = config['RNA_dir']
    output:
        outTrimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trimmomatic"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/trimm.log')
    conda:
        '../envs/trimm.yaml'
    threads: 20
    resources:
        mem_mb = 30000
    params:
        runtime = '40:00:00',
        trimlog = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/trimlog.log'),
        adapters_file = config['adapters_file']
    shell:
        """
        mkdir results/2_braker/align_RNA/out_trimmomatic
        for filepath in {input.RNA_dir}/*_1.fastq.gz
        do
            i=$(basename "$filepath" _1.fastq.gz)
            echo "$i"
            trimmomatic PE -threads {threads} -phred33 -trimlog {params.trimlog} {input.RNA_dir}/$i\_1.fastq.gz {input.RNA_dir}/$i\_2.fastq.gz results/2_braker/align_RNA/out_trimmomatic/$i\_R1_paired.fq.gz results/2_braker/align_RNA/out_trimmomatic/$i\_R1_unpaired.fq.gz results/2_braker/align_RNA/out_trimmomatic/$i\_R2_paired.fq.gz results/2_braker/align_RNA/out_trimmomatic/$i\_R2_unpaired.fq.gz ILLUMINACLIP:/{params.adapters_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
        done
        """



rule fastqc_trimmed:
    """Create a Database for RepeatModeler and execute RepeatModeler"""
    input:
        outTrimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trimmomatic"))
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_fastQC_trimmed"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/fastqc_trimmed.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 20
    resources:
        mem_mb = 30000
    params:
        runtime = '30:00:00'
    shell:
        """
        mkdir {output.outdir}
        cd {input.outTrimm}
        for file in $(ls {input.outTrimm})
        do
            SAMPLE=$(basename $file)
            fastqc -t {threads} $SAMPLE -o {output.outdir}
        done
        """


rule hisat2:
    """ (Be careful, if genome fasta has long header it will cause problems afterward, so cut header with cut -d ' ' -f1 your_file.fa > new_file.fa) """
    input:
        out_trimm = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/out_trimmomatic")),
        asm = config['asm']
    output:
        aln_out = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2")),
        aln_summary = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/splicesite.txt")
    threads: 20
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2/hisat2.log')
    resources:
        mem_mb = 400000
    params:
        index_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/index"),
        runtime = '50:00:00'
    conda:
        '../envs/hisat.yaml'
    shell:
        """
        mkdir {params.index_dir}
        (hisat2-build {input.asm} {params.index_dir}/genome_index) 2> {log}
        for filepath in {input.out_trimm}/*_R1_paired.fq.gz
        do
             i=$(basename "$filepath" _R1_paired.fq.gz)
             echo "$i"
             hisat2 --phred33 --new-summary --novel-splicesite-outfile {output.aln_summary} -p {threads} -x {params.index_dir}/genome_index -1 {input.out_trimm}/$i\_R1_paired.fq.gz -2 {input.out_trimm}/$i\_R2_paired.fq.gz -S {output.aln_out}/$i\_accepted_hits.sam
        done
        """



rule to_bam:
    """ sam into bam file and sort + index for Apollo """
    input:
        aln_out = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2"))
    output:
        aln_sam = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/test.txt")
    threads: 4
    resources:
        mem_mb = 30000
    params:
        index_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/index"),
        runtime = '20:00:00'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        cd {input.aln_out}
        for file in *.sam
        do
            SAMPLE=${{file%%_*}}
            samtools view -bS -o ${{SAMPLE}}_accepted_hits.bam $file
            samtools sort -o ${{SAMPLE}}_accepted_hits.sorted.bam ${{SAMPLE}}_accepted_hits.bam
            samtools index ${{SAMPLE}}_accepted_hits.sorted.bam 
        done
        echo "ok" > {output.aln_sam}
        """



# I should do like a output function to have proper output file not sheeting as I do















