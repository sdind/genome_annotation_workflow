

rule atropos_error:
    #Estimates the sequencing error-rate using the error function of Atropos. This can help to refine subsequent trimming parameters.
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['R1'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['R2'] if config['samples'][wildcards.sample]['type'] == 'paired-end' else None
    output:
        error_report = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/trimming/{sample}/{sample}_error_report.txt'),
        error_rate = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_error_rate.txt')
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/trimming/{sample}_error.log')
    conda:
        '../envs/atropos.yaml'
    threads: 1
    resources:
        mem_mb = 30000
    params:
        runtime = '03:00:00'
    shell:
        """
        if [ -z "{input.R2}" ]; then
            # For single-end
            atropos error -se {input.R1} -a quality -o {output.error_report} 2> {log}
            tail -n 1 {output.error_report} | awk '{{print substr($0,13,4);}}' | awk '{{printf $0/10}}' > {output.error_rate} 2>> {log}
        else
            # For paired-end
            atropos error -pe1 {input.R1} -pe2 {input.R2} -a quality -o {output.error_report} 2> {log}
            tail -n 1 {output.error_report} | awk '{{print substr($0,13,4);}}' | awk '{{printf $0/10}}' > {output.error_rate} 2>> {log}
        fi
        """


rule detect_adapters:
    # This rule detects adapter sequences in raw sequencing reads using Atropos
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['R1'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['R2'] if config['samples'][wildcards.sample]['type'] == 'paired-end' else None
    output:
        adapter_R1 = "results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_adapters_1.fasta"
    params:
        tmp = 'results/2_braker/align_RNA/trimming/{sample}/{sample}_adapters_tmp',
        adapter_R2 = lambda wildcards: f"results/2_braker/align_RNA/trimming/{wildcards.sample}/{wildcards.sample}_atropos_adapters_2.fasta" if config['samples'][wildcards.sample]['type'] == 'paired-end' else None,
        runtime = '30:00:00'
    log:
        'logs/2_braker/align_RNA/trimming/{sample}_atropos_adapters.log'
    conda:
        '../envs/atropos.yaml'
    threads: 1
    resources:
        mem_mb = 30000
    shell:
        """
        echo "Detecting adapters" > {log}
        if [ -z "{input.R2}" ]; then
            # Single-end
            atropos detect -se {input.R1} --no-cache-contaminants -O fasta -o {params.tmp} &>> {log}
            mv {params.tmp}.0.fasta {output.adapter_R1} 2>> {log}
        else
            # Paired-end
            atropos detect -pe1 {input.R1} -pe2 {input.R2} --no-cache-contaminants -O fasta -o {params.tmp} &>> {log}
            mv {params.tmp}.0.fasta {output.adapter_R1} 2>> {log}
            mv {params.tmp}.1.fasta {params.adapter_R2} 2>> {log}
        fi
        """


rule trim_reads:
    # This rule uses Atropos to trim reads by removing identified adapters and error rate
    input:
        R1 = lambda wildcards: config['samples'][wildcards.sample]['R1'],
        R2 = lambda wildcards: config['samples'][wildcards.sample]['R2'] if config['samples'][wildcards.sample]['type'] == 'paired-end' else None,
        adapter_R1 = "results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_adapters_1.fasta",
        error_rate = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_error_rate.txt')
    output:
        trim_R1 = "results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_trimmed_1.fastq",
        report = 'results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_trimming_report.txt'
    log:
        'logs/2_braker/align_RNA/trimming/{sample}_atropos_trimming.log'
    conda:
        '../envs/atropos.yaml'
    threads: 32
    resources:
        mem_mb = 30000
    params:
        adapter_R2 = lambda wildcards: f"results/2_braker/align_RNA/trimming/{wildcards.sample}/{wildcards.sample}_atropos_adapters_2.fasta" if config['samples'][wildcards.sample]['type'] == 'paired-end' else None,
        trim_R2 = lambda wildcards: f"results/2_braker/align_RNA/trimming/{wildcards.sample}/{wildcards.sample}_atropos_trimmed_2.fastq" if config['samples'][wildcards.sample]['type'] == 'paired-end' else None,
        runtime = '10:00:00',
    shell:
        """
        # Extract the error rate from the file
        ERR_RATE=$(cat {input.error_rate})
        ADJUSTED_ERR_RATE=$(echo "if ($ERR_RATE < 0.1) $ERR_RATE else 0.1" | bc)
        if [ -z "{input.R2}" ]; then
            # Single-end
            atropos trim -a "A{{20}}" -a file:{input.adapter_R1} -n 2 -se {input.R1} -o {output.trim_R1} -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 --threads {threads} --stats post -e $ADJUSTED_ERR_RATE --no-cache-adapters --report-file {output.report} &>> {log}
        else
            atropos trim -a "A{{20}}" -a file:{input.adapter_R1} -A file:{params.adapter_R2} -n 2 -pe1 {input.R1} -pe2 {input.R2} -o {output.trim_R1} -p {params.trim_R2} -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 --threads {threads} --stats post -e $ADJUSTED_ERR_RATE --no-cache-adapters --report-file {output.report} &>> {log}
        fi
        """


r1_paths = [sample_data['R1'] for sample_data in config['samples'].values()]
r2_paths = [sample_data.get('R2', None) for sample_data in config['samples'].values() if sample_data.get('R2')]

rule fastqc:
    #RNA-seq quality check using FastQC
    input:
        r1 = r1_paths,
        r2 = r2_paths
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/QC_RNA/out_fastQC"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/QC_RNA/fastqc.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 20
    resources:
        mem_mb = 30000
    params:
        runtime = '30:00:00'
    shell:
        """
        mkdir -p {output.outdir}
        for file in {input.r1} {input.r2}
        do
            if [ -f $file ]; then
                fastqc -t {threads} $file -o {output.outdir}
            fi
        done    
        """


trimmed_r1_paths = ["results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_trimmed_1.fastq".format(sample=sample_name) for sample_name in config['samples'].keys()]
trimmed_r2_paths = ["results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_trimmed_2.fastq".format(sample=sample_name) for sample_name in config['samples'].keys() if config['samples'][sample_name].get('R2')]

rule fastqc_trimmed:
    #FastQC quality check on trimmed reads
    input:
        r1 = trimmed_r1_paths
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/QC_RNA/out_fastQC_trimmed"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/QC_RNA/fastqc_trimmed.log')
    conda:
        '../envs/fastqc.yaml'
    threads: 20
    resources:
        mem_mb = 30000
    params:
        runtime = '30:00:00',
        r2 = trimmed_r2_paths
    shell:
        """
        mkdir -p {output.outdir}
        for file in {input.r1} {params.r2}
        do
            if [[ -f $file ]]; then
                fastqc -t {threads} $file -o {output.outdir}
            fi
        done    
        """


rule multiqc:
    #Aggregate all QC reports into a single MultiQC report
    input:
        fastqc_raw =  directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/QC_RNA/out_fastQC")),
        fastqc_trimmed =  directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/QC_RNA/out_fastQC_trimmed")),
    output:
        report = os.path.join(config['snakemake_dir_path'], "results/2_braker/QC_RNA/multiqc_report.html")
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = '05:00:00' 
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/QC_RNA/multiqc.log')
    conda:
        '../envs/multiqc.yaml'
    shell:
        """
        multiqc {input.fastqc_raw} {input.fastqc_trimmed} --filename {output.report}
        """


rule build_hisat2_index:
    #Build HISAT2 index.
    input:
        asm = config['asm']
    output:
        index_dir = directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/index"))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2_index_build.log')
    conda:
        '../envs/hisat.yaml'
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = '03:00:00'
    shell:
        """
        mkdir -p {output.index_dir}
        hisat2-build {input.asm} {output.index_dir}/genome_index 2> {log}
        """


rule hisat2:
    #Align RNAseq reads to the genome using HISAT2.
    input:
        r1_trimmed = "results/2_braker/align_RNA/trimming/{sample}/{sample}_atropos_trimmed_1.fastq",
        asm = config['asm'],
        index_dir = directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/index"))
    output:
        aln = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sam"),
        aln_summary = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_splicesite.txt")
    threads: 20
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2/{sample}_hisat2.log')
    resources:
        mem_mb = 400000
    params:
        runtime = '50:00:00',
        r2_trimmed = lambda wildcards: f"results/2_braker/align_RNA/trimming/{wildcards.sample}/{wildcards.sample}_atropos_trimmed_2.fastq" if config['samples'][wildcards.sample]['type'] == 'paired-end' else None
    conda:
        '../envs/hisat.yaml'
    shell:
        """
        if [ -f "{params.r2_trimmed}" ]; then
            hisat2 --phred33 --new-summary --novel-splicesite-outfile {output.aln_summary} -p {threads} -x {input.index_dir}/genome_index -1 {input.r1_trimmed} -2 {params.r2_trimmed} -S {output.aln}
        else
            hisat2 --phred33 --new-summary --novel-splicesite-outfile {output.aln_summary} -p {threads} -x {input.index_dir}/genome_index -U {input.r1_trimmed} -S {output.aln}
        fi
        """


rule to_bam:
    #Convert sam into bam file, then sort and index for Apollo and braker
    input:
        aln_sam = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sam")
    output:
        aln_bam = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sorted.bam"),
        bam_index = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sorted.bam.bai")
    threads: 4
    resources:
        mem_mb = 30000
    params:
        runtime = '20:00:00'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        # Convert SAM to BAM
        samtools view -bS -o {output.aln_bam}.temp.bam {input.aln_sam}
        # Sort the BAM file
        samtools sort -o {output.aln_bam} {output.aln_bam}.temp.bam 
        # Index the sorted BAM
        samtools index {output.aln_bam}
        # Remove the temporary BAM
        rm {output.aln_bam}.temp.bam
        """


rule mapping_stats_samtools:
    #This rule produces mapping statistics using samtools.
    input:
        aln_bam = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sorted.bam")
    output:
        samtools_stats = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/hisat2/mapping_stats_samtools/{sample}_mapping_stats_samtools.txt'),
        read_length = temp('results/2_braker/align_RNA/hisat2/mapping_stats_samtools/{sample}_mapping_stats_samtools_RL.txt')
    params:
        wd_samtools = 'results/2_braker/align_RNA/hisat2/mapping_stats_samtools/{sample}/'
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2/{sample}_samtools_stats.log')
    conda:
        '../envs/QC_stat.yaml'
    resources:
        mem_mb = 10000
    threads: 2
    shell:
        """
        samtools stats --threads {threads} {input.aln_bam} > {output.samtools_stats} 2>> {log}
        grep ^RL {output.samtools_stats} | cut -f 2- > {output.read_length} 2>> {log}
        plot-bamstats -p {params.wd_samtools}{wildcards.sample}_mapping_stats_samtools {output.samtools_stats} 2>> {log}
        """


rule mapping_stats_qualimap_bamqc:
    #This rule produces mapping statistics using Qualimap.
    input:
        aln_bam = os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sorted.bam")
    output:
        bamqc_pdf = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/hisat2/mapping_stats_qualimap_bamqc/{sample}_mapping_stats_qualimap_bamqc.pdf'),
        bamqc_html = os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/hisat2/mapping_stats_qualimap_bamqc/{sample}_mapping_stats_qualimap_bamqc.html')
    params:
        wd_bamqc = 'results/2_braker/align_RNA/hisat2/mapping_stats_qualimap_bamqc/{sample}/'
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/align_RNA/hisat2/{sample}_qualimap_stats.log')
    conda:
        '../envs/QC_stat.yaml'
    resources:
        mem_mb = 10000
    threads: 4
    shell:
        """
        qualimap bamqc -bam {input.aln_bam} -hm 3 --collect-overlap-pairs -nr 1000 -nw 500  \
        -outformat PDF:HTML -outdir {params.wd_bamqc} -nt {threads} --java-mem-size={resources.mem_mb}M &>> {log}
        mv {params.wd_bamqc}report.pdf {output.bamqc_pdf} 2>> {log}
        mv {params.wd_bamqc}qualimapReport.html {output.bamqc_html} 2>> {log}
        """






















