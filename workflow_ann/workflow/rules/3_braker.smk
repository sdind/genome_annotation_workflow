
rule braker:
    # Uses braker3 to predict protein-coding genes annotation using previouly align RNA-seq and curated protein data 
    input:
        asm_masked = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatMasker/', os.path.basename(config['asm'])) + '.masked',
        bams = expand(os.path.join(config['snakemake_dir_path'], "results/2_braker/align_RNA/hisat2/{sample}_accepted_hits.sorted.bam"), sample=config['samples'].keys())
    output:
        braker_aa = os.path.join(config['snakemake_dir_path'], "results/2_braker/out_braker/braker/braker.aa")
    threads: get_threads('braker')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('braker', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('braker', attempt)
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/out_braker/out_braker.log')
    params:
        protDB = config['prot'],
        out_dir = directory(os.path.join(config['snakemake_dir_path'], "results/2_braker/out_braker"))
    singularity:
        'docker://teambraker/braker3:latest'
    shell:
        """
        bams_comma_sep="$(echo '{input.bams}' | tr ' ' ',')"
        mkdir -p {params.out_dir}
        cd {params.out_dir}
        (braker.pl --genome={input.asm_masked} --prot_seq={params.protDB} --bam="$bams_comma_sep" --softmasking --threads {threads} --gff3) 2> {log}
        """


rule eval:
    # Uses Busco5 to evaluate the completness of the predicted annotation
    input: 
        braker_aa = os.path.join(config['snakemake_dir_path'], "results/2_braker/out_braker/braker/braker.aa")
    output: 
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_braker/braker_busco'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_braker/busco/busco.log')
    conda:
        '../envs/busco5.yaml'
    threads: get_threads('eval')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('braker', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('braker', attempt)
    params:
        out_path =  os.path.join(config['snakemake_dir_path'], 'results/2_braker'),
        out_name = 'braker_busco',
        lineage = config['busco_phylum']
    shell:
        """
        cd {params.out_path}
        (busco -f -m prot -i {input.braker_aa} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads/
        """

