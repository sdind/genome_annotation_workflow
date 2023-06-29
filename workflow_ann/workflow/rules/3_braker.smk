

rule braker:
    input: 
        asm_masked = os.path.join(config['snakemake_dir_path'], 'results/1_MaskRepeat/RepeatMasker/asm_decontaminated.fasta.masked'),
        aln_sam = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2/test.txt") # to replace by proper output
    output:
        os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker/braker/braker.aa") 
    threads: 20
    resources:
        mem_mb = 100000
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_braker/out_braker/out_braker.log')
    params:
        protDB = os.path.join(config['snakemake_dir_path'], 'workflow/files/proteins.fasta'),
        bams_dir = os.path.join(config['snakemake_dir_path'],"results/2_braker/align_RNA/hisat2"),
        runtime = '40:00:00',
        out_dir = directory(os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker"))

    singularity:
        'docker://teambraker/braker3:latest'
    shell:
        """
        list_aln=(`ls {params.bams_dir}/*.sorted.bam`)
        mkdir {params.out_dir}
        cd {params.out_dir}
        (braker.pl --genome={input.asm_masked} --prot_seq={params.protDB} --bam="${{list_aln[@]}}" --softmasking --threads {threads} --gff3) 2> {log}
        """



rule eval:
    input: 
        braker_aa = os.path.join(config['snakemake_dir_path'],"results/2_braker/out_braker/braker/braker.aa")
    output: 
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_braker/braker_busco'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_braker/busco/busco.log')
    conda:
        '../envs/busco5.yaml'
    threads: 10
    resources:
        mem_mb = 30000
    params:
        out_path =  os.path.join(config['snakemake_dir_path'], 'results/2_braker'),
        out_name = 'braker_busco',
        lineage = config['busco_phylum'],
        runtime = '4:00:00',
    shell:
        """
        cd {params.out_path}
        (busco -f -m prot -i {input.braker_aa} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads/
        """
