# Genome Annotation Workflow
This workflow performs genome annotation by integrating RNA-seq and protein data in a series of rules implemented in Snakemake. It includes steps for repeat identification and filtering, genome masking, RNA-seq data quality control and trimming, alignment, gene prediction, and evaluation.
The workflow includes the following steps:
* **Repeat Identification**: Utilizes [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) to identify and classify repetitive elements in the genome assembly.
* **Repeat Filtering**: filters the RepeatModeler library for TEs resembling proteins using [transposonPSI](https://transposonpsi.sourceforge.net), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [ProtExcluder](https://github.com/NBISweden/ProtExcluder) and scripts from [GAAS](https://github.com/NBISweden/GAAS). This step is adapted from the [workflow developed by Verena Kutschera, July 2020](https://github.com/NBISweden/repeatlib_filtering_workflow)
* **Genome Masking**: Masks repetitive elements in the genome assembly using [RepeatMasker](https://www.repeatmasker.org), resulting in a masked assembly file.
* **RNA-Seq Reads QC**: Performs quality checks on RNA-seq data using [FastQC](https://github.com/s-andrews/FastQC).
* **RNA-Seq Reads trimming**: Trims RNA-seq reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and re-evaluates the quality using FastQC
* **RNA-Seq Alignment**: Aligns the trimmed RNA-seq reads to the genome assembly using [HISAT2](http://daehwankimlab.github.io/hisat2/).
* **Annotation**: Performs gene prediction using [BRAKER3](https://github.com/Gaius-Augustus/BRAKER).
* **Annotation evaluation**: Evaluates the gene predictions using [BUSCO](https://busco.ezlab.org).

## Configuration
**Cluster Configuration**: This workflow is tailored for our specific SLURM cluster setup. Please adjust the profile settings accordingly.
**Local Execution**: If you're interested in running this workflow locally, a customised version suitable for local environments is available at [WorkflowHub](https://workflowhub.eu/workflows/569).

The workflow requires a configuration file (config.yaml) to specify the input and output directories, file paths, and other parameters. Make sure to update the configuration file with your specific paths and settings before running the workflow.
RNA-seq data: The raw RNA-seq data should be organized in a directory specified in the configuration file.
The following parameters have to be customized in the configuration file:
* `asm`: The path to the genome assembly file.
* `snakemake_dir_path`: The directory path where the Snakemake workflow files are located.
* `name`: The short name of your species or assembly run
* `RNA_dir`: Directory containing raw RNA-seq reads. **Important:** files must ends with _1.fastq.gz and _2.fastq.gz suffix to indicate reads 1 and 2 in the paired-end data and need to be compressed (to be changed in future)
* `adapters_file`: text file containing the specific sequences of the adapters used for sequencing setup.
* `busco_phylum`: The BUSCO database identifier for the phylum of the organism being assembled. For example, ‘hymenoptera_odb10’ represents the hymenoptera phylum.
All the required tools mentioned in the workflow, will be automatically installed via conda using the provided YAML file during the workflow execution.

**Currently, the workflow exclusively supports paired-end RNA-seq reads.
New changes coming very soon !! The forthcoming update will introduce significant enhancements, including parallelised RNA-seq alignments, automated adapter detection using atropos, and support for single-end reads. 
Please note that this workflow is still a work in progress.**

### Directory Structure
```
.
├── config.yaml    # Configuration file specifying input data and parameters
├── logs        # Log files for each step
├── results      # Directory containing output files for each step
├── README.md     # This README file
└── workflow      # Output directory (generated during workflow execution)
  ├── Snakefile   # Global Snakemake file
  ├── envs      # Environment YAML files for required tools
  ├── rules     # Snakemake rules for each step of the workflow
    ├── Snakefile_1.smk  # Snakemake file for preprocessing and quality assessment steps
    ├── Snakefile_2.smk  # Snakemake file for genome assembly and evaluation steps
  └── scripts    # Custom scripts (if any)
```
