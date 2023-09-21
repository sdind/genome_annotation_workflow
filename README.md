# Genome Annotation Workflow
This workflow performs genome annotation by integrating RNA-seq and protein data in a series of rules implemented in Snakemake. It includes steps for repeat identification and filtering, genome masking, RNA-seq data quality control and trimming, alignment and mapping statistics, gene prediction and evaluation.
The workflow includes the following steps:
* **Repeat Identification**: Utilizes [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) to identify and classify repetitive elements in the genome assembly.
* **Repeat Filtering**: filters the RepeatModeler library for TEs resembling proteins using [transposonPSI](https://transposonpsi.sourceforge.net), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [ProtExcluder](https://github.com/NBISweden/ProtExcluder) and scripts from [GAAS](https://github.com/NBISweden/GAAS). This step is adapted from the [workflow developed by Verena Kutschera, July 2020](https://github.com/NBISweden/repeatlib_filtering_workflow)
* **Genome Masking**: Masks repetitive elements in the genome assembly using [RepeatMasker](https://www.repeatmasker.org), resulting in a masked assembly file.
* **RNA-Seq Reads Trimming and Adapters Detection**: [Atropos](https://github.com/jdidion/atropos) is employed for automated adapter detection, error rate evaluation, and trimming. 
* **RNA-Seq Reads QC**: Performs quality checks on RNA-seq data before and after trimming using [FastQC](https://github.com/s-andrews/FastQC).A single report is provided using [multiQC](https://github.com/ewels/MultiQC).
* **Parallelized RNA-Seq Alignment**: Utilizing [HISAT2](http://daehwankimlab.github.io/hisat2/), trimmed RNA-seq read alignments to the genome assembly are executed in parallel, optimizing performance.
* **Mapping Statistics Analysis**: Generate detailed mapping statistics, qualitative alignment analysis using [samtools](https://github.com/samtools/samtools) and [qualimap](http://qualimap.conesalab.org/).
* **Annotation**: Use the mapped RNAseq reads and the uniprot sequences to create hints for gene prediction using [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) on the masked genome
* **Annotation evaluation**: Run [BUSCO](https://busco.ezlab.org) to evaluate the completeness of the annotation produced.

<br/>

## Prerequisites

The following programs are required to run the workflow and the listed version were tested. It should be noted that older versions of snakemake are not compatible with newer versions of singularity as is noted here: https://github.com/nextflow-io/nextflow/issues/1659.

* conda v 23.7.3
* singularity v 3.7.3
* snakemake v 7.32.3

You will also need to acquire a licence key for Genemark and place this in your home directory with name ~/.gm_key The key file can be obtained from the following location, where the licence should be read and agreed to: http://topaz.gatech.edu/GeneMark/license_download.cgi


### Required Input data
* Reference genome in fasta format
* RNAseq data in fastq format
* protein database in fasta format

<br/>

## Configuration
**Cluster Configuration**: This workflow is tailored for our specific SLURM cluster setup. Please adjust the profile settings accordingly.
**Local Execution**: If you're interested in running this workflow locally, a customised version suitable for local environments is available at [WorkflowHub](https://workflowhub.eu/workflows/569).

The workflow requires a configuration file (config.yaml) to specify the input and output directories, file paths, and other parameters. Make sure to update the configuration file with your specific paths and settings before running the workflow.
The following parameters have to be customized in the configuration file:
* `asm`: The path to the genome assembly file.
* `snakemake_dir_path`: The directory path where the Snakemake workflow files are located.
* `name`: The short name of your species or assembly run
* `busco_phylum`: The BUSCO database identifier for the phylum of the organism being assembled. For example, ‘hymenoptera_odb10’ represents the hymenoptera phylum.
* `prot`: protein fasta sequences in zipped fasta format (e.g from [Swiss-Prot](https://www.uniprot.org/help/downloads)) or [OrthoDB](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/).
* `samples`: Contains details about the RNA-seq samples including type and paths to the fastq files.
   - Subkeys:
     - `Sample Name` (e.g., `17_L4_001`):
       - `type`: Defines the RNA-seq sample as 'paired-end' or 'single-end'.
       - `R1`: Path to the first/fastq file in the pair or the single-end read file.
       - `R2`: Path to the second/fastq file in the pair (only for paired-end).
   - Example:
     ```
     samples:
       17_L4_001:
         type: 'paired-end'
         R1: '/path/to/sample1.fastq_1.gz'
         R2: '/path/to/sample1.fastq_2.gz'
     ```

**Important**: RNA-seq files must ends with _1.fastq.gz and _2.fastq.gz suffix to indicate reads 1 and 2 in the paired-end data and need to be compressed (to be changed in future).

All the required tools mentioned in the workflow, will be automatically installed via conda using the provided YAML file during the workflow execution.

**Please note that this workflow is still a work in progress.**

### Directory Structure
```
.
├── config.yaml    # Configuration file specifying input data and parameters
├── logs        # Log files for each step
├── results      # Directory containing output files for each step
├── README.md     # This README file
└── workflow      # Workflow Directory 
  ├── Snakefile   # Global Snakemake file
  ├── envs      # Environment YAML files for required tools
  ├── rules     # Snakemake rules for each step of the workflow
    ├── 1_MaskRepeat.smk  # Snakemake file for identifying and masking repeats in the genome
    ├── 2_alignRNA.smk  # Snakemake file to trim and map reads to the input genome
    ├── 3_braker.smk  # Snakemake file to predict protein coding gene using braker3
  └── scripts    # Custom scripts for the workflows
```
