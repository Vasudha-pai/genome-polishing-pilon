# genome-polishing-pilon
polishing the genome of the species Lactobacillus iners using pilon tool
In genome assembly, pivotal for analyzing sequences from diverse populations, high-throughput technologies play a central role. They aid in identifying, monitoring, and distinguishing genome variants, leading to multiple assembly versions based on data production, error corrections, and improvements.

The assembly process, involving data acquisition, assembly, and quality control, necessitates careful consideration of tools and metrics to ensure accuracy. Pre-processing raw data, assembling contigs into scaffolds, and quality assessment using metrics like N50 are crucial steps. Misassemblies can arise due to various factors like repetitive regions, sequencing errors, assembly algorithms, and data quality, highlighting the need for meticulous analysis and appropriate algorithm selection.

For Lactobacillus iners, raw sequence data was acquired from the Sequence Read Archive and assessed using FastQC for quality analysis. SPAdes was utilized for sequence assembly, generating multiple contigs for further evaluation. The assembly was assessed using QUAST and ICARUS, revealing misassemblies that required correction. The polished genome was refined using Pilon, correcting 1,337,021 bases, including insertions concerning the reference file used for comparison. Comparative evaluation with the reference genome using tools like MAUVE helped align and rearrange contigs and gaps, elucidating the improvements made in the genome sequence.

The study's methodology, encompassing data acquisition, assembly, error correction, and comparative assessment, underscores the importance of high-throughput technologies and sophisticated tools in advancing genome assembly and improving sequence accuracy.

We first begin with the installation of conda in our systems
```   conda
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      chmod +x Miniconda3-latest-Linux-x86_64.sh
      ./Miniconda3-latest-Linux-x86_64.sh
```
next we add our required chanels along with the tools to download and analyse our data
``` conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n assembly2 nanoplot flye bandage bwa
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download
sudo install bwa
sudo-apt get install fast-qc
```
We also installed the SRA toolkit 
## Downloading and Installing SRA Toolkit

The SRA Toolkit is a set of utilities for interacting with sequencing data from the NCBI Sequence Read Archive. Follow the steps below to download and install the SRA Toolkit on your system.

### Prerequisites

Before you begin, ensure you have the following:
- An operating system compatible with the SRA Toolkit (Windows, macOS, or Linux)
- Internet connection for downloading the toolkit
- Basic command line proficiency

### Steps to Download and Install SRA Toolkit

1. **Download the SRA Toolkit:**
   - Visit the [SRA Toolkit download page](https://ftp.ncbi.nlm.nih.gov/sra/sdk/current/).
   - Select the appropriate version for your operating system.

2. **Extract the downloaded archive:**
   - For Linux and macOS:
     ```sh
     tar -xvzf sratoolkit.current-<OS>.tar.gz
     ```
     Replace `<OS>` with your specific operating system identifier (e.g., `ubuntu64` for Ubuntu 64-bit).

   - For Windows:
     - Use an archive extraction tool (such as WinRAR or 7-Zip) to extract the contents of the downloaded `.zip` file.

3. **Add the SRA Toolkit to your system's PATH:**
   - For Linux and macOS:
     ```sh
     export PATH=$PATH:/path/to/sratoolkit.current-<OS>/bin
     ```
     Replace `/path/to/sratoolkit.current-<OS>/bin` with the path to the extracted directory.

   - For Windows:
     - Open the Control Panel and navigate to `System and Security` -> `System` -> `Advanced system settings`.
     - Click on `Environment Variables`.
     - In the `System variables` section, find the `Path` variable and click `Edit`.
     - Add the path to the `bin` directory of the extracted SRA Toolkit.

4. **Verify the installation:**
   - Open a terminal or command prompt and run:
     ```sh
     fastq-dump --version
     ```

## Downloading Sequence Data from SRA

This section provides steps to download sequence data from the NCBI Sequence Read Archive (SRA) using the SRA Toolkit. As an example, we'll download the sequence data for the accession number `SRR26772089`.

### Prerequisites

Before you begin, ensure you have the following:
- SRA Toolkit installed on your system (see the previous section for installation instructions)
- Internet connection for downloading the sequence data
- Basic command line proficiency

### Steps to Download Sequence Data from SRA

1. **Verify SRA Toolkit installation:**
   - Open a terminal or command prompt and run:
     ```sh
     fastq-dump --version
     ```
   - Ensure you see output indicating the version of the SRA Toolkit, confirming that it is installed correctly.

2. **Download the sequence data:**
   - Open a terminal or command prompt and run the following command:
     ```sh
     fastq-dump SRR26772089
     ```
   - This command will download the sequence data for `SRR26772089` and save it as a `.fastq` file in your current directory.

3. **Verify the downloaded file:**
   - After the download completes, you should see a file named `SRR26772089.fastq` in your current directory.
   - You can check the first few lines of the file to ensure it has been downloaded correctly:
     ```sh
     head SRR26772089.fastq
     ```
   - This command will display the first ten lines of the `.fastq` file.


now that we have our sequence from the SRA, we need to assess the quality of our sequence, which can be done using the Fast-qc tool
### Installing FastQC

If you haven't installed FastQC yet, you can download and install it from [FastQC's official website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
 ```sh
 fastqc SRR26772089.fastq
```
 This command helps to generate a factqc report 
 After FastQC completes, it will generate a .zip file and an .html file in the same directory.
The .html file contains a detailed quality report. Open this file in a web browser to view the report.
## Running QUAST for Quality Assessment

QUAST (Quality Assessment Tool for Genome Assemblies) is a tool for evaluating genome assemblies by comparing them to a reference genome. This section provides steps to run QUAST on your assembly data using Python.

### Prerequisites

Before you begin, ensure you have the following:
- QUAST installed on your system
- Assembly file (e.g., `contigs.fasta`)
- Reference genome file (e.g., `refseq.fasta`)

### Installing QUAST

If you haven't installed QUAST yet, you can download and install it from [QUAST's official website](http://quast.sourceforge.net/).

- For Linux:
  ```sh
  sudo apt-get install quast
  ```
Alternatively you can also use the pip install command
```sh
pip install quast
```
Quast can be run for the following using 
```sh
Python quast.py –large contigs.fasta -r refseq.fasta
```
here the refseq.fasta was taken from the SRA 
The next step that I followed was to assemble the Sequence into assemblies 
Run SPAdes using Python:
```sh
spades.py -k 21,33,55,77 --careful -1 SRR26772089_1.fastq -2 SRR26772089_2.fastq -m 1000 -t 10 -o SRR26772089out
```
This command will run SPAdes with the following options:

-k 21,33,55,77: Specifies k-mer sizes to be used in the assembly.
--careful: Reduces the number of mismatches and short indels.
-1 SRR26772089_1.fastq: Specifies the file containing the forward reads.
-2 SRR26772089_2.fastq: Specifies the file containing the reverse reads.
-m 1000: Allocates up to 1000 GB of memory for the assembly process.
-t 10: Uses 10 threads for the assembly process.
-o SRR26772089out: Specifies the output directory for the assembly results.
Now remember the process of assembly takes time, Be patient. it depends upon the RAM size of the system. Also ensure that u have enough memory if using a larger sequence.

the next step involves the usage of few tools so that we can use this sequence as out input for pilon
Steps to Align Sequences and Process Alignments
Index the reference genome with BWA:
```sh
bwa index refseq.fasta
```
Align the sequences to the reference genome:
```sh
bwa mem refseq.fasta SRR26772089_1.fastq SRR26772089_2.fastq > seq.sam
```
Convert the SAM file to a BAM file using SAMtools:
```sh
samtools view -@ 8 -Sb -o bamfile.bam seq.sam
```
Sort the BAM file:
```sh
samtools index sortedbam.bam
```

Pilon is a software tool designed to improve genome assemblies by correcting bases, fixing misassemblies, and filling gaps. It is particularly useful for polishing assemblies generated from high-throughput sequencing data.

### Key Features of Pilon

- **Base correction:** Identifies and corrects base errors in the assembly.
- **Gap filling:** Detects and fills gaps in the assembly.
- **Misassembly detection and fixing:** Identifies potential misassemblies and corrects them.
- **Structural variation detection:** Finds and reports structural variations such as insertions and deletions.
- **Annotation improvement:** Enhances existing annotations by providing more accurate sequence data.

### Prerequisites

Before using Pilon, ensure you have the following:
- Java 1.8 or higher installed on your system (Pilon is a Java-based tool).
- A genome assembly to be improved.
- Paired-end reads or other high-throughput sequencing data that can be aligned to the assembly.

### Installing Pilon

1. **Download Pilon:**
   - Visit the [Pilon GitHub releases page](https://github.com/broadinstitute/pilon/releases) and download the latest version of Pilon.

2. **Set up Pilon:**
   - Extract the downloaded archive to a directory of your choice.
   - Ensure that the `pilon.jar` file is accessible.
#run pilon  
```sh
pilon --genome GCF_002884705.1_ASM288470v1_genomic.fna --frags SRR5945998.fastq --output polished_assembly
```
 This command will produce an improved assembly and save the results to the specified output directory (`pilon_output`).
   

