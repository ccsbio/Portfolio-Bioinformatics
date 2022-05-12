# FASTA/FASTQ alignments with the reference genome of *Saccharomyces cerevisiae*  

In this directory can be found two scripts (Bash/Python) that, given an input FASTA/FASTQ file, implement the following pipeline:  

1. Counts the abundance of each possible 3-mers (histogram)
2. Splits the input in FASTA/FASTQ into files of only 1 sequence each
3. Hard-trims (from the right) all sequences from all files 20nt
4. Using the genome mapping tool BWA and the reference genome of the Scaromice Cerevisiae (any strain will do), aligns each of the files producing its corresponding SAM file.
5. Merges all SAM files ignoring headers
6. Sorts the SAM file by chromosome and position
7. Computes how many reads have been aligned  

## Dependencies  
Both scripts need the following dependencies:  

1. They need a Linux environment. They were programmed with the Ubuntu distribution and some commands are used by the terminal while running.
2. **Bwa** package should be installed. It can be done easely by using two command lines sequentially:

```
sudo apt-get update
```  

```  
sudo apt-get install bwa
```  

3. The file called S.cerevisiae_refg.fna which is the reference genome of *Saccharomyces cerevisiae*.
4. The Python script called fast_processing_CCS.py, which is used in case the **Bash** script has been used. This script pre-processes the FASTA/FASTQ files.

## Usage  

- With **Python** you can use the follqing command line input:  

```
python Alignment_CCS.py fasta_test.fasta
```  

- With **Bash** you can use the following command line input:  

```  
bash Alignment_CCS.sh fasta_test.fasta
```

## Expected output  
- With **Python**: a file called sorted_merged.sam containing the alignsments against *Saccharomyces cerevisiae*. In addition, the histogram of possible trimers and the number of alignments are displeayed on the command line.
- With **Bash**: a directory called FASTA_aln_output containing a text file with the histogram and a file called sorted_merged.sam with the alignments.  

In addition:  
The **Bwa** package creates five files for indexing the reference genome of the input. These files have the following extensions:  
- .amb
- .ann
- .bwt
- .pac
- .sa
