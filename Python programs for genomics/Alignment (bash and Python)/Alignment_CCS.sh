FILE=$1
FIRST=$(head -n 1 $FILE)

# Function for creating sam files with bwa and for merging-sorting them
function bwa_aln_SC_sam {

# Creating index for bwa alignments with S.cerevisiae reference genome
bwa index S.cerevisiae_refg.fna

# Generate SAM files for each splitted file
        # Loop counter
        n=0
        # Create a directory for SAM files
	DIR=SAM_FILES
	mkdir $DIR
        # Loop for trimming and aligning sequences
        for file in $PATH_FILES
        do
                OUT_NAME=$file".trimmed"

                # Trimming each generated sequence from right (20 nc)
                python3 fast_processing_CCS.py --input $file --output $OUT_NAME --operation trim --trim-right 20
                rm $file

                # Starting to do bwa alignments
                # Create a SAM file with mem algorithm
                bwa mem S.cerevisiae_refg.fna $OUT_NAME > aln-$n.sam
                # Place alignments in SAM_files directory
                mv aln-$n.sam $DIR
                # Update counter
                n=$(expr $n + 1)
        done

	# Merging SAM files ignoring tags
        cat ./$DIR/*.sam | grep -vE ^"@" > merged_sam_files.sam

        # Sorting merged SAM file by chromosome (field 3 and 4) with Refseq annotation
        cat merged_sam_files.sam | sort -r -k 3 -k 4 > sorted_merged.sam

	# Removing directory with splitted sam files and non-sorted sam file
	rm -r $DIR
	rm merged_sam_files.sam
}

# FASTA files
if [[ $FIRST =~ ^">" ]]
then
	# Possible trimers (histogram)
	cat <(grep -vE "^>" $FILE | fold -w 3) \
	<(grep -vE "^>" $FILE | fold -w 1 | tail +2 | tr -d "\n" | fold -w 3) \
	<(grep -vE "^>" $FILE | fold -w 1 | tail +3 | tr -d "\n" | fold -w 3) \
	| sort | uniq -c | awk '{if (length($2)==3) print $0}' \
	| awk '{ sub(/^[ \t]+/, ""); print}' > histogramF.txt 
	# This was for filtering just trimers: awk '{if (length($2)==3) print $0}'
	# This was for to trim whitspaces: awk '{ sub(/^[ \t]+/, ""); print}'

	# Splitting one record per file
	mkdir output_sequences_fasta
	split -d -l 2 $FILE output_sequences_fasta/sequence.fasta.
	PATH_FILES=output_sequences_fasta/sequence.fasta.*

	# Calling function for bwa alignments in order to generate a merged and sorted sam file
	bwa_aln_SC_sam

	# Gather output files into a directory and remove splitted files
        mkdir FASTA_aln_output
        mv histogramF.txt FASTA_aln_output
        mv sorted_merged.sam FASTA_aln_output
	rm -r output_sequences_fasta

	# Printing on the terminal the number of reads aligned
	NO_READS=$(cut -f3 ./FASTA_aln_output/sorted_merged.sam | grep -vE "\*" | wc -l)
	echo -e "\n \n \n *** THE NUMBER OF ALIGNED READS IS: $NO_READS ***"

# FASTQ file
elif [[ $FIRST =~ ^"@" ]]
then 
	# Possible trimers (histogram)
	cat <(cat $FILE | paste - - - - | cut -f2 | fold -w 3) \
	<(cat $FILE | paste - - - - | cut -f2 | fold -w 1 | tail +2 | tr -d "\n" | fold -w 3) \
	<(cat $FILE | paste - - - - | cut -f2 | fold -w 1 | tail +3 | tr -d "\n" | fold -w 3) \
	| sort | uniq -c | awk '{if (length($2)==3) print $0}' \
        | awk '{ sub(/^[ \t]+/, ""); print}' > histogramFQ.txt

	# Splitting one record per file
        mkdir output_sequences_fastq
        split -d -l 4 $FILE output_sequences_fastq/sequence.fastq.
	PATH_FILES=output_sequences_fastq/sequence.fastq.*

        # Calling function for bwa alignments in order to generate a merged and sorted sam file
        bwa_aln_SC_sam

	# Gather output files into a directory and remove splitted files
	mkdir FASTQ_aln_output
	mv histogramFQ.txt FASTQ_aln_output
	mv sorted_merged.sam FASTQ_aln_output
	rm -r output_sequences_fastq

	# Printing on the terminal the number of reads aligned
	NO_READS=$(cut -f3 ./FASTQ_aln_output/sorted_merged.sam | grep -vE "\*" | wc -l)
        echo -e "\n \n \n *** THE NUMBER OF ALIGNED READS IS: $NO_READS ***"

fi
