import os
import re
import sys

def open_input_output(file):
    f = open(file, "rt") 
    line = f.readline()
    # If FASTA or FASTQ format is checked by looking for the first character of the first line (tag)
    if line[0] == "@":
         input_format = "FASTQ"
    elif line[0] == ">":
        input_format = "FASTA"
    # We've already read the first line of the file to know the format. Let's restart the reading order from the beginning
    f.seek(0)
    return f, input_format

# Histogram
def kmer_profile(sequence,k=3): 
    profile = {}
    for i in range(0,len(sequence)-k+1): 
        # In order to avoid masked bases we use capitalization (same entries in the dictionary)
        kmer = sequence[i:i+k].upper()
        if kmer not in profile:
            profile[kmer] = 0 
        profile[kmer] += 1
    return profile

# Create histogram for each file and split file into one sequence files
def split_file_hist(file):
    # Empty histogram
    file_histogram = {}
    # Loop counter
    n = 0
    # We assign local variables returned from open_input_output function to local variables of the outer function
    (f, input_format) = open_input_output(file)
    # Let's start reading the file. This function read and writes both: FASTA and FASTQ formats
    while True:
        # Tag
        tag = f.readline()[1::].rstrip("\n")
        if not tag: break
        # Sequence
        sequence = f.readline().rstrip("\n")
        # Update histogram with each sequence
        trimers = kmer_profile(sequence)
        for trimer in trimers.items():
            if trimer[0] not in file_histogram: file_histogram[trimer[0]] = trimer[1]    
            else: file_histogram[trimer[0]] += trimer[1]
        # Write FASTA or FASTQ
        if input_format == "FASTQ":
            # Ignoring "+"
            f.readline().rstrip("\n")
            # Qualities' line
            qualities = f.readline().rstrip("\n")
            f_out = open("file%d.fastq" % (n), "wt")
            f_out.write("@%s\n%s\n+\n%s\n" % (tag, sequence, qualities))
            f_out.close()
        else:
            f_out = open("file%d.fasta" % (n), "wt")
            f_out.write(">%s\n%s\n" % (tag, sequence))
            f_out.close()
        # Update the counter
        n += 1
    f.close()
    return file_histogram


def trimmed_sequences():
    # List of files in the directory
    files = os.listdir('.')
    # Create a directory for trimmed sequences
    os.mkdir("./TRIMMED_FILES")
    for i in files:
        # We search on the directory the previously splitted files with Regex
        match = re.search(r'(file)[1234567890]+(.fasta|.fastq)', i)
        if match:
            f = open(i, "rt")
            # Tag
            tag = f.readline().rstrip("\n")
            # Sequence 
            sequence = f.readline().rstrip("\n")
            if tag[0] == "@":
                # +
                f.readline().rstrip("\n")
                # Qualities
                qualities = f.readline().rstrip("\n")
            # When we open again a document to be written, even if there is already content, it starts to write at the first line
                # Writting FASTQ
                writef = open(i, "wt")
                writef.write("%s\n%s\n+\n%s\n" % (tag, sequence[:-20:], qualities[:-20:]))
            # Writting FASTA
            else: 
                writef = open(i, "wt")
                writef.write("%s\n%s\n" % (tag, sequence[:-20:]))
            # Stop reading and writting
            f.close()
            writef.close()
            # Move the trimmed file to a directory
            os.system("mv %s TRIMMED_FILES" % (i))
        
def bwa_aln_sam_SC():
    # Create an index with the reference genome
    os.system("bwa index S.cerevisiae_refg.fna")
    # Create a directory for new SAM files
    os.mkdir("./SAM_FILES")
    # List all the files for alignment
    files = os.listdir('./TRIMMED_FILES')
    # Set a counter
    n = 0
    for file in files:
        os.system("bwa mem S.cerevisiae_refg.fna TRIMMED_FILES/%s > aln-%d.sam" % (file, n))
        os.system("mv aln-%d.sam SAM_FILES" % (n))
        n += 1
    # Remove previous directory with trimmed files
    os.system("rm -r ./TRIMMED_FILES")

def merge_sam_files():
    # Find all the sam files previously created
    files = os.listdir('./SAM_FILES')
    # Openning a file for writting all sequences merged
    mf = open("merged_files.sam", "wt")
    # Sam headears pattern
    pattern = r'^@(SQ|HD|RG|PG|CO){1}'
    for file in files:
        with open("SAM_FILES/%s" % (file), "rt") as f:
            for line in f:
                match = re.search(pattern, line)
                # Ignoring sam headers
                if match:
                    continue
                # Writting reads
                else:
                    mf.write(line)
    mf.close()
    # Remove directory with splitted sam files
    os.system("rm -r SAM_FILES")

# Compares chromosome and position of sam records
def compare_sam(sam1,sam2):
    # Split the strings in order to be able to campare chromosome and position 
    list1 = sam1.split("\t")
    list2 = sam2.split("\t")
    # Chromosome 
    if list1[2] < list2[2]:
        return -1
    elif list1[2] > list2[2]:
        return 1 
    # Position
    else:
        if list1[3] < list2[3]: 
            return -1
        elif list1[3] > list2[3]:
            return 1
        else:
            return 0
    
# Merge lists when they are sorted with the function compare_sam
def merge_sam_lists(list1,list2,merged): 
    p1 = 0
    p2 = 0
    pmerged = 0
    # Merge both list by comparing each element at the head 
    while p1 < len(list1) and p2 < len(list2):
        if compare_sam(list1[p1],list2[p2]) <= 0: 
            merged[pmerged] = list1[p1]
            p1 += 1
        else:
            merged[pmerged] = list2[p2] 
            p2 += 1
        pmerged += 1
    # Append the remaining elements of list1 
    while p1 < len(list1):
        merged[pmerged] = list1[p1] 
        p1 += 1
        pmerged += 1
    # Append the remaining elements of list2
    while p2 < len(list2): 
        merged[pmerged] = list2[p2] 
        p2 += 1
        pmerged += 1
    
# Function with divide and conquer algorithm for sorting and merging a list holding lines from sam file
def merge_sort_sam_list(sam_list):
    # Check list length
    if len(sam_list) == 1:
        return 
    else:
        middle_point = len(sam_list) // 2
        left_list = sam_list[:middle_point] # Left-half of the list 
        right_list = sam_list[middle_point:] # Right-half of the list 
        # Sort recursively the halves
        merge_sort_sam_list(left_list) 
        merge_sort_sam_list(right_list)
        # Merge both halves into the final sorted list 
        merge_sam_lists(left_list,right_list,sam_list)

# This function writes in a file sorted lines from a sam file
def sort_sam_file(file):
    with open(file, "rt") as f:
        # Storing lines of merged sam file
        lines = f.readlines()
        # Calling a recursive function with divide & conquer algorithm for sorting the list by chromosome and position
        merge_sort_sam_list(lines)
        # Let's write the sorted lines in a new file
        f_out = open("sorted_merged.sam", "wt")
        # Reverse the list with sorted lines in order to put chromosome before "*" character (non-aligned sequences)
        for elem in lines[::-1]:
            f_out.write(elem)
        f_out.close()
        # Remove merged but not sorted file
        os.system("rm %s" % (file))
        
def compute_alignments_sam(file):
    with open(file, "rt") as f:
        no_alignments = 0
        for line in f:
            fields = line.split("\t")
            if fields[2] == "*":
                continue
            else:
                no_alignments += 1
    return no_alignments

# Call sequentially each process and store the histogram and the number of alignments
histogram = split_file_hist(str(sys.argv[1]))
trimmed_sequences()
bwa_aln_sam_SC()
merge_sam_files()
sort_sam_file("merged_files.sam")
no_alignments = compute_alignments_sam("sorted_merged.sam")

# Print on the terminal the results 
print("\n")
print("The histogram of the whole file is:")
print(histogram)
print("\n")
print("*** THE NUMBER OF READS ALIGNED IS:  %d ***" % (no_alignments))
