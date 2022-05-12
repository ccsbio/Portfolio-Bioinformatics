# Let's start parsing the arguments with the argparse module

import argparse

# A parser is created and the arguments are defined
# Only --input, --output and --operation are obligatory
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--operation", required=True)
parser.add_argument("--trim-right", type=int, required=False)
parser.add_argument("--trim-left", type=int, required=False)
parser.add_argument("--adaptor", required=False)

# This variable gathers the parsed arguments
args = parser.parse_args()

# In the following lines, first of all the functions and dictonaries needed are written 

# Dictionary for nucleotides conversion
nucleotides = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# Dictionaries for statistics
stats = {"reads": 0, "bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0} 
trimmed_stats = {"bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
adaptor_stats = {"adaptors-found": 0}

# Function for opening the file that is going to be read (--input) and a file to be written (--output)
def open_input_output():
    f = open(args.input, "rt") 
    f_output = open(args.output, "wt")
    line = f.readline()
    # If FASTA or FASTQ format is checked by looking for the first character of the first line (tag)
    if line[0] == "@":
         input_format = "FASTQ"
    elif line[0] == ">":
        input_format = "FASTA"
    # We've already read the first line of the file to know the format. Let's restart the reading order
    # from the beginning
    f.seek(0)
    return f, f_output, input_format


# This function would be called for reading and processing a file and, additionally, for writting the output in the original file format
def process_file(file = args.input):
    # We assign local variables returned from open_input_output function to local variables of the outer function
    (f, f_output, input_format) = open_input_output()
    # Let's start reading the file. This function read and writes both: FASTA and FASTQ formats
    while True:
        # Tag line
        # The first character, which is different in FASTA and FASTQ, is ignored in order to recycle tag line independently of the format 
        tag = f.readline()[1::].rstrip("\n")
        # If this condition is met the while loop ends iteration
        if not tag: break
        # Sequence line 
        sequence = f.readline().rstrip("\n")
        # In the following line, a function that processes the previously stored sequence is called in order to do the operation inputted by the user.
        # A boolean is also returned in order to be a linker between sequences' process and qualities' process 
        (processed_sequence, bool_result) = calling_operation(sequence)
        # Now additional lines for FASTQ format are stored
        if input_format == "FASTQ":
            # Ignoring "+"
            f.readline().rstrip("\n")
            # Qualities' line
            qualities = f.readline().rstrip("\n")
            # We have to modify qualitie's line also in order to suit qualities with the processed DNA sequence
            # Important, I use the same function for processing sequence and qualities. So, the returned boolean in the previous step is important 
            # We store again booleans in order to not return them to the output file
            (processed_qualities, qbool_result) = calling_operation(qualities, boolean = bool_result)
            # Writting FASTQ format by adding "@" and "+" manually
            f_output.write("@%s\n%s\n+\n%s\n" % (tag,processed_sequence,processed_qualities))
        else:
            # Writting FASTA format by adding ">" manually in the tag's line
            f_output.write(">%s\n%s\n" % (tag, processed_sequence))
        # Update stats
        stats_update(sequence)
    print_summary()
    f.close()
    f_output.close()


# Trimming sequence operation. Can be used equally for DNA sequences or qualities
# Stats of trimmed fragments are also computed, just in case of DNA sequences' lines
# If boolean is False, then the line to process is a DNA sequence. But if is True it has been already processed, so it's qualities and stats aren't needed 
def trim_sequence(sequence, boolean = False):
    # Here the sequence is cut from the beginning and from the end
    if args.trim_right is not None and args.trim_left is not None:    
        # Trimmed fragments stats: first trimmed fragments are stored and then they are used as a parameter for the trimmed_stats_update just if needed
        trimmed_fragments = sequence[0:args.trim_left:] + sequence[len(sequence)-args.trim_right::]
        if boolean == False: trimmed_stats_update(trimmed_fragments)
        # Sequence trimmed
        return sequence[args.trim_left:len(sequence)-args.trim_right:], True
    # I cut the sequence from the beggining
    elif args.trim_left is not None:
        # Trimmed fragments stats
        trimmed_fragments = sequence[0:args.trim_left:]
        if boolean == False: trimmed_stats_update(trimmed_fragments)
        # Sequence trimmed
        return sequence[args.trim_left::], True
    # I cut the sequence from the end, so is length of sequence minus the number of nt to cut
    elif args.trim_right is not None:
        # Trimmed fragments stats
        trimmed_fragments = sequence[len(sequence)-args.trim_right::]
        if boolean == False: trimmed_stats_update(trimmed_fragments)
        # Sequence trimmed
        return sequence[:len(sequence)-args.trim_right:], True


# Adaptor removal operation
def adaptor_rem(sequence, adaptor, boolean = False):
    # This step checks if the adaptor is equal to the beginning of the sequence (same length as adaptor to see if it is found here)
    # Additionally, we also check if it is a sequence or a quality line with boolean = True
    if adaptor == sequence[0:len(adaptor)] or boolean == True:
        # This computes adaptor count in adaptor_stats dictionary, just in case of DNA sequences' lines (not qualities)
        if adaptor == sequence[0:len(adaptor)]: adaptor_stats["adaptors-found"] += 1
        # Line is returned without the beginning, just in adaptor-found cases (for both sequence and quality). True is also returned to process also quality's line of the current case 
        return sequence[len(adaptor)::], True
    # The original sequence is returned if the adaptor is not found and the boolean False, therefore quality line for this sequence won't be processed
    else:
        return sequence, False


# Generic stats, independently of the operation called. All these updates are stored in the stats dictionary 
def stats_update(f_line):
    stats["reads"] += 1
    stats["bases"] += len(f_line)
    for nt in f_line:
        stats[nt] += 1


# Stats of trimmed fragments of DNA. Just when --operation in trim. A trimmed fragment has to be used as a function parameter
def trimmed_stats_update(fragments):
    for nt in fragments:
        trimmed_stats[nt] += 1
    if args.trim_right is not None and args.trim_left is not None: trimmed_stats["bases"] += args.trim_right + args.trim_left
    elif args.trim_right is not None: trimmed_stats["bases"] += args.trim_right
    elif args.trim_left is not None: trimmed_stats["bases"] += args.trim_left

    
# This function checks and calls the operations inputted by the user. 
# Some arguments of this function have to be equally introduced in the called function such as: f_line or boolean. 
def calling_operation(f_line, operation = args.operation, boolean = False):
    # Reversing operation
    if operation == "rc":
        # When boolean is True, the DNA sequence has been already processed. 
        # As qualities don't need a complementary chain, by introducing as well the boolean as argument they only will be reversed when boolean = TRUE
        if boolean == True: 
            return f_line[::-1], True
        # This returns complementary and reversed chain. So it is only though for DNA sequences.
        else:
            return "".join(nucleotides[nt] for nt in f_line[::-1]), True
    # Trimming sequence 
    elif operation == "trim":
        return trim_sequence(f_line, boolean = boolean)
    # Adaptor removal
    elif operation == "adaptor-removal":
        return adaptor_rem(sequence = f_line, adaptor = args.adaptor, boolean = boolean)
    # Warning error for operations mistyped or don't defined in the program
    else:
        print("Operation '%s' not recognized" % operation)
        exit(1)

        
def print_summary(): 
    # Translate operation
    operation = args.operation
    if operation == "rc": operation = "reversed-complemented" 
    elif operation == "trim": operation = "hard-trimmed" 
    else: operation = "processed " 
    # Print summary
    print("File '%s' has been successfully %s ('%s')" % (args.input, operation, args.output))
    print("Summary:")
    print("\t%d reads processed" % stats["reads"])
    # The stats for each base is converted into a percentage. Numbers need to be rounded in order to get percentages without decimals. In some cases, the sum of percentages won't be 100%, but a close value
    print("\t%d bases processed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (stats["bases"],round(stats["A"]/stats["bases"]*100,0),round(stats["C"]/stats["bases"]*100,0),round(stats["G"]/stats["bases"]*100,0),round(stats["T"]/stats["bases"]*100,0),round(stats["N"]/stats["bases"]*100,0))) 
    if operation == "hard-trimmed":
        print("\t%d bases trimmed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (trimmed_stats["bases"],round(trimmed_stats["A"]/trimmed_stats["bases"]*100,0),round(trimmed_stats["C"]/trimmed_stats["bases"]*100,0),round(trimmed_stats["G"]/trimmed_stats["bases"]*100,0),round(trimmed_stats["T"]/trimmed_stats["bases"]*100,0),round(trimmed_stats["N"]/trimmed_stats["bases"]*100,0)))
    elif operation == "processed ":
        print("\t%d adaptors-found" % (adaptor_stats["adaptors-found"]))


# Last step is a conditional statement to call the whole process, just in case additional arguments are also inputted if needed.
# Otherwise, warning messages appear on the command line.
if args.operation == "adaptor-removal" and not args.adaptor:
        print("An adaptor is required")
elif args.operation == "trim" and args.trim_right == None and args.trim_left == None:
        print("A number of bases to be trimmed is required")
else:
	process_file()
