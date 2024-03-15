# %% [markdown]
# import all relevant modules

# %%
from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
import Bio.Align
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
import os
import matplotlib.pyplot as plt
import numpy as np 
import scipy as sp
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from pymsaviz import MsaViz
import time
from biotite import sequence as seq


# %% [markdown]
# Setup the directory

# %%
#set current directory to my_project

#complete these variables with the input files directory
dog_breeds = r"/workspaces/Coursework/my_project/data/dog_breeds.fa"
mystery_breed = r"/workspaces/Coursework/my_project/data/mystery.fa"
output = r"/workspaces/Coursework/my_project/results"
ind_breeds = r"/workspaces/Coursework/my_project/results/individual_breed_sequences"
all_alignments = r"/workspaces/Coursework/my_project/results/All_alignments"


# %% [markdown]
# Set up basic functions 

# %%
def read_fasta(filename):
    """function which reads a fasta file and returns a list of sequiences found in the file"""
    sequences = []
    for record in SeqIO.parse(filename, "fasta"): 
        sequences.append(record.seq)
    return(sequences)

def create_output(content, filename:str, filetype:str):
        """creates a file in the results folder with the content provided in the the correct format """
        #create an empty file by openening it in a write format 
        filepath = f"{output}/{filename}"
        with open(filename, "w") as f:
                if filetype == "fasta":
                        SeqIO.write(content, filepath, filetype) 
                elif filetype == "txt":
                        f.write(content)
def simple_alignment(seq1,seq2):
    """align two sequences given as parameters and return alignment score and alignment"""
    aligner = Bio.Align.PairwiseAligner() 
    score = aligner.score(seq1, seq2)
    alignment = aligner.align(seq1, seq2)
    return score, alignment



# %% [markdown]
# Set up the Breed class 

# %%
class Breed():
    """This class stores information about every breed and their sequence"""
    #create a list to store all instances of the class
    all_instances = []
    #define an init function which stores the sequence as a sequence object, the breed as a string and the fasta format of the sequence
    def __init__(self, sequence, breed, fasta):
        self.sequence = sequence
        self.breed = breed
        self.fasta = fasta
        #store all initialised instances in the defined list
        Breed.all_instances.append(self) 
#open the fasta file and save the sequence breed and sequence name in a Breed Class to be accessed later
def initialise_Breed(filename = dog_breeds, format = "fasta"):
    #parse through the dog_breeds file
    for record in SeqIO.parse(filename, format):
        #get the description of each sequence to find out what breed it is 
        #split the description on "[" as these are used to separate each descriptor
        for key in list(record.description.split("[")): 
            #if the keyword breed is in the description but the keyword "isolate" is not in thr description it will define the dog breed 
            if "breed" in key and "isolate" not in key:
                #each sequence in filename gets assigned a name based on the breed identified in the desciption
                breed_name = (key[6:-2]).upper()
                #intialises an object of class Breed which contains the sequence, breed name and the full record 
                record.name = Breed(record.seq, breed_name, record)
#call the function
initialise_Breed()

# %% [markdown]
# Create consensus files for each breed

# %%
def unique_breeds():
    """set up a list that contains all unique breed names in a list"""
    all_breeds = set()
    #itterate through all instances in Breed class
    for key in Breed.all_instances:
        #access and add every breed to the all_breeds set
        all_breeds.add(key.breed)
    return all_breeds
    

def breed_sequences(directory = ind_breeds):
    """write a fasta file containing all sequences that belong to the same breed """
    for breed in unique_breeds():
        #creates a temporary variable corresponding to each individual dog breed 
        temp = breed
        #create a directory for filename
        filename = f"{directory}/{breed}.fa"
        #itterate throough all instances of the class Breed
        sequences = []
        for key in Breed.all_instances: 
            #check if the breedd is the same as the current breed in the loop stored in the temp variable 
            if key.breed == temp:
                sequences.append(key.fasta)
                #adds the sequences to a file with the name of the breed as a filename
        SeqIO.write(sequences, filename, "fasta") 

breed_sequences()
def consensus_seq(filename):
    """Function that takes in a file containing a number of sequences and returns a consensus sequence"""
    os.chdir(ind_breeds) 
    #align all the sequences in each file
    alignments = AlignIO.parse(filename, "fasta") 
    #assign filename to varible recordname which will be used to create a name for this record 
    recordname = f"{filename}"
    for alignment in alignments: 
        #get summary info of each alignment to create a consensus 
        summary = SummaryInfo(alignment) 
        #create a consensus of each alignment 
        consensus = summary.dumb_consensus() 
        #create a fasta format sequence using the consensus sequence and recordname
        seq_record = SeqRecord(Seq(consensus), id=recordname) 
        #add each consensus seq to a list 
    return seq_record


def consensus_file(directory=output):
    """create a consensus sequence for each breed and store it in a list"""
    #create a list to store the sequences
    consensus_sequences = [] 
    for file in unique_breeds():
        #run the consensus_seq function to get the consensus file of each breed in the consensus list 
        consensus_sequences.append(consensus_seq(file))
    return consensus_sequences
        
def add_mystery_to_consensus():
    consensus_sequences = list(consensus_file())
    unknown_sequence = (read_fasta(mystery_breed))[0]
    #add mystery sequence to the consensus file list 
    consensus_sequences.append(SeqRecord(unknown_sequence, id="MYSTERY SEQUENCE"))
    return consensus_sequences

#store consensus sequences in a results folder in a file called consensus_sequences
create_output(consensus_file(), "consensus_sequences", "fasta" )
#create output with the mystery sequence
create_output(add_mystery_to_consensus(), "consensus_sequences_with_mystery", "fasta" )
#hold the list as a value for easy acess later on 
consensus_sequences_with_unknown = add_mystery_to_consensus()

# %% [markdown]
# Create a class to store consensus sequences 

# %%
#get consensus sequences
consensus_file = f"{output}/consensus_sequences"
consensus_sequences = read_fasta(consensus_file)
#get mystery sequences
mystery_sequence = read_fasta(mystery_breed)

class Breed_consensus():
    """This class stores information about every breed and their sequence"""
    #create a list to store all instances of the class
    all_instances = []
    #define an init function which stores the sequence as a sequence object, the breed as a string and the fasta format of the sequence
    def __init__(self, sequence, breed, fasta):
        self.sequence = sequence
        self.breed = breed
        self.fasta = fasta
        #store all initialised instances in the defined list
        Breed_consensus.all_instances.append(self) 
#open the fasta file and save the sequence breed and sequence name in a Breed Class to be accessed later
def initialise_Breed_consensus(filename = consensus_file, format = "fasta"):
    #parse through the dog_breeds file
    for record in SeqIO.parse(filename, format):
            #intialises an object of class Breed which contains the sequence, breed name and the full record 
            breed_name = record.description.replace("<unknown description>", "")
            record.id = Breed_consensus(record.seq, breed_name, record)
#call the function
initialise_Breed_consensus()


# %% [markdown]
# Find the top alignment and store it in a file (this function is the one which takes the longest)

# %%
#align_consensus sequences with the mystery sequence
def align_consensus(sequences_list=consensus_sequences, unknown_sequence = mystery_sequence[0]):
    """uses the simple alignment function to cretae an alignment between unknown sequence and database. Returns the score and alignment of the top scoring alignment."""
    #set up a list for all alignment score 
    current_best = 0
    current_best_alignment = 0
    current_sequence = 0
    top_breed = ""
    #itterate through the sequences in the database
    for sequence in sequences_list:
        #determine the score of each alignment
        score = (simple_alignment(sequence, unknown_sequence)[0])
        if score > current_best:
            current_best = score
            alignment = (simple_alignment(sequence, unknown_sequence)[1])
            current_sequence = sequence
            current_best_alignment = alignment
    for key in Breed_consensus.all_instances:
        if key.sequence == current_sequence:
            top_breed = key.breed
    print(current_best, current_best_alignment, top_breed)
    return current_best, current_best_alignment, top_breed


top_alignment_details = align_consensus()
top_alignment = top_alignment_details[1][0]

#create a directory for the top alignment file 
filename = f"{output}/top_alignment_output"  
#create a file to store the top scoring alignment as clustal file 
with open(filename, "w"):
        Bio.Align.write(top_alignment, filename, "clustal")

# %% [markdown]
# Calculate percentage similarity of the top scoring sequence 

# %%
def percentage_similarity(aln):
    """Given a pairwise alignment calculates the percentage similarity between the two sequences"""
    #create a varieble to store the instances where columns are identical
    identical_columns = float()
    #itterate through every column of the alignment
    for a in range(len(aln[0])): 
        #check if first and second sequence is the same at point a, a describing the column 
        if aln[0,a] == aln[1,a]: 
            #if base at position a is the same in both sequences add 1 to the amount of identical columns 
            identical_columns += 1
    #calculate the percentage based on the identical columns number and the length of the alignment 
    percentage = 100 * identical_columns / float(len(aln[0])) 
    #return percentage to 3dp
    return  round(percentage,3) 

#calculate percentage similarity of the top scoring sequence 
top_percentage = percentage_similarity(top_alignment)

# %% [markdown]
# Save the results to a text file containing an explanation of the values generated 
# -create a pymasviz grpah to better visualise the alignment, stored in results 

# %%
#get results aka the breed, its sequence and percentage similarity 
results = "The breed most similar to the mystery DNA file is the", top_alignment_details[2], "its percent identity is", top_percentage, "%" ,"the alignment of the mystery dog breed and", top_alignment_details[2], "is displayed here\n", top_alignment

#create a string representing thr results that can be written to the results file 
results_str = str()
for key in results:
    results_str += str(key) 

#create a directory and filename for the details of the top alignment
filename2 = f"{output}/top_alignment_details"  
#create a file to store details about the  top alignment as a txt file 
with open(filename2, "w") as file:
        file.write(results_str)

#create a graph for each 100 bases long sequence alignment figure 

read_alignment = AlignIO.read(f"{output}/top_alignment_output", "clustal")  

l = len(top_alignment[0]) 
mv = MsaViz(read_alignment, format="clustal", start=1, end=l, wrap_length=100)
mv.savefig(f"{output}/top_alignment_image")

# %% [markdown]
# Create a phylogenic tree describing the relationship between the breeds 
# - organise this into functions 

# %%
#create a multisequence alignment between unknownn dna and all consensus sequences
# Create a MultipleSeqAlignment object
MSA_alignment = Bio.Align.MultipleSeqAlignment(consensus_sequences_with_unknown)
#create a directory and filename for the MSA alignment
filename3 = f"{output}/MSA_alignment"  
#save the multiple sequence alignment in a clustal file 
AlignIO.write(MSA_alignment, filename3, "clustal")


# Calculate the distance matrix
calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(MSA_alignment)

# Build the tree using the neighbor-joining method
constructor = DistanceTreeConstructor(calculator, method="nj")
breeds_tree = constructor.build_tree(MSA_alignment)
# Save the tree to a new file 

Phylo.write(breeds_tree, "breeds_tree.xml", "phyloxml")
# Convert the tree to a different format
Phylo.convert("breeds_tree.xml", "phyloxml", "breeds_tree.nex", "nexus")

breeds_nex = Phylo.read("breeds_tree.nex", "nexus")
breeds_nex.rooted = True
# Create a custom label function that returns None for inner clade labels
def custom_label_func(node):
    if node.is_terminal():
        return (node.name).replace("'", "")
    else:
        return None
    
import matplotlib
fig = plt.figure(figsize=(20,15), dpi=100, frameon=False)
matplotlib.rc("font", size=12)
ax = plt.gca()
Phylo.draw(breeds_nex, show_confidence=True, axes=ax, label_func=custom_label_func)
plt.savefig(f"{output}/Phylogenetic_tree")

# %% [markdown]
# #find regions of the most diversity 
# #define diversity as more than 10 sequences differing in that region 

# %%
#make a plot to represent the conservation levels of the MSA

MSA_alignment = AlignIO.read(f"{output}/MSA_alignment", "clustal")  

l = len(MSA_alignment[0]) 
mv = MsaViz(MSA_alignment, format="clustal", start=50, end=150, wrap_length=100, show_consensus=True)
mv.savefig(f"{output}/MSA_alignment_image")


# %% [markdown]
# Probability 

# %%
#convert X to " "
def ambiguous_letter_deltion():
    """this function removes all the instances of unknown nucleotides from the sequences in Breed_consensus class. Fasta format remains the same"""
    for instance in  Breed_consensus.all_instances:
            if instance == "AIDI":
                  print(instance.sequence)
            #sequence = (instance.sequence).replace("X", "")
            #instance.self = Breed_consensus(sequence, instance.breed, instance.fasta)

ambiguous_letter_deltion()



# %%

def alignment_scores(unknown_sequence = mystery_sequence[0]):
    """uses the simple alignment function to create an alignment between unknown sequence and consenus sequences. Returns a list of breeds with the corresponding alignment score."""
    #set up a list for all alignment score 
    scores = {}
    #itterate through the sequences in the database
    for instance in Breed_consensus.all_instances:
        #determine the score of each alignment
        score = (simple_alignment(instance.sequence, unknown_sequence)[0])
        scores[instance.breed] = score
    return scores

breed_scores = alignment_scores()

def compute_e_value(scores = breed_scores):
    """this function takes in a dictionary of alignment scores and breed names to return a dictionary of breed names and E_values in decadic logarithm form"""
    e_values = {}
    for breed, score in zip(scores.keys(), scores.values()):
        e_value = biotite.log_evalue(score, l, l)
        e_values[breed] = e_value
    return e_values

print(compute_e_value())
#log_evalue(score, seq1_length, seq2_length)

# %% [markdown]
# Tests 

# %%

query = biotite.sequence.NucleotideSequence("CGACGGCGTCTACGAGTCAACATCATTC")
hit = biotite.sequence.NucleotideSequence("GCTTTATTACGGGTTTACGAGTTCAACATCACGAAAACAA")
example = biotite.sequence.io.fasta.get_sequence(mystery_breed)


matrix = biotite.sequence.align.SubstitutionMatrix.std_nucleotide_matrix()
gap_penalty = (-12, -2)
alignment = biotite.sequence.align.align_optimal(query, hit, matrix, gap_penalty, local=True)[0]
print(alignment)
print(alignment.score)

# Ensure deterministic results
np.random.seed(0)
# Sequences in database have a GC content of 0.6
background = np.array([0.2, 0.3, 0.3, 0.2])
estimator = biotite.sequence.align.EValueEstimator.from_samples(example.alphabet, matrix, gap_penalty, background, sample_length=100)


log_e = estimator.log_evalue(alignment.score, len(query), 100 * len(hit))
print(f"E-value = {10**log_e:.2e}")




# %% [markdown]
# 
