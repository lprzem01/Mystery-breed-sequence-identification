from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
import Bio.Align
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import random 
import pandas as pd
import numpy as np 
import pytest

dog_breeds = r"../data/dog_breeds.fa"
mystery_breed = r"../data/mystery.fa"
output = r"../results" 
ind_breeds = r"../results/individual_breed_sequences"

def read_fasta(filename):
    """function which reads a fasta file and returns a list of sequences found in the file"""
    sequences = []
    #read the file
    for record in SeqIO.parse(filename, "fasta"):
        #add each sequence to the list
        sequences.append(record.seq)
    return(sequences)

def create_output(content, filename:str, filetype:str):
        """creates a file in the results folder with the content provided in the the 
        correct format (either txt or fasta)"""
        #create an empty file by openening it in a write format 
        filepath = f"{output}/{filename}"
        with open(filename, "w") as f:
                #if the filetype is fasta crete a file using SeQIO
                if filetype == "fasta":
                        SeqIO.write(content, filename, filetype)
                #otherwise if the file is a txt write the content to the file 
                elif filetype == "txt":
                        f.write(content)
                #if neither is true return an error message         
                else:
                    return "Provide a valid file type, this function only creates fasta files and txt files."
                
def simple_alignment(seq1,seq2):
    """align two sequences given as parameters and return alignment score and alignment"""
    #set up aligner
    aligner = Bio.Align.PairwiseAligner()
    #get alignemnt score 
    score = aligner.score(seq1, seq2)
    #get alignemnt
    alignment = aligner.align(seq1, seq2)
    return score, alignment

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

def initialise_Breed(filename, format = "fasta"):
    #parse through the dog_breeds file
    for record in SeqIO.parse(filename, format):
        #get the description of each sequence to find out what breed it is 
        #split the description on "[" as these are used to separate each descriptor
        for key in list(record.description.split("[")): 
            #if the keyword breed is in the description but the keyword "isolate" is not in thr description it will define the dog breed 
            if "breed" in key and "isolate" not in key:
                #each sequence in filename gets assigned a name based on the breed identified in the desciption
                breed_name = ((key[6:-2]).upper()).replace(" ", "_")
                #intialises an object of class Breed which contains the sequence, breed name and the full record 
                record.name = Breed(record.seq, breed_name, record)

def unique_breeds():
    """Set up a list that contains all unique breed names from the class Breed in a list"""
    all_breeds = set()
    #itterate through all instances in Breed class
    for key in Breed.all_instances:
        #access and add every breed to the all_breeds set
        all_breeds.add(key.breed)
    return all_breeds

def breed_sequences(directory):
    """Write a fasta file containing all sequences that belong to the same breed """
    for breed in unique_breeds():
        #creates a temporary variable corresponding to each individual dog breed 
        temp = breed
        #create a directory for filename
        filename = f"{directory}/{breed}"
        #itterate throough all instances of the class Breed
        sequences = []
        for key in Breed.all_instances: 
            #check if the breedd is the same as the current breed in the loop stored in the temp variable 
            if key.breed == temp:
                sequences.append(key.fasta)
        #adds the sequences to a file with the name of the breed as a filename
        SeqIO.write(sequences, filename, "fasta") 

def consensus_seq(filename):
    """Function that takes in a file containing a number of sequences and returns a 
    consensus sequence"""
    #align all the sequences in each file
    alignments = AlignIO.parse(f"{ind_breeds}/{filename}", "fasta") 
    for alignment in alignments: 
        #get summary info of each alignment to create a consensus 
        summary = SummaryInfo(alignment) 
        #create a consensus of each alignment 
        consensus = summary.dumb_consensus() 
        #create a fasta format sequence using the consensus sequence and recordname
        seq_record = SeqRecord(Seq(consensus), id=filename) 
        #add each consensus seq to a list 
    return seq_record

def consensus_file():
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
def initialise_Breed_consensus(filename, format = "fasta"):
    #parse through the dog_breeds file
    for record in SeqIO.parse(filename, format):
            #intialises an object of class Breed which contains the sequence, breed name and the full record 
            breed_name = record.description.replace("<unknown description>", "")
            record.id = Breed_consensus(record.seq, breed_name, record)

def align_consensus(sequences_list, unknown_sequence):
    """uses the simple alignment function to cretae an alignment between unknown 
    sequence and database. Returns the score and alignment of the top scoring alignment."""
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
    return current_best, current_best_alignment, top_breed

def percentage_similarity(aln):
    """Given a pairwise alignment calculates the percentage similarity between the 
    two sequences"""
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

def MSA_alignment(all_sequences):
    """Create a multisequence alignment between unknownn dna and all consensus sequence"""
    #read the file with MSA alignment 
    sequences = SeqIO.parse(all_sequences, "fasta")
    # Create a MultipleSeqAlignment object
    MSA_alignment = Bio.Align.MultipleSeqAlignment(sequences)
    #create a directory and filename for the MSA alignment
    filename = f"{output}/MSA_alignment"  
    #save the multiple sequence alignment in a clustal file 
    AlignIO.write(MSA_alignment, filename, "clustal")
    return MSA_alignment

def build_phylo_tree(MSA_alignment):
    """Create a distance matrix from the MSA alignment"""
    # Calculate the distance matrix
    calculator = DistanceCalculator("identity")
    # Build the tree using the neighbor-joining method
    constructor = DistanceTreeConstructor(calculator, method="nj")
    breeds_tree = constructor.build_tree(MSA_alignment)
    # Save the tree to a new file 
    Phylo.write(breeds_tree, f"{output}/breeds_tree.xml", "phyloxml")
    return breeds_tree

def custom_label_func(node):
    """Creates better phylogenetic tree visualisation by removing unwanted labels 
    and label characters"""
    if node.is_terminal():
        #remove the '' marks in names
        return (node.name).replace("'", "")
    #delete non terminal node labels 
    else:
        return None
    
def base_content(mystery_sequence):
    """generate the proportion of each base in a DNA sequence"""
    bases = ["G", "C", "T", "A"]
    bases_count = {}
    total_count = len(mystery_sequence)
    #itterate through each base of the sequence to find the percentage intiger of each base 
    for base in bases:
        bases_count[base] = int(mystery_sequence.count(base) / total_count * 100)
    return list(bases_count.values()) 

def random_DNA(length, base_count):
    """Generate random sequences with the same proportion of base content as the 
    real sequence"""
    bases = ["G", "C", "T", "A"]
    #generate the DNA in form of a list of bases
    random_DNA = random.choices(bases, weights=base_count, k=length)
    #output the sequence with bases joied together
    return (''.join(random_DNA) )

def random_sequences(database_size, seq_len, base_count):
    """generates a list of sequences of the size and length of the original database"""
    rand_seq_database = []
    #reapeats the random DNA generation process to miror the database length
    for number in range(0, database_size):
        rand_seq = random_DNA(seq_len, base_count)
        #add each generated sequence to a list 
        rand_seq_database.append(rand_seq)
    return rand_seq_database

def alignment_scores(target_seq, sequences):
    """Generates alignment scores for all the sequences found in the random 
    sequences list"""
    scores = []
    #itterates through the seequences in a list and generates an alignment score
    for sequence in sequences:
        alignment = simple_alignment(target_seq, sequence)
        #alignment score is added to a list
        scores.append(alignment[0])
    #list of alignment scores is returned
    return scores

def breed_alignment_scores(unknown_sequence):
    """uses the simple alignment function to create an alignment between unknown 
    sequence and consenus sequences. Returns a list of breeds with the corresponding 
    alignment score."""
    #set up a list for all alignment score 
    scores = {}
    #itterate through the sequences in the database
    for instance in Breed_consensus.all_instances:
        #determine the score of each alignment
        score = (simple_alignment(instance.sequence, unknown_sequence)[0])
        #add score and sequence to scores dictionary
        scores[instance.breed] = [score, len(instance.sequence), len(unknown_sequence[0])]
    return scores

# The probability density function of the extreme value distribution
def pdf(x, l, u):  # noqa: E741
    t = np.exp(-l * (x - u))
    return l * t * np.exp(-t)

def probability_df(breeds, scores, l, u, df_name):  # noqa: E741
    """Function takes in a list of breeds and a list of alignmet scores acompanying 
    the breed as well as l and u constants used to calculate pdf and returns a 
    dataframe of breed names and alignment probabilities."""
    probabilities = []
    #itterate through alignemnt scores of the actual dataset 
    for breed, score in zip(breeds, scores):
        #get the breed name and probability of alignment occuring by chance
        line = breed, pdf(score[0], l, u)
        #add the breed and probability to a list
        probabilities.append(line)
    #convert the list of breeds and probabilities to a dataframe
    df = pd.DataFrame(probabilities, columns=["breed", "probability"])
    #sor the dataframe so that the least likelly to occur by chance alignment is at the top 
    df = df.sort_values("probability")
    #save the dataframe to a file  
    df.to_csv(f'{output}/{df_name}')
    #funtion returns the dataframe create for testing and further modificaations purpouses
    return df

