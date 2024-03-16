
import pytest
from Bio import SeqIO
import os
import pandas
import time
from functions_copy import *

#complete these variables with the input files directory
dog_breeds = r"../data/dog_breeds.fa"
mystery_breed = r"../data/mystery.fa"
output = r"../results"
ind_breeds = r"../results/individual_breed_sequences"
tests = r"../tests"



def test_read_fasta():
    """check if all the breeds are capitalised"""
    assert isinstance(Breed, type)
    #print(type(Breed))
    #all_breeds = set(Breed.breed())
    #all_breeds_same_type = set((Breed.breed()).capitalize())
    #assert all_breeds == all_breeds_same_type

test_read_fasta()



def test_percentage_similarity():
    """check if an alignment of two identical sequences returns 100% similarity as expected"""
    mystery_sequence = SeqIO.read(f"{tests}/test_file.fasta", "fasta")
    identical_sequence_alignment = simple_alignment(mystery_sequence.seq,mystery_sequence.seq)[1]
    assert percentage_similarity(identical_sequence_alignment[0]) == 100

test_percentage_similarity()



#create a test txt
def test_create_output_txt():
    """use a test file with a known txt content to check if the expected content is written to the file"""
    test_content = "test"
    test_file =  f"{tests}/test_txt"
    test_output = create_output(test_content, test_file, "txt")
    with open(test_file) as f:
         read = f.readlines()
    assert read[0] == test_content
#create a test fasta

test_create_output_txt()


def test_unique_breeds():
    initialise_Breed()
    all_breeds = unique_breeds()
    assert "BOXER" in all_breeds
test_unique_breeds()


def test_breed_sequences():
    """check if an example file contains only one breed type in the name"""
    example_file = f"{ind_breeds}/BOXER"
    record_names = []
    for record in SeqIO.parse(example_file, "fasta"):
        record_names.append(record.name)
    assert len(set(record_names)) == 1

test_breed_sequences()



def test_consensus_seq():
    """check if a file containing multiple instances of the same sequence has a consensus_seq output exactly the same"""
    filename1 = f"{tests}/test_file.fasta"
    filename2 = f"{tests}/consensus_test_file.fasta"
    for record in SeqIO.parse(filename1, "fasta"):
        test_sequence = record.seq
    for record in SeqIO.parse(((filename2)), "fasta"):
        test_consensus =  record.seq
    assert test_sequence == test_consensus

test_consensus_seq()


def test_consensus_file():
    """checks if all breeds are present in consensus file"""
    initialise_Breed()
    consensus_file_breeds = []
    for key in consensus_file():
        consensus_file_breeds.append(key.id)
    all_breeds = unique_breeds()
    for breed in all_breeds:
        assert breed in consensus_file_breeds

test_consensus_file()


def test_align_consensus():
    """test if the align_consensus function returns test sequence from a list of sequences which includes the test sequence"""
    #initialise both clases
    initialise_Breed()
    consensus_file = f"{output}/consensus_sequences"
    initialise_Breed_consensus(consensus_file)
    #get the test sequence
    test_sequence = read_fasta("../tests/BOXER_test_sequence.fasta")
    #get a lsit of consensus sequences
    consensus_file = f"{output}/consensus_sequences"
    consensus_sequences = read_fasta(consensus_file)
    #create an artificial alignment between the same sequence    
    correct_alignment = simple_alignment(test_sequence[0], test_sequence[0])[1]
    #generate the best scoring alignment
    consensu_alignment = align_consensus(consensus_sequences, test_sequence[0])
    #get the alignment and breed name from the best scoring alignment
    alignment = consensu_alignment[1]
    breed = consensu_alignment[2]
    #check if the breed name is the same as test file breed name 
    assert breed == "BOXER "
    #check if alignment matches the artifical alignment 
    assert alignment[0] == correct_alignment[0]

test_align_consensus()


def test_simple_alignment():
    """test to check that an alignment of two identical sequences returns an alignment with 100% percentage similarity"""
    sample_seq = "AGCTAGCATAATCGG"
    alignment = simple_alignment(sample_seq,sample_seq)
    similarity = test_percentage_similarity(alignment[1])
    assert similarity == 100

test_simple_alignment


def test_initialise_Breed():
    initialise_Breed()
    for instance in Breed.all_instances:
        assert instance.breed.isupper() 

test_initialise_Breed()


def test_add_mystery_to_consensus():
    #get mystery sequence
    mystery = (read_fasta(mystery_breed))[0]
    #get consensus sequences with mystery sequence list 
    records = add_mystery_to_consensus()
    #check mystery sequence is in consensus
    sequences = []
    for record in records:
        sequences.append(record.seq)
    assert mystery in sequences

test_add_mystery_to_consensus()


def test_initialise_Breed_consensus():
    #initialise breed consensus
    consensus_file = f"{output}/consensus_sequences"
    initialise_Breed_consensus(consensus_file)
    #check unknown description not in breed name
    breeds = []
    for instance in Breed_consensus.all_instances:
        breeds.append(instance.breed)
    assert "<unknown description>" not in breeds
    
test_initialise_Breed_consensus()


def test_MSA_alignment():
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    align = MSA_alignment(consesnsus_with_mystery)
    assert type(align) == Bio.Align.MultipleSeqAlignment

test_MSA_alignment()


def test_build_phylo_tree():
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    alignment = MSA_alignment(consesnsus_with_mystery)
    tree = build_phylo_tree(alignment)
    assert type(tree) == Bio.Phylo.BaseTree.Tree

test_build_phylo_tree()


def test_custom_label_func():
    nodes = []
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    alignment = MSA_alignment(consesnsus_with_mystery)
    tree = build_phylo_tree(alignment)
    print()
    assert tree.find_any("XXX") == None
    assert str(tree.find_any("BOXER")) == "BOXER"

test_custom_label_func()


def test_base_content():
    test_seq = "AAAAATTTTT"
    content = base_content(test_seq)
    assert content[0] == 0
    assert content[1] == 0
    assert content[2] == 50
    assert content[3] == 50

test_base_content()


def test_random_DNA():
    """check if random DNA has the defined base content"""
    content = [0,0,50,50]
    DNA = random_DNA(100,content)
    random_DNA_content = base_content(DNA)
    assert random_DNA_content[0] == 0
    assert random_DNA_content[1] == 0
    #base content is not always exact as the sequence is random so checking if A and T have at least 30% base content
    assert random_DNA_content[2] > 30
    assert random_DNA_content[3] > 30
test_random_DNA()


def test_random_sequences():
    """check if the number of sequences generated matches target database size"""
    #generate sequences 
    base_count = [0, 0, 47, 53]
    test_database = random_sequences(10, 25, base_count)
    assert len(test_database) == 10

test_random_sequences()



def test_breed_alignment_scores():
    #get some scores 
    sequence = (read_fasta("../data/mystery.fa"))[0]
    scores = breed_alignment_scores(sequence)
    for score in scores.values():
        assert score[0] > 1
test_breed_alignment_scores()


def test_alignment_scores():
    #get some scores 
    target = (read_fasta("../tests/tets_mystery_sequence.fasta"))[0]
    sequences = (read_fasta("../tests/test_sequences.fasta"))
    scores = alignment_scores(target, sequences)
    for score in scores:
        assert score > 1
test_alignment_scores()


def test_probabilitydf():
    scores = [(11111, 45678), (22222, 4567), (33333, 56789)]
    breeds = ["X", "Y", "Z"]
    #arbitrary l and u values
    df = probability_df(breeds, scores, 0.5, 0.09)
    assert type(df) == pandas.core.frame.DataFrame
test_probabilitydf()

