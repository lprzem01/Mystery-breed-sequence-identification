import pytest
from Bio import SeqIO
import os
from Main_code import *
import Main_code

#complete these variables with the input files directory
dog_breeds = r"C:\Users\User\Downloads\Coursework-main\Coursework-main\my_project\data\dog_breeds.fa"
mystery_breed = r"C:\Users\User\Downloads\Coursework-main\Coursework-main\my_project\data\mystery.fa"
output = r"C:\Users\User\Downloads\Coursework-main\Coursework-main\my_project\results"
ind_breeds = r"C:\Users\User\Downloads\Coursework-main\Coursework-main\my_project\results\individual_breed_sequences"
tests = r"C:\Users\User\Downloads\Coursework-main\Coursework-main\my_project\tests"



def test_read_fasta_test():
    """check if all the breeds are capitalised"""
    assert isinstance(Breed, type)
    #print(type(Breed))
    #all_breeds = set(Breed.breed())
    #all_breeds_same_type = set((Breed.breed()).capitalize())
    #assert all_breeds == all_breeds_same_type

test_read_fasta_test()


def test_percentage_similarity():
    """check if an alignment of two identical sequences returns 100% similarity as expected"""
    identical_sequence_alignment = simple_alignment(mystery_sequence[0],mystery_sequence[0])[1]
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
    print(test_sequence)
    for record in SeqIO.parse((consensus_seq(filename2)), "fasta"):
        test_consensus =  record.seq
    print(test_consensus)
    assert test_sequence == test_consensus

test_consensus_seq()

def test_consensus_file():
    """checks if all breeds are present in consensus file"""
    all_breeds = unique_breeds()
    for breed in all_breeds:
        assert breed in consensus_file()

test_consensus_file()

def test_align_consensus():
    """test the alignment of two identical sequences returns the correct output"""
    for key in Breed.all_instnces:
        if key.breed == "BOXER":
            test_sequence = key.sequence
    sequences_list = test_sequence
    unknown_sequence = test_sequence
    correct_alignment = simple_alignment(sequences_list,unknown_sequence)
    alignment = (align_consensus(sequences_list, unknown_sequence))[1]
    breed = (align_consensus(sequences_list, unknown_sequence))[2]
    assert breed == "BOXER"
    assert alignment == correct_alignment

test_align_consensus()


def test_runtime():
    """check the program runs for less than 10 min in total"""
    start_time = time.time()
    Main_code()
    total_time = ("--- %s seconds ---" % (time.time() - start_time)) 
    assert total_time  < 10*60
