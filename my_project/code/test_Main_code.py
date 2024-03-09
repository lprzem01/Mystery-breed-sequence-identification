import pytest
from Bio import SeqIO
import os

import Main_code


#complete these variables with the input files directory
dog_breeds = r"/workspaces/Coursework/my_project/data/dog_breeds.fa"
mystery_breed = r"/workspaces/Coursework/my_project/data/mystery.fa"
output = r"/workspaces/Coursework/my_project/results"
ind_breeds = r"/workspaces/Coursework/my_project/results/individual_breed_sequences"
tests = r"/workspaces/Coursework/my_project/tests"


def test_read_fasta():
    """check if all the breeds are capitalised"""
    all_breeds = set(Main_code.Breed.breed())
    all_breeds_same_type = set((Main_code.Breed.breed()).capitalize())
    assert all_breeds == all_breeds_same_type


def test_percentage_similarity():
    """check if an alignment of two identical sequences returns 100% similarity as expected"""
    identical_sequence_alignment = Main_code.simple_alignment(Main_code.mystery_sequence,Main_code.mystery_sequence)
    assert Main_code.percentage_similarity(identical_sequence_alignment) == 100

#create a test txt
def test_create_output_txt():
    """use a test file with a known txt content to check if the expected content is written to the file"""
    test_content = "test"
    test_file =  f"{tests}/test_txt"
    test_output = Main_code.create_output(test_content, test_file, "txt")
    with open(test_file) as f:
         read = f.readlines()
    assert read == test_content
#create a test fasta

def test_unique_breeds():
    all_breeds = Main_code.all_breeds()
    assert "BOXER" in all_breeds

def test_breed_sequences():
    """check if an example file contains only one breed type in the name"""
    example_file = f"{ind_breeds}/BOXER"
    record_names = []
    for record in SeqIO.parse(example_file, "fasta"):
        record_names.append(record.name)
    assert len(set(record_names)) == 1

def test_concensus_seq():
    """check if a file containing multiple instances of the same sequence has a concensus_seq output exactly the same"""
    test_sequence = SeqIO.read("./my_project/tests/test_file.fasta", "fasta")
    test_concensus = Main_code.concensus_seq("./my_project/tests/consensus_test_file.fasta")
    assert test_sequence == test_concensus