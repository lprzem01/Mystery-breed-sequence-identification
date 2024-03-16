
from Bio import SeqIO
import pandas

from functions_copy import *  # noqa: F403

#complete these variables with the input files directory
dog_breeds = r"../data/dog_breeds.fa"
mystery_breed = r"../data/mystery.fa"
output = r"../results"
ind_breeds = r"../results/individual_breed_sequences"


def test_read_fasta():
    """check if the sequences are returned as sequence objects"""
    test_read = read_fasta(mystery_breed)
    assert type(test_read[0]) == Bio.Seq.Seq


def test_percentage_similarity():
    """check if an alignment of two identical sequences returns 100% similarity as expected"""
    #load test sequence 
    mystery_sequence = SeqIO.read(r"test_file.fasta", "fasta")
    #create an alignment between the test sequence and itself
    identical_sequence_alignment = simple_alignment(mystery_sequence.seq,mystery_sequence.seq)[1]  
    #check if the alignment has 100% similarity
    assert percentage_similarity(identical_sequence_alignment[0]) == 100  

#create a test txt
def test_create_output_txt():
    """use a test file with a known txt content to check if the expected content is written to the file"""
    #create arbitrary test content
    test_content = "test"
    #create test content filename
    test_file =  "test_txt"
    #create the file using create_output function
    create_output(test_content, test_file, "txt")  
    with open(test_file) as f:
         read = f.readlines()
    #read the file and assert it is the same as the content provided
    assert read[0] == test_content

def test_unique_breeds():
    """tests if expected breed is found in all_breeds"""
    #initialise the Breeds class
    initialise_Breed(dog_breeds)  
    all_breeds = unique_breeds()  
    #check BOXER is found in the output of calling unique_breeds
    assert "BOXER" in all_breeds

def test_breed_sequences():
    """check if an example file contains only one breed type in the name"""
    #find example file
    example_file = f"{ind_breeds}/BOXER"
    record_names = []
    #find record names
    for record in SeqIO.parse(example_file, "fasta"):
        record_names.append(record.name)
    #check there is only one name used 
    assert len(set(record_names)) == 1


def test_consensus_seq():
    """check if a file containing multiple instances of the same sequence has a consensus_seq output exactly the same"""
    filename1 = r"test_file.fasta"
    filename2 = r"consensus_test_file.fasta"
    #get a sequence of the test_sequence and test_sequence consensus
    for record in SeqIO.parse(filename1, "fasta"):
        test_sequence = record.seq
    for record in SeqIO.parse(((filename2)), "fasta"):
        test_consensus =  record.seq
    #check that the sequence and a consensus of the same sequence but copied multiple times results in the same output 
    assert test_sequence == test_consensus

def test_consensus_file():
    """checks if all breeds are present in consensus file"""
    #initialise the breed class
    initialise_Breed(dog_breeds) 
    consensus_file_breeds = []
    #get all breed names from the consenus file and add them to a list
    for key in consensus_file():  
        consensus_file_breeds.append(key.id)
    all_breeds = unique_breeds()
    #check if all of the breeds are dound in the consensus file   
    for breed in all_breeds:
        assert breed in consensus_file_breeds



def test_align_consensus():
    """test if the align_consensus function returns test sequence from a list of sequences which includes the test sequence"""
    #initialise both clases
    initialise_Breed(dog_breeds)  
    consensus_file = f"{output}/consensus_sequences"
    initialise_Breed_consensus(consensus_file)  
    #get the test sequence
    test_sequence = read_fasta(r"BOXER_test_sequence.fasta")  
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


def test_simple_alignment():
    """test to check that an alignment of two identical sequences returns an alignment with 100% percentage similarity"""
    #get a random DNA sequence
    sample_seq = "AGCTAGCATAATCGG"
    #get the alignment and score of the sequence with itself
    alignment_score = simple_alignment(sample_seq,sample_seq) 
    #access the alignment only
    alignment = alignment_score[1]
    #check percentage similarity of the alignment  
    similarity = percentage_similarity(alignment[0])
    #assert the percentage similarity is 100  
    assert similarity == 100


def test_initialise_Breed():
    """check if all breeds in the class Breed are capittalised"""
    #initialise the class
    initialise_Breed(dog_breeds)
    #find all instances of the class  
    for instance in Breed.all_instances:
        #check the breed for each instance is capitalised  
        assert instance.breed.isupper() 

def test_add_mystery_to_consensus():
    """Checks if the mystery file sequence hass succesfully been added to consensus sequences"""
    #get mystery sequence
    mystery = (read_fasta(mystery_breed))[0]  
    #get consensus sequences with mystery sequence list 
    records = add_mystery_to_consensus()  
    #check mystery sequence is in consensus
    sequences = []
    #find all the sequences in the consensus
    for record in records:
        sequences.append(record.seq)
    #check that mustery sequence is in that list 
    assert mystery in sequences


def test_initialise_Breed_consensus():
    """Test to check that <unknown description> has been sucessfully removed from breeds name during initalisation of the Breed_consensus class"""
    #initialise breed consensus
    consensus_file = f"{output}/consensus_sequences"
    initialise_Breed_consensus(consensus_file)  
    #check unknown description not in breed name
    breeds = []
    for instance in Breed_consensus.all_instances:  
        breeds.append(instance.breed)
    assert "<unknown description>" not in breeds
    
def test_MSA_alignment():
    """check MSA alignment is returned in the correct format"""
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    #create a MSA alignment
    align = MSA_alignment(consesnsus_with_mystery)  
    #check that the alignnt is a MultipleSeqAlignment
    assert type(align) == Bio.Align.MultipleSeqAlignment  

def test_build_phylo_tree():
    """check if the function returns the tree in correct format"""
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    alignment = MSA_alignment(consesnsus_with_mystery)  
    tree = build_phylo_tree(alignment)
    #check if the tree is a Phylogentic Tree
    assert type(tree) == Bio.Phylo.BaseTree.Tree  

def test_custom_label_func():
    """check expected labels are found but unexpected arent"""
    #create an MSA alignment
    consesnsus_with_mystery = f"{output}/consensus_sequences_with_mystery"
    alignment = MSA_alignment(consesnsus_with_mystery) 
    #use the alignment to make a tree
    tree = build_phylo_tree(alignment)
    #check the unexpected labels are not found   
    assert tree.find_any("XXX") is None
    #check expected labels are found 
    assert str(tree.find_any("BOXER")) == "BOXER"

def test_base_content():
    """check is a random sequence is created with the correct wights of bases"""
    #create an arbitrary tets sequence
    test_seq = "AAAAATTTTT"
    #check the sequence base content 
    content = base_content(test_seq) 
    #check there is no G 
    assert content[0] == 0
    #check there is no C
    assert content[1] == 0
    #check if half of the content is A
    assert content[2] == 50
    #check if the other half is T
    assert content[3] == 50

def test_random_DNA():
    """check if random DNA has base content close to the one which was defined"""
    #define base content 
    content = [0,0,50,50]
    #create a random DNa sequence with predifned base content
    DNA = random_DNA(100,content)  
    random_DNA_content = base_content(DNA)  
    #check bases with defined weight of 0 are not present
    assert random_DNA_content[0] == 0
    assert random_DNA_content[1] == 0
    #base content is not always exact as the sequence is random so checking if A and T have at least 30% base content
    assert random_DNA_content[2] > 30
    assert random_DNA_content[3] > 30

def test_random_sequences():
    """check if the number of sequences generated matches target database size"""
    #generate sequences 
    base_count = [0, 0, 47, 53]
    test_database = random_sequences(10, 25, base_count)
    #check the correct number of sequences  was generated  
    assert len(test_database) == 10

def test_breed_alignment_scores():
    """checks if the scores generated are a number greater than 1"""
    #get some scores 
    sequence = (read_fasta(mystery_breed))[0]  
    scores = breed_alignment_scores(sequence)
    #check if the scores are more than 1  
    for score in scores.values():
        assert score[0] > 1

def test_alignment_scores():
    """checks if the scores generated are a number greater than 1"""
    #get some scores 
    target = (read_fasta("tets_mystery_sequence.fasta"))[0]  
    sequences = (read_fasta("test_sequences.fasta"))  
    scores = alignment_scores(target, sequences) 
    #check if the scores are greater than 1 
    for score in scores:
        assert score > 1

def test_probabilitydf():
    """check that the probability_df function returns a pandas dataframe"""
    scores = [(11111, 45678), (22222, 4567), (33333, 56789)]
    breeds = ["X", "Y", "Z"]
    #arbitrary l and u values
    df = probability_df(breeds, scores, 0.5, 0.09, "test_df")  
    assert type(df) == pandas.core.frame.DataFrame



