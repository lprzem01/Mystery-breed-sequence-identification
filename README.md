![Process map](https://github.com/lprzem01/Coursework/blob/main/my_project/doc/Process%20map.png)
# Mystery Breed Sequence Identification

## Overview

This project provides a Python-based bioinformatics pipeline for identifying an unknown dog breed sequence by comparing it against a database of known breed sequences. The workflow performs sequence alignment, statistical analysis of alignment significance, and phylogenetic tree construction to determine the most likely breed relationship.

The project demonstrates the application of sequence analysis, statistical modelling, and phylogenetics using BioPython and related scientific Python libraries.

### Key Features

- Identification of the closest matching breed sequence
- Pairwise sequence alignment and visualisation
- Statistical evaluation of alignment significance
- Probability estimation for observed alignments
- Phylogenetic tree generation
- Automated result export and visual outputs

---

## Analysis Workflow

The pipeline performs three primary analyses:

### 1. Breed Identification

The unknown sequence is aligned against all reference sequences in the database to identify the closest matching breed.

**Outputs:**

- Alignment report (.txt)
- Alignment visualisation image
- Best matching breed identification

### 2. Alignment Significance Analysis

The probability of each alignment occurring by chance is estimated using simulated alignment score distributions.

**Outputs:**

- Alignment probability table
- Statistical summary of alignment significance

### 3. Phylogenetic Analysis

A phylogenetic tree is constructed to visualise the evolutionary relationships between all known breeds and the unknown sequence.

**Outputs:**

- Phylogenetic tree image
- Comparative breed relationship analysis

---

## Project Structure

```text
my_project/
├── code/
│   ├── Main_code.py
│   ├── functions.py
│   └── __init__.py
├── data/
│   ├── dog_breeds.fa
│   └── mystery.fa
├── results/
├── docs/
│   └── requirements.txt
└── README.md
```

---

## Requirements

### Python

- Python 3.8+

### Packages

- BioPython
- NumPy
- Matplotlib
- pyMSAviz

Install dependencies using:

```bash
pip install -r docs/requirements.txt
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/lprzem01/Mystery-breed-sequence-identification.git
cd Mystery-breed-sequence-identification
```

Install the required dependencies:

```bash
pip install -r docs/requirements.txt
```

---

## Usage

### Input Files

Replace the provided FASTA files with your own sequences if desired:

- `dog_breeds.fa` – Reference breed database
- `mystery.fa` – Unknown sequence to be identified

Alternatively, use the example datasets included in the repository.

### Run the Analysis

Execute:

```bash
python Main_code.py
```

The analysis will automatically generate outputs in the `results/` directory.

---

## Outputs

The pipeline generates:

- Best matching breed identification
- Pairwise alignment reports
- Alignment visualisations
- Alignment probability statistics
- Phylogenetic tree visualisations

---

## Future Improvements

### Visualisation

- Multiple sequence alignment (MSA) heatmap
- Interactive alignment viewer
- Improved result dashboards

### Statistical Analysis

- E-value calculation
- Larger simulation datasets
- Background nucleotide frequency estimation using the full reference database

### Performance

- Optimised alignment algorithms
- Improved runtime efficiency
- Expanded automated testing coverage

---

## Technologies Used

- Python
- BioPython
- Sequence Alignment
- FASTA Processing
- Statistical Analysis
- Phylogenetics
- Data Visualisation

---

## Acknowledgements

This project was developed using concepts and resources from:

- https://www.biotite-python.org/examples/gallery/sequence/local_alignment_statistics.html
- https://pypi.org/project/pyMSAviz/
- https://biopython.org/wiki/Phylo



