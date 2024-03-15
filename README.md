![Process map](https://github.com/lprzem01/Coursework/blob/main/my_project/doc/Process%20map.png)
## Project's Title
<p> Mystery breed sequence identification </p>

## Project Description
<p> The code found here provides a relatively simple approach analyse an unknown sequence against a database of known sequences </p>
<p> It analyses a database provided to carry out three primary functions</p>
    <p> 1. Find the most closely aligned dog breed to the unknown dog breed sequence -> outputs a txt file with alignment details and an alignment image </p>
    <p> 2. Finds the proabbility of each alignment occuring by chance -> outputs a database of alignment probability of each breed </p>
    <p> 3. Creates a phylogentic tree of all the breeds in the database including the unknown breed -> saved as an image file </p>

### Limitations/next steps  
#### Visual output improvements
    <p>Create a heatmap representing the MSA for better visualisation of the alignment</p>
    <p>Output a scrollable best scoring sequence alignment</p>
#### Statistical significance improvements
    <p>Output e scores not just probability scores</p>
    <p>Bigger sample size to generate better estimates for pdf</p>
    <p>Calculate base content based on whole database rather than one sequence</p>
#### General improvements
    <p>Improve speed of the code</p>
    <p>More comprehensive testing </p>

## How to Install and Run the Project
<p> To run this program download the my project folder from this repository and configure your enviroment using the requiremnts.txt file found in the doc folder. Adapt the dog_breeds.fa and mystery.fa files to include your database and unknown sequence respectivelly. Alternatively use the data provided. Run the Main code to populate the results folder with all graphs and txt files outlining the results of the comparassion. </p>

## Credits
### Most significant contributions to produce this project came from the following sources 
<p> https://www.biotite-python.org/examples/gallery/sequence/local_alignment_statistics.html </p>
<p> https://pypi.org/project/pyMSAviz/ </p>
<p> https://biopython.org/wiki/Phylo </p>
