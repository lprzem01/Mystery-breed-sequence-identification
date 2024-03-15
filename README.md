![Process map](https://github.com/lprzem01/Coursework/blob/main/my_project/doc/Process%20map.png)
## Project's Title
<p>DNA identification service of the closest sequence in the database to the provided sequence.</p>

## Project Description
<p> The code in this file analyses a database to carry out 3 functions</p>
    <p> 1. Find the most closely aligned dog breed to the unknown dog breed sequence -> outputs a txt file with alignment details and an alignment image </p>
    <p> 2. Finds the proabbility of each alignment occuring by chance -> outputs a database of alignment probability of each breed </p>
    <p> 3. Creates a phylogentic tree of all the breeds in the database including the unknown breed -> saved as an image file </p>
### Limitations/next steps  
•	Visual output improvements
o	Create a heatmap representing the MSA for better visualisation of the alignment
o	Output a scrollable best scoring sequence alignment
•	Statistical significance improvements
o	Output e scores not just probability scores
o	Bigger sample size to generate better estimates for pdf
o	Calculate base content based on whole database rather than one sequence
•	General improvements
o	Improve speed of the code
o	More comprehensive testing 

  The code in this file can be used to identify a dog’s breed based on DNA sequence.
It takes a sequence as input and compares it to known breeds in a database to identify the most probable breed and the percentage of confidence based on alignment.<
•	Why you used the technologies you used,
•	Some of the challenges you faced and features you hope to implement in the future.</p>

## How to Install and Run the Project
<p> To run this program download the my project folder from this repository and configure your enviroment using the requiremnts.txt file found in the doc folder. Adapt the dog_breeds.fa and mystery.fa files to include your database and unknown sequence respectivelly. Alternatively use the data provided. Run the Main code to populate the results folder with all graphs and txt files outlining the results of the comparassion. </p>

## Credits
<p> If you worked on the project as a team or an organization, list your collaborators/team members. You should also include links to their GitHub profiles and social media too.
Also, if you followed tutorials or referenced a certain material that might help the user to build that particular project, include links to those here as well.
This is just a way to show your appreciation and also to help others get a first hand copy of the project. </p>
