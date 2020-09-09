# Soenmezer_2020_SMF
Scripts for analysing SMF data from Krebs Lab @ EMBL.
These scripts depend on functions (included) for calling methylation and single molecule TF analysis
Exemplary input files and objects are provided for all analyses

# To Do
When done (or starting directly from alignments), use an input file for QuasR to point at alignments. 
Run context_methylation_call to call and analyse methylation genome-wide with context information. 
Use sort_TF_states and sort_TF_states_clusters to generate single molecule methylation vectors based on your analysis type. We recommend starting with single motifs and using sort TF_states.

# Working / Output Directories 
Please set your working directory as specified in the script, or use our downloaded github folder.
If you create 2 subdirectories in our downloaded GitHub folder, titled “rds” and “tmp”, the scripts will generate the output in these folders.

# Exemplary Data 
We provide an exemplary input file for alignment via QuasR.
We provide an exemplary input file for analysing BAMs via QuasR.

We provide exemplary R objects for TF motif lists
