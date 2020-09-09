# Soenmezer_2020_SMF
Scripts for analysing SMF data from Krebs Lab @ EMBL.
These scripts depend on functions (included) for calling methylation and single molecule TF analysis.
Exemplary input files and objects are provided for all analyses.

# Workflow
When starting from raw sequencing reads, use an input file to point at reads (exemplified), run SMF_align script and align reads via QuasR.
When done (or starting directly from alignments), use an input file (exemplified) for QuasR to point at alignments. 
Run context_methylation_call to call and analyse methylation genome-wide with context information. 
Use sort_TF_states and sort_TF_states_clusters to generate single molecule methylation vectors based on your analysis type. We recommend starting with single motifs and using sort TF_states.
Run plot_single_molecule_examples script to plot the data, including average methylation and single molecule stacks.

# Working / Output Directories 
Please set your working directory as specified in the script, or use our downloaded github folder.
If you create 2 subdirectories in our downloaded GitHub folder, titled “rds” and “tmp”, the scripts will generate the output in these folders.

# Exemplary Data 
We provide an exemplary input file for alignment via QuasR. (QuasR_aln_input)
We provide an exemplary input file for analysing BAMs via QuasR. (QuasR_input)

We provide an exemplary R objects for TF motif lists, needed in sort_TF_states. (mapped_jaspar_ChIP_bound_motifs.rds)
This object is also used to generate an object for TF clusters, needed in sort_TF_states_clusters.

We provide an exemplary R object for amplicon coordinates, which is needed if analysing PCR-based data.
This can be used as an example, if you are interested in plotting information for your regions of interest.
