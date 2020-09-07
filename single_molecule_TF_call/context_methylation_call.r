#author: arnaud.krebs@embl.de
#created 02.07.2020
#script to call methylation average in GC and CG context
#splits the call for individual chromosomes to avoid memory overload (only needed for large genomes)

## Define your working directory to be SMF github directory
WD = "./github/Soenmezer_2020_SMF"
setwd(WD)
##Define output folder
out_path=paste(WD,"/rds/",sep = "")
 
	library(BSgenome)
	library(BSgenome.Mmusculus.UCSC.mm10)
	library(QuasR)


	source('./single_molecule_TF_call/functions/context_methylation_functions.r') 	
			
	########
	#Arguments
	########
	Qinput='./examples/QuasR_input.txt' #QuasR alignement file containing a sample
	cO=19 #minimal coverage
	out.dir='/g/krebs/krebs/Rscripts/Github/methCall/tmp/' #output directory
	nb.cores=20 #nb of cores to be used

	cluObj=makeCluster(nnodes=nb.cores)

#############################
# Load Data
#############################
	#create new file with non unique IDs (for merging with QuasR)
	NOMEproj=qAlign(sampleFile=Qinput, 
            genome="BSgenome.Mmusculus.UCSC.mm10", 
            paired="fr", 
            bisulfite="undir",
            clObj=cluObj)
	
	NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])
	NOMEproj@aligner = "Rbowtie"

#############################
# Partition a genome by chromosome ("natural partitioning")
#############################
	musmus_length=seqlengths(Mmusculus) 														 # Length of all chromosome of MusMusculus
	tiles <- tileGenome(musmus_length, tilewidth=max(musmus_length),cut.last.tile.in.chrom=TRUE) # Get a range of all chro (start,end ...1 to seqlengths(Mmusculus))
	cO=19                  

#############################
# Call the methylation genome wide for all Cs, loop/chromosome
#############################
	lapply(1:21,function(i){																		# Loop from the chromosome 1 to 21
		print(i)
		
		meth_gr <- qMeth(NOMEproj,mode="allC",  tiles[1:21][i])#,clObj=cluObj) 			
		#call context methylation with coverage cutOff ( % of methylation)
		#cO=19 																				# cutoff 20...20 depht min of each base...otherwise is NA
		contextMet=call_context_methylation_v2(meth_gr ,cO,Mmusculus) 							# contextMet : list of 2 objects...CG and GC object..First part GRange, secoond part elementMetadata
		saveRDS(contextMet,paste(out_path,'Context_met_call',string.split(input_file,'input',2),'_',as.character(seqnames( tiles[1:21][i])),'_Co',as.character(cO),'.rds',sep=''))
		})


#############################
# Call the methylation genome wide for all Cs, loop/chromosome
#############################
		
	AllCf=mclapply(1:21,function(i){ 															# AllC final: list of all C of all chro...each chro an objetc
		contextMet=readRDS(paste(out_path,'Context_met_call',string.split(input_file,'input',2),'_',as.character(seqnames( tiles[1:21][i])),'_Co',as.character(cO),'.rds',sep=''))
		CG=contextMet[[1]]																		# extract CG object
		GC=contextMet[[2]] 																		# extract GC object	
		AllC=c(CG,GC)																			# Bind CG and GC data
		met=elementMetadata(AllC) 																# Extract elementMetadata of AllC (delete GRanges column)
		met2=met[,1:(ncol(met)-1)] 																# delete last colonne to have just the muneric column ( % of methylation)
	#	cov.inx=!rowSums(is.na(met2))==ncol(met2) 												# True if NA, False if not
		AllCf=AllC#[cov.inx] 																	#delete all line with NA for each sample...leaves lines with min 1 numeric data
		AllCf
	},mc.cores=nb.cores)

	AllC=unlist(GRangesList(AllCf)) 															# Fusionne les donnees de tt les chro...qui etaient en objetc
	AllC=sort(AllC)																				# Range par chro ( car AllC GRange)

	saveRDS(AllC, paste(out_path,'../Context_methylation_call',string.split(input_file,'input',2),'.rds',sep=''))

