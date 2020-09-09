###SMF data aliqnment & QC


library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(QuasR)
library(viridis)
library(pheatmap)
## Important: Define WD and set working directory to our github  directory
## If working with different directories mmake sure your WD contains a subdirectory called rds
WD = "./github/Soenmezer_2020_SMF"
setwd(WD)

##Define output folder
out_path=paste(WD,"/rds/",sep = "")

## Reads are to be trimmed Pre-alignment 
#Set a working directory 

#Set the desired location of your alignments 
path=out_path
my.alignmentsDir=out_path
##Set a tmp folder for reports
tmp=paste(WD,'/tmp/',sep='')

#load the experiments
## An input file is required, specifiying sample names and matching read files, ideally in fq.gz format.
seqExp=read.table('./examples/QuasR_aln_input.txt',sep='\t',header=T)



#Generate forward and reverse reads for Paired End data:
F1=string.split(as.character(seqExp$FileName1),',',1)
F2=string.split(as.character(seqExp$FileName2),',',1)
sampleName=paste(as.character(seqExp$SampleName))


samples=as.data.frame(cbind(FileName1=F1,
                            FileName2=F2,
                            SampleName=sampleName))


#create logs to detect potential issue in trimming
con <- file("trimming.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


#start trimming 

for (i in sl(samples[,1])){
  
  spID=as.character(samples$SampleName  [i])
  
  #clip the low quality bases #remove adapters
  system(paste(
    'java -jar /g/krebs/utilities/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 5 ',  
    samples$FileName1[i],' ', samples$FileName2[i], ' ',
    './tmp/',samples$SampleName[i],'_forward_paired.fq.gz ', 
    './tmp/',samples$SampleName[i],'_forward_unpaired.fq.gz ',
    './tmp/',samples$SampleName[i],'_reverse_paired.fq.gz ', 
    './tmp/',samples$SampleName[i],'_reverse_unpaired.fq.gz ',
    'ILLUMINACLIP:/g/krebs/utilities/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18',
    sep='')
  )
}

sink() 
sink(type="message")


## Brind data into proper format for alignment 
#Paired end format

AlnInput=as.data.frame(cbind(
  FileName1=paste(tmp,samples$SampleName,'_forward_paired.fq.gz',sep=''),
  FileName2=paste(tmp,samples$SampleName,'_reverse_paired.fq.gz',sep=''),
  SampleName=as.character(samples$SampleName)))


## Create a new .txt file for the freshly trimmed samples
## This will contain the name of the samples and the temporary files for alignment

write.table(AlnInput,'QuasR_aln_input_CS18.txt',quote=F,row.names=F,sep='\t')


## STARTING ALIGNMENT 

## Can be ran as a single instance as follows:



#######
# Define Parameters

input_file <- './examples/QuasR_aln_input.txt'

QuasRdef='-k 2 --best --strata' 
print('starting alignment')

#  ALIGN with Rbowtie via QuasR
NOMEproj=qAlign(sampleFile=input_file,
                genome="BSgenome.Mmusculus.UCSC.mm10", 
                aligner="Rbowtie", 
                paired="fr", 
                bisulfite="undir", 
                projectName="MMNomeSeq", 
                alignmentsDir=my.alignmentsDir,
                alignmentParameter=paste('-e 70 -X 1000 ',QuasRdef,sep=''),
                cacheDir = tmp)




sampleNames=unique(alignments(NOMEproj)[[1]][,2])
NOMEproj@aligner = "RBowtie"


## Generate a QC file
qQCReport(NOMEproj,paste(path,'/QC_files/QC_CS_Amplicons_18_can.pdf',sep=''))
dev.off()



## Check for coverage of targets 
## Analyze methylation 

