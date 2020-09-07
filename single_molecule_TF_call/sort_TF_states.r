#author: arnaud.krebs@embl.de
#created 30.06.2020
#Script to classify reads around binding sites of TFs genome-wide
#designed for TFs binding in isolation (no binding sites within 35bp of the center of the motif)
#output #1: a list of read IDs sorted by states that can be used for state quantification and single molecule plotting
#output #2: count and frequency matrices describing the binding behaviour at individual motifs


#Dependancies:
#QuasR
#BSGenome
# setwd("/g/krebs/sonmezer/Rscripts/github/Soenmezer_2020_SMF")
## Important: Define WD and set working directory to our github  directory
## If working with different directories mmake sure your WD contains a subdirectory called rds
WD = "/g/krebs/sonmezer/Rscripts/github/Soenmezer_2020_SMF"
setwd(WD)
##Define output folder
out_path=paste(WD,"/rds/",sep = "")

	#load single molecule functions

	source('./single_molecule_TF_call/functions/single_molecule_manipulation_functions.r') #load the ranges

	#####################
	#example call
	###################

	library(QuasR) 
	library("BSgenome.Mmusculus.UCSC.mm10")
		
	
#use a set of reference genomic regions (i.e TF binding motifs mapped on the genome)
	wmMapped_SOu=readRDS('./single_molecule_TF_call/example_sets/mapped_jaspar_ChIP_bound_motifs.rds')
	
#load the reference alignment files created by QuasR
	sampleSheet='./examples/QuasR_input.txt'
	
	NOMEproj=qAlign(sampleFile=sampleSheet, 
			genome="BSgenome.Mmusculus.UCSC.mm10", 
			paired="fr", 
			bisulfite="undir")
	NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])

	## Important must: "all necessary alignment files found" 

#################
#Sort reads around TFs
##################
	#############
	#step1 
	#############
	
	#split the genome in chunks to avoid memory overusage
	#create genomic chunks (based on TF positions)
	#create tiles
	#sort TFobject
	
	wms=(wmMapped_SOu[order(start(wmMapped_SOu))])
	#split by chr
	wms.c=split(wms,seqnames(wms))
	regs=unlist(GRangesList(lapply(sl(wms.c), function(i){
		x=wms.c[[i]]
		l=length(x)
		if(l>0){
		#cut by TF position
		seq.i=seq(1,l,round(l/25))
		st.seq=start(x[seq.i])
		lowB=st.seq[1:(length(st.seq)-1)]-1000
		lowB[lowB<0]<-1
		highB=st.seq[2:(length(st.seq))]+1000
		highB[highB>seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]]=seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]
		gr=GRanges(seqnames(x[1]),IRanges(lowB,highB))
		gr
		}else(GRanges())
		})))
	
	#############
	#step2 
	#############
	
	##Important: Set Output path to your desired directory
	## This is where we will save the single molecule vectors
	getwd()
	out_path=paste(WD,"/rds/",sep = "")

## Define which samples and how many cores
		sampleNames=unique(alignments(NOMEproj)[[1]][,2])
		nb_cores=10
	## 	
	
		for(sp in seq(length(sampleNames))){
			regsD <- cbind(as.data.frame(regs)[,1:3],sample=sampleNames[sp])
			regsD=regsD[regsD[,1] %in% paste('chr', 1:19,sep=''),]
			regsL <- split(regsD,1:nrow(regsD))
  
			sRs <- mclapply( sl(regsL),function(i){
				sR=sorTFbyReg_MM_m2(regsL[[i]],sampleSheet,wmMapped_SOu,inMs=c(-7,7),upMs=c(-35,-25),doMs=c(25,35))  	
				sR
				},mc.cores=nb_cores)
			
			 sRs.u <- do.call(c, sRs)
			 names(sRs.u)=unlist(lapply(sRs,names))
	
			 saveRDS(sRs.u,paste(out_path,'tmp_',sampleNames[sp],'_sorted_reads_TFmats.rds',sep=''))
		}


 	########################################
	#step3 Compute binding frequencies for each TF
	#########################################
			#create a verctor describing all 2^3 states 
			allPos=expand.grid(c(0,1),c(0,1),c(0,1))
			patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
			
			#create a individual count matrices per sample
			#only covered TFs will be represented
		
				countMats=lapply(seq(length(sampleNames)), function(sp){
					sRs.u=readRDS(paste(out_path,'tmp_',sampleNames[sp],'_sorted_reads_TFmats.rds',sep=''))

					Rcounts=mclapply(sl(sRs.u),function(i){
						counts=unlist(lapply(sRs.u[[i]],length))
						oV=rep(0,length(patternStrings))
						names(oV)=patternStrings	
						oV[string.split(names(counts),'_',3)]=counts
						oV
						},mc.cores=10)
					names(Rcounts)=names(sRs.u)
					countMat=do.call(rbind,Rcounts)
					countMat		
					})
	
			#create a individual count matrices per sample
			#only covered TFs will be reprensented
			#apply coverage cutoff
				cutoff=20 #minimal coverage of the three bins
				
				freqMat=mclapply(seq(length(sampleNames)),function(sp){
						totC=1+apply(countMats[[sp]],1,function(x){sum(x,na.rm=T)})
	
						apply((countMats[[sp]]),2,function(x){
							fid=totC>cutoff
							y=(x/totC)*100
							y[!fid]<-NA
							y
						})},mc.cores=20)
		
			#Expand the matrix to all the TFs in the input object
			#to enable comparisons between samples
		
			
			freqMat.e=mclapply(seq(length(sampleNames)),function(sp){
					m=matrix(NA,ncol=length(patternStrings),nrow=length(wmMapped_SOu))
					colnames(m)=patternStrings;rownames(m)=names(wmMapped_SOu)
						m[rownames(freqMat[[sp]]),]=freqMat[[sp]]
					m},mc.cores=10)
	
				names(freqMat.e)=sampleNames

 	
				########################################
				#step4 Group states 
				#########################################
				#Several states relate to nucleosome occupancy (various postions of the nucleosome) 
				#and can be grouped for the quantitative analysis	
				########################################
				#generate grouping factors
				#########################################
				allPos=expand.grid(c(0,1),c(0,1),c(0,1))
				patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
				#using only 'pure' states
				states=list(
				  nucleosome=patternStrings[c(1,sl(patternStrings)[!sl(patternStrings)%in%c(1,6,8)])],
				  unbound=patternStrings[8],
				  bound=patternStrings[6]
				)
				statesF=as.factor(unlist(lapply(sl(states),function(i){rep(names(states[i]),length(states[[i]]))}))[order(unlist(states))])
				statesF=factor(statesF,levels=names(states))
				#create grouped counts object 
				groupedCounts.l=mclapply(seq(length(sampleNames)),function(i){
				  countMats=countMats
				  t(apply(countMats[[i]],1,function(x){tapply(x,statesF,sum )} ))
				},mc.cores=10)
				names(groupedCounts.l)=(sampleNames)
			
					#get a grouped frequency matrix
				cutoff=20
				groupedFreq.l=lapply(sl(groupedCounts.l),function(j){
				  groupedCounts=groupedCounts.l[[j]]
				  groupedFreq=mclapply(sl(sampleNames),function(spi){
				    totC=1+apply(groupedCounts.l[[spi]],1,function(x){sum(x,na.rm=T)})
				    apply((groupedCounts.l[[spi]]),2,function(x){
				      fid=totC>cutoff
				      y=(x/totC)*100
				      y[!fid]<-NA
				      y
				    })},mc.cores=10)
				  names(groupedFreq)=sampleNames
				  groupedFreq
				})
				saveRDS(groupedFreq.l[[1]],paste(out_path,'tmp_grouped_freqencies_TFmats.rds',sep=''))

				
				
				
				
				
				
				
				
	