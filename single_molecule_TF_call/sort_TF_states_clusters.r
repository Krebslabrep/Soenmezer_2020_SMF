#author: arnaud.krebs@embl.de
#created 30.06.2020
#Script to classify reads around binding sites of TF clusters genome-wide
#designed for analysing cluster of binding sites and determine frequency of co-occupancy of TFs
#outputs a list of read IDs sorted by states that can be used for state quantification and single molecule plotting

#Dependancies:
#QuasR
#BSGenome

## Important: Define WD and set working directory to our github  directory
## If working with different directories mmake sure your WD contains a subdirectory called rds
WD = "./github/Soenmezer_2020_SMF"
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
		
	
#use a set of reference genomic regions
## Here we exemplify TF clusters 

## Cluster Object 
## This is a list of all TF clusters we perform co-occupancy analysis on

wmMapped_SOu=readRDS('./single_molecule_TF_call/example_sets/mapped_jaspar_ChIP_bound_motifs.rds')


	

	#load the reference alignment files created by QuasR
	sampleSheet='./examples/QuasR_input.txt'
	
	NOMEproj=qAlign(sampleFile=sampleSheet, 
			genome="BSgenome.Mmusculus.UCSC.mm10", 
			paired="fr", 
			bisulfite="undir")
	NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])

	

#################
#Sort reads around TFs
##################
	#############
	#step1 
	#############
	
	#################
	#Define TF clusters to perform the sorts
			#identify TFBS clusters
			#annotate the clusters
			#number of TFBS (to set the number of states
			#identity and postion of the TFBS.
	#################
		
	#find sites with pairs of TF bound within 150bp
	#75b corespond to the size of the distance from TFBS center to largest collection bin
		ov=as.matrix(findOverlaps(resize(wmMapped_SOu,1,fix='center'),resize(wmMapped_SOu,150,fix='center'),ignore.strand=T))
	#assign distance between partners
		pair_dist=start(resize(wmMapped_SOu[ov[,2],],1,fix='center'))-start(resize(wmMapped_SOu[ov[,1],],1,fix='center'))

	#remove self & consider only positive distances (object is symetrical) > important to avoit negative distances
	#make sure that distance >25bp top avoid overlapping TFBS
		ov2=ov[ pair_dist>25 ,]	

	#assign  name of partners
		pair_name=paste(wmMapped_SOu[ov2[,1],]$name,wmMapped_SOu[ov2[,2],]$name,sep='_')
		
	#create a new object with all pairs summarised
	#start from first motif center
	#end from second motif center
		
		MotDB_pair=GRanges(
		  seqnames(wmMapped_SOu[ov2[,1]]),
		  IRanges(start(resize(wmMapped_SOu[ov2[,1]],1,fix='center')),start(resize(wmMapped_SOu[ov2[,2]],1,fix='center'))),
		  TF1=wmMapped_SOu[ov2[,1]]$name,
		  TF2=wmMapped_SOu[ov2[,2]]$name,
		  within_baits=ifelse(!(is.null(wmMapped_SOu$within_baits)), wmMapped_SOu[ov2[,1]]$within_baits&wmMapped_SOu[ov2[,2]]$within_baits, NA),
		  TF1_isBound=ifelse(!(is.null(wmMapped_SOu$isBound)), wmMapped_SOu[ov2[,1]]$isBound, NA),
		  TF2_isBound=ifelse(!(is.null(wmMapped_SOu$isBound)), wmMapped_SOu[ov2[,2]]$isBound, NA)
		)
		
	
	#reduce the object to create cluster
		TF_cluster=reduce(MotDB_pair)
		
	#numer of TFBS per cluster
		TFBS_number=countOverlaps(TF_cluster,wmMapped_SOu)
		
	#defne TF cluster as the GRanges covering the whole cluster
		TF_cluster$numer_of_TF=TFBS_number
		names(TF_cluster)=paste('TFBS_cluster_',sl(TF_cluster),sep='')
	
	#overlap cluster with TFBS
		ov3=as.matrix(findOverlaps(TF_cluster,wmMapped_SOu))
	
	
	
	#create a GRAngesList with the TFBS composing each cluster
		TF_list=split(wmMapped_SOu[ov3[,2]],ov3[,1])
		names(TF_list)=paste('TFBS_cluster_',sl(TF_list),sep='')

	#################
	#Load the TFBS clusters
	#################
	
	
	#split the genome in chunks to avoid memory over-usage
	#create genomic chunks (based on TF positions)
	#makes sure the chunks do not span a TF motif
	wmMapped_SOu=TF_cluster

	#sort TFobject
		wms=(TF_cluster[order(start(TF_cluster))])
		#split by chr
		wms.c=split(wms,seqnames(wms))
		regs=unlist(GRangesList(lapply(sl(wms.c), function(i){
			#i=21	
			print(i)
			x=wms.c[[i]]
			l=length(x)
			if(l>0){
			#cut by TF position
			seq.i=seq(1,l,ceiling(l/25))
			st.seq=start(x[seq.i])
			lowB=st.seq[1:(length(st.seq)-1)]-1000
			lowB[lowB<0]<-1
			highB=st.seq[2:(length(st.seq))]+1000
			highB[highB>seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]]=seqlengths(Mmusculus)[as.character(seqnames(wms.c[[i]])[1])]
			gr=GRanges(seqnames(x[1]),IRanges(lowB,highB))
			gr
			}else(GRanges())
			})))
			
	
	########################################
	#step2 perform vectorial read extraction
	#########################################
	
		sampleNames=unique(alignments(NOMEproj)[[1]][,2])

	for(sp in seq(length(sampleNames))){

		regsD <- cbind(as.data.frame(regs)[,1:3],sample=sampleNames[sp])
		regsD=regsD[regsD[,1] %in% paste('chr', 1:19,sep=''),]
		regsL <- split(regsD,1:nrow(regsD))
  
		print(as.character(sampleNames[sp]))
		sRs <- mclapply( sl(regsL),function(i){
			sR=sorTF_clusters_MM_m2(regsL[[i]],	sampleSheet,TF_list,inMs=c(-7,7),upMs=c(-35,-25),doMs=c(25,35))  	
			sR
			},mc.cores=20)
			 sRs.u <- do.call(c, sRs)
		 names(sRs.u)=unlist(lapply(sRs,names))
		
		
		saveRDS(sRs.u,paste(out_path,'tmp_',sampleNames[sp],'_sorted_reads_TFmats_cluster.rds',sep=''))


		}

 	########################################
	#step3 Compute binding counts and frequencies for pairs of TFs within the cluster
	#########################################
				########################################
				#Define temporary state classification for pairs of TFs
				#iterates over TFs within a cluster and considers all possible pairs of TFs
				#calculates binding frequency by the respective fator and co-vbinding frequencies
				########################################
	
				#define states for each factor separately
					allPos=expand.grid(c(0,1),c(0,1))
					TF1=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
					names(TF1)=c('nucleosome','unassigned','bound','accessible')
	
					TF2=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
					names(TF2)=c('nucleosome','bound','unassigned','accessible')
	
				#create a combined state factor			
					combined_states=unlist(lapply(sl(TF1),function(i){
						lapply(sl(TF2),function(j){
							paste(TF1[i],TF2[j],sep='')
							})}))
				
					names(combined_states)=unlist(lapply(sl(TF1),function(i){
						lapply(sl(TF2),function(j){
							paste(names(TF1)[i],names(TF2)[j],sep='_')
							})}))
	
				#factorise
				combined_statesF=as.factor(unlist(lapply(sl(combined_states),function(i){rep(names(combined_states[i]),length(combined_states[[i]]))}))[order(unlist(combined_states))])
				combined_statesF=factor(combined_statesF,levels=names(combined_states))
	
				stateM=cbind(string.split(as.character(combined_statesF),'_',1),string.split(as.character(combined_statesF),'_',2))

				#group states (important for the nucleosome occupied state)
					grouped_states=	list(
												combined_states[stateM[,1]=='bound'& stateM[,2]=='bound'],
												combined_states[stateM[,1]=='bound'& !stateM[,2]=='bound' & !stateM[,2]=='nucleosome'],
												combined_states[!stateM[,1]=='bound'& !stateM[,1]=='nucleosome' &stateM[,2]=='bound'],
												combined_states[(!stateM[,1]=='bound'& !stateM[,2]=='bound')|(stateM[,1]=='bound'& stateM[,2]=='nucleosome')|(stateM[,2]=='bound'& stateM[,1]=='nucleosome')])

					grouped_states=rev(grouped_states)
	
	
			########################################
			#create the count matrices
			#########################################
		
			nb_TF=lengths(TF_list)
			TFcat=sort(unique(nb_TF))

			#remove cluster with more than 10TFBS (optional)
			TFcat=TFcat[TFcat<=10]

		
			
			
			countMats=lapply(sl(sampleNames),function(sp){
							#load sorted reads
							 sRs.u=readRDS(paste(out_path,'tmp_',sampleNames[sp],'_sorted_reads_TFmats_cluster.rds',sep=''))
							#create summary count table 
							#calculate count matrices for each cluster type separately
							countMats=lapply(TFcat,function(nbc){
								#nbc=2
								allPos=expand.grid(lapply(seq(nbc+2),function(i){c(0,1)}))
								patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
								sRs.u.c=sRs.u[ names(sRs.u)%in%  names(nb_TF)[nb_TF==nbc] ]
								if(length(sRs.u.c)>0){
									Rcounts=mclapply(sl(sRs.u.c),function(i){
										#i=1	
										counts=unlist(lapply(sRs.u.c[[i]],length))
										oV=rep(0,length(patternStrings))
										names(oV)=patternStrings
							
										oV[string.split(names(counts),'_',4)]=counts
										oV
										},mc.cores=10)
									names(Rcounts)=names(sRs.u.c)
									countMat=do.call(rbind,Rcounts)
									countMat
							
									}else(NA)})
								names(countMats)=TFcat
								countMats
								})
					names(countMats)=sampleNames
			
			########################################
			#expand the matrices to include prositions that are not covered to ease sample comparison
			#########################################
				countMat.es=lapply(sl(sampleNames),function(sp){		
							mats=countMats[[sp]]
							
							freqMat.e=mclapply(TFcat,function(nbc){	
								allPos=expand.grid(lapply(seq(nbc+2),function(i){c(0,1)}))
						
								patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
								TFs=names(nb_TF)[nb_TF==nbc]
								if(!is.na(mats[[as.character(nbc)]])){
								m=matrix(NA,ncol=ncol(mats[[as.character(nbc)]]),nrow=length( TFs))
								colnames(m)=patternStrings;rownames(m)=TFs
								m[rownames(mats[[as.character(nbc)]]),]=mats[[as.character(nbc)]]
								m
								}else(NA)
								},mc.cores=10)
								names(freqMat.e)=names(mats)
								freqMat.e
								})
				names(countMat.es)=sampleNames
				
			
			
					########################################
					#Calculate frequency matrices
					#########################################	
						cutoff=20
						freqMats=lapply(sl(sampleNames),function(sp){
								mats=countMats[[sp]]
								fmat=mclapply(sl(mats),function(i){
									mat2=mats[[i]]
									if(sum(!is.na(mat2))>0){				
									totC=apply(mat2,1,function(x){sum(x,na.rm=T)})
									mat3=apply((mats[[i]]),2,function(x){
										fid=totC>cutoff
										y=(x/totC)*100
										y[!fid]<-NA
										y
									})		
								if(!is.matrix(mat3)){
									mat3=matrix(mat3,nrow=1,ncol=length(mat3))
									rownames(mat3)=rownames(mat2)
									colnames(mat3)=colnames(mat2)
									}
									mat3		
								}else(NA)
								},mc.cores=10)
								names(fmat)=names(mats)
								fmat
								})			
						names(freqMats)=sampleNames

				########################################
				#expand the matrices to include prositions that are not covered
				#########################################
	
						freqMat.es=lapply(sl(sampleNames),function(sp){
							#sp=10
							print(sp)		
							mats=freqMats[[sp]]
									freqMat.e=mclapply(TFcat,function(nbc){
									#nbc=5
										allPos=expand.grid(lapply(seq(nbc+2),function(i){c(0,1)}))
						
										patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
										TFs=names(nb_TF)[nb_TF==nbc]
										mat2=mats[[as.character(nbc)]]
											if( sum(!is.na(mat2))>0 ){
											
											if(!is.matrix(mats[[as.character(nbc)]])){
											mat2=matrix(mat2,nrow=1,ncol=length(mat2))
											}
											m=matrix(NA,ncol=ncol(mat2),nrow=length( TFs))
											colnames(m)=patternStrings;rownames(m)=TFs
											m[rownames(mat2),]=mat2
											m
											}else(NA)
										},mc.cores=10)
										names(freqMat.e)=names(mats)
										freqMat.e
										})
						names(freqMat.es)=sampleNames

 			