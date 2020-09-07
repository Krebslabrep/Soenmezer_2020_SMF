## Plottinf of single molecule information 
## Requires single molecule vectors inn your out_path

WD = "./github/Soenmezer_2020_SMF"
setwd(WD)
##Define output folder
out_path=paste(WD,"/rds/",sep = "")

library(QuasR) 
library("BSgenome.Mmusculus.UCSC.mm10")

library(QuasR)
library(RColorBrewer)
library(rtracklayer)
library(caTools)
library(gplots)
library(data.table)


#load single molecule functions

source('./single_molecule_TF_call/functions/single_molecule_manipulation_functions.r') #load the ranges
source('./single_molecule_TF_call/functions/context_methylation_functions.r') 	

#####################
#example call
###################



#use a set of reference genomic regions (i.e TF binding motifs mapped on the genome)
wmMapped_SOu=readRDS('./single_molecule_TF_call/example_sets/jaspar2018_mm10_mapped_subset.rds')

## load the clusters 
## These are targeted by the second part of the analysis 
TF_cluster=readRDS('./single_molecule_TF_call/example_sets/MotDB_clusters.rds')[[1]]


#load amplicon definition
ampliconfile <- read.table("./single_molecule_TF_call/example_sets/amplicon_coordinates.txt", header = T)
amplicons <- GRanges(ampliconfile)


		#####################
		#Load the data
		#####################

		#load the alignments files
	
#load the reference alignment files created by QuasR
sampleSheet='./examples/QuasR_input.txt'

NOMEproj=qAlign(sampleFile=sampleSheet, 
                genome="BSgenome.Mmusculus.UCSC.mm10", 
                paired="fr", 
                bisulfite="undir")
NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])

		sp.sbs=sl(sampleNames)
	

	
		
		#definition of the states
		allPos=expand.grid(c(0,1),c(0,1),c(0,1))
		patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
				states=list(
							unassigned=patternStrings[!sl(patternStrings)%in%c(1,2,5,6,8)],
				
						nucleososme=patternStrings[c(1,2,5)],
						unbound=patternStrings[8],
						bound=patternStrings[6]
						)


		statesF=as.factor(unlist(lapply(sl(states),function(i){rep(names(states[i]),length(states[[i]]))}))[order(unlist(states))])
		statesF=factor(statesF,levels=names(states))


		
			#####################################################
			#Plot example locus
			#####################################################
			#definition of the states
				allPos=expand.grid(c(0,1),c(0,1),c(0,1))
				patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
				#using only 'pure' states
				states=list(
						closed=patternStrings[1],
						accessible=patternStrings[8],
						bound=patternStrings[6],
						unassigned=patternStrings[!sl(patternStrings)%in%c(1,6,8)]
						)


				statesF=as.factor(unlist(lapply(sl(states),function(i){rep(names(states[i]),length(states[[i]]))}))[order(unlist(states))])
				statesF=factor(statesF,levels=names(states))
		
		### Here, Ensure that single molecule vectors can be read succesfully. 
			
				for (sp in sp.sbs){
			
				sRs.u=readRDS(paste(out_path,sampleNames[sp.sbs][sp],'_sR_TFmats_vect_inM14.rds',sep=''))
			
				#separate different categories
			
				hits=resize(amplicons,1,fix='center')
				strand(hits)='*'
				#i=10
				
				#####################################################
				#using single TF classification
				####################################################
	
				lapply(sl(hits),function(i){
						print(i)
						TFid=hits[i]					
						TF=TFid$TFBSname
						st=resize(TFid,600,fix='center')
						mots=subsetByOverlaps(MotDBf,st,ignore.strand=T)
					
						sm.mat=getCMethMatrix(NOMEproj,st,sampleNames[sp.sbs][sp])
					
					
						if (length(sm.mat)>0 & length(unlist(sRs.u[[as.character(TF)]]))>9){
						#extractMultipleMatrices(st,NOMEproj,sampleNames[sp],Mmusculus,1)
						GC_CGmat=getGCMatrix(sm.mat,chr=as.character(seqnames(st)[1]),genome=Mmusculus,destrand=TRUE)
					
						#separate context
						st.seq=getSeq(Mmusculus,st)
				
				
						GCpos=start(st)+start(matchPattern(DNAString("DGCHN"),st.seq[[1]],fixed="subject"))+1
						CGpos=start(st)+start(matchPattern(DNAString("NWCGW"),st.seq[[1]],fixed="subject"))+1
					
					
						GC_mat=GC_CGmat[['matGC']][,colnames(GC_CGmat[['matGC']])%in% as.character(GCpos)]
						CG_mat=GC_CGmat[['matCG']][,colnames(GC_CGmat[['matCG']])%in% as.character(CGpos)]
	

						#sort the reads
						read_sort=sRs.u[[as.character(TF)]]
						read_sort=read_sort[!is.na(names(read_sort))	]
				
						stN=string.split(names(read_sort),'_',3)	
						names(read_sort)=stN
						read_sort.e=vector("list", length(unlist(states))) 
						names(read_sort.e)=unlist(states)
						read_sort.e[stN]=read_sort
						read_sort=read_sort.e
				
						read_sort=lapply(sl(read_sort),function(i){read_sort[[i]][read_sort[[i]]%in%rownames(GC_mat)]})
						names(read_sort)=unlist(states)
						read_sort=read_sort[(unlist(states))]


						#plot 
						 png(paste(out_path,TF,'_',i,'_',sampleNames[sp.sbs][sp],'.png',sep=''),width = 480, height = 480)
						# postscript(paste('./plots/example_locus/example_',TF,'_',TFid,'_',sampleNames[sp],'.ps',sep=''),width = 480, height = 480)
					
						layout(matrix(c(1,2,3,6,4,5), 3, 2, byrow = F),	widths=c(3,1), heights=c(1,1))

						#plot the average
							plot(NA,xlim=c(start(st),end(st)),ylim=c(-0.2,1),xlab='',ylab='SMF (%)',main=TFid)	
								#GCs
							points(colnames(GC_mat),1-colMeans(GC_mat,na.rm='T'),type='l')
							points(colnames(GC_mat),1-colMeans(GC_mat,na.rm='T'),pch=20)
								#CGs

							abline(h=0)
							rect(start(mots),-0.2,end(mots),-0.15)
							text(start(resize(mots,1,fix='center')),rep(-0.1,length(mots)),mots$name,cex=0.8)
				
	# 					#Plot the single molecules
							readIDs=unlist(read_sort,recursive=T,use.names =F)#[unlist(read_sort,recursive=T,use.names =F)%in%rownames(sm.mat[[1]])]

							###################
							#plot the single molecules
							##############
							#using hierarchical clustering
								if(length(readIDs)>500){ #subset to 500 molecules to avoid problem with Hc
									GC_mat_sbs=GC_mat[sample(readIDs,500),]
								}else{GC_mat_sbs=GC_mat[readIDs,]}
								hc=hclust(dist(GC_mat_sbs))

						
							#using the classification
								M=GC_mat[readIDs,]			
								vR1=VectorizeReads(st,M)
								BR=c(col=c('black','grey'))					
								colors=BR
								plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)),ylim=c(-1,length(unique(vR1[[2]]))))#,main=paste(mots$name[1],'only'))
								points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)))

							#using the hc
								M_hc=GC_mat_sbs[hc$order,]
								vR1=VectorizeReads(st,M_hc)
								BR=c(col=c('black','grey'))					
								colors=BR
								plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)),ylim=c(-1,length(unique(vR1[[2]]))))#,main=paste(mots$name[1],'only'))
								points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)))
					
				
						#add the classification plot
							#define counts (for side bar)
							counts=unlist(lapply(read_sort,length))
							grouped.counts=unlist(lapply(sl(states),function(i){sum(counts[states[[i]]])}))
							names(grouped.counts)=names(states)
				
							TF1c=colorRampPalette(brewer.pal(9,"Set1"))(9)[c(2,3,4,9 )]
							names(TF1c)=names(states)
							boundaries=cumsum(grouped.counts)
							boundaries=boundaries
							TF1colv=lapply(seq(length(boundaries)),function(j){
								bd=grouped.counts
								rep(TF1c[names(grouped.counts)][j],bd[j])#length(seq(bd[j],bd[j+1]-1,1)))
					
							})
						
							plot(NA,xlim=c(0,3),ylim=c(0,sum(lengths(read_sort))))
							points(rep(1, sum(grouped.counts)-1),seq(sum(grouped.counts)-1),col=unlist(TF1colv),pch='_',cex=2)
							text(rep(2,length(grouped.counts)),boundaries,round(grouped.counts/sum(grouped.counts)*100))		
				

						dev.off()
					}})
					}
	
				#####################################################
				#using single cluster TF classification
				## This is the analysis where we measure TF co-occupancy
				####################################################
	
			########################################
 			#Define temporary state classification for pairs
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
			
		
			combined_statesF=as.factor(unlist(lapply(sl(combined_states),function(i){rep(names(combined_states[i]),length(combined_states[[i]]))}))[order(unlist(combined_states))])
			combined_statesF=factor(combined_statesF,levels=names(combined_states))
			
			stateM=cbind(string.split(as.character(combined_statesF),'_',1),string.split(as.character(combined_statesF),'_',2))


			grouped_states=	list(
										combined_states[stateM[,1]=='bound'& stateM[,2]=='bound'],
										combined_states[stateM[,1]=='bound'& !stateM[,2]=='bound' & !stateM[,2]=='nucleosome'],
										combined_states[!stateM[,1]=='bound'& !stateM[,1]=='nucleosome' &stateM[,2]=='bound'],
										combined_states[(!stateM[,1]=='bound'& !stateM[,2]=='bound')|(stateM[,1]=='bound'& stateM[,2]=='nucleosome')|(stateM[,2]=='bound'& stateM[,1]=='nucleosome')])
	
	

				# sp=sp.sbs[1]
				for (sp in sp.sbs){
			
				sRs.u=readRDS(paste(out_path,'tmp_',sampleNames[sp],'_sR_TFcluster_vect_inM14.rds',sep=''))
			
				#separate different categories
			
				hits=resize(amplicons,1,fix='center')
				strand(hits)='*'
			
				hits=subsetByOverlaps(resize(hits,600,fix='center'),TF_cluster[names(TF_cluster)%in%names(sRs.u) & TF_cluster$numer_of_TF==2])
				#i=1
			
				lapply(sl(hits),function(i){
						print(i)
						TFid=hits[i]					
						TF=TFid$TFBSname
						st=resize(TFid,600,fix='center')
						mots=subsetByOverlaps(MotDBf,st,ignore.strand=T)
					
						sm.mat=getCMethMatrix(NOMEproj,st,sampleNames[sp])
					
						TFBS_cluster=names(subsetByOverlaps(TF_cluster[names(sRs.u)],st))
						
						if (length(sm.mat)>0 & length(unlist(sRs.u[[as.character(TFBS_cluster)]]))>9){
						#extractMultipleMatrices(st,NOMEproj,sampleNames[sp],Mmusculus,1)
						GC_CGmat=getGCMatrix(sm.mat,chr=as.character(seqnames(st)[1]),genome=Mmusculus,destrand=TRUE)
					
						#separate context
						st.seq=getSeq(Mmusculus,st)
				
				
						GCpos=start(st)+start(matchPattern(DNAString("DGCHN"),st.seq[[1]],fixed="subject"))+1
						CGpos=start(st)+start(matchPattern(DNAString("NWCGW"),st.seq[[1]],fixed="subject"))+1
					
					
						GC_mat=GC_CGmat[['matGC']][,colnames(GC_CGmat[['matGC']])%in% as.character(GCpos)]
						CG_mat=GC_CGmat[['matCG']][,colnames(GC_CGmat[['matCG']])%in% as.character(CGpos)]
	

						#sort the reads

						read_sort=sRs.u[[TFBS_cluster]]
						read_sort=read_sort[!is.na(names(read_sort))	]
				
						stN=string.split(names(read_sort),'_',4)
						names(read_sort)=stN
					
					
						#create a read count vector
						rC=rep(0,length(unlist(grouped_states)))
						names(rC)=unlist(grouped_states) 
						rC[names(read_sort)]=lengths(read_sort)
# 						read.sort=rC

						read_sort=read_sort[unlist(rev(grouped_states))]
						read_sort=read_sort[!is.na(names(read_sort))	]
# 						
						groupF=lapply(sl(rev(grouped_states)),function(i){
							rep(i, length(rev(grouped_states)[[i]]))
						})
						
						


						#plot 
						 png(paste(out_path,TF,'_',i,'_',sampleNames[sp],'.png',sep=''),width = 480, height = 480)
						# postscript(paste('./plots/example_locus/example_',TF,'_',TFid,'_',sampleNames[sp],'.ps',sep=''),width = 480, height = 480)
					
						layout(matrix(c(1,2,3,6,4,5), 3, 2, byrow = F),	widths=c(3,1), heights=c(1,1))

						#plot the average
							plot(NA,xlim=c(start(st),end(st)),ylim=c(-0.2,1),xlab='',ylab='SMF (%)',main=TFid)	
								#GCs
							points(colnames(GC_mat),1-colMeans(GC_mat,na.rm='T'),type='l')
							points(colnames(GC_mat),1-colMeans(GC_mat,na.rm='T'),pch=20)
								#CGs

							abline(h=0)
							rect(start(mots),-0.2,end(mots),-0.15)
							text(start(resize(mots,1,fix='center')),rep(-0.1,length(mots)),mots$name,cex=0.8)
				
	# 					#Plot the single molecules
							readIDs=unlist(read_sort,recursive=T,use.names =F)#[unlist(read_sort,recursive=T,use.names =F)%in%rownames(sm.mat[[1]])]

							###################
							#plot the single molecules
							##############
							#using hierarchical clustering
								if(length(readIDs)>500){ #subset to 500 molecules to avoid problem with Hc
									GC_mat_sbs=GC_mat[sample(readIDs,500),]
								}else{GC_mat_sbs=GC_mat[readIDs,]}
								hc=hclust(dist(GC_mat_sbs))

						
							#using the classification
								M=GC_mat[readIDs,]			
								vR1=VectorizeReads(st,M)
								BR=c(col=c('black','grey'))					
								colors=BR
								plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)),ylim=c(-1,length(unique(vR1[[2]]))))#,main=paste(mots$name[1],'only'))
								points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)))

							#using the hc
								M_hc=GC_mat_sbs[hc$order,]
								vR1=VectorizeReads(st,M_hc)
								BR=c(col=c('black','grey'))					
								colors=BR
								plot(NA,pch='_',col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)),ylim=c(-1,length(unique(vR1[[2]]))))#,main=paste(mots$name[1],'only'))
								points(vR1[[1]],vR1[[2]],pch='_',cex=1.2,col=colors[as.factor(vR1[[3]])],xlim=c(start(st),end(st)))
					
				
						#add the classification plot
							TF1c=colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,9,4,9)]
							TF2c=colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,4,9,9)]
							names(TF1c)=TF1
							names(TF2c)=TF2
				
							boundaries=cumsum(lengths(read_sort))
							TF1colv=lapply(seq(length(boundaries)),function(i){
							print(i)
					
								bd=c(0,boundaries)
								co=rep(TF1c[substr(names(read_sort),1,2)][i],length(seq(bd[i],(bd[i+1]-1),1)))
								co
								#length(co)
							
							})
					
							TF2colv=lapply(seq(length(boundaries)),function(i){
								bd=c(0,boundaries)
								rep(TF2c[substr(names(read_sort),3,4)][i],length(seq(bd[i],bd[i+1]-1,1)))
							})
					
							plot(NA,xlim=c(0,3),ylim=c(0,sum(lengths(read_sort))))
							points(rep(1, sum(lengths(read_sort))-1),seq(sum(lengths(read_sort))-1),col=(unlist(TF1colv)),pch='_',cex=2)
							points(rep(1.5, sum(lengths(read_sort))-1),seq(sum(lengths(read_sort))-1),col=(unlist(TF2colv)),pch='_',cex=2)
						

						dev.off()
					}})
	
	}

		
		
	
			