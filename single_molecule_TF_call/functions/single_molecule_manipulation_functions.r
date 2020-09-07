


	#basic dependancies

	sl=function(x){seq(length(x))}
	string.split=function(string,sep,pos){unlist(lapply(string,function(x){lapply(strsplit(x,sep),'[',pos)}))}

	##################################################
	#sort reads around TF binding sites
	#to be used for factors binding in isolation
	#requires presence of GC or CG in each of the three collection bins
	#creates 3 bins 
	#output: list of sorted read IDs for each of the covered instance in the target_range (TFBS)
	##################################################

	sorTFbyReg_MM_m2<- function(regDF,sampleSheet,target_range,inMs=c(-15,15),upMs=c(-35,-25),doMs=c(25,35)){
				#regDF - genomic chunk where single molecule methylation will be extracted (keep it ~1e8bp for optimal memory usage)
				#sampleSheet - QuasR input file containing sample information. 
					#SampleNames have to be tagged with the type of treatment _NO_: only GCs will be used, _SS_: only CGs will be considered, _DE_: CG and GCs will be considered
				#target_range - GRanges of the TF binding motifs to be analysed
				#inMs - borders of the bin covering the TF binding site
				#upM -  borders of the bin upstream of the binding site(relative to the center of the motif)
				#doMs - borders of the bin downstream of the binding site  (relative to the center of the motif)
			 
				#reload libraries (required for multicore)
					library(QuasR)
					library("BSgenome.Mmusculus.UCSC.mm10")
					
					#load path to the alignments
					projAll <- qAlign(sampleSheet,"BSgenome.Mmusculus.UCSC.mm10",paired="fr",bisulfite="dir")
					projAll@aligner = "Rbowtie"
					proj <- projAll[alignments(projAll)$genome$SampleName==regDF[1,4]]
				
					#define regions to extact methylation information
					reg <- GRanges(seqnames = Rle(regDF[1,1]),ranges = IRanges(regDF[1,2], end = regDF[1,3]),seqlengths=seqlengths(Mmusculus))
					regExp <- reg # expand the region at the end to catch all the pairs within the window
					end(regExp) <- pmin(end(regExp)+2000,seqlengths(regExp)[as.character(seqnames(regExp))])


					#identify Cs positions in various contexts
					#for DE experiments
					#GC All			  
					dm3_GCs_chr_XSV <- matchPattern("GC",Mmusculus[[as.character(seqnames(reg))]])
					dm3_GCs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_GCs_chr_XSV), end = end(dm3_GCs_chr_XSV)),strand="+")
					start(dm3_GCs_chr) <- start(dm3_GCs_chr)+1
					dm3_GCs_chr_coord <- start(dm3_GCs_chr)
				
					#CG All
					dm3_CGs_chr_XSV <- matchPattern("CG",Mmusculus[[as.character(seqnames(reg))]])
					dm3_CGs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_CGs_chr_XSV), end = end(dm3_CGs_chr_XSV)),strand="+")
					end(dm3_CGs_chr) <- end(dm3_CGs_chr)-1
					dm3_CGs_chr_coord <- start(dm3_CGs_chr)
				
				
					#for NO experiments
					#GC excluding all ambiguous contexts			  
					dm3_DGCHNs_chr_XSV <- matchPattern(DNAString("DGCHN"),Mmusculus[[as.character(seqnames(reg))]],fixed=F)
					dm3_DGCHNs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_DGCHNs_chr_XSV), end = end(dm3_DGCHNs_chr_XSV)),strand="+")
					dm3_DGCHNs_chr <-resize(dm3_DGCHNs_chr,2,fix='center')
					start(dm3_DGCHNs_chr) <- start(dm3_DGCHNs_chr)+1
					dm3_DGCHNs_chr_coord <- start(dm3_DGCHNs_chr)
				
					#call methylation vectors
					methAln <- qMeth(proj,regExp,mode="allC",reportLevel="alignment")[[1]]
					
					if(length(grep("_NO_",regDF[1,4])) > 0){
						#GC
						methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]+1
						selCs <- methAln$Cid %in% dm3_DGCHNs_chr_coord
					}else if(length(grep("_SS_",regDF[1,4])) > 0){
						#CG
						methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]-1
						selCs <- methAln$Cid %in% dm3_CGs_chr_coord
					}else if(length(grep("_DE_",regDF[1,4])) > 0){
						#GC&CG simultaneously. strand collapse is more complicated. prioritize GC in GCG minus strand collapse
						plusStrandSel <- methAln$strand=="+"
						SelOvGC_plus <- ((methAln$Cid %in% dm3_GCs_chr_coord) & plusStrandSel)
						SelOvGC_minus <- ((methAln$Cid+1) %in% dm3_GCs_chr_coord) & !plusStrandSel
						SelOvCG_plus <- ((methAln$Cid %in% dm3_CGs_chr_coord) & plusStrandSel)
						SelOvCG_minus <- ((methAln$Cid-1) %in% dm3_CGs_chr_coord) & !plusStrandSel
						methAln$Cid[SelOvGC_minus] <- methAln$Cid[SelOvGC_minus]+1
						methAln$Cid[!SelOvGC_minus & SelOvCG_minus] <- methAln$Cid[!SelOvGC_minus & SelOvCG_minus]-1
						selCs <- SelOvGC_plus | SelOvGC_minus | SelOvCG_plus | SelOvCG_minus
					}else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}
	
					regT=GRanges( regDF[1,1],IRanges(regDF[1,2],regDF[1,3]))
					TFs=subsetByOverlaps(target_range,regT,type='within')
					#filter SelCs
					if(length(methAln[[1]])>0){ #control that at leas one read is found in the region
						sRs=	sortTFreads_vect_m2(methAln,selCs,TFs,inMs,upMs,doMs)
					sRs}else{list()}
					}
 
 
	 
	sortTFreads_vect_m2<-function(methAln,selCs,st,inMs=c(-15,15),upMs=c(-35,-25),doMs=c(25,35)){
			#daughter function to sorTFbyReg_MM_m2 
			Crange=IRanges(methAln[[2]][selCs],methAln[[2]][selCs])
			#create collecting intervals
			st=resize(st,1,fix='center')
			midP=start(st)
			inMP=GRanges(seqnames(st),IRanges(ifelse(as.logical(strand(st) =='+'),midP+inMs[1],midP-inMs[2]),ifelse(strand(st)=='+',midP+inMs[2],midP-inMs[1])))
			upMP=GRanges(seqnames(st),IRanges(ifelse(as.logical(strand(st) =='+'),midP+upMs[1],midP-upMs[2]),ifelse(strand(st)=='+',midP+upMs[2],midP-upMs[1])))
			doMP=GRanges(seqnames(st),IRanges(ifelse(as.logical(strand(st) =='+'),midP+doMs[1],midP-doMs[2]),ifelse(strand(st)=='+',midP+doMs[2],midP-doMs[1])))
			inO=as.matrix(findOverlaps(ranges(inMP),Crange)  )
			upO=as.matrix(findOverlaps(ranges(upMP),Crange)  )
			doO=as.matrix(findOverlaps(ranges(doMP),Crange)  )
   
		   #compute the methylation vectors
			inM=methAln[[4]][selCs][inO[,2]]
			upM=methAln[[4]][selCs][upO[,2]]
			doM=methAln[[4]][selCs][doO[,2]]
		
			doID=paste(names(st)[doO[,1]],methAln[[1]][selCs][doO[,2]],sep='_')
			upID=paste(names(st)[upO[,1]],methAln[[1]][selCs][upO[,2]],sep='_')
			inID=paste(names(st)[inO[,1]],methAln[[1]][selCs][inO[,2]],sep='_')
	
			#group Cs per region/per_read
			g.inM=round(tapply(  inM,inID,mean))
			g.upM=round(tapply(  upM,upID,mean))
			g.doM=round(tapply(  doM,doID,mean))
		
		
			u.inID=sort(unique(inID))
			u.upID=sort(unique(upID))
			u.doID=sort(unique(doID))
	
			#intersect IDS
			ids=sl( u.inID)[u.inID %in% intersect(u.inID , u.upID)& u.inID %in% intersect(u.inID , u.doID)] #makes sure the three bins are covered

			#binary met vectors
			sID=u.inID[ids]
			s.inM=g.inM[sID]
			s.upM=g.upM[sID]
			s.doM=g.doM[sID]
	
			pattern=paste(as.character(s.upM),as.character(s.inM),as.character(s.doM),sep='')
			st.id=paste(string.split(sID,'_',1),string.split(sID,'_',2),sep='_')
			read.id=string.split(sID,'_',3)
	
			#SORTED READ LISTS
			if(length(ids)>0){
			sR=split(read.id,paste(st.id,pattern,sep='_'))
			sRl=split(sR,paste(string.split(names(sR),'_',1),string.split(names(sR),'_',2),sep='_'))
	# 		sRl=sRl[order(as.numeric(names(sRl[[1]])))]
			}else{sRl=list()}
			sRl
			}
			



	##################################################
	#sort reads around clusters of TF binding sites
	#to be used for factors binding in clusters
	#creates n+2 bins (n=number of TFs in the cluster)
	#output: list of sorted read IDs for each of the covered instance in the target_range (TFBS)
	##################################################



	sorTF_clusters_MM_m2<- function(regDF,sampleSheet,target_range,inMs=c(-7,7),upMs=c(-35,-25),doMs=c(25,35),conv.rate=0,remove.duplicates=F){
		 	 	library(QuasR)
				#regDF - genomic chunk where single molecule methylation will be extracted (keep it ~1e8bp for optimal memory usage)
				#sampleSheet - QuasR input file containing sample information
							#SampleNames have to be tagged with the type of treatment _NO_: only GCs will be used, _SS_: only CGs will be considered, _DE_: CG and GCs will be considered
				#target_range - GRanges of the TF binding motifs to be analysed
				#inMs - borders of the bin covering the TF binding site
				#upM -  borders of the bin upstream of the binding site(relative to the center of the motif)
				#doMs - borders of the bin downstream of the binding site  (relative to the center of the motif)
			 	#conv.rate - filter for reads with low conversion - integer 0-100 - minimal percentage of conversion based on non CG non GC Cytosines - in development
			 	#remove.duplicates - remove identical reads - in development
			 	
			 	#reload libraries
				library("BSgenome.Mmusculus.UCSC.mm10")
				
				projAll <- qAlign(sampleSheet,"BSgenome.Mmusculus.UCSC.mm10",paired="fr",bisulfite="dir")
				projAll@aligner = "Rbowtie"
				proj <- projAll[alignments(projAll)$genome$SampleName==regDF[1,4]]
				
				#define regions to extact methylation information
				reg <- GRanges(seqnames = Rle(regDF[1,1]),ranges = IRanges(regDF[1,2], end = regDF[1,3]),seqlengths=seqlengths(Mmusculus))
				regExp <- reg # expand the region at the end to catch all the pairs within the window
				end(regExp) <- pmin(end(regExp)+2000,seqlengths(regExp)[as.character(seqnames(regExp))])


				#identify Cs positions in various contexts
				#for DE experiments
				#GC All			  
				dm3_GCs_chr_XSV <- matchPattern("GC",Mmusculus[[as.character(seqnames(reg))]])
				dm3_GCs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_GCs_chr_XSV), end = end(dm3_GCs_chr_XSV)),strand="+")
				start(dm3_GCs_chr) <- start(dm3_GCs_chr)+1
				dm3_GCs_chr_coord <- start(dm3_GCs_chr)
				
				#CG All
				dm3_CGs_chr_XSV <- matchPattern("CG",Mmusculus[[as.character(seqnames(reg))]])
				dm3_CGs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_CGs_chr_XSV), end = end(dm3_CGs_chr_XSV)),strand="+")
				end(dm3_CGs_chr) <- end(dm3_CGs_chr)-1
				dm3_CGs_chr_coord <- start(dm3_CGs_chr)
				
				#for NO experiments
				#GC excluding all ambiguous contexts			  
				dm3_DGCHNs_chr_XSV <- matchPattern(DNAString("DGCHN"),Mmusculus[[as.character(seqnames(reg))]],fixed=F)
				dm3_DGCHNs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_DGCHNs_chr_XSV), end = end(dm3_DGCHNs_chr_XSV)),strand="+")
				dm3_DGCHNs_chr <-resize(dm3_DGCHNs_chr,2,fix='center')
				start(dm3_DGCHNs_chr) <- start(dm3_DGCHNs_chr)+1
				dm3_DGCHNs_chr_coord <- start(dm3_DGCHNs_chr)
				
				#call methylation vectors
			  	methAln <- qMeth(proj,regExp,mode="allC",reportLevel="alignment")[[1]]
				if(length(methAln[[1]])>0){ #control that at leas one read is found in the region
				
					#filter PCR duplicates and low conversion reads
					if(!is.na(conv.rate)| remove.duplicates ){
					
						#WCW for conversion and de-duplication
						dm3_WCWs_chr_XSV_p <- matchPattern("HCH",Mmusculus[[as.character(seqnames(reg))]],fixed="subject")
						dm3_WCWs_chr_XSV_n <- matchPattern(reverseComplement(DNAString("HCH")),Mmusculus[[as.character(seqnames(reg))]],fixed="subject")
	
						dm3_WCWs_chr_p <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_WCWs_chr_XSV_p), end = end(dm3_WCWs_chr_XSV_p)),strand="+")
						dm3_WCWs_chr_n <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_WCWs_chr_XSV_n), end = end(dm3_WCWs_chr_XSV_n)),strand="-")
		
						#start(dm3_WCWs_chr) <- start(dm3_WCWs_chr)+1
						dm3_WCWs_chr_coord <- c(start(dm3_WCWs_chr_p)+1,start(dm3_WCWs_chr_n)+1)
	
						###################
						#identify molecules with low conversion rates
						###################
						methAln_WCW<-methAln
						selCs_WCW <- methAln_WCW$Cid %in% dm3_WCWs_chr_coord 
						readIDs=methAln_WCW[[1]][selCs_WCW]
						metV=methAln_WCW[[4]][selCs_WCW]
						molConvRates=tapply(metV,readIDs,function(x){
								mean(x,na.rm=T)
							})
						conversion.id=(1-molConvRates)*100>conv.rate 
						# 	
						###################
						#remove duplicates based on conversion errors (useful for amplicon data)
						###################
							#identify duplicates
						if(remove.duplicates==T){
							sm=split(metV,readIDs)
							sm2=lapply(sl(sm),function(i){paste(sm[[i]],sep='',collapse ='')})
							unique.id=!duplicated(sm2)
							selReads=names(conversion.id)[conversion.id&unique.id]
						methAlnF=lapply(sl(methAln),function(i){methAln[[i]][methAln[['aid']]%in% selReads]})	
							
						}else{
						selReads=names(conversion.id)[conversion.id]
						methAlnF=lapply(sl(methAln),function(i){methAln[[i]][methAln[['aid']]%in% selReads]})	
						}
				
						names(methAlnF)=names(methAln)
						methAln=methAlnF
						}
				
					
				if(length(grep("_NO_",regDF[1,4])) > 0){
					#GC
					methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]+1
					selCs <- methAln$Cid %in% dm3_DGCHNs_chr_coord
				}else if(length(grep("_SS_",regDF[1,4])) > 0){
					#CG
					methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]-1
					selCs <- methAln$Cid %in% dm3_CGs_chr_coord
				}else if(length(grep("_DE_",regDF[1,4])) > 0){
					#GC&CG simultaneously. strand collapse is more complicated. prioritize GC in GCG minus strand collapse
					plusStrandSel <- methAln$strand=="+"
					SelOvGC_plus <- ((methAln$Cid %in% dm3_GCs_chr_coord) & plusStrandSel)
					SelOvGC_minus <- ((methAln$Cid+1) %in% dm3_GCs_chr_coord) & !plusStrandSel
					SelOvCG_plus <- ((methAln$Cid %in% dm3_CGs_chr_coord) & plusStrandSel)
					SelOvCG_minus <- ((methAln$Cid-1) %in% dm3_CGs_chr_coord) & !plusStrandSel
					methAln$Cid[SelOvGC_minus] <- methAln$Cid[SelOvGC_minus]+1
					methAln$Cid[!SelOvGC_minus & SelOvCG_minus] <- methAln$Cid[!SelOvGC_minus & SelOvCG_minus]-1
					selCs <- SelOvGC_plus | SelOvGC_minus | SelOvCG_plus | SelOvCG_minus
				}else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}
	
				regT=GRanges( regDF[1,1],IRanges(regDF[1,2],regDF[1,3]))
				
				#subset the target_range for ranges within the collection bin 
					nb_ov=countOverlaps(unlist(target_range),regT,type='within')
					element_length= lengths(target_range)
					fact=string.split(names(unlist(target_range)),'_',3)
					fact=as.factor(fact)
					cluster.i=tapply(nb_ov,fact,sum)
					ids=paste('TFBS_cluster_',string.split(names(which(cluster.i>0)),'\\.',1),sep='')
					TFs=target_range[ids]
				
				#TFs=subsetByOverlaps(unlist(target_range),regT,type='within')
			
				#filter SelCs
					sRs=sortTF_clustersreads_vect_m2(methAln,selCs,TFs,inMs,upMs,doMs)
				sRs}else{list()}
				}
 

	 
	sortTF_clustersreads_vect_m2<-function(methAln,selCs,st,inMs=c(-15,15),upMs=c(-35,-25),doMs=c(25,35)){
		 #Daugther function of sorTF_clusters_MM_m2
		 #create  Crange
			Crange=IRanges(methAln[[2]][selCs],methAln[[2]][selCs])
 		
 		
 			#create collecting intervals
				#get center of each TF binding site
				TFcenter=resize(unlist(st),1,fix='start')
				inMP=GRanges(seqnames(TFcenter),IRanges(start(TFcenter)+inMs[1],start(TFcenter)+inMs[2]))
				names(inMP)=string.split(names(TFcenter),'\\.',2)
				inMPs=split(inMP,string.split(names(TFcenter),'\\.',1))
				#get borders of the clusters
				upMPs=unlist(GRangesList(lapply(sl(st),function(i){
					#calculate borders
					TFcluster=st[[i]]
					minV=order(start(TFcluster),decreasing=F)[1]
					upMP=GRanges(seqnames(TFcluster[minV]),IRanges(start(TFcluster[minV])+upMs[1],start(TFcluster[minV])+upMs[2]))
					upMP
					})))	
				names(upMPs)=names(st)
				
				doMPs=unlist(GRangesList(lapply(sl(st),function(i){
					#calculate borders
					TFcluster=st[[i]]
					maxV=order(start(TFcluster),decreasing=T)[1]
					doMP=GRanges(seqnames(TFcluster[maxV]),IRanges(start(TFcluster[maxV])+doMs[1],start(TFcluster[maxV])+doMs[2]))
					doMP
					})))
				names(doMPs)=names(st)
				
			#generate dynamic bins depending on the conposition of the TF cluster
			nb_TF=lengths(st)
		
		
			#group clusters by number of TFs
			#create dynamic bins
			#d_bins are grouped by binning type (number of TFs in the cluster)
			d_bins=lapply(sl(sort(unique(nb_TF))),function(i){
					#i=2
					nbTF=sort(unique(nb_TF))[i]
					TFrange=inMPs[lengths(inMPs)==nbTF]
					
					se=lapply(sl(TFrange),function(cl){
						lapply(seq(nbTF),function(pos){TFrange[[cl]][pos]})
					})
					
					inMPsd=lapply(seq(nbTF),function(pos){
						do.call(c,lapply(se, `[[`, pos))
					})
					
					cluster_id=names(st)[i]					
					bins=GRangesList(c(list(upMPs[names(TFrange)]),inMPsd,list(doMPs[names(TFrange)])	))
					bins
				})
				names(d_bins)=paste('clusterOf_',sort(unique(nb_TF)),'_TF',sep='')
 		
 	
				#iterate ober bin types
				sRls=lapply(sl(d_bins),function(bin_t){
					#print(bin_t)
					bins=d_bins[[bin_t]]
					#find overlaping Cs
					ovs=lapply(sl(bins),function(i){
						as.matrix(findOverlaps(ranges(bins[[i]]),Crange)  )
						})
				   #compute the methylation vectors
					Mvs=lapply(sl(bins),function(i){
						methAln[[4]][selCs][ovs[[i]][,2]]
						})
				  #compute the readIDs
					IDvs=lapply(sl(bins),function(i){
						paste(names(bins[[1]])[ovs[[i]][,1]],methAln[[1]][selCs][ovs[[i]][,2]],sep='_')
						})	
					#group Cs per region/per_read
					g.Ms=lapply(sl(bins),function(i){
						round(tapply(  Mvs[[i]],IDvs[[i]],mean))
						})
					#find unique IDs
					u.IDs=lapply(sl(bins),function(i){
						sort(unique(IDvs[[i]]))
						})
		
		
					intMat=do.call(cbind,lapply(sl(bins),function(i){
						u.IDs[[1]]%in% u.IDs[[i]]
						}))
	
					#intersect IDS
					ids=sl( u.IDs[[1]])[rowSums(intMat)==ncol(intMat)] #makes sure the three bins are covered
		
				#  	s.outID=sl( u.outID)[u.outID %in% intersect( u.inID ,u.outID)]
					#binary met vectors
					sID=u.IDs[[1]][ids]
		
		
					s.M=lapply(sl(bins),function(i){
						g.Ms[[i]][sID]
						})
			
					patternMat=do.call(cbind,s.M)
					#print(patternMat)
					pattern=apply(patternMat,1,function(x){paste(as.character(x),sep='',collapse='')})
		
					st.id=paste(string.split(sID,'_',1),string.split(sID,'_',2),string.split(sID,'_',3),sep='_')
					read.id=string.split(sID,'_',4)
					#TABULAR COUNTS
				# 	p.counts=table(st.id,pattern)
		
					#SORTED READ LISTS
					if(length(ids)>0){
						sR=split(read.id,paste(st.id,pattern,sep='_'))
						sRl=split(sR,paste(string.split(names(sR),'_',1),string.split(names(sR),'_',2),string.split(names(sR),'_',3),sep='_'))
			
			
					}else{sRl=list()}
	
					sRl
					})
					#merge the results
					
					sRlsm=do.call(c, sRls)			
		}
  		 
