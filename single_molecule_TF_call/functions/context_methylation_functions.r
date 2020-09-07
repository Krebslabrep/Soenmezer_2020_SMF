

#basic dependancies

sl=function(x){seq(length(x))}
string.split=function(string,sep,pos){unlist(lapply(string,function(x){lapply(strsplit(x,sep),'[',pos)}))}



getCMethMatrix<-function(proj,range,samp){
	Cs=qMeth(proj, query=range,mode="allC",reportLevel="alignment")
	# use data.table to get a 1,0 matrix of methylation profiles
	all.cids=unique(Cs[[samp]]$Cid) # get all possible C locations
	# make the data.table object
	dt=data.table(meth=Cs[[samp]]$meth ,aid=Cs[[samp]]$aid ,cid=Cs[[samp]]$Cid)
	# this function converts cids to columns
	myfun2<-function(x,all.cids){
	vec=rep(-1,length(all.cids))
	names(vec)=as.character(all.cids)
	b=as.list((vec))
	b[ as.character(x$cid)]=as.double(x$meth)
	return(b)
}
	dtm=dt[,myfun2(.SD,all.cids), by=aid]
	ronames=dtm$aid
	dtm[,aid:=NULL] # remove unwanted row
	CpGm=as.matrix(dtm)
	CpGm[CpGm == -1]=NA # put NAs
	rownames(CpGm)=ronames
	return(CpGm)
} 



	mergeGC_CGmat=function(GC_CG_mat){
			CGmat=GC_CG_mat$matCG	
			GCmat=GC_CG_mat$matGC		
			uReads=unique(c(rownames(CGmat),rownames(GCmat)))
			uCs=sort(unique(c(colnames(CGmat),colnames(GCmat))))
			matCGo=matrix(nrow=length(uReads),ncol=length(uCs))
			colnames(matCGo)=uCs
			rownames(matCGo)=uReads
			matCGo[rownames(CGmat),colnames(CGmat)]=CGmat#get CGs in
			matCGo[rownames(GCmat),colnames(GCmat)]=GCmat#get GCs in
	 		matCGo

		}




	call_context_methylation_v2=function(meth_gr,cO,genome=Mmusculus){
			dm3_CGs <- vmatchPattern("CG",genome)
			dm3_CGs <- dm3_CGs[strand(dm3_CGs)=="+",]
			strand(dm3_CGs) <- "*"

			dm3_GCs <- vmatchPattern("GC",genome)
			dm3_GCs <- dm3_GCs[strand(dm3_GCs)=="+",]
			strand(dm3_GCs) <- "*"

			sel_CGs <- meth_gr %over% dm3_CGs
			sel_GCs <- meth_gr %over% dm3_GCs

			meth_CGs_gr <- meth_gr[sel_CGs]
			meth_GCs_gr <- meth_gr[sel_GCs]
		##################
		# collapse strands
		##################
			meth_CGsCol_gr <- meth_CGs_gr[seq(1,length(meth_CGs_gr),by=2)]
			end(meth_CGsCol_gr) <- end(meth_CGsCol_gr)+1
			values(meth_CGsCol_gr) <- as.matrix(values(meth_CGs_gr[seq(1,length(meth_CGs_gr),by=2)]))+as.matrix(values(meth_CGs_gr[seq(2,length(meth_CGs_gr),by=2)]))

			meth_GCsCol_gr <- meth_GCs_gr[seq(2,length(meth_GCs_gr),by=2)]
			start(meth_GCsCol_gr) <- start(meth_GCsCol_gr)-1
			values(meth_GCsCol_gr) <- as.matrix(values(meth_GCs_gr[seq(1,length(meth_GCs_gr),by=2)]))+as.matrix(values(meth_GCs_gr[seq(2,length(meth_GCs_gr),by=2)]))

		#####################
		#filter for coverage
		####################
			Tcounts=grep('_T\\>',colnames(elementMetadata(meth_CGsCol_gr)))
			Mcounts=grep('_M\\>',colnames(elementMetadata(meth_CGsCol_gr)))
	
			######
			#CGs
			######

			CG.met.mat=as.matrix(elementMetadata(meth_CGsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])
			#filter for coverage
			CovFilter=as.matrix(elementMetadata(meth_CGsCol_gr)[,Tcounts])>cO
			for (i in sl(CG.met.mat[1,])){CG.met.mat[!CovFilter[,i],i]=NA}
				#bind the GRanges with the scores
			CG.met=resize(meth_CGsCol_gr,1,fix='start')
			elementMetadata(CG.met)=CG.met.mat


			######
			#GCs
			######

			GC.met.mat=as.matrix(elementMetadata(meth_GCsCol_gr)[,Mcounts])/as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])
			#filter for coverage
			CovFilter=as.matrix(elementMetadata(meth_GCsCol_gr)[,Tcounts])>cO
			for (i in sl(GC.met.mat[1,])){GC.met.mat[!CovFilter[,i],i]=NA}
				#bind the GRanges with the scores

			GC.met=resize(meth_GCsCol_gr,1,fix='end')
			elementMetadata(GC.met)=GC.met.mat

			###########################
			#disclose context
			###########################
			#resize to single C
	#		getSeq(Mmusculus,resize(CG.met,3,fix='center')[1:10])
		#	GCcontext=getSeq(Mmusculus,resize(GC.met,3,fix='center'))
			oGCG=as.matrix(findOverlaps(GC.met,CG.met,type='equal'))
			GC.met$type='GCH'
			CG.met$type='CGH'
			GC.met$type[oGCG[,1]]='GCG'
			CG.met$type[oGCG[,2]]='GCG'

			umet=list(CG.met,GC.met)
			names(umet)=c('CG','GC')
			return(umet)
	
		}


