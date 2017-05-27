#visualizations of SVU features and simularity metrics
library('ggplot2');
library('reshape2');

get_sample_name<-function(full_path){
	s <- strsplit(full_path,'/')[[1]]
	s[length(s)-1]
}

get_matrix<-function(data,c){
	types <- unique(data[,'type'])
	L <- vector('list',length(types)) #store matricies by type
	names(L) <- types
	c_l <- sort(union(unique(data[,'C1']),unique(data[,'C2'])))
	n <- length(c_l)
	if(c=='S'){
		for(t in types){
			D <- data[data[,'type']==t,]
			m <- dim(D)[1]
			M <- matrix(0.0,nrow=n,ncol=n); #square matrix here
			rownames(M) <- c_l;
			colnames(M) <- c_l;
			#populate this M
			for(i in 1:m){
				M[D[i,'C1'],D[i,'C2']] <- D[i,'I_p'];
				M[D[i,'C2'],D[i,'C1']] <- D[i,'I_p'];
			}
			L[[toString(t)]] <- M;
		} 
	}
	if(c=='N'){
		for(t in types){
			D <- data[data[,'type']==t,]
			m <- dim(D)[1]
			M <- matrix(0.0,nrow=n,ncol=n); #square matrix here
			rownames(M) <- c_l;
			colnames(M) <- c_l;
			#populate this M
			for(i in 1:m){
				M[D[i,'C1'],D[i,'C2']] <- D[i,'n_C1']/(1+D[i,'n_C2']+D[i,'n_C1'])
				M[D[i,'C2'],D[i,'C1']] <- D[i,'n_C2']/(1+D[i,'n_C1']+D[i,'n_C2'])
			}
			L[[toString(t)]] <- M;
		}	
	}
	if(c=='D'){
		for(t in types){
			D <- data[data[,'type']==t,]
			m <- dim(D)[1]
			M <- matrix(0.0,nrow=n,ncol=n); #square matrix here
			rownames(M) <- c_l;
			colnames(M) <- c_l;
			#populate this M
			for(i in 1:m){
				M[D[i,'C1'],D[i,'C2']] <- D[i,'D1_p']/(1+D[i,'D2_p'])
				M[D[i,'C2'],D[i,'C1']] <- D[i,'D2_p']/(1+D[i,'D1_p'])
			}
			L[[toString(t)]] <- M;
		}	
	}	
	L
}

#read in the svu_features.csv file for each sample in the base directory: svu
svu_base_dir <- '~/Documents/CourseWork/15_2015_Fall/svu'
files   <- list.files(path=svu_base_dir,pattern='svu_features.csv',full.names=T,recursive=T);
s <- length(files)
samples <- as.vector(matrix('',nrow=s,ncol=1));
for(i in 1:length(files)){
	samples[i]<-get_sample_name(files[i]);
}
A <- vector('list',s) #store all the samples for averaging
names(A)<-samples

#28-1 samples here, as one went over the CPU walltime and can be rerun later
#for each sample, read the pairwise features for comparision
#C1,C2, t, C1_n, C2_n, i_prop, d1_prop, d2_prop
col_names <- c('C1','C2','type','n_C1','n_C2','I_p','D1_p','D2_p')
for(i in 1:s){
	#read the rows for one sample into the data frame
	#these are not uniform, so the plot has to do each seperately or we can bin them...
	data <- read.table(files[i],header=F,sep=",",col.names=col_names,stringsAsFactors=F);
	A[[samples[i]]]<- get_matrix(data,'D')
}
#average them up here
types <- unique(data[,'type'])
L <- vector('list',length(types)) #store matricies by type
names(L) <- types
c_l <- sort(union(unique(data[,'C1']),unique(data[,'C2'])))
n <- length(c_l)
for(t in types){
	L[[toString(t)]] <- matrix(0.0,nrow=n,ncol=n); #square matrix here
	rownames(L[[toString(t)]]) <- c_l;             #emtpy values
	colnames(L[[toString(t)]]) <- c_l;
	j <- s
	for(i in samples){
		X <- L[[toString(t)]]
		Y <- A[[i]][[toString(t)]]
		if(!is.null(Y) && dim(X)[1]==dim(Y)[1]){
			L[[toString(t)]] <- X+Y
		} else { j <- j - 1; }
	}
	L[[toString(t)]] <- L[[toString(t)]]/j
}
#do this for each metric and each type
melted_mat <- melt(L[['5']])
colnames(melted_mat)<-c('C1','C2','value')
ggplot(data = melted_mat, aes(x=C2,y=C1,fill=value))+
	geom_tile(color = "black")+
 	scale_fill_gradient2(low ="white",high="green",limit = c(0.0,1.0), name="Unique Area\nProportion") +
  	theme_minimal()+ 
 	theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1))+
 	coord_fixed()


