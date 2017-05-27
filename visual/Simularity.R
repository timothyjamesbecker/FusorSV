#simulairity heat map plot for each type: DEL, INV, DUP

#visualizations of SVU features and simularity metrics
library('ggplot2');
library('grid');
library('gridExtra');
library('plyr');
library('reshape2');

get_sample_name<-function(full_path){
	s <- strsplit(full_path,'/')[[1]]
	s[length(s)-1]
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

change_caller_name<-function(data,n1,n2){
	data[which(data[,'C1']==n1),'C1'] <- n2;
	data[which(data[,'C2']==n1),'C2'] <- n2;
	data
}

sim_matrix_path <- '/Users/tbecker/Desktop/meta_caller_NP4_result/visual/sim_matrix.0123456789ABCDEF101112131415161718191A.tsv';
#sim_matrix_path <- '~/Desktop/TEMP/varsim_R1_big/visual/sim_matrix.tsv';
col_names <- c('C1','C2','type','bin','j');
data <- read.table(sim_matrix_path,header=T,sep='	',col.names=col_names,stringsAsFactors=F);
types <- c('DEL','DUP','INV'); #dig out the interesting ones here
t_size <- 18;

data <- change_caller_name(data,'G1K-P3','True');

#dig out averages accross the bins
DEL <- data[which(data[,'type']=='DEL'),c('C1','C2','bin','j')];
DEL <- ddply(DEL, .(C1,C2),summarize,j=mean(j)); #take averages
DUP <- data[which(data[,'type']=='DUP'),c('C1','C2','bin','j')];
DUP <- ddply(DUP, .(C1,C2),summarize,j=mean(j)); #take averages
INV <- data[which(data[,'type']=='INV'),c('C1','C2','bin','j')];
INV <- ddply(INV, .(C1,C2),summarize,j=mean(j)); #take averages

del_callers <- rev(c('True','cnMOPS','CNVnator','Delly','GenomeSTRiP','Hydra','Lumpy','BreakDancer'));
p1 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(DEL,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='red',limit=c(0.0,1.0),name='Jaccard\nSimularity') +
	ylim(del_callers) + xlim(rev(del_callers)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle=90,vjust=1,size=t_size,hjust=1),axis.title.x = element_blank()) +
	theme(axis.text.y = element_text(angle=0,vjust=1,size=t_size,hjust=1),axis.title.y = element_blank()) +
	#theme(legend.position="none") +
	ggtitle("DEL");

dup_callers <- rev(c('True','cnMOPS','CNVnator','Delly','GenomeSTRiP','Hydra','Lumpy'));
p2 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(DUP,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='red',limit=c(0.0,1.0),name='Jaccard\nSimularity') +
	ylim(dup_callers) + xlim(rev(dup_callers)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle=90,vjust=1,size=t_size,hjust=1),axis.title.x = element_blank()) +
	theme(axis.text.y = element_text(angle=0,vjust=1,size=t_size,hjust=1),axis.title.y = element_blank()) +
	#theme(legend.position="none") +
	ggtitle("DUP");

inv_callers <- rev(c('True','Delly','Hydra','Lumpy','BreakDancer'));
p3 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(INV,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='red',limit=c(0.0,1.0),name='Average\nBase Pair\nSimularity') +
	ylim(inv_callers) + xlim(rev(inv_callers)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle=90,vjust=1,size=t_size,hjust=1),axis.title.x = element_blank()) +
	theme(axis.text.y = element_text(angle=0,vjust=1,size=t_size,hjust=1),axis.title.y = element_blank()) +
	#theme(legend.position="none") +
	ggtitle("INV");


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lgnd<-g_legend(p1 + theme(legend.key.size=unit(1.0,'in')) + 
               theme(legend.position="right",legend.direction="horizontal",legend.text=element_text(angle=0,size=t_size)));
	
#svg(filename = '~/Desktop/TEMP/fusionSVU/visual/sim_matrix.svg',width=12,height=5)
pdf('~/Desktop/sim_matrix_g1k_50X_27.pdf',width=28,height=12)
g <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"),
                              p3 + theme(legend.position="none"),
                              nrow=1),
                   lgnd, nrow=2,heights=c(8, 1.5))
dev.off()
       
       
