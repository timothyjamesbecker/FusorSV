#simulairity heat map plot for each type: DEL, INV, DUP

#Rscript cmd_parser.R commands...
options(warn=-1);
cmd_args <- commandArgs();
scriptpath <- gsub('Simularity.R','',gsub('--file=','',cmd_args[4]));
source(paste(scriptpath,'cmd_parser.R',sep=''));
source(paste(scriptpath,'plot_utils.R',sep=''));

#sim_matrix_path <- '/fusorsv_result/visual/sim_matrix.HEX.tsv';
in_path    <- args['in_path'][[1]];
in_file_split <- strsplit(in_path,'/')[[1]];
in_file    <- in_file_split[length(in_file_split)];
out_path   <- args['out_path'][[1]];
out_file<-strsplit(in_file,'.tsv')[[1]]
out_width  <- args['out_width'][[1]];
out_height <- args['out_height'][[1]];
t_size   <- args['t_size'][[1]]; #t_size <- 18;

sim_matrix_path <- in_path;
col_names <- c('C1','C2','type','bin','j');
data <- read.table(sim_matrix_path,header=T,sep='	',col.names=col_names,stringsAsFactors=F);
types <- c('DEL','DUP','INV'); #dig out the interesting ones here

data <- change_caller_name(data,'G1K-P3','True');

#dig out averages accross the bins
DEL <- data[which(data[,'type']=='DEL'),c('C1','C2','bin','j')];
DEL <- ddply(DEL, .(C1,C2),summarize,j=mean(j)); #take averages
DUP <- data[which(data[,'type']=='DUP'),c('C1','C2','bin','j')];
DUP <- ddply(DUP, .(C1,C2),summarize,j=mean(j)); #take averages
INV <- data[which(data[,'type']=='INV'),c('C1','C2','bin','j')];
INV <- ddply(INV, .(C1,C2),summarize,j=mean(j)); #take averages

del_callers <- rev(c('True','cnMOPS','CNVnator','Delly','Hydra','Lumpy','BreakDancer'));
p1 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(DEL,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='blue',limit=c(0.0,1.0),name='Jaccard\nSimularity') +
	ylim(del_callers) + xlim(rev(del_callers)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle=90,vjust=1,size=t_size,hjust=1),axis.title.x = element_blank()) +
	theme(axis.text.y = element_text(angle=0,vjust=1,size=t_size,hjust=1),axis.title.y = element_blank()) +
	#theme(legend.position="none") +
	ggtitle("DEL");

dup_callers <- rev(c('True','cnMOPS','CNVnator','Delly','Hydra','Lumpy'));
p2 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(DUP,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='blue',limit=c(0.0,1.0),name='Jaccard\nSimularity') +
	ylim(dup_callers) + xlim(rev(dup_callers)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle=90,vjust=1,size=t_size,hjust=1),axis.title.x = element_blank()) +
	theme(axis.text.y = element_text(angle=0,vjust=1,size=t_size,hjust=1),axis.title.y = element_blank()) +
	#theme(legend.position="none") +
	ggtitle("DUP");

inv_callers <- rev(c('True','Delly','Hydra','Lumpy','BreakDancer'));
p3 <- ggplot(aes(x=C2,y=C1,fill=value),data=melt(INV,id.var=c('C1','C2'))) +
	geom_tile(color='black') +
	scale_fill_gradient2(low='white',high='blue',limit=c(0.0,1.0),name='Average\nBase Pair\nSimularity') +
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

lgnd<-g_legend(p1 + theme(legend.key.size=unit(0.5,'in')) + 
               theme(legend.position="right",legend.direction="horizontal",legend.text=element_text(angle=0,size=t_size)));
	
#svg(filename = '~/Desktop/TEMP/fusionSVU/visual/sim_matrix.svg',width=12,height=5)
pdf(paste(out_path,out_file,'.pdf',sep=''),width=out_width,height=out_height)
g <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"),
                              p3 + theme(legend.position="none"),
                              nrow=1),
                   lgnd, nrow=2,heights=c(8, 1.5))
dev.off()
       
       
