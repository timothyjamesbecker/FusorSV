#Caller Jaccard Simularity to G1K.P3 by Type and Bin

#Rscript cmd_parser.R commands...
options(warn=-1);
cmd_args <- commandArgs();
scriptpath <- gsub('BinnedCallerPerformance.R','',gsub('--file=','',cmd_args[4]));
source(paste(scriptpath,'cmd_parser.R',sep=''));
source(paste(scriptpath,'plot_utils.R',sep=''));

#p_size  <- 3;
#p_alpha <- 0.7;
#m_size  <- 5;
#t_size  <- 12;

#'/fusorsv_result/visual/callers_tbj.tsv'
in_path    <- args['in_path'][[1]];
in_file_split <- strsplit(in_path,'/')[[1]];
in_file    <- in_file_split[length(in_file_split)];
out_path   <- args['out_path'][[1]];
out_file<-strsplit(in_file,'.tsv')[[1]]
#grap specfic listings here-----------------------
out_width  <- args['out_width'][[1]];
out_height <- args['out_height'][[1]];
p_size     <- args['p_size'][[1]];
p_alpha    <- args['p_alpha'][[1]];
m_size     <- args['m_size'][[1]];
t_size     <- args['t_size'][[1]];
cat('input commands are',in_path,out_path);

#red=1,royalblue=2,orange=3,lightgreen=4,green=5,pink=6,lightorange=7,lightblue=8,
#light_purple=9,dark_purple=10,yellow=11,black=12,brown=13
cs <- c('#ff6666','#1f78b4','#ff7f00','#b2df8a',
        '#33a02c','#e7298a','#999999','#a6cee3',
        '#cab2d6','#6a3d9a','#e6ab02','#000000','#b15928');
colors <- c(cs[3],cs[10],cs[4],cs[5],cs[11],
            cs[2],cs[1],cs[6],cs[8],cs[13],cs[9],cs[7]);
col_names <- c('caller','type','bin','j');
data <- read.table(in_path,header=T,sep='	',col.names=col_names,stringsAsFactors=F);
types <- c('DEL','DUP','INV');

#subset callers if needed for plots
data <- data[which(data[,'caller']!='GATK' & data[,'caller']!='Tigra' & data[,'caller']!='G1K-P3' & data[,'caller']!='fusorSV' & data[,'caller']!='MetaSV' & data[,'caller']!='alpha' & data[,'caller']!='Pindel'),] ;

#clean up formatting for ggplot
DEL <- data[which(data[,'type']=='DEL'),col_names];
DEL[,'caller'] <- as.factor(DEL[,'caller']);
DEL[,'svlen'] <- as.factor(DEL[,'bin']);
DEL[,'svlen'] <- factor(DEL[,'svlen'],levels=bin_order(levels(DEL[,'svlen'])))

DUP <- data[which(data[,'type']=='DUP'),col_names];
DUP[,'caller'] <- as.factor(DUP[,'caller']);
DUP[,'svlen'] <- as.factor(DUP[,'bin']);
DUP[,'svlen'] <- factor(DUP[,'svlen'],levels=bin_order(levels(DUP[,'svlen'])))

INV <- data[which(data[,'type']=='INV'),col_names];
INV[,'caller'] <- as.factor(INV[,'caller']);
INV[,'svlen'] <- as.factor(INV[,'bin']);
INV[,'svlen'] <- factor(INV[,'svlen'],levels=bin_order(levels(INV[,'svlen'])))

p1 <- ggplot(aes(x=svlen,y=j,color=caller,group=caller,fill=caller),data=DEL) +
	geom_line(size=p_size,alpha=p_alpha) +
	geom_point(size=2.0*p_size,alpha=p_alpha/2.0) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.6)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
	theme(axis.text.x = element_text(angle=45,vjust=1,size=t_size,hjust=1)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("DEL")

p2 <- ggplot(aes(x=svlen,y=j,color=caller,group=caller,fill=caller),data=DUP) +
	geom_line(size=p_size,alpha=p_alpha) +
	geom_point(size=2.0*p_size,alpha=p_alpha/2.0) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.3)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
	theme(axis.text.x = element_text(angle=45,vjust=1,size=t_size,hjust=1)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("DUP")

p3 <- ggplot(aes(x=svlen,y=j,color=caller,group=caller,fill=caller),data=INV) +
	geom_line(size=p_size,alpha=p_alpha) +
	geom_point(size=2.0*p_size,alpha=p_alpha/2.0) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.55)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
	theme(axis.text.x = element_text(angle=45,vjust=1,size=t_size,hjust=1)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("INV")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lgnd<-g_legend(p1 + theme(legend.position="bottom",legend.direction="horizontal"));

pdf(paste(out_path,out_file,'.pdf',sep=''),width=out_width,height=out_height)
g <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"),
                              p3 + theme(legend.position="none"),
                              nrow=3),
                   lgnd, ncol=1,heights=c(8,1.5))
dev.off()


