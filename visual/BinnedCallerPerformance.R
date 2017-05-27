#Caller Jaccard Simularity to G1K.P3 by Type and Bin

library('ggplot2');
library('grid');
library('gridExtra');
library('plyr');
library('reshape2');

bin_order <- function(bin){
	#looks like c('1-50','50-100','1.0K-1.5K',etc)
	m <- as.data.frame(matrix(0,nrow=length(bin),ncol=1)); #extract value
	rownames(m) <- bin;
	colnames(m) <- 'bin';
	for(b in 1:length(bin)){
		upper <- strsplit(bin[b],'-')[[1]][2];
		if(grepl('K',upper)){ m[b,1] <- as.numeric(strsplit(upper,'K')[[1]])*1E3;  }
		if(grepl('M',upper)){ m[b,1] <- as.numeric(strsplit(upper,'M')[[1]])*1E6;  }
		if(grepl('G',upper)){ m[b,1] <- as.numeric(strsplit(upper,'G')[[1]])*1E9;  }
		if(grepl('T',upper)){ m[b,1] <- as.numeric(strsplit(upper,'T')[[1]])*1E12; }
		if(!(grepl('K',upper) || grepl('M',upper) || grepl('G',upper) || grepl('T',upper))){
			m[b,1] <- as.numeric(upper);
		}
	}
	bin[order(m)]
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

p_size  <- 3;
p_alpha <- 0.7;
m_size  <- 5;
t_size  <- 12;

#red=1,royalblue=2,orange=3,lightgreen=4,green=5,pink=6,lightorange=7,lightblue=8,
#light_purple=9,dark_purple=10,yellow=11,black=12,brown=13
cs <- c('#ff6666','#1f78b4','#ff7f00','#b2df8a',
        '#33a02c','#e7298a','#999999','#a6cee3',
        '#cab2d6','#6a3d9a','#e6ab02','#000000','#b15928');
colors <- c(cs[3],cs[10],cs[4],cs[5],cs[11],
            cs[2],cs[1],cs[6],cs[8],cs[13],cs[9],cs[7]);
in_path <- '/Users/tbecker/Desktop/meta_caller_NP4_result/visual/callers_tbj.0123456789ABCDEF101112131415161718191A.tsv'
#in_path <- '~/Desktop/TEMP/varsim_R1_big/visual/callers_tbj.tsv'
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
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.6)) +
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
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.3)) +
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
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.55)) +
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

pdf('~/Desktop/caller_perform_T_B_g1k_50X_27.pdf',width=14,height=14)
#(filename = '~/Documents/CourseWork/16_2016_Spring/meta_caller_overlap/visual/callers_by_type_and_bin.png',width=800,height=600)
g <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"),
                              p3 + theme(legend.position="none"),
                              nrow=3),
                   lgnd, ncol=1,heights=c(8,1.5))
dev.off()


