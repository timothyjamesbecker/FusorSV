#Caller To Target Perfomance by Sample
#take in the cross_fold_stat matrix as a .tsv and do the dot plot

#Rscript cmd_parser.R commands...
options(warn=-1);
cmd_args <- commandArgs();
scriptpath <- gsub('PerformanceBySample.R','',gsub('--file=','',cmd_args[4]));
source(paste(scriptpath,'cmd_parser.R',sep=''));
source(paste(scriptpath,'plot_utils.R',sep=''));

#p_size  <- 1;
#p_alpha <- 0.6;
#m_size  <- 5;
#t_size  <- 16;

#'/fusorsv_result/visual/cross_fold_stats.HEX.tsv'
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
t_size     <- args['t_size'][[1]];
m_size     <- args['m_size'][[1]];
cross_fold <- args['cross_fold'][[1]];
cat('input commands are',in_path,out_path);

#red=1,royalblue=2,orange=3,lightgreen=4,green=5,pink=6,lightorange=7,lightblue=8,
#light_purple=9,dark_purple=10,yellow=11,black=12,brown=13
cs <- c('#ff6666','#1f78b4','#ff7f00','#b2df8a',
        '#33a02c','#e7298a','#999999','#a6cee3',
        '#cab2d6','#6a3d9a','#e6ab02','#000000','#b15928');
colors <- c(cs[3],cs[10],cs[4],cs[5],cs[11],cs[12],
            cs[2],cs[1],cs[6],cs[13],cs[8],cs[9],cs[7]);
#loop over all the possible cross_fold leave one out validations and load into the same plot...         
#cross_fold <- F;
#cross_fold_stats_folder <- '~/Desktop/TEMP/kfold_plots/';
cross_fold_stats_path <- in_path;
col_names <- c('sample','caller','type','prec','rec','f1','j','n','m','l_mu','r_mu','l_sd','r_sd');

if(cross_fold){
	cat('\nreading input directory: ',cross_fold_stats_folder,'\n')
	loocv <- list.files(path=cross_fold_stats_folder,pattern='*.tsv');
	data <- data.frame();
	for(i in 1:length(loocv)){
		data <- rbind(data,read.table(paste(cross_fold_stats_folder,loocv[i],sep=''),
	                                  header=T,sep='	',col.names=col_names,stringsAsFactors=F))
	}
	cat('finished reading directory: ',cross_fold_stats_folder, '\n')
}else{
	cat('\nreading input file: ', cross_fold_stats_path,'\n')
	data <- read.table(cross_fold_stats_path,header=T,sep='	',col.names=col_names,stringsAsFactors=F);
	cat('\nfinished reading file: ', cross_fold_stats_path,'\n')	
}
types <- c('DEL','DUP','INV'); #dig out the interesting ones here
#subset callers if needed for plots
data <- data[which(data[,'caller']!='GATK' & data[,'caller']!='fusorSV' & data[,'caller']!='MetaSV' & data[,'caller']!='Tigra' & data[,'caller']!='Pindel'),];

DEL <- data[which(data[,'type']=='DEL'),c('caller','prec','rec','f1','j','n','m')];
DEL[,'nIm'] <- DEL[,'n']*DEL[,'prec'];
DEL[,'dn'] <- DEL[,'n']-DEL[,'nIm'];
DEL[,'dm'] <- DEL[,'m']-DEL[,'nIm'];
DEL[,'caller'] <- as.factor(DEL[,'caller']);

DUP <- data[which(data[,'type']=='DUP'),c('caller','prec','rec','f1','j','n','m')];
DUP[,'nIm'] <- DUP[,'n']*DUP[,'prec'];
DUP[,'dn'] <- DUP[,'n']-DUP[,'nIm'];
DUP[,'dm'] <- DUP[,'m']-DUP[,'nIm'];
DUP[,'caller'] <- as.factor(DUP[,'caller']);

INV <- data[which(data[,'type']=='INV'),c('caller','prec','rec','f1','j','n','m')];
INV[,'nIm'] <- INV[,'n']*INV[,'prec'];
INV[,'dn'] <- INV[,'n']-INV[,'nIm'];
INV[,'dm'] <- INV[,'m']-INV[,'nIm'];
INV[,'caller'] <- as.factor(INV[,'caller']);

DEL_N_M <- get_unique(DEL,c('caller','dn','dm'));
DEL_N_M <- merge(DEL_N_M,aggregate(cbind(mean.dn=dn,mean.dm=dm)~caller,DEL_N_M,mean),by="caller");
q1 <- ggplot(aes(x=dn,y=dm,color=caller,fill=caller),data=DEL_N_M) + 
	geom_point(size=p_size+2,alpha=0.9) +
	geom_vline(xintercept = mean(DEL[,'n'])) +
	geom_point(aes(x=mean.dn,y=mean.dm),size=m_size*2,alpha=p_alpha/length(levels(DEL[,'caller']))) +
	scale_color_manual(values=colors) +
	scale_fill_manual(values=colors) +
	#annotate("text", label = "average P3", x = mean(DEL[,'n'])*0.98, 
	#         y = mean(DEL[,'m'])/4, size = p_size+3, colour = "black",angle=90) +
	#scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(DEL[,'n']))) + 
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(DEL[,'m']/2))) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=90,vjust=0,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=0,vjust=0,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("DEL EVENTS");

DEL_P_R <- get_unique(DEL,c('caller','prec','rec'));
DEL_P_R <- merge(DEL_P_R,aggregate(cbind(mean.rec=rec,mean.prec=prec)~caller,DEL_P_R,mean),by="caller");
p1 <- ggplot(aes(x=rec,y=prec,color=caller,fill=caller),data=DEL_P_R) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(DEL_P_R),"caller",find_hull,'rec','prec'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.rec,y=mean.prec),size=m_size,alpha=p_alpha/length(levels(DEL_P_R[,'caller']))) +
  	#geom_segment(aes(x=mean.rec, y=mean.prec, xend=rec, yend=prec),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.9)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.8)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("DEL");

DUP_N_M <- get_unique(DUP,c('caller','dn','dm'));
DUP_N_M <- merge(DUP_N_M,aggregate(cbind(mean.dn=dn,mean.dm=dm)~caller,DUP_N_M,mean),by="caller");
q2 <- ggplot(aes(x=dn,y=dm,color=caller,fill=caller),data=DUP_N_M) + 
	geom_point(size=p_size+2,alpha=0.9) +
	geom_vline(xintercept = mean(DUP[,'n'])) +
	geom_point(aes(x=mean.dn,y=mean.dm),size=m_size*2,alpha=p_alpha/length(levels(DUP[,'caller']))) +
	scale_color_manual(values=colors) +
	scale_fill_manual(values=colors) +
	#annotate("text", label = "average P3", x = mean(DUP[,'n'])*0.98, 
	#         y = mean(DUP[,'m'])/4, size = p_size+3, colour = "black",angle=90) +
	#scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(DUP[,'n']))) + 
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(DUP[,'m'])/2)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=90,vjust=0,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=0,vjust=0,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("DUP EVENTS");

DUP_P_R <- get_unique(DUP,c('caller','prec','rec'));
DUP_P_R <- merge(DUP_P_R,aggregate(cbind(mean.rec=rec,mean.prec=prec)~caller,DUP_P_R,mean),by="caller");	
p2 <- ggplot(aes(x=rec,y=prec,color=caller,fill=caller),data=DUP_P_R) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(DUP_P_R),"caller",find_hull,'rec','prec'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.rec,y=mean.prec),size=m_size,alpha=p_alpha/length(levels(DUP_P_R[,'caller']))) +
  	#geom_segment(aes(x=mean.rec, y=mean.prec, xend=rec, yend=prec),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) + 
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.4)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.35)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="none") +
	ggtitle("DUP");

INV_N_M <- get_unique(INV,c('caller','dn','dm'));
INV_N_M <- merge(INV_N_M,aggregate(cbind(mean.dn=dn,mean.dm=dm)~caller,INV_N_M,mean),by="caller");
q3 <- ggplot(aes(x=dn,y=dm,color=caller,fill=caller),data=INV_N_M) + 
	geom_point(size=p_size+2,alpha=0.9) +
	geom_vline(xintercept = mean(INV[,'n'])) +
	#annotate("text", label = "average P3", x = mean(INV[,'n'])*0.95, 
	#         y = mean(INV[,'m'])/4, size = p_size+3, colour = "black",angle=90) +
	geom_point(aes(x=mean.dn,y=mean.dm),size=m_size*2,alpha=p_alpha/length(levels(INV[,'caller']))) +
	scale_color_manual(values=colors) +
	scale_fill_manual(values=colors) +
	#scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(INV[,'n']))) + 
    scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,mean(DUP[,'m'])/2)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=90,vjust=0,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=0,vjust=0,size=t_size,hjust=0)) +
	#theme(legend.position="bottom") +
	ggtitle("INV EVENTS");

INV_P_R <- get_unique(INV,c('caller','prec','rec'));
INV_P_R <- merge(INV_P_R,aggregate(cbind(mean.rec=rec,mean.prec=prec)~caller,INV_P_R,mean),by="caller");		
p3 <- ggplot(aes(x=rec,y=prec,color=caller,fill=caller),data=INV_P_R) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(INV_P_R),"caller",find_hull,'rec','prec'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.rec,y=mean.prec),size=m_size,alpha=p_alpha/length(levels(INV_P_R[,'caller']))) +
  	#geom_segment(aes(x=mean.rec, y=mean.prec, xend=rec, yend=prec),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) + 
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.65)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,1.0)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0)) +
	#theme(legend.position="none") +
	ggtitle("INV");

DEL_J_F <- get_unique(DEL,c('caller','j','f1'));
DEL_J_F <- merge(DEL_J_F,aggregate(cbind(mean.j=j,mean.f1=f1)~caller,DEL_J_F,mean),by="caller");
p4 <- ggplot(aes(x=f1,y=j,color=caller,fill=caller),data=DEL_J_F) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(DEL_J_F),"caller",find_hull,'f1','j'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.f1,y=mean.j),size=m_size,alpha=p_alpha/length(levels(DEL_J_F[,'caller']))) +
  	#geom_segment(aes(x=mean.f1, y=mean.j, xend=f1, yend=j),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.7)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.6)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0));
	#theme(legend.position="none") #+
	#ggtitle("DEL");

DUP_J_F <- get_unique(DUP,c('caller','j','f1'));
DUP_J_F <- merge(DUP_J_F,aggregate(cbind(mean.j=j,mean.f1=f1)~caller,DUP_J_F,mean),by="caller");
p5 <- ggplot(aes(x=f1,y=j,color=caller,fill=caller),data=DUP_J_F) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(DUP_J_F),"caller",find_hull,'f1','j'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.f1,y=mean.j),size=m_size,alpha=p_alpha/length(levels(DUP_J_F[,'caller']))) +
  	#geom_segment(aes(x=mean.f1, y=mean.j, xend=f1, yend=j),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.3)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.3)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0));
	#theme(legend.position="none") #+
	#ggtitle("DUP");

INV_J_F <- get_unique(INV,c('caller','j','f1'));
INV_J_F <- merge(INV_J_F,aggregate(cbind(mean.j=j,mean.f1=f1)~caller,INV_J_F,mean),by="caller");	
p6 <- ggplot(aes(x=f1,y=j,color=caller,fill=caller),data=INV_J_F) +
	geom_point(size=p_size,alpha=0.9) +
	geom_polygon(data = ddply(na.omit(INV_J_F),"caller",find_hull,'f1','j'),alpha=p_alpha/1.5) +
	geom_point(aes(x=mean.f1,y=mean.j),size=m_size,alpha=p_alpha/length(levels(INV_J_F[,'caller']))) +
	#geom_segment(aes(x=mean.f1, y=mean.j, xend=f1, yend=j),alpha=p_alpha) +
	scale_color_manual(values=colors) + 
	scale_fill_manual(values=colors) +
    #scale_x_continuous(expand=c(0.0,0.0),limits=c(0.0,0.65)) + 
    #scale_y_continuous(expand=c(0.0,0.0),limits=c(0.0,0.95)) +
    theme_minimal() +
    theme(panel.background = element_rect(colour = "black")) +
	theme(axis.title = element_text(size=t_size)) +
	theme(axis.text.x = element_text(angle=0,vjust=1,size=t_size,hjust=0)) +
	theme(axis.text.y = element_text(angle=90,vjust=1,size=t_size,hjust=0));
	#theme(legend.position="none") #+
	#ggtitle("INV");

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lgnd<-g_legend(p1 + theme(legend.position="right",legend.direction="vertical")+
                    theme(legend.key.size=unit(0.8,'in'),legend.text=element_text(angle=0,size=t_size)));

pdf(paste(out_path,in_file,'.pdf',sep=''),width=out_width,height=out_height)
# #svg(filename = '~/Desktop/TEMP/k3_fold_performance_by_sample.svg',width=6,height=8)
g1 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               p3 + theme(legend.position="none"),
                               p4 + theme(legend.position="none"),
                               p5 + theme(legend.position="none"),
                               p6 + theme(legend.position="none"),
                               q1 + theme(legend.position="none"),
                               q2 + theme(legend.position="none"),
                               q3 + theme(legend.position="none"),
                               nrow=3),
                   lgnd, ncol=2,widths=c(8, 1.5))
dev.off()


# lgnd<-g_legend(q1 + theme(legend.position="right",legend.direction="vertical")+
                    # theme(legend.key.size=unit(0.8,'in'),legend.text=element_text(angle=0,size=t_size)));
# #png('~/Desktop/fusorSV_paper/counts_g1k_P3_27.png',1400,800)
# g2 <- grid.arrange(arrangeGrob(q1 + theme(legend.position="none"),
                               # q2 + theme(legend.position="none"),
                               # q3 + theme(legend.position="none"),
                              # nrow=1),
                   # lgnd, ncol=2,widths=c(8, 1.5))
# #dev.off()








