library('ggplot2');
library('grid');
library('gridExtra');
library('plyr');
library('reshape2');

get_sample_name<-function(full_path){
	s <- strsplit(full_path,'/')[[1]]
	s[length(s)-1]
}

get_unique<-function(data,cols){
  data[!duplicated(data[,cols]),]
}

find_hull <- function(df,x,y) { df[chull(df[,x],df[,y]),] }

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