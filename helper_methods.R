### calculate confusion matrix ###
getConfusionMatrix <- function(rf) {
  
  tbl = table(predict(rf), rf$y)
  class.error = vector()
  
  for (i in 1:nrow(tbl)) {
    rowSum = sum(tbl[i,])
    accurate = diag(tbl)[i]
    error = rowSum - accurate
    
    class.error[i] = error / rowSum
  }   
  return(cbind(tbl, class.error))
}

twoclassAUC <- function(rf, data, labels, group)
{
  require(pROC)
  roc(labels, predict(rf, data, type="prob")[,2])
}

### calculate a pseudo AUC value ###
multiclassAUC <- function(rf, labels){
  library(HandTill2001)
  library(randomForest)
  auc(multcap(
    response = labels
    , predicted = as.matrix(predict(rf,type="prob"))
  ))
}

### plot the most significant features using a random forest result ###
plotTopFeatures <- function(rf, threshold=1, top=NULL, lookup=NA, return.df=T, alias2sym=T, main=NULL){
  library(ggplot2)
  
  if(is.null(top)){
    data.rrf.top <- rf$importance[order(rf$importance, decreasing=T),]
    data.rrf.top <- data.rrf.top[data.rrf.top > threshold]
  }
  else{
    data.rrf.top <- rf$importance[order(rf$importance, decreasing=T),][1:top]
  }

  if(!is.na(lookup)){
    gene.names <- lookup[names(data.rrf.top)]
  }
  else if(alias2sym)
  {
    library(limma)
    gene.names <- alias2Symbol(names(data.rrf.top))
  }
  else{
    gene.names <- names(data.rrf.top)
  }
  
  data.rrf.top <- data.frame(Gene=gene.names, MeanDecreaseGini=data.rrf.top)
  data.rrf.top <- transform(data.rrf.top, Gene = reorder(as.character(Gene), -MeanDecreaseGini))
  
  p <- qplot(x=Gene, y=MeanDecreaseGini, width=.8, data=data.rrf.top, geom="bar", stat="summary", fun.y="mean") + theme(axis.text.x=element_text(angle=45, hjust=1))
  if(!is.null(main)){
    p <- p + ggtitle(main)
  }
  
  if(return.df){
    print(p)
    return(data.rrf.top[with(data.rrf.top, order(-MeanDecreaseGini)),])
  }
  else return(p)
}

get.PAM50 <- function(){
  library(genefu) #bioconductor http://www.bioconductor.org/packages/2.13/bioc/manuals/genefu/man/genefu.pdf
  
  data(pam50)
  PAM50.genes <- pam50$centroids.map$probe
  #replace two gene names
  setdiff(PAM50.genes, colnames(brca.gene.expr)) == 2
  PAM50.genes[PAM50.genes%in%c("CDCA1")] <- "NUF2"
  PAM50.genes[PAM50.genes%in%c("KNTC2")] <- "NDC80"
  rm(pam50)
  return(PAM50.genes)
}

### multiple plots with with ggplots ###
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
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
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}