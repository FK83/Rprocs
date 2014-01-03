##############################################
# fcstplot - make time series plot of realizations and Gaussian density forecasts
# Input: 
# - dat, data frame of dim. T x 3, rows represent time periods 
# - dat$m: mean forecasts
# - dat$s: predictive standard deviations
# - dat$rlz: (ex-post) realizations
# - ax, T x 1 vector of (numeric) forecast target dates (default: integers 1 to T)
# - ci, scalar between 0.5 and 1, sign. level for confidence bands (shown as shaded areas)
# - lx, character, label for X axis (defaults to "time")
# - ly, character, label for Y axis (defaults to "y")
# Output:
# ggplot object, time series plot 
##############################################
fcstplot <- function(dat,ax=NULL,ci=0.95,lx="time",ly="y"){
  require(ggplot2)
  # Create numeric time variable (quarterly)
  if (is.null(ax)) ax <- 1:dim(dat)[1]
  dat$tdate <- ax
  n <- dim(dat)[1]
  # Compute prediction intervals of Gaussian forecasts
  lb <- qnorm(rep((1-ci)*0.5,n),dat$m,dat$sd)
  ub <- qnorm(rep(1-(1-ci)*0.5,n),dat$m,dat$sd)
  # Auxiliary data frame
  auxf <- data.frame(lb=lb,ub=ub,rlz=dat$m,tdate=dat$tdate)
  # Plot
  ggplot(dat,aes(tdate,rlz)) + geom_point() + geom_line(data=auxf) +
    geom_ribbon(data=auxf,aes(ymin=lb,ymax=ub),alpha=0.3) + xlab(lx) + ylab(ly) 
}

##############################################
# ggmts - make gg plot of multiple time series
# Input: 
# - ax, numeric time axis (defaults to integers starting from 1)
# - dat, T x m matrix of data (rows represent time, cols represent variables)
# - legs, character vector of length m, labels (defaults to NULL)
# - type, either "lines" (default) or "symb"
# - lsz, size of lines/symbols (can be left unspecified)
# - lx, character, label for X axis (defaults to "time") 
# - ly, character, label for Y axis (defaults to "y")
# - ltit, character, title (defaults to empty)
# - savename, character, location to save graph (defaults to null, graph not saved)
# Output:
# ggplot object, time series plot
##############################################
ggmts <- function(dat,ax=NULL,legs=NULL,type="lines",lsz=NULL,lx="time",ly="y",ltit="",savename=NULL){
  # Packages
  require(ggplot2)
  require(reshape)
  # Axis
  if (is.null(ax)) ax <- 1:dim(dat)[1]
  aux <- data.frame(date=ax,dat)
  # Legends
  if (is.null(legs)){
    names(aux)[2:ncol(aux)] <- paste("V",1:(ncol(aux)-1),sep="")
  } else {
    names(aux)[2:ncol(aux)] <- legs
  }
  # Convert data to long format
  aux_long <- melt(aux, id="date")  
  # Size of lines/symbols
  if (is.null(lsz)){
    if (type == "lines"){
      lsz <- 1.2
    } else {
      lsz <- 3
    }
  }
  # Either lines or symbols
  if (type == "lines"){
  g <- ggplot(data=aux_long,aes(x=date, y=value, colour=variable)) +
     xlab(lx) + ylab(ly) + ggtitle(ltit) + scale_colour_brewer(palette="Set1") + theme(legend.title=element_blank()) + geom_line(size=lsz)  
  # Remove legend?
  if (is.null(legs)) g <- g + guides(colour=FALSE) 
  } else {
  g <- ggplot(data=aux_long,aes(x=date, y=value, shape=variable)) +
    xlab(lx) + ylab(ly) + ggtitle(ltit) + scale_colour_brewer(palette="Set1") + theme(legend.title=element_blank()) + geom_point(size=lsz)
  # Remove legend?
  if (is.null(legs)) g <- g + guides(shape=FALSE) 
  }

  # Save?
  if (!is.null(savename)) ggsave(savename,g)  
  g
}