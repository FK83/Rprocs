##############################################
# qdiff - Compute difference (in quarters) between two quarterly dates
# Input:
# - d1, quarterly date, STRING in the format "YYYYQx", where x is between 1 and 4
# - d2, quarterly date, same format as d1
# Output: 
# - difference (in quarters) between d1 and d2, positive if d2 is later in time than d1
##############################################
qdiff <- function(d1,d2){
  y1 <- substr(d1,1,4)
  y2 <- substr(d2,1,4)
  q1 <- substr(d1,6,6) 
  q2 <- substr(d2,6,6)
  return((as.numeric(y2)-as.numeric(y1))*4+as.numeric(q2)-as.numeric(q1))
}

##############################################
# plusq - Find the date 'q' quarters after d1
# Input:
# - d1, quarterly date, STRING in the format "YYYYQx", where x is between 1 and 4
# - q, integer (possibly negative)
# Output: 
# - date `q' quarters after d1, same format as d1
##############################################
plusq <- function(d1,q){
  y1 <- as.numeric(substr(d1,1,4))
  q1 <- as.numeric(substr(d1,6,6))
  aux <- y1 + (q1-1)/4 + q/4
  y2 <- floor(aux)
  q2 <- (aux-y2)*4+1
  return(paste(y2,"Q",q2,sep=""))
}

##############################################
# sel.complete - Select matrix rows with complete entries 
# Input:
# - dat, matrix/data frame
# Output:
# - matrix containing complete rows of dat (if any)
##############################################
sel.complete <- function(dat){
  sel <- apply(dat,1,function(z) !any(is.na(z)))
  dat[sel,] 
}

##############################################
# agr - Compute growth rate from x to y 
# Input:
# - y, numeric scalar or vector
# - x, same format as y
# - lg, indicator, TRUE for log growth rate, else geometric growth (defaults to TRUE)
# Output:
# - growth rate from x to y
##############################################
agr <- function(y,x,lg=TRUE){
  if (lg == TRUE){
    100*(((y/x)^4)-1)
  } else {
    400*(log(y)-log(x))
  }
}

##############################################
# merge.ts - Merge up to three time series 
# Input:
# - d1, data frame with two columns, 1st col holds date ("2009Q4"), 2nd col holds numeric value of time series
# - d2, same format as d1
# - Optional: Third data frame d3, character vector of labels for merged data frame ("labs")
# Output:
# - data frame with three (four) columns, containing values at dates that are contained in ALL input data frames
##############################################
merge.ts <- function(d1,d2,d3=NULL,labs=NULL){
aux <- intersect(d1[,1],d2[,1])  
dat <- sel.complete(cbind(d1[d1[,1]%in%aux,],d2[d2[,1]%in%aux,2]))
if (!is.null(d3)){
  aux <- intersect(dat[,1],d3[,1])
  dat <- sel.complete(cbind(dat,d3[d3[,1]%in%aux,2]))
}  
if (is.null(labs)){
  names(dat) <- paste0("V",1:dim(dat)[2])
} else {
  names(dat) <- labs
}
dat  
}

##############################################
# read.table2, write.table2 - Wrappers for read.table, write.table, fixing my usual settings
##############################################
read.table2 <- function(s) read.table(s,sep=",",header=TRUE,stringsAsFactors=FALSE)  
write.table2 <- function(dat,s) write.table(dat,s,sep=",",row.names=FALSE,col.names=TRUE)