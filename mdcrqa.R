mdcrqa <- function(ts1,ts2,emb,del,norm,rad,dline,vline,zscore, metric = "euclidean") {
  # This function performs Multidimensional Cross-Recurrence Quantification Analysis (MdCRQA)
  # 
  # Inputs:
  #
  #  ts1 = a M x N matrix, where N is the number of indivdual time series that are to be correlated and M is the number of data points in each time series
  #  ts2 = a M x N matrix, where N is the number of indivdual time series that are to be correlated and M is the number of data points in each time series
  #  emb = embedding dimension parameter, where emb = 1 mean no embedding, emb = 2 mean one creation of a surrogate copy of the original data etc.
  #  del = delay parameter, where del is the number of data points to be skipped for each additional embedding
  #  norm = norm parameter, where phase-space can be normalized by Euclidean distance 'euc', Maximum distance 'max', Minimum distance 'min', nor not 'non'
  #  rad = radius parameter, were rad is the distance that defines the neighborhood around a coordinate in phase-space
  #  dline = selection of minimum diagonal line length for the computation of DET, ADL and DENTR. The default value is 2 (i.e., including all lines that consist of at least two diagonally adjacent recurrence points)
  #  dline = selection of minimum vertical line length for the computation of LAM, AVL and VENTR. The default value (and minimum) value is 2 (i.e., including all lines that consist of at least two vertically adjacent recurrence points)
  #  zscore = indicats, whether the data (i.e., the different columns of ts1 and ts2, being the different signals or dimensions of a signal) should be z-scored before performing mdcrqa, where 0 = no z-scoring, and 1 = z-score columns of ts1 and ts2. The default value is zscore = 0.
  #  metric = indicates the type of distance measure to apply, default euclidean but see help(cdist) for more options
  #
  # Outputs:
  #
  #  SizeCRP = size of the cross-recurrence plot
  #  REC = percent recurrence
  #  DET = percent determinism
  #  ADL = average diagonal line length
  #  MDL = maximum diagonal line length
  #  DENTR = diagonal line entropy
  #  LAM = percent laminarity
  #  AVL = average vertical line length
  #  MVL = maximum vertical line length
  #  VENTR = vertical line entropy
  #  DIM = dimensions of input data (i.e., number of time series)
  #  EMB = embedding dimension
  #  DEL = delay
  #  NORM = phase-space normalization
  #  RAD = radius
  #  CRP = recurrence plot as sparse matrix
  #
  # Reference:
  #
  #  Wallot, S. (2018). Multidimensional Cross-Recurrence Quantification
  #  Analysis (MdCRQA) - a method for quantifying correlation between
  #  multivariate time-series. ???
  #
  # Version:
  #
  # v1.0, 13. April 2018
  # by Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
  #
  # v1.1, 09 April 2025
  # - fixed typoes in the description of the function
  # 
  # The author does not give any warranty whatsoever on the reliability of
  # the results obtained by this software.
  
  # Load cdist from rdist and Matrix form Matrix
  library(rdist)
  library(Matrix)
  
  # check input variables
  if (exists("ts1")) {
  } else {
    stop("No data has been specified for ts1.")
  }
  
  # check input variables
  if (exists("ts2")) {
  } else {
    stop("No data has been specified for ts2.")
  }
  
  # check machting of input data 
  if (dim(as.matrix(ts1))[1] == dim(as.matrix(ts2))[1]) {
  } else {
    stop("Input data (ts1 and ts2) do not match in terms of number of data points (number of rows).")
  }
  
  if (dim(as.matrix(ts1))[2] == dim(as.matrix(ts2))[2]) {
  } else {
    stop("Input data (ts1 and ts2) do not match in terms of dimensionality (number of columns).")
  }
  
  if (exists("emb")) {
  } else {
    emb <- 1
  }
  
  if (exists("del")) {
  } else {
    del <- 1
  }
  
  if (dim(as.matrix(ts1))[1] > (emb-1)*del) {
  } else {
    stop("Insufficient number of data points to embedd time-series (number of data points < (emb-1)*del).")
  }  
  
  if (exists("norm")) {
  } else {
    norm <- "non"
  }
  
  if (exists("rad")) {
  } else {
    rad <- 1
  }
  
  if (exists("dline")) {
    if (dline < 2) {
      dline <- 2
    }
  } else {
    dline <- 2
  }
  
  if (exists("vline")) {
    if (vline < 2) {
      vline <- 2
    }
  } else {
    vline <- 2
  }
  
  if (exists("zscore")) {
  } else {
    zscore <- 0
  }
  
  if (zscore == 0)  {
  } else {
    ts1 <- scale(ts1)
    ts2 <- scale(ts2)
  }
  
  dims <- dim(ts1)[2]
  
  # store parameter settings
  PARAMETERS <- c(dims, emb, del, rad, norm, dline, vline)
  
  # check whether embedding needs to be performed
  if (emb > 1) {
    newLength <- dim(ts1)[1] - (emb-1)*del
    tempTs1 <- ts1[1:newLength,]
    for (i in seq(2,emb)) {
      tempTs1 <- cbind(tempTs1,ts1[(1+(del*(i-1))):(newLength+del*(i-1)),])
    }
    ts1 <- tempTs1
    rm(tempTs1)
    tempTs2 <- ts2[1:newLength,]
    for (i in seq(2,emb)) {
      tempTs2 <- cbind(tempTs2,ts2[(1+(del*(i-1))):(newLength+del*(i-1)),])
    }
    ts2 <- tempTs2
    rm(tempTs2)
  }
  
  # create distance matrix
  CRP <- as.matrix(cdist(ts1,ts2, metric = metric))
  
  CRP <- abs(CRP)*-1
  
  # store dimensions of the matrix
  nrw = nrow(CRP)
  ncl = ncol(CRP)
  CRP[is.na(CRP)] <- 0
  
  # apply norm, radius and threshold matrix
  if (norm == "euc") {
    CRP <- CRP/abs(sum(CRP)/(dim(CRP)[1]^2-dim(CRP)[1]))
  } else if (norm == "min") {
    CRP <- CRP/abs(min(CRP))
  } else if (norm == "max") {
    CRP <- CRP/abs(max(CRP))
  } else {
  }
  
  CRP <- CRP + rad
  CRP <- ifelse(CRP < 0, 0, 1) ## the ifelse statement takes some time
  
  # Computing the line counts
  numrecurs = length(which(CRP == TRUE));
  
  # check here if there any recurrences at all
  if (numrecurs > 0){ ## there is nothing
    
    RR = (numrecurs/((nrw*ncl)))*100;
    
    # calculate diagonal and vertical line distributions
    diagLine <- split(CRP, row(CRP) - col(CRP)) 
    diaglines = sort(as.numeric(unlist(lapply(diagLine, 
                                              function(x) {
                                                runs <-  rle(x); 
                                                lx   = runs$lengths[runs$values == 1]
                                                return (lx)}))),
                     decreasing = T)
    
    # delete line counts less than the minimum diagonal.
    dcrit = which(diaglines < dline)
    
    if (length(dcrit) > 0){ diaglines = diaglines[-dcrit]}
    
    if(length(diaglines) != 0){
      
      DNRLINES = length(diaglines) 
      MDL      = max(diaglines) # extract the max length of diag
      ADL      = mean(diaglines)
      
      tabled = as.data.frame(table(diaglines))
      
      total = sum(tabled$Freq)       
      p = tabled$Freq/total
      
      # remove zero probability
      del = which(p == 0 )
      if (length(del) > 0) {
        p = p[-del]
      }
      
      DENTR  = - sum(p*log(p))    
      DET = sum(diaglines)/numrecurs*100
      
    } else {
    	DNRLINES = 0
    	MDL = 0
    	ADL = 0
    	DENTR = 0
    	DET = 0
    }
    
    
    
    vertLine <- split(CRP,col(CRP)) 
    vertlines = sort(as.numeric(unlist(lapply(vertLine, 
                                              function(x) {
                                                runs <-  rle(x); 
                                                lx   = runs$lengths[runs$values == 1]
                                                ## we could return also the black lines by just changing the above parameter to 0
                                                return (lx)}))),
                     decreasing = T)
    
    # delete line counts less than the minimum vertical.
    vcrit = which(vertlines < vline)
    
    if (length(vcrit) > 0){ vertlines = vertlines[-vcrit]}
    
    if(length(vertlines) != 0){
      
      VNRLINES = length(vertlines) # extract the number of vertical lines
      MVL      = max(vertlines)  # longest vertical line
      AVL      = mean(vertlines) # this is the TT
      
      tabled   = as.data.frame(table(vertlines))
      
      total    = sum(tabled$Freq)       
      p        = tabled$Freq/total
      
      #remove zero probability
      del = which(p == 0 )
      if (length(del) > 0) {
        p = p[-del]
      }
      
      VENTR  = - sum(p*log(p))    
      LAM = sum(vertlines)/numrecurs*100 ## laminarity
      
    } else {

    	VNRLINES = 0
    	MVL = 0
    	AVL = 0
    	VENTR = 0
    	LAM = 0
    }

    
  } else { # in case you find no recurrence
    DNRLINES = 0; MDL = 0; ADL = 0; DET = NA;
    DENTR    = NA;
    VNRLINES = 0; MVL = 0; AVL = 0; LAM = NA;
    VENTR    = NA;
    RP = NA
  }
  
  # make CRP a sparse matrix
  CRP = Matrix(CRP, sparse = TRUE)
  
  results = list(RR = RR, DET = DET, 
                 DNRLINES = DNRLINES, MDL = MDL,
                 ADL = ADL, DENTR = DENTR, 
                 LAM = LAM, AVL = AVL, MVL = MVL, VNRLINES = VNRLINES, 
                 VENTR = VENTR,
                 RP = CRP)
  
  return (results)
  
}
