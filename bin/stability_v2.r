#### Stability script ###############

## Notations are those in 
## "Algebraic Comparison of Partial Lists in Bioinformatics" by
## Giuseppe Jurman,  Samantha Riccadonna,  Roberto Visintainer,  Cesare Furlanello
## PLoS ONE 7(5): e36540. https://doi.org/10.1371/journal.pone.0036540
##
## edited by Davide Chicco <davidechicco@davidechicco.it>
##

# function that prints two decimals of a number
dec_two <- function(x) {
  return (format(round(x, 2), nsmall = 2));
}

# Function that checks if a list contain all the same elements
allSameElements <- function(x) length(unique(x)) == 1

# fromFileToLists function
fromFileToLists <- function(listCsvFileName, verboseFlag=FALSE){

   if (verboseFlag) cat("[INPUT FILE] Reading the list in the ", listCsvFileName, " file (each row in the csv file is a list)\n", sep="")

  theseLists <- lapply(strsplit(readLines(listCsvFileName),","),as.character)  

  return(theseLists)
}



# harmonic function
harmonic <- function(n){
  h=0
  for(i in 1:n)
    h= h+1/i
  return(h)
}

# delta function
Delta <- function(a,b,c){
  if(c<a)
  d=b-a+1-2*c*(harmonic(b+c)-harmonic(a+c-1))
  
  if((c>=a)&(c<=b))  
  d=2*c*(2*harmonic(2*c)-harmonic(a+c-1)-harmonic(b+c)-1)+b+a-1  
    
  if(c>b)  
  d=2*c*(harmonic(b+c)-harmonic(a+c-1))-b+a-1
  
  return(d)
}

# nd function
Nd<- function(a,b,c){
  t <- 0
  for(i in a:b){
    t <- t+abs(c-i)/(c+i)
  }
    return(t)
}

# epsilon
epsilon <- function(k,s){
  return(0.5*(s-k)*(s+k+1)*harmonic(s+k+1)+0.5*k*(k+1)*harmonic(k+1)+0.25*s*(2*k-s-1))
}

# xi
xi <- function(s){
  return(harmonic(2*s+1)*(s+0.5)**2-0.125*harmonic(s)-0.25*(2*s*s+s+1))
}

# expval
expval <- function(p){
  return((2*p+2+1/(2*p))*harmonic(2*p)-(2*p+2+1/(4*p))*harmonic(p)-(p+1.5))
}

# Canberra list
canberra <- function(ranked_list_left, ranked_list_right, universe){
  p=length(universe)
  maxValueLeftList=max(ranked_list_left, na.rm=TRUE)
  maxValueRightList=max(ranked_list_right, na.rm=TRUE)
  
  cbr=0
  if((maxValueLeftList==p) & (maxValueRightList==p)){
    cbr = sum(abs(ranked_list_left-ranked_list_right)/(ranked_list_left+ranked_list_right))/expval(p)
    }else{
    if(maxValueLeftList>maxValueRightList){
          tmp=ranked_list_left
          ranked_list_left =  ranked_list_right 
          ranked_list_right = tmp
    }
    maxValueLeftList=max(ranked_list_left,na.rm=TRUE)
    maxValueRightList=max(ranked_list_right,na.rm=TRUE) 
    
    A = sum(is.na(ranked_list_left | ranked_list_right))/((p-maxValueLeftList)*(p-maxValueRightList))
    idx=which(as.numeric(ranked_list_left & ranked_list_right)==1, arr.ind = TRUE)
    
    T1=0
    if(length(idx)>0){
      for(i in idx){
        Z <- abs(ranked_list_left[i]-ranked_list_right[i])/(ranked_list_left[i]+ranked_list_right[i])
        Z1 <- Delta(maxValueRightList+1,p,ranked_list_left[i])/(p-maxValueRightList)
        Z2 <- Delta(maxValueLeftList+1,p,ranked_list_right[i])/(p-maxValueLeftList)
        T1=T1+Z-Z1-Z2
      }
    }
 

    T2 <- 1/(p-maxValueRightList)* (maxValueLeftList*(p-maxValueRightList)-2*epsilon(p,maxValueLeftList) +2*epsilon(maxValueRightList,maxValueLeftList)) + 1/(p-maxValueLeftList)* (maxValueLeftList*(p-maxValueLeftList)+4*epsilon(maxValueLeftList,maxValueLeftList)+2*xi(maxValueRightList)-2*xi(maxValueLeftList)-2*epsilon(maxValueLeftList,maxValueRightList)-2*epsilon(p,maxValueRightList)+(p+maxValueLeftList)*(maxValueRightList-maxValueLeftList)+maxValueLeftList*(maxValueLeftList+1)-maxValueRightList*(maxValueRightList+1))
    T3 <- A*(2*xi(p)-2*xi(maxValueRightList)-2*epsilon(maxValueLeftList,p) + 2*epsilon(maxValueLeftList,maxValueRightList)-2*epsilon(p,p) +2*epsilon(p,maxValueRightList) + (p+maxValueLeftList)*(p-maxValueRightList) +maxValueRightList*(maxValueRightList+1) - p*(p+1))
    
      
      cbr=(T1+T2+T3)/expval(p)
        
    }
  return(cbr)
}


# stability function
stability <- function(the_user_lists, verboseFlag=FALSE, the_universe=unique(unlist(the_user_lists))){

    lengthIndex <- 1
    listOfLengths <- as.numeric(array((summary(the_user_lists))[,lengthIndex]))
    
    if (verboseFlag)cat("[NUMBER OF LISTS] The program has read ", length(the_user_lists), " lists\n", sep="")
    
    if ( allSameElements(listOfLengths) == TRUE ) {
    
        if (verboseFlag) cat("[LENGTH OF EACH LIST] All the lists have the same length, which is ", listOfLengths[1], "\n", sep="")
    
    } else {
        if (verboseFlag) cat("[LENGTH OF EACH LIST] The lists have different lengths, which are the following: ")
        if (verboseFlag) cat(listOfLengths, "\n", sep=" ")     
    }

  
  LISTS_START_INDEX <- 1
  
  N <- length(the_user_lists)
  my_idx_lists <- list()
  for(i in LISTS_START_INDEX:N)
    my_idx_lists[[i]] = match(the_universe,the_user_lists[[i]])
  
  MATRIX_COL_NUM <- 3

  dist_mat <- matrix(NA, nrow =N*(N-1)/2, ncol = MATRIX_COL_NUM)

  LEFT_LIST_INDEX <-  1
  RIGHT_LIST_INDEX <-  2
  CANBERRA_LIST_INDEX <- 3
  
  k <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      k <- k+1
      dist_mat[k, LEFT_LIST_INDEX] <- i
      dist_mat[k, RIGHT_LIST_INDEX] <- j
      dist_mat[k, CANBERRA_LIST_INDEX] <- canberra(my_idx_lists[[i]],my_idx_lists[[j]],the_universe)
    }
  }
  
  distances <- as.data.frame(dist_mat)
  names(distances) <- c("L_left","L_right","Canberra_distance")
  
  INTERVAL <- 2
  # normalization: minimum 0.00, and maximum 1.00
  resultStability <- mean(distances$Canberra_distance) / INTERVAL
  
  if (verboseFlag) { 
        cat("[STABILITY] The stability of the lists is ", dec_two(resultStability), " in the [0, 1] interval\n", sep="") 
        cat("[STABILITY] where 0.00 means identical lists, and 1.00 means maximum instability\n")
  }
    
  return(list(stability=resultStability, all_distances=distances))
  
  }

#### Instructions:
## Stability(LISTS, ALL_FEATURES)
## LISTS: a list of vectors, each vector consisting of a ranked set of features; these ranked sets can be "complete" if they rank 
##        all the available features or "partial", if they include only a subset of all possible features
## ALL_FEATURES: the complete set of all features of the experiment; as the default, the set of all elements of all components of LISTS are
##               taken as the full set unique(unlist(LISTS))


########### Examples
##
## LISTS are read from a csv file; each row of the file becomes an element of LISTS
##
##
###########################  Example 1: complete lists
listCsvFileName <- "../data/rflists.csv"
# my_lists <- lapply(strsplit(readLines(listCsvFileName),","),as.character)
thisVerboseFlag <- TRUE
if(thisVerboseFlag) cat("\n")
my_lists <- fromFileToLists(listCsvFileName, thisVerboseFlag)
thisStabilityCompleteLists <- stability(my_lists, thisVerboseFlag)$stability

# ###########################  Example 2: partial lists
listCsvFileName <- "../data/fake_lists.csv"
if(thisVerboseFlag) cat("\n")
my_lists <- fromFileToLists(listCsvFileName, thisVerboseFlag)
thisStabilityPartialLists <- stability(my_lists, thisVerboseFlag)$stability

# ###########################  Example 3: simple lists
listCsvFileName <- "../data/simple_lists.csv"
if(thisVerboseFlag) cat("\n")
my_lists <- fromFileToLists(listCsvFileName, thisVerboseFlag)
thisStabilitySimplelLists <- stability(my_lists, thisVerboseFlag)$stability

# ###########################  Example 4: simple lists
tryList <- list("aaa", "fff")
if(thisVerboseFlag) cat("\n")
thisStabilitySimplelLists <- stability(tryList, thisVerboseFlag)$stability

