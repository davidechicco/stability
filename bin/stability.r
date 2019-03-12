#### Stability script ###############

## Notations are those in 
## "Algebraic Comparison of Partial Lists in Bioinformatics" by
## Giuseppe Jurman,  Samantha Riccadonna,  Roberto Visintainer,  Cesare Furlanello
## PLoS ONE 7(5): e36540. https://doi.org/10.1371/journal.pone.0036540

### Accessory functions #############
harmonic <- function(n){
  h=0
  for(i in 1:n)
    h= h+1/i
  return(h)
}

Delta <- function(a,b,c){
  if(c<a)
  d=b-a+1-2*c*(harmonic(b+c)-harmonic(a+c-1))
  
  if((c>=a)&(c<=b))  
  d=2*c*(2*harmonic(2*c)-harmonic(a+c-1)-harmonic(b+c)-1)+b+a-1  
    
  if(c>b)  
  d=2*c*(harmonic(b+c)-harmonic(a+c-1))-b+a-1
  
  return(d)
}

Nd<- function(a,b,c){
  t <- 0
  for(i in a:b){
    t <- t+abs(c-i)/(c+i)
  }
    return(t)
}

epsilon <- function(k,s){
  return(0.5*(s-k)*(s+k+1)*harmonic(s+k+1)+0.5*k*(k+1)*harmonic(k+1)+0.25*s*(2*k-s-1))
}

xi <- function(s){
  return(harmonic(2*s+1)*(s+0.5)**2-0.125*harmonic(s)-0.25*(2*s*s+s+1))
}

expval <- function(p){
  return((2*p+2+1/(2*p))*harmonic(2*p)-(2*p+2+1/(4*p))*harmonic(p)-(p+1.5))
}

canberra <- function(ranked_list_a,ranked_list_b,universe){
  p=length(universe)
  l1=max(ranked_list_a,na.rm=TRUE)
  l2=max(ranked_list_b,na.rm=TRUE)
  
  cbr=0
  if((l1==p) & (l2==p)){
    cbr = sum(abs(ranked_list_a-ranked_list_b)/(ranked_list_a+ranked_list_b))/expval(p)
    }else{
    if(l1>l2){
          tmp=ranked_list_a
          ranked_list_a =  ranked_list_b 
          ranked_list_b = tmp
    }
    l1=max(ranked_list_a,na.rm=TRUE)
    l2=max(ranked_list_b,na.rm=TRUE) 
    
    A = sum(is.na(ranked_list_a | ranked_list_b))/((p-l1)*(p-l2))
    idx=which(as.numeric(ranked_list_a & ranked_list_b)==1, arr.ind = TRUE)
    
    T1=0
    if(length(idx)>0){
      for(i in idx){
        Z <- abs(ranked_list_a[i]-ranked_list_b[i])/(ranked_list_a[i]+ranked_list_b[i])
        Z1 <- Delta(l2+1,p,ranked_list_a[i])/(p-l2)
        Z2 <- Delta(l1+1,p,ranked_list_b[i])/(p-l1)
        T1=T1+Z-Z1-Z2
      }
    }
 

    T2 <- 1/(p-l2)* (l1*(p-l2)-2*epsilon(p,l1) +2*epsilon(l2,l1)) + 1/(p-l1)* (l1*(p-l1)+4*epsilon(l1,l1)+2*xi(l2)-2*xi(l1)-2*epsilon(l1,l2)-2*epsilon(p,l2)+(p+l1)*(l2-l1)+l1*(l1+1)-l2*(l2+1))
    T3 <- A*(2*xi(p)-2*xi(l2)-2*epsilon(l1,p) + 2*epsilon(l1,l2)-2*epsilon(p,p) +2*epsilon(p,l2) + (p+l1)*(p-l2) +l2*(l2+1) - p*(p+1))
    
      
      cbr=(T1+T2+T3)/expval(p)
        
    }
  return(cbr)
}


stability <- function(the_user_lists,the_universe=unique(unlist(the_user_lists))){
  
  N <- length(the_user_lists)
  my_idx_lists <- list()
  for(i in 1:N)
    my_idx_lists[[i]] = match(the_universe,the_user_lists[[i]])
  
  
  dist_mat <- matrix(NA,nrow =N*(N-1)/2,ncol = 3 )
  k <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      k <- k+1
      dist_mat[k,1] <- i
      dist_mat[k,2] <- j
      dist_mat[k,3] <- canberra(my_idx_lists[[i]],my_idx_lists[[j]],the_universe)
    }
  }
  distances <- as.data.frame(dist_mat)
  names(distances) <- c("L1","L2","distance")
  return(list(stability=mean(distances$distance), all_distances=distances))
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
filename <- "./rflists.csv"
my_lists <- lapply(strsplit(readLines(filename),","),as.character)
stability(my_lists)

# $stability
# [1] 0.1799231
# 
# $all_distances
# L1 L2   distance
# 1   1  2 0.11609515
# 2   1  3 0.06012716
# 3   1  4 0.09005133
# 4   1  5 0.34282180
# 5   1  6 0.10582204
# 6   1  7 0.08553748
# 7   1  8 0.21007874
# 8   1  9 0.06729111
# 9   1 10 0.09566066
# 10  2  3 0.17615460
# 11  2  4 0.18028733
# 12  2  5 0.38999676
# 13  2  6 0.17947767
# 14  2  7 0.15040042
# 15  2  8 0.13737950
# 16  2  9 0.09005133
# 17  2 10 0.14408213
# 18  3  4 0.14957443
# 19  3  5 0.37400057
# 20  3  6 0.13706585
# 21  3  7 0.08556271
# 22  3  8 0.26897322
# 23  3  9 0.12735055
# 24  3 10 0.12709964
# 25  4  5 0.31696264
# 26  4  6 0.13128537
# 27  4  7 0.12145579
# 28  4  8 0.18215802
# 29  4  9 0.13148328
# 30  4 10 0.10582204
# 31  5  6 0.23761444
# 32  5  7 0.40227719
# 33  5  8 0.46157602
# 34  5  9 0.30155348
# 35  5 10 0.24777360
# 36  6  7 0.19093867
# 37  6  8 0.25076354
# 38  6  9 0.09011904
# 39  6 10 0.06174949
# 40  7  8 0.18350344
# 41  7  9 0.10159638
# 42  7 10 0.15562718
# 43  8  9 0.22582277
# 44  8 10 0.25148402
# 45  9 10 0.05403080

###########################  Example 2: partial lists
filename <- "./fake_lists.csv"
my_lists <- lapply(strsplit(readLines(filename),","),as.character)
stability(my_lists)

# $stability
# [1] 0.9983511
# 
# $all_distances
# L1 L2  distance
# 1    1  2 1.0184316
# 2    1  3 1.0378070
# 3    1  4 0.9769800
# 4    1  5 1.0816564
# 5    1  6 0.9117459
# 6    1  7 0.9791354
# 7    1  8 0.8296671
# 8    1  9 1.0492437
# 9    1 10 0.9059253
# 10   1 11 0.9395743
# 11   1 12 1.0542130
# 12   1 13 1.0941113
# 13   1 14 0.8617471
# 14   1 15 0.9605580
# 15   1 16 1.0550804
# 16   1 17 1.0362476
# 17   1 18 0.9926427
# 18   1 19 1.1007954
# 19   2  3 1.0143234
# 20   2  4 1.0262921
# 21   2  5 1.0309787
# 22   2  6 1.0287651
# 23   2  7 1.0214982
# 24   2  8 0.9693663
# 25   2  9 1.0262921
# 26   2 10 1.0086182
# 27   2 11 0.8362533
# 28   2 12 0.9693663
# 29   2 13 1.0268640
# 30   2 14 1.0184316
# 31   2 15 1.0285785
# 32   2 16 1.0184316
# 33   2 17 1.0255667
# 34   2 18 1.0294536
# 35   2 19 1.0287651
# 36   3  4 1.0506949
# 37   3  5 1.0345456
# 38   3  6 0.9750006
# 39   3  7 0.9197170
# 40   3  8 0.8631661
# 41   3  9 1.0589035
# 42   3 10 1.0143234
# 43   3 11 0.8389065
# 44   3 12 1.0478524
# 45   3 13 0.9007180
# 46   3 14 1.0378070
# 47   3 15 1.0665282
# 48   3 16 0.9585673
# 49   3 17 0.9690157
# 50   3 18 1.0523689
# 51   3 19 0.9133385
# 52   4  5 1.0851157
# 53   4  6 1.0105590
# 54   4  7 0.8796280
# 55   4  8 0.9387961
# 56   4  9 0.8428740
# 57   4 10 0.9867117
# 58   4 11 1.1572517
# 59   4 12 0.8586039
# 60   4 13 0.9806473
# 61   4 14 0.9676339
# 62   4 15 1.0805984
# 63   4 16 0.8656107
# 64   4 17 0.8474488
# 65   4 18 0.8885011
# 66   4 19 0.9242019
# 67   5  6 0.9605556
# 68   5  7 1.0912194
# 69   5  8 1.1187368
# 70   5  9 1.0047399
# 71   5 10 1.0248576
# 72   5 11 1.0502555
# 73   5 12 1.2545774
# 74   5 13 1.1885382
# 75   5 14 0.9680460
# 76   5 15 0.8267826
# 77   5 16 1.0283395
# 78   5 17 0.9680569
# 79   5 18 0.9821988
# 80   5 19 0.8587484
# 81   6  7 0.9591221
# 82   6  8 1.0328301
# 83   6  9 0.9768867
# 84   6 10 0.8362533
# 85   6 11 0.8678936
# 86   6 12 1.2158316
# 87   6 13 1.1251266
# 88   6 14 1.0709102
# 89   6 15 0.9328820
# 90   6 16 0.8993616
# 91   6 17 1.0034930
# 92   6 18 0.9742515
# 93   6 19 1.0558436
# 94   7  8 0.9626307
# 95   7  9 0.9817158
# 96   7 10 1.0214982
# 97   7 11 0.8583765
# 98   7 12 0.9951447
# 99   7 13 0.9020791
# 100  7 14 0.9791354
# 101  7 15 1.1332564
# 102  7 16 0.9117548
# 103  7 17 1.0248831
# 104  7 18 1.1518271
# 105  7 19 0.7091310
# 106  8  9 1.0609242
# 107  8 10 1.0288233
# 108  8 11 0.9480299
# 109  8 12 0.9311424
# 110  8 13 1.0148976
# 111  8 14 0.9092785
# 112  8 15 1.1334400
# 113  8 16 0.9652077
# 114  8 17 0.9070335
# 115  8 18 1.0029303
# 116  8 19 1.1633228
# 117  9 10 1.0293079
# 118  9 11 1.1325622
# 119  9 12 1.0725884
# 120  9 13 0.8932665
# 121  9 14 1.0898166
# 122  9 15 0.7391595
# 123  9 16 1.0679078
# 124  9 17 0.9997513
# 125  9 18 0.8860675
# 126  9 19 1.0176096
# 127 10 11 1.0087019
# 128 10 12 1.0087019
# 129 10 13 1.0268640
# 130 10 14 1.0184316
# 131 10 15 0.8362533
# 132 10 16 1.0184316
# 133 10 17 1.0255667
# 134 10 18 1.0087019
# 135 10 19 1.0287651
# 136 11 12 1.0367198
# 137 11 13 1.0649161
# 138 11 14 1.0079355
# 139 11 15 1.1996242
# 140 11 16 1.0014693
# 141 11 17 1.0601590
# 142 11 18 0.9465324
# 143 11 19 0.9674092
# 144 12 13 0.9349497
# 145 12 14 1.0476176
# 146 12 15 1.0368663
# 147 12 16 0.9175461
# 148 12 17 1.1513162
# 149 12 18 1.0287289
# 150 12 19 1.0958675
# 151 13 14 1.0553161
# 152 13 15 0.9549269
# 153 13 16 1.0722263
# 154 13 17 1.0818015
# 155 13 18 1.0468178
# 156 13 19 1.0132065
# 157 14 15 1.0564496
# 158 14 16 0.9692848
# 159 14 17 0.9036568
# 160 14 18 1.0496623
# 161 14 19 0.9310582
# 162 15 16 1.0937325
# 163 15 17 1.0545682
# 164 15 18 0.9400843
# 165 15 19 1.2363651
# 166 16 17 0.7549037
# 167 16 18 1.0399917
# 168 16 19 0.9822070
# 169 17 18 1.0387952
# 170 17 19 0.9383377
# 171 18 19 1.1166949