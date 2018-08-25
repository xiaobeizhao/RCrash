
##' ************************************************************************
##' These are chunks of codes for the R Crash lecture in
##' Bioinformatics of high throughput analyses (BioHTA2012)
##' 
##'
##' Xiaobei Zhao (xiaobei@binf.ku.dk)
##' Bioinformatics Centre, University of Copenhagen
##' April 23, 2012
##'
##' options(width=64)
##' ************************************************************************


##' ========================================================================
##' Part I
##' Introduction to R
##' ========================================================================


##' ------------------------------------------------------------------------
##' Getting help!
##' ------------------------------------------------------------------------

help(mean)
help("mean")                      
?mean                             # an alternative 

help.search("heatmap")
??heatmap                         # an alternative 

RSiteSearch("heatmap")


##' ------------------------------------------------------------------------
##' Installing package
##' ------------------------------------------------------------------------

install.packages("gplots")
## --- Please select a CRAN mirror for use in this session ---

source("http://bioconductor.org/biocLite.R")
biocLite("ChromHeatMap") # DO NOT RUN; it might take some time for downloading  


##' ------------------------------------------------------------------------
##' Loading package
##' ------------------------------------------------------------------------

library(gplots)
require(gplots)



##' ========================================================================
##' Part II
##' R basics
##' ========================================================================

##' ------------------------------------------------------------------------
##' Working directory
##' ------------------------------------------------------------------------

getwd()
## list directory contents
dir()

## store current directory
oldDir <- getwd()
## set wd
setwd("TMP")
## check whether the specified wd is set
getwd()

## switch back
setwd(oldDir)
getwd()


##' ------------------------------------------------------------------------
##' Expressions and assignments
##' ------------------------------------------------------------------------

5+6

x <- 5+6
x
y <- x+7
y


##' ------------------------------------------------------------------------
##' Assignment Operators
##' ------------------------------------------------------------------------

help(assignOps,package="base")


##' ------------------------------------------------------------------------
##' Types of objects
##' ------------------------------------------------------------------------

?integer
?numeric
?logical
?character
?factor

##' ------------------------------------------------------------------------
##' Types of objects - integer and numeric
##' ------------------------------------------------------------------------

typeof(1)
is.numeric(1)

typeof(as.integer(1))
is.numeric(as.integer(1))

typeof(1:3)
typeof(matrix(1:12,ncol=4))


as.integer(12)
## as.single(12)
as.double(12)



##' ------------------------------------------------------------------------
##' Types of objects - logical
##' ------------------------------------------------------------------------

as.logical("TRUE")
as.logical("True")
as.logical("true")
as.logical("T")
as.logical("t")                   # wont't work
as.logical("r")                   # wont't work

as.logical(1)
as.logical(10)
as.logical(0)
as.logical(-1)
as.logical(-1.5)



##' ------------------------------------------------------------------------
##' Types of objects - character
##' ------------------------------------------------------------------------

## ## create a character vecter with empty strings
## character(length=10)    

## Or with anything within a pair of quotation marks " and '
c("abc",'123', "Hello, R!", "")

## a invalid string
""Thank you", she said."

"\"Thank you\", she said."

'"Thank you", she said.'
"Breakfast at Tiffany's"

as.character(TRUE)      # convert from a logical
as.character(12)        # convert from a numeric



##' ------------------------------------------------------------------------
##' Types of objects - factor
##' ------------------------------------------------------------------------

## x0 <- c("low", "medium", "high")
## set.seed(123)
## x <- sample(x0,6,replace=TRUE)

x <- c("low","high","medium","high","high","low")
x

y <- factor(x)
y
sort(y)

y <- factor(x, levels=c("low", "medium", "high"))
y
sort(y)


##' ------------------------------------------------------------------------
##' Constants
##' ------------------------------------------------------------------------

help(Constants,package="base")
pi
print(pi, digits=16)

letters


##' ------------------------------------------------------------------------
##' Special values
##' ------------------------------------------------------------------------

## Get help
?NA
?NULL
?Inf
?NaN

##
x <- NA
is.na(x)
length(x)

y <- NULL
is.null(y)
length(y)

1/0
-1/0
0/0



##' ========================================================================
##' Part II (Cont'd)
##' R basics - Data structures
##' ========================================================================

##' ------------------------------------------------------------------------
##' Vector
##' ------------------------------------------------------------------------

letters
is.vector(letters)
length(letters)

## a vector of length one
pi
is.vector(pi)
length(pi)


## Create a vector
1:10
10:1
seq(1,10)
seq(1,10,2)
seq(1,10,length.out=3)


c(2,4,10,"z","a")
x <- c(21:25,letters[1:5],NA)
x
typeof(x)
mode(x)


## Indexing
x <- seq(21,30)
x
x[1]                  # the first element
x[length(x)]          # the last element
x[2:4]                # the 2nd, 3rd and 4th elements
x[c(1,3,5)]           # the 1st, 3rd and 5th elements

## replacement
x
x[2] <- 99
x
x[3] <- NA
x

## add element(s)
x <- 1:5
x[10] <- 99
x
x <- letters[1:3]
x[10] <- "hi!"
x


##' ------------------------------------------------------------------------
##' List
##' ------------------------------------------------------------------------

## Create a list
x <- list(21:25, letters[1:6])                # unnamed list
x
is.list(x)                                    # check if it is a `list'
is.vector(x)                                  # check if it is a `vector'
length(x)                                     # get the length
str(x)                                        # display the structure

x <- list(number=21:25, letter=letters[1:6])  # named list
x
names(x)                                      # get the names

## Create a new list by combining elements using c
x <- list(21:25, letters[1:6])
x
y <- c(x,pi)
y
str(y)

## try this?
y2 <- c(x,pi,1:3)
str(y2)

## Indexing
x <- list(number=21:25, letter=letters[1:6], unit=c("C","F")) 
x
## Indexing: to access an element
x[[1]]
x$unit

## Indexing: to obtain a list of element(s)
x[c(1,3)]


## replacement
x <- list(letter=letters[1:6], unit=c("C","F")) 
x
x$letter <- letters[7:1]
x
x$letter[1] <- "z"
x

## Add element(s) by [[]]
x <- list(letter=letters[1:6], unit=c("C","F"))
x
length(x)
x[[5]] <- 21:25

## Add element(s) by $
x <- list(letter=letters[1:6], unit=c("C","F"))
x
x$id <- 1:3
x
x[[3]]


##' ------------------------------------------------------------------------
##' Matrix
##' ------------------------------------------------------------------------

## matrix
matrix(1:12,ncol=4)               # fill in elements by column
matrix(1:12,ncol=4,byrow=TRUE)    # fill in elements by row

## Creat a matrix and save to a variable
x <- matrix(1:12,ncol=4)
x

## Get size
dim(x)                            # retrieve the dimension
nrow(x)                           # return the number of rows
ncol(x)                           # return the number of columns
length(x)                         # total number of elements

## Set size
dim(x) <- c(4,3)
x

## Get names
rownames(x)                       # retrieve row names 
colnames(x)                       # retrieve column names

## Set names
rownames(x) <- 1:nrow(x)          # set the row names
colnames(x) <- letters[1:ncol(x)] # set the column names
rownames(x)                       
colnames(x)                       
x

## Indexing by row/column numbers
x[c(1,3),]                        # get rows by row numbers
x[,2:3]                           # get columns by column numbers

## Indexing by row/column names
x["2",]                           # by row names
x[,c("a","c")]                    # by column names

## Indexing by both row and column
x[1:2,2:3] 

## Indexing by conditions,
## select rows by values in the 1st ("a") column
x
x[x[,"a"]>=2,]

## select columns by values in the 2nd row
x[,x[2,]>5]

## select given combined conditions on both a row and a column
x[x[,"a"]>2,x[1,]>=5]


## Replacement
## Creat a matrix with missing values
y <- matrix(ncol=4,nrow=4)
y

## Replace the diagonal with zeros, by indices
y[1,1] <- y[2,2] <- y[3,3] <- y[4,4] <- 0
## diag(y) <- 0                     # an alternative
y

## Replace the lower triangle with numbers, by conditions
## require(gdata)                   # load an add-on package for `lowerTriangle()' 
## lowerTriangle(y) <- 1:length(lowerTriangle(y))
row(y)>col(y)
y[row(y)>col(y)] <- 1:sum(row(y)>col(y))
y

## Add row(s) using \Sobj{rbind}
x <- matrix(1:8,ncol=2)
x

rbind(x,11:12)
rbind(x,11:12,21:22)              # add multiple rows at once


## Add column(s) using \Sobj{cbind}
x <- matrix(1:8,ncol=2)
x

cbind(x,11:14)
cbind(x,11:14,21:24)              # add multiple columns at once


## Add a column when the vector is shorter
x <- matrix(1:8,ncol=2)
x

cbind(x,11)
cbind(x,c(21,22))
cbind(x,c(31,32,33))


## Combine two matrices
x <- matrix(1:12,ncol=3)
x
y4 <- matrix(16:21,ncol=3)
y4
rbind(x,y4)

## ... if corresponding dimensions do not match?
cbind(x,y4)

## add cells by new indices?
x <- matrix(1:12,ncol=3)
x[5,6] <- 99
x[5,] <- 99
x[,6] <- 99

x <- matrix(1:12,ncol=3)
x[20] <- 99
x


## multi-dimentional matrix/array (Not required for the course.)
a <- array(1:600,c(2,3,10,10))



##' ------------------------------------------------------------------------
##' data.frame
##' ------------------------------------------------------------------------

## Create a data.frame
x <- data.frame(1:3,letters[1:3],seq(1,10,length.out=3))
x
x2 <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3))
x2
str(x2)
x3 <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                 stringsAsFactors=FALSE)
x3
str(x3)

as.list(x3)


## Size and names
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
dim(x)                            # retrieve the dimension
length(x)                         # number of columns
rownames(x)
colnames(x)

## Indexing
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x[1:2,]
x[,3]
x[,c("id","measure")]

x$measure

## Replacement
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x$measure <- x$measure*2
x
x[2,2:3] <- c("f",50)
x

## Add row(s) by rbind
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x0 <- c(0,"control",0)
y <- rbind(x0, x)
y
class(y)

## Add column(s) by cbind
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x0 <- c("BRIC","BINF","BINF")
y <- cbind(x,where=x0)
y
str(y)

## Add column(s) by $
y$who <- c("xb","xb","as")
y

## Add row(s) by new indices?
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)

x[4,] <- c(4,NA,NA)
x[6,] <- c(6,NA,NA)
x


## Add column(s) by new indices?
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x[,4] <- NA
x
x[,6] <- NA


## Add cell(s) by new indices?
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3),
                stringsAsFactors=FALSE)
x
x[4,1] <- 4
x[6,1] <- 6
x[1,4] <- NA
x[1,6] <- NA



## clean the workspace by removing all objects
save.image(file="workspace01.RData") # save current workspace
rm(list=ls())                        # clean the workspace
ls()


## create a data.frame and attach it
x <- data.frame(id=1:3,label=letters[1:3],measure=seq(1,10,length.out=3))
x
head(ls(pos=2))
label                                # columns are not accessible directly.
attach(x)
label                                # columns are accessible now!
id
head(ls(pos=2))

##
detach(x)
head(ls(pos=2))


##' ========================================================================
##' Part II (Cont'd)
##' R basics - Oprations
##' ========================================================================


##' ------------------------------------------------------------------------
##' Arithmetic
##' ------------------------------------------------------------------------

help(Arithmetic,package="base")


##' ------------------------------------------------------------------------
##' Conditions
##' ------------------------------------------------------------------------

help(Logic,package="base")


##' ------------------------------------------------------------------------
##' Vector operation
##' ------------------------------------------------------------------------

v <- 1:5
v
v**2                              
v%%2
v + 1
v + c(1,3)                             # c(1,3) is recycled
v*c(1,3) 


toupper(c("a","b","c"))
sqrt(c(1,2,3))

##options(width=72)
paste("chr",c(1:22,"X","Y"),sep="")    # Concatenate Strings



##' ------------------------------------------------------------------------
##' applys
##' ------------------------------------------------------------------------

## Apply over vector by \Sobj{vapply}
x <- 1:10
sapply(x,Fibonacci)

## Apply over list by \Sobj{lapply}
x <- list(a=1:5,b=letters[1:10])
x
lapply(x,length)
length(x)


## Apply over matrix by \Sobj{apply}
x <- matrix(1:12,ncol=4)
x
apply(x,1,max)
apply(x,2,max)
apply(x,c(1,2),max)



##' ------------------------------------------------------------------------
##' indexing
##' ------------------------------------------------------------------------

## Indexing by indixes/positions
x <- data.frame(a=letters[1:6],b=7:12)
x
x[c(1,3,5),]                      # select 1st, 3rd and 5th rows

## Indexing by conditions
x[(1:nrow(x))%%2==0,]             # select even rows
x[x[,2]%%3 == 0,]                 # select rows if the elements in the 2nd column
                                  # is a multiple of 3.



##' ------------------------------------------------------------------------
##' Sorting
##' ------------------------------------------------------------------------

##
mydata <- data.frame(label=letters[1:6],score=c(7, 13, 19, 11, 17, 5))
mydata
attach(mydata)
sort(score)
sort(score, decreasing=TRUE)          # sort in decreasing order
detach(mydata)

mydata
attach(mydata)
sort.list(score)                      # original indices of rearranged "score" 
label[sort.list(score)]               # rearrange "label" by "score"

mydata[sort.list(score),]
detach(mydata)

## A graphical representation of sorting-related functions
## x0 <- c(5,7,11,13,17,19)
## set.seed(2012)
## x <- sample(x0)
## x

v.ori <- c(7, 13, 19, 11, 17, 5)
i.ori1 <- 1:length(v.ori)
i.ori1
i.new1 <- rank(v.ori)
i.new1

sort(v.ori)
v.new <- sort(v.ori) 
v.new
v.ori[sort.list(v.ori)]

i.ori2 <- sort.list(v.ori)
i.ori2
i.new2 <- 1:length(v.new)
i.new2


v.ori[sort.list(v.ori)] == v.new
v.new[rank(v.ori)] == v.ori



##' ------------------------------------------------------------------------
##' Tabulation
##' ------------------------------------------------------------------------

## set.seed(2012)
## n <- 8
## x <- data.frame(id=1:n,
##                 gender=sample(c("F","M"),n,TRUE),
##                 grade=sample(1:3,n,TRUE),
##                 specialty=sample(c("bio","math"),n,TRUE))


x <- data.frame(id=1:8,
                gender=c("F", "M", "F", "M", "M", "M", "M", "M"),
                grade=c(1, 2, 2, 3, 2, 1, 2, 1)
                )
x

table(x$gender)
table(x$grade,x$gender)

mytable <- table(x$grade,x$gender)
mytable
dimnames(mytable)         # a 2D table!
mytable["2","M"]          # retrieve the male at grade 2

as.matrix(mytable)
as.data.frame(mytable)


x$specialty <- c("math", "bio", "bio", "math", "math", "math", "math", "bio")
mytable2 <- table(x$grade,x$gender,x$specialty)
mytable2
dimnames(mytable2)        # a 3D table!
mytable2[,,"math"]         # retrieve "math" specialists
mytable2[,"F","math"]      # retrieve female "math" specialists
mytable2["3","M","math"]     # retrieve male "math" specialists at grade 3




##' ========================================================================
##' Part II (Cont'd)
##' R basics - Control flows
##' ========================================================================


##' ------------------------------------------------------------------------
##' Logic-control flow
##' ------------------------------------------------------------------------

## x <- readline("Please enter an integer: ")
## x <- as.numeric(x)
## if (x<0) {
##   print("A negtive number!")
## } else if (x==0) {
##   print("A zero!")
## } else if (x%%2) {
##   print("A positive odd number!")
## } else {
##   print("A positive even number!")
## }

##
mycoin <- function(x){
  if (x<0 | x>1) {
    print("Not a probability!")
  } else if (x>=0 & x<=0.5) {
    print("A head!")
  } else {
    print("A tail!")
  } 
}

mycoin(-3)
mycoin(0)
mycoin(0.2)
mycoin(0.5)
mycoin(0.7)
mycoin(1)
mycoin(1.5)


##' ------------------------------------------------------------------------
##' Loop-control flow
##' ------------------------------------------------------------------------

n <- 10                         # set the number of elements required
v <- numeric(n)                 # initialize a vector with length "n"
v 
v[1] <- 1                       # initialize the 1st element in "v"
for (i in 2:length(v)){         
  v[i] <- v[i-1]*2
}
print(v)

## There are at least two ways to looping
## variant 1 - loop over the elements in "v"
v <- 1:3
for (e in v){
  print(e*2)
}

## variant 2 - loop over the index of the elements in "v"
for (i in 1:length(v)){
  print(v[i]*2)
}



##' ========================================================================
##' Part II (Cont'd)
##' R basics - Function
##' ========================================================================


Fibonacci <- function(n,x1=0,x2=1){
  if(n<=0){
    stop("`n' must be positive.")
  }
  
  v <- numeric(n)                 # initialize a vector with length "n"
  v[1] <- x1                      # set the initial seeds
  v[2] <- x2
  
  if(n<=2){
    return(v[1:n])
  } else {
    for (i in 3:n){
      v[i] <- v[i-1]+v[i-2]
    }
  }
  return(v)
}


Fibonacci(10)                                 # using default
Fibonacci(10,5,6)                             # by positional arguments
Fibonacci(x1=5,x2=6,n=10)                     # by named arguments



##' ========================================================================
##' Part III
##' A case study using biological data
##' ========================================================================


##' ------------------------------------------------------------------------
##' Part III (Cont'd)
##' Data handling
##' ------------------------------------------------------------------------

## x <- iris
## colnames(x)
## colnames(x) <- tolower(colnames(x))


x <- read.table("/Users/xiaobei/Project/RCrash/iris.txt",header=TRUE,sep="\t")
str(x)
head(x)

write.table(x,file="iris.txt",sep="\t",
            row.names=FALSE,col.names=TRUE,
            quote=FALSE)


## x <- read.table("iris.txt",header=TRUE,sep="\t",
##                 stringsAsFactors=FALSE)
## str(x)


##' ------------------------------------------------------------------------
##' Part III (Cont'd)
##' Case studies - Visualization
##'
##' fv=`ls *.pdf`; for f in ${fv[*]}; do pdf2eps $f; done
##' ------------------------------------------------------------------------

attach(x)
plot(sepal.width,sepal.length,main="Fisher's Iris Data")
dev.copy2pdf(file="iris-plot.pdf",width=6,height=6)

mycolor <- c("red","green3","blue")
plot(sepal.width,sepal.length,main="Fisher's Iris Data",
     col=mycolor[species])
dev.copy2pdf(file="iris-color.pdf",width=6,height=6)

mypch <- 0:2
plot(sepal.width,sepal.length,main="Fisher's Iris Data",
     col=mycolor[species],
     pch=mypch[species])
dev.copy2pdf(file="iris-pch.pdf",width=6,height=6)


legend("topright",legend=levels(species),
       col=mycolor,pch=mypch)

dev.copy2pdf(file="iris-legend.pdf",width=6,height=6)

##
require(plyr)
mymeans <- ddply(.data=x,
                  .variables=c("species"),
                  .fun=summarise,
                  sepal.width=mean(sepal.width),
                  sepal.length=mean(sepal.length)
                  )
mymeans
mypch2 <- 15:17
points(x=mymeans$sepal.width,y=mymeans$sepal.length,
       col="black",pch=mypch2)
dev.copy2pdf(file="iris-points.pdf",width=6,height=6)

##

par(mfrow=c(2,2))
for (e in levels(species)){
  plot(sepal.width[species==e],sepal.length[species==e],main=e,
       xlab="sepal.width",ylab="sepal.length",
       xlim=range(sepal.width),ylim=range(sepal.length))
}

dev.copy2pdf(file="iris-subplots.pdf",width=8,height=8)


## Does petal length have a normal distribution? Let's see ...
dev.new()
qqnorm(petal.length)

dev.copy2pdf(file="iris-qqnorm.pdf",width=6,height=6)


##
barplot(mymeans$sepal.length,names=mymeans$species,
        ylab="Mean sepal length", main="Fisher's Iris Data")
dev.copy2pdf(file="iris-barplot.pdf",width=6,height=6)


## Make a histogram
hist(sepal.length)
dev.copy2pdf(file="iris-hist.pdf",width=6,height=6)

binCount <- 10
mybreaks <- seq(min(sepal.length),max(sepal.length),length.out=binCount+1)
hist(sepal.length,breaks=mybreaks)
dev.copy2pdf(file="iris-hist2.pdf",width=6,height=6)

## add a density line
d <- density(sepal.length)        # determine densities
myhist <- hist(sepal.length,breaks=mybreaks)
hist(sepal.length,breaks=mybreaks,
     xlim=range(d$x),ylim=range(myhist$density,range(d$y)),
     probability=TRUE)
lines(d, col="red")
dev.copy2pdf(file="iris-lines.pdf",width=6,height=6)


## boxplot
boxplot(sepal.length,sepal.width,petal.length,petal.width,
        xlab="Flower characteristics",ylab="Centimeter",
        names = c("Sepal length","Sepal width","Petal length","Petal width"),
        main="A descriptive statistics of flower characteristics")
dev.copy2pdf(file="iris-boxplot.pdf",width=8,height=8)

detach(x)

## lattice

library(lattice)
xyplot(sepal.length~sepal.width|species,data=x,main="Fisher's Iris Data")
dev.copy2pdf(file="iris-xyplot.pdf",width=8,height=8)


## ggplot2
library(ggplot2)
ggplot(data=x,aes(x=sepal.width,y=sepal.length)) +
  geom_point() +
  facet_wrap(~species,ncol=2) +
  opts(title="Fisher's Iris Data")
dev.copy2pdf(file="iris-ggplot.pdf",width=8,height=8)





##' ========================================================================
##' Part III (Cont'd)
##' Case studies - Statistics
##' ========================================================================

##
cor(x$petal.length,x$petal.width)
cor(x$sepal.length,x$sepal.width)

cor.test(x$petal.length,x$petal.width)


##
t.test(petal.length~species, data=x)

##So, subset ...
t.test(petal.length~species, data=x[x$species%in%c("virginica","versicolor"),])

##
wilcox.test(petal.length~species,
            data=x[x$species%in%c("virginica","versicolor"),])

##
oneway.test(petal.length~species,data=x)
kruskal.test(petal.length~species,data=x)

##
mylm <- lm(petal.length~petal.width,data=x)
mylm
plot(x$petal.width,x$petal.length,main="Fisher's Iris Data")
abline(mylm,col="red")
dev.copy2pdf(file="iris-lm.pdf",width=6,height=6)


anova(lm(petal.length~species,data=x))


##
require(datasets)
HairEyeColor
##dimnames(HairEyeColor)

maleHairEyeColor <- HairEyeColor[,,"Male"]

maleHairEyeBrBlBu <- maleHairEyeColor[c("Brown","Blond"),c("Brown","Blue")]
maleHairEyeBrBlBu

chisq.test(maleHairEyeBrBlBu)

fisher.test(maleHairEyeBrBlBu)



##' ========================================================================
##' Appendix A. Exercise on R basics and solutions
##' ========================================================================


v <- c(7, 13, 19, 11, 17, 5)
v

for (element in v){
  if (element > 10){
    print(element)
  }
}


myfun1 <- function(x){
  res <- c()
  for (i in 2:length(x)){
    if (x[i]>x[i-1]){
      res <- c(res,x[i])
    }
  }
  return(res)
}
myfun1(v)


myfun2 <- function(x){
  the.element <- x[2:length(x)]
  the.neighbour <- x[1:(length(x)-1)]  
  res <- the.element[the.element>the.neighbour]
  return(res)
}
myfun2(v)




