##' Cost function for a logistic classifier
##'
##' calculates the cost and derivatives for a simple
##' logistic classifier
##
##' @param theta numeric vector of parameters
##' @param X design matrix: an array of numeric variables or features
##' including intercepts
##' @param y binary classification vector
##' @param lambda regularization coefficient
##' @author Marco D. Visser
##' @rdname logisCost
##' @examples
##' X<-array(rnorm(200),dim=c(100,2))
##' y<-sqrt(X[,2]^2+X[,1]^2)<.5
##' designX<-mapFeat(X)
##' designX<-cbind(rep(1,100),designX)
##' theta<-rep(0,ncol(designX))
##' lambda<-1
##' logisCost(theta,designX,y,lambda)
##' logisGrad(theta,designX,y,lambda)
##' optim(theta,logisCost,logisGrad,X=designX,y=y,lambda=1,method="BFGS")
##' @export
logisCost <- function(theta,X,y,lambda){

  m <- length(y)
  K <- length(theta)
  grad <- numeric(K)

  pred <- plogis(X%*%theta)

  ## Cost function
  J <- (sum(-y*log(pred)-(1-y)*log(1-pred))*(m)^-1)+
  (lambda/(2*m)*sum(theta[2:K]^2))


  return(J)
}


##' Gradient function for a logistic classifier
##' @rdname logisCost
##' @export
logisGrad <- function(theta,X,y,lambda){

  m <- length(y)
  K <- length(theta)
  grad <- numeric(K)

  pred <- plogis(X%*%theta)

 grad <- (1/m)*t(pred-y)%*%X +(lambda/m)*theta
 grad[1] <- (1/m)*sum((pred-y)*X[,1])

  return(grad)
}

##' Prediction function for a fitted logistic classifier
##'
##' Classifies from variables X 
##
##' @param theta numeric vector of parameters
##' @param X design matrix: an array of numeric variables or features
##' including intercepts
##' @param thres decision boundry
##' @author Marco D. Visser
##' @param type type of predictions on logit scale ("logit") or
##' "reponse" or as classification ("class")? Defaults to response
##' @rdname logisPred
##' @examples
##' X<-array(runif(2000),dim=c(1000,2))
##' y<-sqrt(2*X[,2]^2+2.4*X[,1]^2)<.5
##' designX<-mapFeat(X)
##' designX<-cbind(rep(1,100),designX)
##' theta<-rep(0,ncol(designX))
##' lambda<-1
##' logisCost(theta,designX,y,lambda)
##' logisGrad(theta,designX,y,lambda)
##' par<-optim(theta,logisCost,logisGrad,X=designX,y=y,lambda=.01,method="BFGS")$par
##' decPlot(par,designX,X,y,thres=.5)
##'@export
logisPred <- function(theta,X,thres=0.5,type="response"){

  if(type=="response") {return(plogis(X%*%theta))}
  if(type=="class") {return(ifelse(plogis(X%*%theta)>=thres,1,0))} else {
    X%*%theta}
}

##' Decision boundry plot function for a fitted logistic classifier
##'
##' Classifies from variables X 
##
##' @param theta numeric vector of parameters
##' @param dX design matrix: an array of numeric variables or features
##' @param X original variables or features
##' @param y original classification
##' including intercepts
##' @param thres decision boundry
##' @author Marco D. Visser
##' @rdname logisPred
##' @export
decPlot <- function(theta,dX,X,y,thres=0.5){

  py <- logisPred(theta=theta,X=dX,thres=thres,type='class')
  sizes <- apply(X,2,range)
  xc <- seq(sizes[1,1],sizes[2,1],length.out=100)
  yc <- seq(sizes[1,2],sizes[2,2],length.out=100)
 
  grid <- outer(xc,yc,function(x,y)
                logisPred(theta=theta,X=cbind(rep(1,length(x)),
                                        mapFeat(cbind(x,y))),
                          thres=thres,type='class')    )
  plot(X,pch=c(19,3)[as.numeric(y)+1],
       col=c("red","blue")[as.numeric(y)+1])
  contour(xc,yc,grid,add=TRUE,col='green')
  
}

##' m- Class classifier
##'
##' Trains a m-class classifer from dataset trainingdata
##
##' @param trainingdata numeric matrix with training data
##' @param y classes corresponding to the rows of trainingdata
##' @param lambda regularization coefficient
##' @author Marco D. Visser
##' @rdname mClass
##' @examples
##' data(digits)
##' X<-digits[[1]]
##' y<-digits[[2]]
##' trained<-mClass(X,y)
##' 
##' 
##' 
##' 
##' 
##' 
##'@export
mClass <- function(trainingdata,y,lambda=1) {

  classes <- unique(y)
  m <- length(classes)
  n <- nrow(trainingdata)
  k <- ncol(trainingdata)
  allTheta <- array(dim=c(m,k+1))
  inits <- rep(0,k+1)
  dX <- cbind(rep(1,n),trainingdata)

  pb <- txtProgressBar(min=1,max=m)
  
  for(i in 1:m) {
    id <- y==classes[i]
    allTheta[i,]<-
      optim(inits,logisCost,logisGrad,X=dX,y=id,
            lambda=lambda,method="CG")$par
  setTxtProgressBar(pb, i)
  }

  close(pb)
return(list(allTheta,classes))
}


##' m-Class predictor
##'
##' uses a trained classifier to predict classes
##
##' @param Data numeric matrix with data to classify
##' @param Theta a previously trained classifier
##' @author Marco D. Visser
##' @rdname mClass
##' @examples
##' data(digits)
##' X<-digits[[1]]
##' y<-digits[[2]]
##' trained<-mClass(X,y)
##' pred<-mPred(X,trained)
##' image1200<-t(as.array(X[1200,]))
##' pred1200<-mPred(image1200,trained)
##' ## force correct rotation of matrix image
##' f <- function(m) t(m)[,nrow(m):1]
##' image(f(matrix(X[1200,]),ncol=20,nrow=20),main=pred1200)
##'@export
mPred <- function(Data,Theta) {

  n<-nrow(Data)
  classes <- Theta[[2]]
  dX <- cbind(rep(1,n),Data)

  p=apply(plogis(dX%*%t(Theta[[1]])),1,which.max)

return(classes[p])
}


##' Digit dataset animation 
##'
##' uses a trained classifier to predict classes
##' and displayes random images with predictions on
##' screeen with a delay
##
##' @param Data numeric matrix with data to classify
##' @param Theta a previously trained classifier
##' @param maxit maximum number of random displayes
##' @param y True classes used in training
##' @author Marco D. Visser
##' @rdname mClass
##' @examples
##' data(digits)
##' X<-digits[[1]]
##' y<-digits[[2]]
##' trained<-mClass(X,y)
##'@export
##'@export
 digitAnim <- function(Data,Theta,y,maxit=10) {

   oldpar <- par("mar","xaxt","yaxt")
   par(mar=c(0,0,2,0),xaxt="n",yaxt="n")
   
 preds <- mPred(Data,Theta)
 f <- function(m) t(m)[,nrow(m):1]


 for( i in 1:maxit){
 n <- sample(1:nrow(Data),1)

 image(f(matrix(X[n,],ncol=20,nrow=20)),main=preds[n],
       xlab='',ylab='',col=grey.colors(400))

 Sys.sleep(0.75)
 }
   print(paste("Accucary of classifier",mean(y==preds)*100,"%"))
   par(oldpar)
}


#' 5000 Handwritten digits 
#' 
#' A dataset containing 5000 20x20 pixel images of handwritten digits
#' The dataset comes from the Standford ML course 2014. It is a list of tow items:
#' 1) a 5000 x 400 matrix containing the grey scale dn values for each image
#' in each row and 2) the classification of each digit as 0 - 9 as c(10,1-9).
#' 
#' \itemize{
#'   \item listitem1. a 5000 x 400 matrix
#'   \item listitem1. a 5000 x 1 vector of classifications
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name digits
#' @usage data(digits)
#' @format A list with 5000 images and their classification
NULL
