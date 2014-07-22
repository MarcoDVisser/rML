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
##' theta<-rep(1,ncol(designX))
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


 grad[1] = (1/m)*sum((pred-y)*X[,1])


 grad[2:K] = sapply(2:K,function(i)
        (1/m)*sum((pred-y)*X[,i])+(lambda/m)*theta[i])

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
##' X<-array(rnorm(200),dim=c(100,2))
##' y<-sqrt(X[,2]^2+X[,1]^2)<.5
##' designX<-mapFeat(X)
##' designX<-cbind(rep(1,100),designX)
##' theta<-rep(1,ncol(designX))
##' lambda<-1
##' logisCost(theta,designX,y,lambda)
##' logisGrad(theta,designX,y,lambda)
##' par<-optim(theta,logisCost,logisGrad,X=designX,y=y,lambda=1,method="BFGS")$par
##' decPlot(par,designX,X,y,thres=.5)
##' plot(logisPred(par,designX))
##' @export
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


