##' Feature mapping function to polynomial feautures
##'
##' mapFeat maps the input variable X into quadratic feautures.
##' for instance, if the input is a two dimensional array the function
##' returns a new array with more features, comprising of 
##'   X1, X2, X1.^2, X2.^2, X1*X2, X1*X2.^2
##' 
##' @param X an array to be mapped. Must be numeric!
##' @param order the final polynomial order to map
##' (defaults to 2), may take a long time if order is large.
##' at this time also only accepts a value of two. 
##' @author Marco D. Visser
##' @examples
##' mapFeat(array(rnorm(20),dim=c(10,2)))
##' @export


mapFeat <- function(X,order=2){

  
  obsize <- dim(X)
  out <- matrix(0,nrow=obsize[1],ncol=obsize[2]+(2*2))
  
  out[,1:2] <- X
  out[,3:4] <- X^2
  out[,5:6] <- c(apply(X,1,prod),apply(X^2,1,prod))
 out 
}



