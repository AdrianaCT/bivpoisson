
dbipois<-function(variables,params){
    if (all(variables == floor(variables))==FALSE) stop('Both variable values must be positive integers')
    if (all(variables>-1)==FALSE) stop('Both variable values must be positive integers')
    if (all(params>0)==FALSE) stop('All parameter values must be positive numbers')
    if (length(variables)>2) warning('More than two variable values were given. Evaluating first two values only',call.=FALSE)
    x=variables[1] ;     y=variables[2];
    lambda1=params[1];    lambda2=params[2];
    xi=params[3];
    density1=exp(-(lambda1+lambda2+xi))
    summation=0
    for (z in 0:min(x,y)){
    summation=summation+( (lambda1^x)/factorial(x-z) *
    (lambda2^y)/factorial(y-z) *(1/factorial(z))*
    (xi/(lambda1*lambda2))^z )
    }
    density=density1*summation
    return(density)
 }


rbipois<-function(n, target.params, envelope.params=list(mean=c(3,3),var=
    5*matrix(c(2,1,1,2),ncol=2)),k=1){
    if (require(mvtnorm, quietly=TRUE)==FALSE) stop('This function needs the mvtnorm package to run. Please install it and try again')
    #SETUP
    mean=envelope.params[[1]]; vcov=envelope.params[[2]]
    #STORAGE
    the.sample=matrix(NA,ncol=2,nrow=n)
    for (j in 1:n){
        calc=1000
        U=calc*k
        while(calc <= k*U){
            X=c(-1,-1)
            while( X[1] < 0 | X[2]< 0 ){
                X=rmvnorm(1,mean=mean,sigma=vcov)
            }
            X=floor(X)
            U=runif(1)
            fofx=dbipois(X,target.params)
            gofx=dmvnorm(X,mean=mean,sigma=vcov)        
        calc=fofx/gofx
        }
        the.sample[j,]=X
    }
    return(the.sample)
} 

