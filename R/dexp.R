# numerically stable evaluation of
# 1 - exp(-x)^2
dexp2 <- function(x,Exp=exp(-x)) {if(Exp<0.7071068) return(1-Exp^2) else return(2*Exp*sinh(x))}

# 1 - exp(-x)^1
dexp1 <- function(x,Exp=exp(-x)) {if(Exp<0.5) return(1-Exp) else return(2*sqrt(Exp)*sinh(x/2))}