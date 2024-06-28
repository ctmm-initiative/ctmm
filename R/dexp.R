# numerically stable evaluation of
# 1 - exp(-x)^2
dexp2 <- function(x,Exp=exp(-x)) { ifelse(Exp<0.7071068,1-Exp^2,2*Exp*sinh(x)) }
# 1 - exp(-x)^1
dexp1 <- function(x,Exp=exp(-x)) { ifelse(Exp<0.5,1-Exp,2*sqrt(Exp)*sinh(x/2)) }