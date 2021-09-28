# get the weights similar as the bernstein_dp
# get the knots
# add the weights and knots propority to the prior curve
set.seed(42)
x=seq(0,1,length=1000)
K=10
r=3 #the degree
L=10000
M=1
G0=runif

# the Stick-breaking process to have Zl
# Sample V
V1 = rbeta(L, 1, M)

# Construct stick-breaking probabilities (p from V)
p = rep(NA, L)
p[1] = V1[1]
for (i in 2:L) {
  p[i] = prod(1 - V1[1:(i - 1)]) * V1[i]
}
p = c(1 - sum(p), p)  # Attach p0 as p[1]

# Sample U from base measure G0
Z = G0(L + 1)

# the Stick-breaking process to have Zl
# Sample V
V2 = rbeta(L, 1, M)

# Construct stick-breaking probabilities (p from V)
q = rep(NA, L)
q[1] = V2[1]
for (i in 2:L) {
  q[i] = prod(1 - V2[1:(i - 1)]) * V2[i]
}
q = c(1 - sum(q), q)  # Attach p0 as p[1]

# Sample U from base measure G0
X = G0(L + 1)

# comeout the knot_diff
K=10
r=3

knot_diff=function(j,K,r,q,X){
  # get Xbound
  boundd=(j-1)/(K-r)
  boundu=j/(K-r)
  # get according q_sum
  out=0
  for (i in 1:length(q)){
    if (X[i]<=boundu&X[i]>boundd){out=out+q[i]}
  }
  return(out)
}

#sequence of knots
interal_knots=rep(0,(K-r+1))
for (i in 2:length(interal_knots)){
  interal_knots[i]=knot_diff(i,K,r,q,X)+interal_knots[i-1]
}
knots=c(rep(0,r),interal_knots,rep(1,r))

# come out the weights
weights=function(j,K,p,Z){
  # get Xbound
  boundd=(j-1)/K
  boundu=j/K
  # get according q_sum
  out=0
    for (i in 1:length(p)){
  if (Z[i]<=boundu&Z[i]>boundd){out=out+p[i]}
  }
  return(out)
}

# get the mixture density
library(splines)
bb=splineDesign(knots, x = x, outer.ok = TRUE)
splinem=function(x,knots){
  out=matrix(0,length(x),K)
  bb=splineDesign(knots, x = x, outer.ok = TRUE)
  for(i in 1:K){
    out[,i]=bb[,i]*weights(i,K,p,Z)
  }
  return(out)
}
splinesmixture=splinem(x,knots)

#plot
par(mfrow = c(1, 2))
plot(range(x), c(0,1), type = "n", xlab = "x", ylab = "",
     main =  "B-splines prior")
abline(v = knots, lty = 3, col = "light gray")
abline(v = knots[c(4,length(knots)-3)], lty = 3, col = "gray10")
lines(x, rowSums(splinesmixture), col = "gray", lwd = 2)
matlines(x, splinesmixture, ylim = c(0,1), lty = 1)
matplot(x, bb, ylim = c(0,1),lty=1,type = 'l', main =  "B-splines mixture")