#Dice's temperature
#install.packages("rootSolve")
library(rootSolve)

Beta = function(E){ #回傳Beta
  domain = c(1:13)
  heatfun = function(B){
    z = 0
    for(i in domain){
      z = z+exp(-B*i)
    }
    s = 0
    for(i in domain){
      s = s+i*exp(-B*i)
    }
    return(s/z-E)
  }
  root = uniroot(heatfun, c(-4, 4))$root #Beta的解
  print(root)
  return(root)
}
S = function(B){ #Entropy of Beta
  domain = c(1:13)
  p = rep(0,length(domain))
  z = 0
  S = 0
  for(i in domain){
    p[i] = exp(-B*i)
    z = z+exp(-B*i)
  }
  p = p/z
  for(i in c(1:length(p))){
    S = S-p[i]*log(p[i])
  }
  plot(p)
  return(S)
}
S(0.5)
S(Beta(8))

#Beta-E
x = vector();y = vector()
for(i in seq(1.1,5.9,by = 0.1)){
  x = append(x,i)
  y = append(y,Beta(i))
}
Tem = 1/y
plot(x,y,xlab = "fixed mean",ylab = "Beta(E)")
plot(x,Tem,xlab = "fixed mean",ylab = "1/Beta = Tem")

#S-Beta
x = vector();y = vector()
for(i in seq(-4,4,by = 0.01)){
  x = append(x,i)
  y = append(y,S(i))
}
plot(x,y,xlab = "Beta",ylab = "Entropy(Beta)")

#S-T
x = vector();y = vector()
for(i in seq(-20,20,by = 0.1)){
  x = append(x,i)
  y = append(y,S(1/i))
}
plot(x,y,xlab = "T",ylab = "Entropy(T)")


#S-E
x = vector();y = vector()
for(i in seq(1.1,12.9,by = 0.1)){
  x = append(x,i)
  y = append(y,S(Beta(i)))
}
plot(x,y,xlab = "fixed mean",ylab = "Entropy(mean)")




#curve(fun(x), 1.5, 5)
#uni <- uniroot(heatfun, c(-8, 8))$root
#points(All, y = rep(0, length(All)), pch = 16, cex = 2)





