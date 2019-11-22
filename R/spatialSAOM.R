library(statnet)
library(RSiena)
library(rio)


spatialSAOM <- function(formula, diffusion=list(), data, net, rateFix=20, maxRounds=10, method="avAlt", projname="saom",...){
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "net"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  net<-as.matrix(net)
  diag(net) <- 0
  net.vals<-attr(table(net),"dimnames")$net
  df=data
  if (!identical(net.vals, c("0","1"))){
    stop("net must have values 0 and 1 only")
  }
  if (dim(net)[1]!=dim(net)[2]){
    stop("net must be N x N matrix")
  }
  if (dim(net)[1]!= dim(df)[1]){
    stop("data must have as many rows as net")
  }
  
  
  # Make all ties structureal zeros and structural ones
  net<-net + 10
  
  # Generate Siena network object
  netData<- sienaDependent(array(c( net, net),dim = c( dim(net)[1], dim(net)[1],2)))

  # Add check for existence of variables in formula in data.frame
  # ...TODO
  
  # Generate Y variable
  dv<-cbind(data[,toString(formula[[2]])], data[,toString(formula[[2]])])
  dv[1,]<-c(0,1)
  dv[2,]<-c(1,0)
  dv.name<-toString(formula[[2]])
  assign(dv.name,sienaDependent(dv, type = "behavior",allowOnly = FALSE))
  
  # Generate X variables
  iv.l<-length(attr(terms(formula),"variables"))
  iv.name<-list()
  for(i in 3:iv.l){
    iv.name[i-2]<-toString(attr(terms(formula),"variables")[[i]])
    assign(iv.name[[i-2]],eval(parse(text = "coCovar(df[,toString(iv.name[[i-2]])])")))
  }
  print(paste0("sienaDataCreate(netData,", toString(dv.name), ",", toString(iv.name), ")"))
  
  
  saom.data = eval(parse(text = paste0("sienaDataCreate(netData,", toString(dv.name), ",", toString(iv.name), ")")))
  saom.data = sienaDataCreate(netData, y)
  # Assign diffusion effect for Y
  saom.eff <- getEffects(saom.data)
  if (method == "avSim"){
    saom.eff <- includeEffects(saom.eff, name = dv.name, avSim, interaction1 = "netData")
  } else {
    saom.eff <- includeEffects(saom.eff, name = dv.name, avAlt, interaction1 = "netData" )
  }

  for(i in 3:iv.l){
     print(i)
     print(toString(attr(terms(formula),"variables")[[i]]))
     saom.eff <- includeEffects(saom.eff, name = dv.name, effFrom, interaction1 = toString(attr(terms(formula),"variables")[[i]]))
  }
  
  # Assign diffusion effect for other X on Y
  if (exists("diffusion")){
     diffusion<-diffusion[diffusion != dv.name]
     diff.l<-length(diffusion)
     if (diff.l != 0) {
       for(i in 1:diff.l){
         diff.name<-diffusion[[i]]
         saom.eff <- setEffect(saom.eff, name = toString(dv.name), shortName=avXAlt, interaction1 = diff.name, interaction2 = "netData")
       }
     }
  }
  saom.eff
  
  # Fix network rate parameter to zero
  # saom.eff[saom.eff$effectName == "basic rate parameter netData", "fix"] <- TRUE
  # saom.eff[saom.eff$effectName == "basic rate parameter netData", "initialValue"] <- 0

  # Exclude linear shape because there is no t1 and t2.
  saom.eff$include[saom.eff$effectName == paste(dv.name, "linear shape")] <- FALSE
  
  # Fix rate parameter for behavior at large enough value.
  # Otherwise the algorithm immediately stops and thinks that
  # it found the perfect solution already to get from fake-t1 to fake-t2.
  saom.eff[saom.eff$effectName == paste("rate",dv.name,"period 1"), "fix"] <- TRUE
  saom.eff[saom.eff$effectName == paste("rate",dv.name,"period 1"), "initialValue"] <- rateFix
  
  # Create algorithm end estimate SAOM
  saom.algorithm <- sienaAlgorithmCreate(useStdInits = FALSE, projname = projname)
  saom.ans <- siena07(saom.algorithm, data = saom.data, eff = saom.eff)
  saom.tm <- saom.ans$tconv.max
  
  # Repeat until convergence is reached
  estimationRound <- 1
  # this should not happen, but it does ...
  if (!is.na(saom.tm)){ 
    while (estimationRound < maxRounds & saom.tm > 0.25 & saom.tm < 10) {
         saom.ans <- siena07(saom.algorithm, data = saom.data, eff = saom.eff, prevAns = saom.ans)
         saom.tm <- saom.ans$tconv.max     
         estimationRound <- estimationRound + 1
         if (is.na(saom.tm)) break
    }
  }
}






# 
# 
# 
# 
# 
# n <- matrix(c(0,1,1,1,0,0,1,0,0), nrow=3)
# df<-data.frame(cbind(c(0,0,1),c(50,10,20), c(0,1,1)))
# names(df)<-c("ssm", "gdp", "religion")
# s<-saom.static(ssm ~ gdp + religion, data=df, net=n)
# 
# 
# # generate Y variable
# dv<-cbind(data[,toString(formula[[2]])], data[,toString(formula[[2]])])
# dv[1,]<-c(0,1)
# dv.name<-toString(formula[[2]])
# assign(dv.name,sienaDependent(dv, type = "behavior",allowOnly = FALSE))
# 
# # generate X variables
# l<-length(attr(terms(formula),"variables"))
# iv.name<-list()
# for(i in 3:l){
#   iv.name[i-2]<-toString(attr(terms(formula),"variables")[[i]])
#   assign(iv.name[[i-2]],eval(parse(text = "coCovar(df[,toString(iv.name[[i-2]])])")))
# }
# 
# 
# 
# 
#   