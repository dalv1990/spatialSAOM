
# EXAMPLE USE

setwd("G:/My Drive/ssm-diffusion/analysis")
longSSM <- import("../data/generated/countries_longitudinal.Rdata")
contiguity <- read.csv("../data/generated/matrices/contiguity-2000.csv", header = FALSE)
idsGeography <- unname(t(contiguity[1,]))
matrixGeography <- as.matrix(contiguity[-1,])
n<-matrixGeography[1:30,1:30]
rownames(matrixGeography) <- colnames(matrixGeography) <- idsGeography
matrixGeography <- matrixGeography[order(idsGeography), order(idsGeography)]
df<-as.data.frame(cbind(floor(runif(30, min=0, max=1)), floor(runif(30, min=0, max=4)),floor(runif(30, min=0, max=4))))
names(df)<-c("y","x1","x2")
ss<-saom.static(y ~ x1, data=df, net=n, maxRound=2, method="avSim")

