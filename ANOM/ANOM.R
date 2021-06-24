mls_data <- read.csv("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/all_mls_data.csv")[ ,-1]
mls_temp = mls_data[ ,c('X15.C','X20.C','X25.C')]
mls_fudr = mls_data[ ,c("no.FUdR.at.20.C","with.FUdR.at.20.C")]
mls_cntry = mls_data[ ,c("USA","China","Germany","Rep..of.Korea","Japan")]
mls_state = mls_data[ ,c("CA","MA","NY","WA","TX","MI")]

##################### Create ANOM-type data.frames ########################

###########################################################################
############################# Temperature #################################
###########################################################################
vec15 <- na.omit(mls_temp[,1])
vec20 <- na.omit(mls_temp[,2])
vec25 <- na.omit(mls_temp[,3])

## Create histogram-style column for temperature
temps <- rep('15C',length(vec15))
temps <- append(temps,rep('20C',length(vec20)), after=length(temps))
temps <- append(temps,rep('25C',length(vec25)), after=length(temps))
temps <- factor(temps)

## Create data column for temperature
mlst <- vec15
mlst <- append(mlst, vec20, after = length(mlst))
mlst <- append(mlst, vec25, after = length(mlst))

## Create dataframe
dftemp <- data.frame(temps,mlst)

###########################################################################
################################ FUDR #####################################
###########################################################################
vecnoF <- na.omit(mls_fudr[,1])
vecwF <- na.omit(mls_fudr[,2])

## Create histogram-style column for FUDR
fudrs <- rep('(-)FUDR',length(vecnoF))
fudrs <- append(fudrs, rep('FUDR',length(vecwF)), after=length(fudrs))
fudrs <- factor(fudrs)

## Create data column for FUDR
mlsf <- vecnoF
mlsf <- append(mlsf, vecwF, after = length(mlsf))

## Create dataframe
dffudr <- data.frame(fudrs,mlsf)

###########################################################################
############################# Countries ###################################
###########################################################################
vecUS <- na.omit(mls_cntry[,1])
vecCh <- na.omit(mls_cntry[,2])
vecGe <- na.omit(mls_cntry[,3])
vecRK <- na.omit(mls_cntry[,4])
vecJa <- na.omit(mls_cntry[,5])

## Create histogram-style column for country data
cntrys <- rep("USA",length(vecUS))
cntrys <- append(cntrys, rep("China",length(vecCh)), after=length(cntrys))
cntrys <- append(cntrys, rep("Germany",length(vecGe)), after=length(cntrys))
cntrys <- append(cntrys, rep("Rep..of.Korea",length(vecRK)), after=length(cntrys))
cntrys <- append(cntrys, rep("Japan",length(vecJa)), after=length(cntrys))
cntrys <- factor(cntrys)

## Create data column for countries
mlsc <- vecUS
mlsc <- append(mlsc, vecCh, after = length(mlsc))
mlsc <- append(mlsc, vecGe, after = length(mlsc))
mlsc <- append(mlsc, vecRK, after = length(mlsc))
mlsc <- append(mlsc, vecJa, after = length(mlsc))

## Create dataframe
dfcntry <- data.frame(cntrys,mlsc)

###########################################################################
############################## States #####################################
###########################################################################
vecCA <- na.omit(mls_state[,1])
vecMA <- na.omit(mls_state[,2])
vecNY <- na.omit(mls_state[,3])
vecWA <- na.omit(mls_state[,4])
vecTX <- na.omit(mls_state[,5])
vecMI <- na.omit(mls_state[,6])

## Create histogram-style column for country data
states <- rep("CA",length(vecCA))
states <- append(states, rep("MA",length(vecMA)), after=length(states))
states <- append(states, rep("NY",length(vecNY)), after=length(states))
states <- append(states, rep("WA",length(vecWA)), after=length(states))
states <- append(states, rep("TX",length(vecTX)), after=length(states))
states <- append(states, rep("MI",length(vecMI)), after=length(states))
states <- factor(states)

## Create data column for states
mlss <- vecCA
mlss <- append(mlss, vecMA, after=length(mlss))
mlss <- append(mlss, vecNY, after=length(mlss))
mlss <- append(mlss, vecWA, after=length(mlss))
mlss <- append(mlss, vecTX, after=length(mlss))
mlss <- append(mlss, vecMI, after=length(mlss))

## Create dataframe
dfstate <- data.frame(states,mlss)

###########################################################################
########################## PERFORM ANOM ###################################
###########################################################################
library(ANOM)
library(multcomp)
library(sandwich)
library(SimComp)

tempmodel <- lm(mlst ~ temps,data=dftemp)
fudrmodel <- lm(mlsf ~ fudrs,data=dffudr)
cntrymodel <- lm(mlsc ~ cntrys,data=dfcntry)
statemodel <- lm(mlss ~ states,data=dfstate)

## Perform standard ANOM 
t <- glht(tempmodel,mcp(temps="GrandMean"),alternative='two.sided')
f <- glht(fudrmodel,mcp(fudrs="GrandMean"),alternative='two.sided')
c <- glht(cntrymodel,mcp(cntrys="GrandMean"),alternative='two.sided')
s <- glht(statemodel,mcp(states="GrandMean"),alternative='two.sided')

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_temps_standard.jpg")
ANOM(t)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_fudrs_standard.jpg")
ANOM(f)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_cntrys_standard.jpg")
ANOM(c)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_states_standard.jpg")
ANOM(s)
dev.off()

## Perform ANOM with a heteroscedasticity-consistent sandwich estimate
t1 <- glht(tempmodel,mcp(temps="GrandMean"),alternative='two.sided',vcov=vcovHC)
f1 <- glht(fudrmodel,mcp(fudrs="GrandMean"),alternative='two.sided',vcov=vcovHC)
c1 <- glht(cntrymodel,mcp(cntrys="GrandMean"),alternative='two.sided',vcov=vcovHC)
s1 <- glht(statemodel,mcp(states="GrandMean"),alternative='two.sided',vcov=vcovHC)

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_temps_sandwich.jpg")
ANOM(t1)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_fudrs_sandwich.jpg")
ANOM(f1)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_cntrys_sandwich.jpg")
ANOM(c1)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_states_sandwich.jpg")
ANOM(s1)
dev.off()

## Perform ANOM considering heteroscedasticity by computing separate degrees of freedom and critical values for each of the contrasts
t2 <- SimCiDiff(data=dftemp,grp="temps",resp="mlst", type="GrandMean",covar.equal=F)
t2p <- SimTestDiff(data=dftemp,grp="temps",resp="mlst", type="GrandMean",covar.equal=F)
f2 <- SimCiDiff(data=dffudr,grp="fudrs",resp="mlsf", type="GrandMean",covar.equal=F)
f2p <- SimTestDiff(data=dffudr,grp="fudrs",resp="mlsf", type="GrandMean",covar.equal=F)
c2 <- SimCiDiff(data=dfcntry,grp="cntrys",resp="mlsc", type="GrandMean",covar.equal=F)
c2p <- SimTestDiff(data=dfcntry,grp="cntrys",resp="mlsc", type="GrandMean",covar.equal=F)
s2 <- SimCiDiff(data=dfstate,grp="states",resp="mlss", type="GrandMean",covar.equal=F)
s2p <- SimTestDiff(data=dfstate,grp="states",resp="mlss", type="GrandMean",covar.equal=F)

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_temps_dof.jpg")
ANOM(t2,stdep=dftemp$mlst,stind=dftemp$temps,pst=t2p)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_fudrs_dof.jpg")
ANOM(f2,stdep=dffudr$mlsf,stind=dffudr$fudrs,pst=f2p)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_cntrys_dof.jpg")
ANOM(c2,stdep=dfcntry$mlsc,stind= dfcntry$cntrys,pst=c2p)
dev.off()

jpeg("/Users/josephcavataio/Documents/Schnell_Lab/N2-Lifespans/ANOM_states_dof.jpg")
ANOM(s2,stdep=dfstate$mlss,stind=dfstate$states,pst=s2p)
dev.off()