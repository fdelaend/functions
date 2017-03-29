
#Weighted Euclidian distance matrix
#positions = matrix with:
#nr of rows= nr of locations; 
#nr of columns = nr of dimensions

#weights= weights attributed to the two dimensions
WEuclid <- function(positions, weights)
{
  #First empty distance matrix
  distances     <- diag(nrow(positions))*0
  #Calculate weighted Euclidian distance between all points
  for (p in 1:nrow(positions))
  {
    position <- as.numeric(positions[p,])
    position <- t(matrix(position,ncol=n, nrow=traits))
    distance <- sqrt(rowSums(((weights)*(positions - position))^2))
    distances[p,]    <- distance
  }
  return(distances)
}




#Partition according to Fox 2005
#M = vector with monocultures
#Yo = vector with observed yields
#Initial = vector with initial densities

Fox2005 <- function(M, Yo, Initial, Check=T)
{
  Richness <- length(M)
  Ye <- Initial/sum(Initial) * M
  RYe <- Initial/sum(Initial)
  RYo <- Yo/M
  DeltaRY <- RYo - RYe
  RYo_RYot__Rye <- RYo/sum(RYo)-RYe
  RYo_RYo__RYot <- RYo-RYo/sum(RYo)
  TIC <- Richness * mean(M) * mean(DeltaRY)
  DOM <- Richness * sum((RYo_RYot__Rye-mean(RYo_RYot__Rye))*(M-mean(M)))/Richness
  TDC <- Richness * sum((RYo_RYo__RYot-mean(RYo_RYo__RYot))*(M-mean(M)))/Richness
  #you can verify that TIC+DOM+TDC = sum(Yo-Ye)
  if (Check==T)
  {
    if (abs((TIC+DOM+TDC) - sum(Yo-Ye)) > abs(sum(Yo-Ye))*1e-5) stop("!")
  }
  return(list(TIC=TIC, DOM=DOM, TDC=TDC))
}

#Calculates 
#...connectances based on an abundance file and a food-web matrix
#...and also nr of links for all pred and all prey
#FwMatrix = food web matrix filename
#AbName = abundance filename
#pred = vector of integers, tells us the places where predator abundances can be found in the abundance files
#prey = same, but for pred

ConnecAbFoodweb <- function(FoodwebName, AbName, Ab=NULL,
                            pred, prey)
{
  FwMatrix <- read.delim(FoodwebName, header=FALSE, sep=" ") #gets the FW matrix with the specified name
  FwMatrix <- data.matrix(FwMatrix)
  if (is.null(Ab))
    {
    Abs <- read.delim(AbName,header=FALSE, sep=" ") #gets the abundances with the specified name
    } else
    {Abs <- Ab}
  PreyAbs <- Abs#creates a copy of the abudnaces with a diff name
  PreyAbs[pred,] <- 0 #sets pred abundances in copy to zero (cause we only want prey in this one)
  PredAbs <- Abs #creates a copy of the abudnaces with diff name
  PredAbs[prey,] <- 0#sets prey abundances in copy to zero (cause we only want pred in this one)
  
  SpPresent <- as.array((Abs>0)*1) #all species present get a value of 1
  PreyPresent <- (PreyAbs>0)*1 #all prey sp present get a value of 1
  PredPresent <- (PredAbs>0)*1 #all pred sp present get a value of 1
  NrPreyPresent <- colSums(PreyPresent) #nr of prey present per iteration = nr of '1's per column
  NrPredPresent <- colSums(PredPresent)#nr of pred present per iteration = nr of '1's per column
  
  NrOfLinksPrey <- SpPresent*(FwMatrix%*%SpPresent) 
    #nr of links for every species and for every iteration, for predator species this will be zero 
    #why? --> tricky to explain, matrix calculus
  NrOfLinksPrey <- NrOfLinksPrey[prey,] #throws away nr of links for the predators (will be zero) 

  NrOfLinksPred <- SpPresent*(t(FwMatrix)%*%SpPresent) 
    #nr of links for every species and for every iteration, for prey species this will be zero 
    #why? --> tricky to explain, matrix calculus
  NrOfLinksPred <- NrOfLinksPred[pred,] #throws away nr of links for the prey (will be zero) 

  TotalNrOfLinks <- colSums(NrOfLinksPrey) #nr of links per iteration; equivalent to TotalNrOfLinks <- colSums(NrOfLinksPred)
  Connectance <- TotalNrOfLinks/(NrPreyPresent*NrPredPresent) #the 'nr of iterations' connectance values  

  return(list(Connectance=Connectance,Prey=NrOfLinksPrey,Pred=NrOfLinksPred))

}


#Plots FW diagram for MI6 2 trophic level food web:
#only 2 TLs: producers (top) and consumers (bottom)
#Food web matrices need to have rows (and cols) 1:x as the x producers
# and rows (and cols) x+1:y as the y-x consumers
#Food web matrices need to have the name FW_[connectance value]
#Will give pdf file with web as output

PlotMI6Web <- function(connectances=c("05"),
                       pathwebs="dynamic/parameters/VarParts/Output/FW_",
                       nrProd=50, nrCons=50, 
                       colProd="green", colCons="black")
{
  require('diagram')
  for (connectance in connectances)
  {
    FW <- read.delim(paste(pathwebs,
                           connectance, sep=""),
                     sep=" ", header=FALSE)
    colnames(FW) <- rownames(FW) <- c(paste("p",prey), paste("z",pred))
    quartz("",4,4,type="pdf",file=paste(connectance,"Fw.pdf"))
    par(mar=c(1.5,1.5,1,1), mgp=c(3,0.5,0))            #2 3
    
    plotmat(t(FW), pos=c(nrProd,nrCons), box.type="round", box.size=0.005,
            curve=0, lwd=0.1,arr.pos=1000,
            box.col=c(rep(colProd, 50), rep(colCons, 50)),
            txt.col=rep("transparent", 100),shadow.col="transparent",
            cex.txt=0, main="")
    dev.off()    
  }
}


#Converts vector into monotonically increasing or decreasing vector.
#Uses weighted averaging (based on sample size of elements of x) for this conversion
#Used in Williams.t.test.WUR function
#x=vector, often of means
#n=sample size, i.e. the number of data points on which the elements of x are based
#order=expected order (decr or incr)
#nr of decimals of the increments (decreases) along x... 
#...in case order=incr(decr)
#As a consequence, 'x' will be rounded to that nr-1
Monoton <- function(x,n,
                    order="decr",
                    decims=4)
{
  delta <- (-x[-1]+x[-length(x)])*(order=="decr") + (x[-1]-x[-length(x)])*(order=="incr")
  
  while(sum(delta<0))
  {
    for (i in c(1:(length(x)-1)))
    {
      delta_i <- (x[i]-x[i+1])*(order=="decr") + (-x[i]+x[i+1])*(order=="incr")
      if (delta_i<0)
      {
        x[i] <- x[i+1] <- (x[i]*n[i] + x[i+1]*n[i+1])/sum(n[i],n[i+1])
      }
    }  
    delta <- (-x[-1]+x[-length(x)])*(order=="decr") + (x[-1]-x[-length(x)])*(order=="incr")
    delta <- round(delta,decims)
  }
  
  return(round(x,decims-1))
}

#Williams t test: to be used when the response variable is assumed to monotonically decrease with increasing dose
#data = data frame in WUR format --> Note that dose levels need to be coded as integers, with '1' as the control
#order= is Index expected to increase (incr) or decrease (decr) with treatment?
#BetasFile (and TsFile) = text file holding the beta coefficients (t's) for different d.f (rows)
#             and treatments (columns) (Williams 1972, table 1)
#alpha = the alpha level: Note that this only changes the critical t at the lowest trt. So best not to change ever.
williams.t.test.WUR <- function(data,
                                TimeName, 
                                TreatmentName,
                                Index="E",
                                order="decr",
                                BetasFile="WilliamsBeta.txt",
                                TsFile="WilliamsT.txt",
                                alpha=0.05)
{
  times <- as.numeric(as.character(unique(data[,TimeName])))
  trts <- as.numeric(as.character(unique(data[,TreatmentName])))
  Sigs <- NULL
  BetasTable <- read.delim(BetasFile, header=FALSE)
  ColNamesBetas <- BetasTable[1,-1]
  RowNamesBetas <- BetasTable[-1,1]
  BetasTable <- BetasTable[-1,-1]
  colnames(BetasTable) <- ColNamesBetas
  rownames(BetasTable) <- RowNamesBetas
  
  TsTable <- read.delim(TsFile, header=FALSE)
  ColNamesTs <- TsTable[1,-1]
  RowNamesTs <- TsTable[-1,1]
  TsTable <- TsTable[-1,-1]
  colnames(TsTable) <- ColNamesTs
  rownames(TsTable) <- RowNamesTs

  for (time in times[which(times>0)])#only include post-exposure!!
  {
    Inds <- which(data[,TimeName]==time)
    Means <- aggregate(data[Inds,"E"],
                       by=list(data[Inds,TreatmentName]),
                       FUN=mean)$x
    Ns <- aggregate(data[Inds,"E"],
                    by=list(data[Inds,TreatmentName]),
                    FUN=length)$x
    Sums <- aggregate(data[Inds,"E"],
                      by=list(data[Inds,TreatmentName]),
                      FUN=sum)$x
    SumSqs <- aggregate(data[Inds,"E"],
                      by=list(data[Inds,TreatmentName]),
                      FUN=function(x){return(sum(x^2))})$x
    SumSumNs <- Sums*Sums/Ns
    MLEs <- Monoton(x=Means,n=Ns,order=order)
    Dfs <- sum(Ns) - length(Ns)
    TotalVar <- (sum(SumSqs) - sum(SumSumNs))/Dfs
    ts <- abs((MLEs[-1]-Means[1])/sqrt(TotalVar/Ns[-1]+TotalVar/Ns[1]))
    ws <- Ns[1]/Ns[-c(1:2)]
    ts_w <- as.numeric(TsTable[as.character(Dfs),as.character(trts[-c(1,2)])])
    Betas <- as.numeric(BetasTable[as.character(Dfs),as.character(trts[-c(1,2)])])
    ts_c <- qt(1-alpha, Dfs)#for treatment 2 (so lowest treatment), critical t=student's t
    ts_c <- c(ts_c, ts_w-Betas*(1-1/ws)/100)
    Sigs <- rbind(Sigs,cbind(time,trts[-1],(ts>ts_c)))
  }
  
  colnames(Sigs) <- c(TimeName,TreatmentName,"Significance")
  Sigs <- rbind(c(1,1,0),Sigs)
  Sigs <- as.data.frame(Sigs)
  Sigs[,TreatmentName] <- as.factor(Sigs[,TreatmentName])
  return(Sigs) 
}
#####################################################################
#####################################################################
#####################################################################

t.test.WUR <- function(data,TimeName, TreatmentName,Index="E")
{
  times <- as.numeric(as.character(unique(data[,TimeName])))
  trts <- as.numeric(as.character(unique(data[,TreatmentName])))
  Ps <- NULL

  for (trt in trts[-1])
  {
    for (time in times[which(times>0)])#only include post-exposure!!
    {
      Inds <- which((data[,TreatmentName]==trt)&(data[,TimeName]==time))
      Inds0 <- which((data[,TreatmentName]==1)&(data[,TimeName]==time))
      if ((length(unique(c(data[Inds,Index], data[Inds0,Index])))>2)&(length(Inds)>1)&(length(Inds0)>1))
      {
        Test <- t.test(data[Inds,Index], 
                       data[Inds0,Index],
                       alternative="l")$p.value
      } else 
        {
        Test <- 1  
        }
      Ps <- rbind(Ps,c(trt, time, Test))
    }
  }
  colnames(Ps) <- c(TreatmentName,TimeName,"p")
  Ps <- rbind(c(1,1,1e10),Ps)
  Ps <- as.data.frame(Ps)
  Ps[,TreatmentName] <- as.factor(Ps[,TreatmentName])
  return(Ps)    
}

#tests if time trends are different between polluted and clean biodiversity (BD)
#It works by substracting the control time series of BD from the treated BD time series and testing if the resulting trend differs from zero
#data = data frame in WUR format
#TimeName and TreatmentName: cf below
#the procedure assumes treatment only occurs when time>0 (as is the case in the WUR format)
#Index= BD index to be calculated
TrendTest <- function(data, TimeName, TreatmentName, Index="E")
{
  require(mgcv)
  #substraction and putting results in new data frame
  dataDelta <- NULL
  times <- as.numeric(as.character(unique(data[,TimeName])))
  trts <- as.numeric(as.character(unique(data[,TreatmentName])))
  
  for (trt in trts[-1])
  {
    for (time in times[which(times>0)])#only include post-exposure!!
    {
      Inds <- which((data[,TreatmentName]==trt)&(data[,TimeName]==time))
      Inds0 <- which((data[,TreatmentName]==1)&(data[,TimeName]==time))
      Delta <- mean(data[Inds,Index]) - mean(data[Inds0,Index])
      dataDelta <- rbind(dataDelta,c(trt, time, Delta))
    }
  }
  colnames(dataDelta) <- c("Treatment", "Time", "BD")
  dataDelta <- as.data.frame(dataDelta)
  dataDelta$Treatment <- as.factor(dataDelta$Treatment)
  
  #the gam fitting
  SmoothersBasis <- "s(Time, by=as.numeric(Treatment=="
  SmootherTail <- "), k=3)"
  Form <- "BD~"
  
  for (trt in unique(dataDelta$Treatment))
  {
    Form <- paste(Form, SmoothersBasis, trt, SmootherTail)
    if (as.numeric(trt)<max(as.numeric(as.character(unique(dataDelta$Treatment))),na.rm=TRUE))
    {
      Form <- paste(Form, "+")
    }
  }
  Model <- gam(as.formula(Form), data=dataDelta)
  PvaluesSmoothers <- c(1,summary(Model)$s.table[,"p-value"])
  names(PvaluesSmoothers) <- c(1:length(PvaluesSmoothers))
  return(p=list(PvaluesSmoothers))
}

#BD calculation based on WUR format
#data = WUR format data
#index = E (Simpson) or R (richness)
#CountCols = columns with count data
BDCalc <- function(data,CountCols,Index)

{
  if (Index=="E") { 
                  BD <- 1/rowSums((data[,CountCols] / rowSums(data[,CountCols]))^2)
                  BD[which(is.na(BD))] <- 0
                  data$E <- BD}
  if (Index=="R") { 
                  BD <- rowSums(data[,CountCols]>0)
                  data$R <- BD} 
  
  return(data)
}

#plots evenness and richness over time using cosm data in the standard WUR file format
#Filename = name of txt file with count data
#TreatmentName = name in the count data file used to indicate treatment level
#TimeName = name in the count data file used to indicate time
#CountCols = columns where species counts are present: a vector with length = 2. If length = 1, it is assumed species counts start at CountCols and continue until the last column 
#FigName = name of pdf file with plots that will be created
#If gamfit=TRUE, the TrendTest function is used to tell if the trend is different from the control trend
#ttest= the type of test that will be 
#       performed to see if BD is different from control at various points in time
#       Defaults to Williams (71 and 72 papers). 
#       Alternatives are t (t test) and None (no test).
PlotEvennessRichness <- function(Filename, CountCols, TimeName, 
                                 TreatmentName, FigName="Plot", 
                                 gamfit=FALSE, ttest="Williams",
                                 ExpectedTrend="decr")
{
  quartz("",6,3,type="pdf",file=paste(FigName, ".pdf", sep=""))
  data <- read.delim(Filename)
  data[,TreatmentName] <- as.factor(data[,TreatmentName])
  #replace NAs by zeros
  nas <- which(is.na(data))#,arr.ind=TRUE
  if (length(nas)>0) 
  {
    data[is.na(data)] <- 0
    print("NAs in data; replaced by zeros") 
  }
  
  if (length(CountCols)==1) {CountCols <- c(CountCols: ncol(data))}
  
  for (Index in c("E", "R"))
  {
    data <- BDCalc(data,CountCols,Index)
    MinY <- floor(min(data[,Index], na.rm=TRUE))
    MaxY <- ceiling(max(data[,Index], na.rm=TRUE))
    dataOriginalTreatmentColumn <- data[,TreatmentName]
    
    if ((length(unique(data[,TimeName]))>4)&(gamfit==TRUE))
    {
      data$Star <- ""
      p <- TrendTest(data=data, TimeName, TreatmentName, Index=Index)[[1]]
      data$Star[which(p[as.numeric(as.character(data$Treatment))]<0.05)]<- "*"
      data[,TreatmentName] <- paste(data[,TreatmentName],data$Star)    
    }  
    form <- paste(Index," ~ ", TimeName, "|", TreatmentName)
    plot1<-xyplot(as.formula(form), data=data, type="p")

    data[,TreatmentName] <- dataOriginalTreatmentColumn
    
    #now make plot with the control for all levels (so as to compare to the control at all times)
    indCtrl <- which(data[,TreatmentName]==1)
    ControlData <- data[indCtrl,]
    AllControlData <- ControlData
    trts <- unique(data[,TreatmentName])
    for (trt in trts[-1])
    {
      NewData <- ControlData
      NewData[,TreatmentName] <- trt
      AllControlData <- rbind(AllControlData, 
                              NewData)
    }
    plotC<-xyplot(as.formula(form), data=AllControlData, type="p",
                  pch=16,col="grey", ylim=c(MinY,MaxY))
    

    if (ttest!="None")
    {
      if (ttest=="t")
      {
        sigs <- t.test.WUR(data=data,TimeName=TimeName, 
                           TreatmentName=TreatmentName,Index=Index)
        sigs[,Index] <- MaxY-1
        sigs[,Index][which(sigs$p>0.05)] <- 1e10
      }
      if (ttest=="Williams") 
      {
        sigs <- williams.t.test.WUR(data=data,TimeName=TimeName, 
                                    TreatmentName=TreatmentName,
                                    Index=Index,order=ExpectedTrend)
        sigs[,Index] <- MaxY-1
        sigs[,Index][which(sigs$Significance==0)] <- 1e10
      }
      
      form <- paste(Index," ~ ", TimeName, "|", TreatmentName)
      plot11<-xyplot(as.formula(form), data=sigs, type="p", 
                     pch="*", cex=3, ylim=c(MinY,MaxY))
      
      print(plotC + as.layer(plot1) + as.layer(plot11))
    } 
    else
    {
      print(plotC + as.layer(plot1))
    }
  }
  dev.off()
}

#gets BD and EF using cosm data in the standard WUR file format
#also: gets EF in systems that are treated but without species loss
#Filename = name of txt file with count data
#TreatmentName = name in the count data file used to indicate treatment level. 
#...Always needs to be 'Treatment'. Control always needs to be '1'.
#TimeName = name in the count data file used to indicate time
#Only time points that are >Affected and <NoAffected are considered
#CountCols = columns where species counts are present: 
#...a vector with length = 2. If length = 1, 
#...it is assumed species counts start at CountCols and 
#...continue until the last column (ncol()). 
#CountCols should be set to NA when there are no counts in the dataframe
#endpoints= a list of endpoints. 
#...Note that only Richness or Evenness can be given
#...Similarity is calculated by default if counts are given
#...(should not be specified in endpoints)
#systemTag = either NULL (no tag) or a string pointing to the 
#...column that contains info on the system identity.
#...If this is included, similarities will only 
#...be calculated between sites that carry the same system tag.
#Binary = passed on to the vegdist function
#x = a fraction; if a species is present less than that fraction 
#...times the nr of observations it is thrown away
BDEF <- function(data, CountCols, TimeName, TreatmentName, 
                 Affected, NoAffected, 
                 endpoints = c("Richness", "Evenness"),
                 systemTag = NULL,
                 Binary=FALSE, 
                 x=0)
{
  data[,TreatmentName] <- as.factor(data[,TreatmentName])
  #do not replace NAs by zeros
  nas <- which(is.na(data))#,arr.ind=TRUE
  if (length(nas)>0) 
  {
    #data[is.na(data)] <- 0
    print("NAs in data; not replaced by zeros") 
  }
  
  #get only data between 'Affected' and 'NoAffected'
  dataPost <- data[which(data[,TimeName]>Affected),]
  dataPost <- dataPost[which(dataPost[,TimeName]<NoAffected),]

  #0/Test if the dataframe contains counts. 
  #If so, make an object 'counts' 
  #...and calculate similarity
  test <- is.na(CountCols)[1]
  if (test==FALSE)
  {
    if (length(CountCols)==1) {CountCols <- c(CountCols: ncol(dataPost))}
    #Take out real count data 
    counts <- dataPost[,CountCols]
    #remove species that occur less than x% of the time
    rare <- which(colSums(counts>0, na.rm=TRUE)<x*nrow(counts))
    if (length(rare) > 0) {counts <- counts[,-rare]}
    #remove empty sites
    empty <- which(as.numeric(rowSums(counts))==0)
    if (length(empty) > 0) 
    {
      counts <- counts[-empty,]
      dataPost <- dataPost[-empty,]
    }
    Sims <- data.matrix(vegdist(counts))
    if (Binary) {Sims <- data.matrix(vegdist(decostand(counts,
                                                       method="pa")))}
    diag(Sims) <- NA #To see why this is needed, check below
    Sims <- 1-Sims #cause vegdist calculates dissimilarity
    #now loop over all rows 
    #to get Sim with control at same time point
    dataPost$Sim <- NA
    for (row in c(1:nrow(dataPost)))
    {
      #look for row with same time point but treatment of 1
      Ind <- which((dataPost$Treatment==1)&(dataPost[,TimeName]==dataPost[row,TimeName]))
      #When 'row' is a control row (i.e. when dataPost$Treatment[row]==1),
      #...this 'row' will also be contained in Ind
      #...and you don't want this since you're comparing a row to itself 
      #...(always has similarity of 1),
      #...so that's why the diagonal of Sims need to be 'NA'.
      if (is.null(systemTag)==FALSE) 
      {
        Ind <- which((dataPost$Treatment==1)&(dataPost[,TimeName]==dataPost[row,TimeName])&(dataPost[,systemTag]==dataPost[row,systemTag]))
      }
      dataPost$Sim[row] <- mean(Sims[row,Ind], na.rm=TRUE)
    }
  }

  #1/Get or calculate richness
  if ("Richness" %in% endpoints)
  {
    #chck if already in data frame; if not do it yrself
    test <- ("Richness" %in% colnames(dataPost)) 
    if (test==TRUE) {richness <- dataPost[,"Richness"]} else {dataPost$Richness <- rowSums(counts>0)}
  }

  #2/Get or calculate evenness
  if ("Evenness" %in% endpoints)
  {
    #chck if already in data frame; if not do it yrself
    test <- ("Evenness" %in% colnames(dataPost)) 
    if (test==TRUE) {evenness <- dataPost[,"Evenness"]}
    else {dataPost$Evenness <- diversity(counts)}
  }
  return(dataPost)
}

RecoveryTime <- function(Effect, EffectCutOff=0.05, DaysBeforeEnd=31) 
  #Effect=vector with effect sizes (BiomassExposed/BiomassCtrl)
  #Cut-off value for effect size: any deviation from 1 that exceeds this value indicates (momentarily) recovery
  #DaysBeforeEnd = nr of days before simulation ends at which recovery should have occurred; if this is not the case, recovery is concluded not to have happenned
{
  EffectExceeded <- which((Effect<1-EffectCutOff)|(Effect>1+EffectCutOff))
  if (length(EffectExceeded)>0)
  {
    Onset <- min(EffectExceeded)
    Recovered <- Effect[c(Onset:length(Effect))]
    Recovered <- ((Recovered>1-EffectCutOff) + (Recovered<1+EffectCutOff)) - 1
    RecoveredTime <- cumsum(rev(Recovered))
    test2<-duplicated(RecoveredTime)
    ind<-min(which(test2==1))
    if (ind*timestep >= DaysBeforeEnd)
    {
      RecoveryTime <- length(Recovered)-ind
      RecoveryTime <- RecoveryTime*timestep         
    } else {
      RecoveryTime <- NA
    }
  }
  else {
    RecoveryTime <- NA
  }
}

Effects <- function(itrs, groups, affected, ref)
{
  effects <- 0
  
  for (i in c(1:itrs))
  {
    from <- (i-1)*groups + 1
    to <- from + groups - 1
    affected_i <- affected[,from:to]
    
    for (j in c(1:itrs))
    {
      from <- (j-1)*groups + 1
      to <- from + groups - 1
      effect <- affected_i / ref[,from:to]
      effects <- cbind(effects, effect)
    }
  }
  return(effects)
}

#function to return plotting positions in a window 
#with a margin left ("ShiftRight" is a nr between 0 and 1) 
#and up ("ShiftDown" is a nr between 0 and 1) 
#to insert titles etc
#ynr = nr of rows in the plot
#xnr = nr of cols in the plot
#ynr*xnr = total nr of plots to be drawn

PlottingPositions <- function(ShiftRight=0.1, 
                              ShiftDown=0, 
                              ShiftUp=0,
                              ynr=2, 
                              xnr=4)
{
  PlotX <- 1-ShiftRight
  fractsx1 <- c(0, rep(PlotX/xnr,(xnr-1))) * c(0:(xnr-1))
  fractsx1 <- ShiftRight + rep(fractsx1, ynr)
  fractsx2 <- rep(PlotX/xnr,xnr) * c(1:xnr)
  fractsx2 <- ShiftRight + rep(fractsx2, ynr)
  
  PlotY <- 1-ShiftDown-ShiftUp
  fractsy1 <- rep(PlotY/ynr,xnr)
  fractsy2 <- rep(PlotY,xnr)
  fractsy1 <- ShiftUp+c(fractsy1,rep(0,xnr))
  fractsy2 <- ShiftUp+c(fractsy2,rep(PlotY/ynr,xnr))
  
  return(FRACTS=list(fractsx1=fractsx1,
                     fractsx2=fractsx2,
                     fractsy1=fractsy1,
                     fractsy2=fractsy2))
}

######################################################
# Plots the predictions by a BD model, as returned by MI6
######################################################
#indsR= rows that will be considered
#preds = data frame with predictions
#COL = color to plot with
#LTY = line type tyo plot with
#Rel = Should relative abundances and ranks be plotted or just absolute values?
#conc_sel = which concentration to plot the SADs for?
#concs_all = all concentrations we have SADs for in the model output
#iterations = nr of itrs in the model output

PlotPreds <- function(indsR, preds, COL, LTY, Rel=FALSE, 
                      conc_sel, concs_all, iterations)
{
  meta.p <- NULL #some meta-data on the predictions: which column represents which concentration
  for (conc in concs_all)
  {
    meta.p    <- c(meta.p,rep(conc, iterations))
  }
  inds <- which(meta.p%in%conc_sel)   
  
  preds <- preds[indsR,]
  
  preds.i <- NULL
  
  for (i in inds)
  {
    pred.i <- as.numeric(preds[,i])
    if (length(unique(pred.i))==1) 
    {
      pred.i <- rep(0,nrow(preds))
    }
    pred.i <- sort(pred.i,decreasing=TRUE)+1
    if (Rel)
    {
      pred.i <- pred.i-1
      if (sum(pred.i)==0)
      {
        pred.i <- pred.i
      } else {
        pred.i <- pred.i / sum(pred.i)
      }
      
    }  
    preds.i <- rbind(preds.i,pred.i)
  }
  
  qpreds.i.05 <- apply(preds.i,2,quantile,0.025)
  qpreds.i.5 <- apply(preds.i,2,quantile,0.5)
  qpreds.i.95 <- apply(preds.i,2,quantile,0.975)
  
  if (Rel)
  {
    points(c(1:length(qpreds.i.05))/length(qpreds.i.05),qpreds.i.05+1e-3, type="l",
           lty=LTY, lwd=2, col=COL)
    points(c(1:length(qpreds.i.95))/length(qpreds.i.95),qpreds.i.95+1e-3, type="l",
           lty=LTY, lwd=2, col=COL) 
    #points(c(1:length(qpreds.i.5))/length(qpreds.i.5),qpreds.i.5+1e-3, type="l",
    #       lty="dotted", lwd=2, col=COL)
  } else {
    points(qpreds.i.05, type="l",
           lty=LTY, lwd=2, col=COL)
    points(qpreds.i.95, type="l",
           lty=LTY, lwd=2, col=COL)    
    #points(qpreds.i.5, type="l",
    #       lty="dotted", lwd=2, col=COL)
  }
  return(list(mean=(rowSums(preds.i>1))))
}

######################################################
# Plots the extinction risk predicted by a BD model, as returned by MI6
######################################################
#indsR= rows that will be considered
#preds = data frame with predictions
#COL = color to plot with
#LTY = line type to plot with
#conc_sel = which concentrationS to plot the extinction risk for?
#concs_all = all concentrations we have SADs for in the model output
#iterations = nr of itrs in the model output
#plotted = do the results need to be added to an existing plot?
#returned = do the results need to be passed on as a function output?
PlotExtinct <- function(indsR, preds, COL,  
                      conc_sel, concs_all, iterations, 
                      plotted=TRUE, returned=FALSE)
{
  preds <- preds[indsR,]
  PresAbs <- colSums((preds==0)*1)
  ExtinctionRisk <- PresAbs/nrow(preds)
  ExtinctionRisk[which(ExtinctionRisk==0)] <- 1 #perfectly even comms result from 100% effect
  
  meta.p <- NULL #some meta-data on the predictions: which column represents which concentration
  for (conc in concs_all)
  {
    meta.p    <- c(meta.p,rep(conc, iterations))
  }
  
  Qs <- NULL
  for (conc in conc_sel)
  {
    inds <- which(meta.p==conc)
    Qs <- rbind(Qs,quantile(ExtinctionRisk[inds], probs=c(0.05, 0.5, 0.95)))    
  }
  
  if (plotted)
  {
    points(conc_sel, Qs[,1], pch=21, cex=2)
    points(conc_sel, Qs[,2], pch=19, cex=2)
    points(conc_sel, Qs[,3], pch=21, cex=2)
  }
  
  if (returned)
  {
    return(Qs)
  }
}

######################################################
# Plots diversity (R or Simpson's E) across a gradient of toxicity and time, 
######################################################

BDDynamics <- function(Control=c("bu"),
                       Subs=c("herb"),
                       Connectance=c("005"),
                       TLs=c("Prey","Pred"),
                       pathPred, #where predictions are at
                       Initial, #data frame with initial E and R with columns Time, E, R, and TL. Should have interations *2 rows 
                       Flag=NA, #vector of size 1xiterations holding different integers, to be used as grouping (or conditional) factors in lattice plot; will not be used unless in formula! 
                       formAgg="E ~ Time + TL", #formula to aggregate the data (later on used for plotting)
                       formPlot="E ~ Time|TL", #formula used in xyplot of lattice
                       ymin=0,
                       concentration="",
                       ymax=25,
                       LWD=1,
                       col="black",
                       plotToo=FALSE) #if TRUE, the dynamics will be plotted as a pdf
{
  combs <- expand.grid(Control,Subs,Connectance)
  
  for (i in c(1:nrow(combs)))
  {
    for (TL in TLs)
    {
      Filename <- paste(pathPred,"Dyns",TL,concentration, 
                        combs[i,1],combs[i,2],
                        combs[i,3],sep="")
      Sims <- read.delim(Filename,header=TRUE, sep=" ")
      colnames(Sims) <- c("Time", "E", "R")
      Sims <- rbind(Initial[which(Initial$TL==TL),c("Time", "E", "R")],Sims)
      Sims$TL <- TL
      if (is.na(Flag[1]))
      {
        Flag <- rep(1,nrow(Sims))
      }
      Sims$Flag <- Flag
      if (match(TL,TLs)==1)
      {
        AllSims <- Sims
      } else
      {
        AllSims <- rbind(AllSims,Sims)
      }
    }
    
    colnames(AllSims) <- c("Time","E","R","TL","Flag")
    AllSims$Time <- AllSims$Time/1000
    
    Q50 <- aggregate(as.formula(formAgg), 
                     data=AllSims, FUN=quantile, 0.5)
    Q50$quantile <- "50"
    
    Q5 <- aggregate(as.formula(formAgg), 
                    data=AllSims, FUN=quantile, 0.05)
    Q5$quantile <- "5"
    
    Q95 <- aggregate(as.formula(formAgg), 
                     data=AllSims, FUN=quantile, 0.95)
    Q95$quantile <- "95"
    
    Qs <- rbind(Q50, Q5, Q95)
    Qs$quantile <- as.factor(Qs$quantile)
    
    p1 <- xyplot(as.formula(formPlot), 
                 type="l", col = c(col, col, col),
                 lwd = rep(LWD,3),
                 lty = c("dashed", "solid", "dashed"),
                 groups=quantile, 
                 data=Qs, xlab="Timesteps/1000",
                 horizontal=FALSE,ylim=c(ymin,ymax))
    
    if (plotToo)
    {
      quartz("",4,3,type="pdf",file=paste(combs[i,1],combs[i,2],
                                          combs[i,3],".pdf"))      
      print(p1)
      dev.off()      
    }
    return(p1)
  }
}

#Will look for a relationship between the nr of links to/from a species and its trait
#Output: a vector of size 1:iterations, holding integers that will be used as grouping (or conditional) factors while plotting with BDDynamics

FlagSearch <- function(Trait, #data frame with iterations columns and nr of rows= nr of species
                       FW, #fw matrix (nr of species x nr of species) 
                       Inds, #rows to focus on (ie the species)
                       Direction="from", #from or to: should the nr of links leaving or entering a species be accd for?
                       alpha=0.01,
                       trend=1, flag="Flag") #if 1, only negative relations will be considered. if -1, the opposite 
{
  iterations <- ncol(Trait)
  NrLinks <- rowSums(FW)[Inds]
  if (Direction=="to")
  {
    NrLinks <- colSums(FW)[Inds]
  }
  ps <- NULL
  estimates <- NULL
  
  for (i in c(1:iterations))
  {
    test <- lm(Trait~NrLinks, data=data.frame(Trait=Trait[Inds,i],NrLinks=NrLinks))
    ps <- c(ps,summary(test)$coefficients["NrLinks","Pr(>|t|)"])
    estimates <- c(estimates,summary(test)$coefficients["NrLinks","Estimate"])
  }
  
  Indices <- which((ps<alpha)&(trend*estimates<0))
  
  Flag <- rep("NoFlag",iterations)
  Flag[Indices] <- flag
  
  return(Flag)
  
}

CConnec <- function(Control=c("bu"),
                    Subs=c("herb"),
                    Connectance=c("005"),
                    Concentration="",
                    Times=seq(1000,10000,1000),
                    pathPred, #where predictions are at
                    Initial, #initial connectances in a data frame with Co and Time as colnames
                    formAgg="Co ~ Time") #formula to aggregate the data (later on used for plotting)
{
  combs <- expand.grid(Control,Subs,Connectance,Times) 
  AllSims <- NULL
  for (i in c(1:nrow(combs)))
  {
    
    Filename <- paste(pathPred,"ConnectanceJ-test",
                      Concentration,
                      combs[i,1],combs[i,2],
                      combs[i,3],combs[i,4],".out",sep="")
    Sims <- read.delim(Filename, header=FALSE, sep=" ")
    Sims <- cbind(combs[i,4], Sims)
    colnames(Sims) <- c("","")
    AllSims <- rbind(AllSims,Sims)
  }
  print(AllSims[c(1:2),])
  AllSims <- rbind(AllSims,Initial)
  colnames(AllSims) <- c("Time","Co")
  AllSims$Time <- AllSims$Time/1000
  Q50 <- aggregate(as.formula(formAgg), 
                   data=AllSims, FUN=quantile, 0.5)
  Q50$quantile <- "50"
  
  Q5 <- aggregate(as.formula(formAgg), 
                  data=AllSims, FUN=quantile, 0.05)
  Q5$quantile <- "5"
  
  Q95 <- aggregate(as.formula(formAgg), 
                   data=AllSims, FUN=quantile, 0.95)
  Q95$quantile <- "95"
  
  Qs <- rbind(Q50, Q5, Q95)
  Qs$quantile <- as.factor(Qs$quantile)
  
  return(Qs)
}

