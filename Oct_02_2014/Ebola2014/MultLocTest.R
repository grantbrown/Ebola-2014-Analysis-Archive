
library(knitr)
library(coda) # Load the coda library for MCMC convergence diagnosis
library(spatialSEIR) # Load the spatialSEIR library to perform the modeling. 
library(XML) # Load the XML library to read in data from Wikipedia
library(parallel) # Load the parallel library to enable multiple chains to be run simultaneously. 

## Define Document Compilation Parameters

#documentCompilationMode = "release"
documentCompilationMode = "debug"
modelDF = 6
pred.days = 14

## Compute number of samples/batches
numBurnInBatches =      ifelse(documentCompilationMode == "release", 1000,  1)
numConvergenceBatches = ifelse(documentCompilationMode == "release", 1000,  10)
convergenceBatchSize =  ifelse(documentCompilationMode == "release", 1000, 50) 
extraR0Iterations =     ifelse(documentCompilationMode == "release", 500,   10)
iterationStride =       ifelse(documentCompilationMode == "release", 1000,   50)

## Read in the data
url = 'http://en.wikipedia.org/wiki/West_Africa_Ebola_virus_outbreak'
tbls = readHTMLTable(url)
dat = tbls[[5]] # This line changes depending on the page formatting.

# One date is now duplicated on the Wikipedia page, due to using different sources. 
# Clean that up first. 

dup.indices = which(as.Date(dat[,1], "%d %b %Y") == as.Date("2014-06-05"))
dat[dup.indices[1],]$V2 = dat[dup.indices[2],]$V2
dat = rbind(dat[1:dup.indices[1],], dat[(dup.indices[2]+1):nrow(dat),])

charDate = as.character(dat[2:nrow(dat),1])
for (i in 1:length(charDate))
{
  charDate[i] = gsub("Sept", "Sep", charDate[i])
}

rptDate = as.Date(charDate, "%d %b %Y")
numDays = max(rptDate) - min(rptDate) + 1
numDays.pred = numDays + pred.days

original.rptDate = rptDate
ascendingOrder = order(rptDate)
rptDate = rptDate[ascendingOrder][2:length(rptDate)]
original.rptDate = original.rptDate[ascendingOrder]


cleanData = function(dataColumn, ascendingOrder)
{
    # Remove commas
    dataColumn = gsub(",", "", dataColumn, fixed = TRUE)
    # Remove +1 -1 stuff
    charCol = as.character(
      lapply(
        strsplit(
          as.character(
            dataColumn)[ascendingOrder], "\n"), function(x){ifelse(length(x) == 0, "—", x[[1]])}
        )
      )
    if (is.na(charCol[1]) || charCol[1] == "—")
    {
      charCol[1] = "0"
    }
    charCol = as.numeric(ifelse(charCol == "—", "", charCol))
    for (i in 2:length(charCol))
    {
      if (is.na(charCol[i]))
      {
        charCol[i] = charCol[i-1]  
      }
    }
    charCol
}

Guinea = cleanData(dat$V4[2:nrow(dat)], ascendingOrder)
Liberia = cleanData(dat$V6[2:nrow(dat)], ascendingOrder)
Sierra.Leone = cleanData(dat$V8[2:nrow(dat)], ascendingOrder)
Nigeria = cleanData(dat$V10[2:nrow(dat)], ascendingOrder)

# Define the plot for the next section
ylim = c(min(c(Guinea, Sierra.Leone, Liberia, Nigeria)), 
             max(c(Guinea, Sierra.Leone, Liberia, Nigeria)))
figure1 = function()
{
      plot(original.rptDate, Guinea, type = "l", 
           main = "Raw Data: Case Counts From Wikipedia", 
           xlab = "Date", 
           ylab = "Total Cases",
           ylim = ylim, lwd = 3)
      abline(h = seq(0,100000, 100), lty = 2, col = "lightgrey")
      lines(original.rptDate, Liberia, lwd = 3, col = "blue", lty = 2)
      lines(original.rptDate, Sierra.Leone, lwd = 3, col = "red", lty = 3)
      lines(original.rptDate, Nigeria, lwd = 3, col = "green", lty = 4)
      legend(x = original.rptDate[1], y = max(ylim), legend = 
               c("Guinea", "Liberia", "Sierra Leone", "Nigeria"), 
             lty = 1:3, col = c("black", "blue","red", "green"), bg="white", cex = 1.1)
}

uncumulate = function(x)
{
    out = c(x[2:length(x)]-x[1:(length(x)-1)])
    ifelse(out >= 0, out, 0)
}
# The "I_star" name will make more sense in a bit
I_star = cbind(uncumulate(Guinea), 
               uncumulate(Liberia), 
               uncumulate(Sierra.Leone),
               uncumulate(Nigeria))

maxIdx = nrow(I_star)

# Define the temporal offset vector to be the number of days reflected in each 
# aggregated record (time between reports).
offsets = uncumulate(original.rptDate)
if (any(offsets <= 0))
{
    cat("Invalid Date Information. The data source has likely changed.\n")
    stop(-1)
}

InfectionPerDay = I_star/(cbind(offsets, offsets, offsets, offsets))

# Define figure 2 for next section
ylim = c(0,max(InfectionPerDay)*1.2)
figure2 = function()
{
    layout(matrix(c(1,2), nrow = 1),
        widths = c(8,4), heights = c(4,4))
    plot(rptDate, InfectionPerDay[,1], main = "Crude Guess at New Case Counts Per Day", 
         xlab = "Date", 
         ylab = "New Cases",
         lty=1, lwd = 2,
         ylim = ylim, type = "l"
         )
    abline(h = seq(0, 1000, 5), lty = 2, col = "lightgrey")
    lines(rptDate, InfectionPerDay[,2], col = "blue",lty=2, lwd = 2)
    lines(rptDate, InfectionPerDay[,3], col = "red", lty = 3, lwd = 2)
    lines(rptDate, InfectionPerDay[,4], col = "green", lty = 4, lwd = 2)
    par(xaxt="n")
    par(yaxt="n")
    par(bty="n")
    par(xpd=TRUE)
    plot(c(0,10),c(0,10), type = "n", main  ="",xlab="",ylab="")
    legend(x=-2,y=10, legend = c("Guinea", "Liberia", "Sierra Leone", "Nigeria"), lty = 1:4,lwd=2, 
           col = c("black", "blue", "red", "green"))
    par(xpd=FALSE)
    par(xaxt="s")
    par(yaxt="s")
    par(bty="o")
}

library(rCharts)
library(stats)

x = rptDate - min(rptDate)
guinea.interp = approx(x,InfectionPerDay[,1],xout = 0:max(x))
liberia.interp = approx(x,InfectionPerDay[,2],xout = 0:max(x))
sierraleone.interp = approx(x,InfectionPerDay[,3],xout = 0:max(x))
nigeria.interp = approx(x,InfectionPerDay[,4],xout = 0:max(x))


interpMatrix = cbind(guinea.interp$y, liberia.interp$y,sierraleone.interp$y, nigeria.interp$y)
cutvals = cut(interpMatrix, breaks = 9)
interpMatrix.cut = matrix(as.numeric(cutvals), nrow = nrow(interpMatrix))


upperVals = as.numeric(lapply(strsplit(c(gsub("[(]", "", gsub("]", "", unique(as.character(cutvals))))), ","), function(x){return(x[2])}))
upperVals = round(upperVals[order(upperVals)],0)


hcol = c("#ffffef", "#fff7bf", "#fee39f", "#fec45f", "#fe993f", "#ec702f", "#cc4c1f", "#993402", "#662520")
color.palette = c(hcol[1],hcol)
fills = setNames(color.palette, c("defaultFill", paste("lt", upperVals, sep = "")))


# GIN, LBR, SLE, 
outList = list()
for (tpt in min(x):max(x))
{
    outList[[as.character(tpt+1)]] = list("GIN" = list("fillKey"=factor(paste("lt", upperVals[interpMatrix.cut[tpt+1,1]], sep =""), 
                                                        levels = names(fills))),
                                          "LBR" = list("fillKey"=factor(paste("lt", upperVals[interpMatrix.cut[tpt+1,2]], sep = ""), 
                                                        levels = names(fills))),
                                          "SLE" = list("fillKey"=factor(paste("lt",upperVals[interpMatrix.cut[tpt+1,3]], sep = ""), 
                                                        levels = names(fills))),
                                          "NGA" = list("fillKey"=factor(paste("lt",upperVals[interpMatrix.cut[tpt+1,4]], sep = ""), 
                                                        levels = names(fills))))
}

cat(numDays)
library(splines)
daysSinceJan = as.numeric(rptDate - as.Date("2014-01-01"))
daysSinceJan.predict = c(max(daysSinceJan) + 1, max(daysSinceJan) 
                         + seq(2,pred.days-2,2))
splineBasis = ns(daysSinceJan, df = modelDF)
splineBasis.predict = predict(splineBasis, daysSinceJan.predict)

# Guinea, Liberia, Sierra Leone, Nigeria
N = matrix(c(10057975, 4128572, 6190280,174507539), nrow = nrow(I_star),ncol = 4, 
           byrow=TRUE)
X = diag(ncol(N))
X.predict = X
Z = splineBasis
Z.predict = splineBasis.predict

# These co-variates are the same for each spatial location, 
# so duplicate them row-wise. 
Z = Z[rep(1:nrow(Z), nrow(X)),]
Z.predict = Z.predict[rep(1:nrow(Z.predict), nrow(X)),]

# For convenience, let's combine X and Z for prediction.
X.pred = cbind(X.predict[rep(1:nrow(X.predict), 
                             each = nrow(Z.predict)/nrow(X)),], Z.predict)


# Define the simple "distance" matrices. There are 4 countries, three of which 
# share borders. Let all of the nations share an overall correlation term, and
# let the three neighboring nations share an additional parameter
DM1 = (1-diag(4))
DM2 = rbind(cbind((1-diag(3)), rep(0, 3)), rep(0,4))
# Make row stochastic:
DM1 = DM1/matrix(apply(DM1, 1, sum), nrow = 4, ncol = 4)
DM2 = DM2/matrix(apply(DM2, 1, sum), nrow = 4, ncol = 4)
DM2 = ifelse(is.na(DM2), 0, DM2)

# Define population sizes for the three countries of interest. This data also 
# from Wikipedia. 

# Define prediction offsets. 
offset.pred = c(1,seq(2,pred.days-2,2))

# There's no reinfection process for Ebola, but we still need to provide dummy
# values for the reinfection terms. This will be changed (along with most of 
# the R level API) Dummy covariate matrix:
X_p_rs = matrix(0)
# Dummy covariate matrix dimension. Why, exactly, am I not just grabbing this 
# kind of thing from Rcpp? No good reason at all: this will be fixed. 
xPrsDim = dim(X_p_rs)
# Dummy value for reinfection params
beta_p_rs = rep(0, ncol(X_p_rs))
# Dummy value for reinfection params prior precision
betaPrsPriorPrecision = 0.5

# Get object dimensions. Again, this will be done automatically in the future
compMatDim = dim(I_star)
xDim = dim(X)
zDim = dim(Z)
    
# Declare prior parameters for the E to I and I to R probabilities. 
priorAlpha_gammaEI = 25;
priorBeta_gammaEI = 100;
priorAlpha_gammaIR = 14;
priorBeta_gammaIR = 100;
# Declare prior parameters for the overdispersion precision
priorAlpha_phi = 1
priorBeta_phi = 1

# Declare prior precision for exposure model paramters
betaPriorPrecision = 0.1

# Declare a function which can come up with several different starting values 
# for the model parameters. This will allow us to assess convergence. 
proposeParameters = function(seedVal, chainNumber)
{
    set.seed(seedVal) 
    
    # 2 to 21 day incubation period according to who
    p_ei = 0.25 + rnorm(1, 0, 0.02) 
    # Up to 7 weeks even after recovery
    p_ir = 0.14 + rnorm(1, 0, 0.01) 
    gamma_ei=-log(1-p_ei)
    gamma_ir=-log(1-p_ir)
      
    # Starting value for exposure regression parameters
    beta = rep(0, ncol(X) + ncol(Z))
    beta[1] = 2.5 + rnorm(1,0,0.5)
    
    rho = 0.1 + rnorm(1,0,0.01) # spatial dependence parameter
    phi = 0.01 # Overdispersion precision
    
    outFileName = paste("./chain_output_ebola_", chainNumber ,".txt", sep = "")
    
    # Make a crude guess as to the true compartments:
    # S_star, E_star, R_star, and thus S,E,I and R
    DataModel = buildDataModel(I_star, type = "overdispersion", 
                               params=c(priorAlpha_phi,priorBeta_phi))
    ExposureModel = buildExposureModel(X, Z, beta, betaPriorPrecision)
    ReinfectionModel = buildReinfectionModel("SEIR")
    SamplingControl = buildSamplingControl(iterationStride=iterationStride,
                                           sliceWidths=c(1, # S_star
                                                         1, # E_star
                                                         1,  # R_star
                                                         1,  # S_0
                                                         1,  # I_0
                                                         0.05,  # beta
                                                         0.0,  # beta_p_rs, fixed in this case
                                                         0.01, # rho
                                                         0.01, # gamma_ei
                                                         0.01,  # gamma_ir
                                                         0.01)) # phi)
    
    InitContainer = buildInitialValueContainer(I_star, N, 
                                               S0 = N[1,]-I_star[1,] - c(86,0,0,0),
                                               I0 = c(86,0,0,0), 
                                               p_ir = 0.5, 
                                               p_rs = 0.00)
    DistanceModel = buildDistanceModel(list(DM1,DM2))
    TransitionPriors = buildTransitionPriorsManually(priorAlpha_gammaEI,priorBeta_gammaEI,
                                                     priorAlpha_gammaIR,priorBeta_gammaIR)
    return(list(DataModel=DataModel,
                ExposureModel=ExposureModel,
                ReinfectionModel=ReinfectionModel,
                SamplingControl=SamplingControl,
                InitContainer=InitContainer,
                DistanceModel=DistanceModel,
                TransitionPriors=TransitionPriors,
                outFileName=outFileName))
}

params =list("estimateR0"=TRUE, "traceCompartments"=TRUE, "seedVal"=12312334,"chainNumber"=4)
buildAndRunModel = function(params)
{
  library(spatialSEIR)
  proposal = proposeParameters(params[["seedVal"]], params[["chainNumber"]])
  SEIRmodel =  buildSEIRModel(proposal$outFileName,
                              proposal$DataModel,
                              proposal$ExposureModel,
                              proposal$ReinfectionModel,
                              proposal$DistanceModel,
                              proposal$TransitionPriors,
                              proposal$InitContainer,
                              proposal$SamplingControl)
  
  SEIRmodel$setRandomSeed(params[["seedVal"]])
  
  # Do we need to keep track of compartment values for prediction? 
  # No sense doing this for all of the chains.
  if (params[["traceCompartments"]])
  {
    SEIRmodel$setTrace(0) #Guinea 
    SEIRmodel$setTrace(1) #Liberia
    SEIRmodel$setTrace(2) #Sierra Leone
    SEIRmodel$setTrace(3) #Nigeria
  }
      
  # Make a helper function to run each chain, as well as update the metropolis 
  # tuning parameters. 
  runSimulation = function(modelObject,
                           numBatches=500, 
                           batchSize=20, 
                           targetAcceptanceRatio=0.2,
                           tolerance=0.05,
                           proportionChange = 0.1
                          )
  {
      for (batch in 1:numBatches)
      {
          modelObject$simulate(batchSize)
          modelObject$updateSamplingParameters(targetAcceptanceRatio, 
                                               tolerance, 
                                               proportionChange)
      }
  }

  # Burn in tuning parameters
  runSimulation(SEIRmodel, numBatches = numBurnInBatches)
  SEIRmodel$compartmentSamplingMode = 14
  SEIRmodel$useDecorrelation = 25
  # Run Simulation
  cat(paste("Running chain ", params[["chainNumber"]], "\n", sep =""))

  tm = 0


  tm = tm + system.time(runSimulation(SEIRmodel, 
	    numBatches=numConvergenceBatches, 
	    batchSize=convergenceBatchSize, 
	    targetAcceptanceRatio=0.2,
	    tolerance=0.025,
	    proportionChange = 0.05))


  
  cat(paste("Time elapsed: ", round(tm[3]/60,3), 
              " minutes\n", sep = ""))
  dat = read.csv(proposal$outFileName)
  
  ## Do we need to estimate R0 for this chain?
  if (params[["estimateR0"]])
  {  
    R0 = array(0, dim = c(nrow(I_star), ncol(I_star), extraR0Iterations))
    for (i in 1:extraR0Iterations)
    {
        SEIRmodel$simulate(iterationStride)
        for (j in 0:(nrow(I_star)-1))
        {
            R0[j,,i] = apply(SEIRmodel$getIntegratedGenerationMatrix(j), 2, sum)
        }
    }
    
    R0Mean = apply(R0, 1:2, mean)
    R0LB = apply(R0, 1:2, quantile, probs = 0.05)
    R0UB = apply(R0, 1:2, quantile, probs = 0.95)
    orig.R0 = R0
    R0 = list("mean"=R0Mean, "LB" = R0LB, "UB" = R0UB)
  } else
  {
     R0 = NULL
     orig.R0 = NULL
  }  
  
  return(list("chainOutput" = dat, "R0" = R0, "rawSamples" = orig.R0))
}



mod = buildAndRunModel(params)
