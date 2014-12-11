# PART 1: Raw Data Processing
pred.days = 120
targetDaysPerRecord = 3
dat = read.csv("../Data/ebola/country_timeseries.csv")[,c(1, 3,4,5)]

charDate = as.character(dat[2:nrow(dat),1])

rptDate = as.Date(charDate, "%m/%d/%Y")
numDays = max(rptDate) - min(rptDate) + 1
numDays.pred = numDays + pred.days

original.rptDate = rptDate
ascendingOrder = order(rptDate)
rptDate = rptDate[ascendingOrder][2:length(rptDate)]
original.rptDate = original.rptDate[ascendingOrder]


cleanData = function(dataColumn, ascendingOrder)
{
  # Remove commas
  charCol = as.character(dataColumn)[ascendingOrder]
  
  if (is.na(charCol[1]))
  {
    charCol[1] = "0"
  }
  charCol = as.numeric(charCol)
  
  for (i in 2:length(charCol))
  {
    if (is.na(charCol[i]))
    {
      charCol[i] = charCol[i-1]  
    }
  }
  # charCol
  # Correct for undercounts
  for (i in seq(length(charCol), 2))
  {
    if (charCol[i-1] > charCol[i])
    {
      charCol[i-1] = charCol[i]
    }
  }
  charCol
}

Guinea = cleanData(dat[,2], ascendingOrder)
Liberia = cleanData(dat[,3], ascendingOrder)
Sierra.Leone = cleanData(dat[,4], ascendingOrder)
rawData = cbind(Guinea, Sierra.Leone, Liberia)
rownames(rawData) = as.character(original.rptDate)
colnames(rawData) = paste(paste("&nbsp;&nbsp;", c("Guinea", "Liberia", "Sierra Leone")), "&nbsp;&nbsp;")

# The data needs to be aggregated: there's some error in the measurements, 
# and we're not actually observing infection times. The data is therefore
# recorded at an artifically high time scale.
uncumulate = function(x)
{
  out = c(x[2:length(x)]-x[1:(length(x)-1)])
  ifelse(out >= 0, out, 0)
}
nDays = uncumulate(original.rptDate)

thinIndices = function(minDays, weights)
{
  keepIdx = c(length(weights))
  currentWeight = 0
  lastIdx = -1
  for (i in seq(length(weights)-1, 1))
  {
    currentWeight = currentWeight + weights[i]
    if (currentWeight >= minDays)
    {
      currentWeight = 0
      keepIdx = c(keepIdx, i)
      lastIdx = i
    }
  }
  if (currentWeight != 0)
  {
    keepIdx = c(keepIdx, lastIdx-1)
  }
  keepIdx
}

keepIdx = thinIndices(targetDaysPerRecord, c(1,nDays))
keepIdx = keepIdx[order(keepIdx)]



# Define the plot for the next section
ylim = c(min(c(Guinea, Sierra.Leone, Liberia)), 
         max(c(Guinea, Sierra.Leone, Liberia)))
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
  
  legend(x = original.rptDate[1], y = max(ylim), legend = 
           c("Guinea", "Liberia", "Sierra Leone"), 
         lty = 1:3, col = c("black", "blue","red"), bg="white", cex = 1.1,
         lwd=3)
}


Guinea = Guinea[keepIdx]
Sierra.Leone = Sierra.Leone[keepIdx]
Liberia = Liberia[keepIdx]

original.rptDate = original.rptDate[keepIdx]
rptDate = original.rptDate[2:length(original.rptDate)]
rptDate = rptDate[2:(length(rptDate)-1)]


# PART 2: Load Chain Output

load("./DF0/chainOutput.Robj")
DF0 = chains
load("./DF3/chainOutput.Robj")
DF3 = chains
load("./DF5/chainOutput.Robj")
DF5 = chains

load("./DF0_CAR/chainOutput.Robj")
DF0_CAR = chains
load("./DF3_CAR/chainOutput.Robj")
DF3_CAR = chains
load("./DF5_CAR/chainOutput.Robj")
DF5_CAR = chains


nTpt = nrow(DF0[[2]]$R0$empiricalR0$mean) -1

EA_R_DF0 = DF0[[2]]$R0$empiricalR0
EA_R_DF3 = DF3[[2]]$R0$empiricalR0
EA_R_DF5 = DF5[[2]]$R0$empiricalR0

R0_DF0 = DF0[[2]]$R0$R0
R0_DF3 = DF3[[2]]$R0$R0
R0_DF5 = DF5[[2]]$R0$R0

EA_R_DF0_CAR = DF0_CAR[[2]]$R0$empiricalR0
EA_R_DF3_CAR = DF3_CAR[[2]]$R0$empiricalR0
EA_R_DF5_CAR = DF5_CAR[[2]]$R0$empiricalR0

R0_DF0_CAR = DF0_CAR[[2]]$R0$R0
R0_DF3_CAR = DF3_CAR[[2]]$R0$R0
R0_DF5_CAR = DF5_CAR[[2]]$R0$R0



color1 = rgb(0,0,0,0.9)
color2 = rgb(0,0,1,0.9)
lwdMean = 3
lwdCI = 1
ltyCI=2
lty1 = 1
lty2 = 4
pch1=1
pch2=3
ltyThreshold = 5
pointCex = 0
colThreshold = rgb(0.5,0.5,0.5,0.5)
fig1Ylim=c(0,2)

R0ComparisonPlot = function(EARList,R0List,colIdx,main, ...)
{  
  EARMean = EARList$mean[2:nTpt,colIdx]
  EARUB = EARList$UB[2:nTpt,colIdx]
  EARLB = EARList$LB[2:nTpt,colIdx]
  
  R0Mean = R0List$mean[2:nTpt,colIdx]
  R0UB = R0List$UB[2:nTpt,colIdx]
  R0LB = R0List$LB[2:nTpt,colIdx]
  
  plot(rptDate,EARMean,
        type = "n", ...)
  title(main=main, line=-2.5)
  abline(h = 1, lty = ltyThreshold, col = colThreshold)
  lines(rptDate, EARMean, col = color1, lwd = lwdMean, lty = lty1,)
  
  points(rptDate, EARMean, col = color1, pch=pch1, cex = pointCex)
  lines(rptDate,EARUB, col = color1, lwd = lwdCI, lty = ltyCI)
  lines(rptDate,EARLB, col = color1, lwd = lwdCI, lty = ltyCI)
  
  lines(rptDate,R0Mean,  col = color2, lwd = lwdMean, lty = lty2)
  points(rptDate, R0Mean, col = color2, pch = pch2, cex = pointCex)
  lines(rptDate,R0UB, col = color2, lwd = lwdCI, lty = ltyCI)
  lines(rptDate,R0LB, col = color2, lwd = lwdCI, lty = ltyCI)
 
}

figureFunc = function(EARList, R0List, titleVal, ylim)
{
  par("oma" = c(1,1,1,1))
  layout(matrix(c(1,2,3,
                  4,4,4), 2, 3, byrow = TRUE),
         widths=c(3,3,3),
         heights = c(5,1))
  
  R0ComparisonPlot(
                    EARList,
                    R0List,
                    1,
                    main = "Guinea",
                    xlab = "Date",
                    ylab = "Reproductive Number",
                    ylim=ylim)
  
  R0ComparisonPlot( EARList,
                     R0List,
                     2,
                    main = "Liberia",
                    xlab = "Date",
                    ylab = "Reproductive Number",
                    ylim=ylim)
  
  R0ComparisonPlot( EARList,
                     R0List,
                     3,
                     main = "Sierra Leone",
                     xlab = "Date",
                     ylab = "Reproductive Number",
                     ylim = ylim)
  par(xaxt="n")
  par(yaxt="n")
  par(xpd=TRUE)
  par(bty="n")
  plot(c(0,10),c(0,10), type = "n",xlab="",ylab="")
  legend(x = 5, 
         y = 5, 
         lty = c(lty1, lty2), 
         col = c(color1, color2), 
         legend=c(expression('R'^{(EA)}), expression('R'[0](t))),
         lwd = c(lwdMean, lwdMean), 
         pt.cex = c(pointCex,pointCex),
         pt.lwd= c(pointCex, pointCex),
         pch = c(pch1, pch2),
         cex=3,
         xjust=0.5,
         yjust=0.5, 
         horiz = TRUE,
         text.width=c(1))
  par(xaxt="s")
  par(yaxt="s")
  par(xpd=FALSE)
  par(bty="o")
  par(cex.main=1.5)
  title(titleVal, outer = TRUE, line = -2)
  par(cex.main=1)
} 



pdf(file="./0DF_Pairwise.pdf", width = 12, height =8)
  figureFunc(EA_R_DF0, R0_DF0, "Underspecified Model\n 0 DF Temporal Basis, Pairwise Spatial Structure", c(0, 3))
dev.off()

pdf(file="./3DF_Pairwise.pdf", width = 12, height =8)
figureFunc(EA_R_DF3, R0_DF3,"Three DF Temporal Basis Model\n Pairwise Spatial Structure", c(0,3))
dev.off()

pdf(file="./5DF_Pairwise.pdf", width = 12, height =8)
figureFunc(EA_R_DF5, R0_DF5, "Five DF Temporal Basis Model\n Pairwise Spatial Structure", c(0,3))
dev.off()





pdf(file="./0DF_CAR.pdf", width = 12, height =8)
figureFunc(EA_R_DF0_CAR, R0_DF0_CAR, "Underspecified Model\n 0 DF Temporal Basis, CAR Spatial Structure", c(0, 3))
dev.off()

pdf(file="./3DF_CAR.pdf", width = 12, height =8)
figureFunc(EA_R_DF3_CAR, R0_DF3_CAR,"Three DF Temporal Basis Model\n Pairwise CAR Structure", c(0,3))
dev.off()

pdf(file="./5DF_CAR.pdf", width = 12, height =8)
figureFunc(EA_R_DF5_CAR, R0_DF5_CAR, "Five DF Temporal Basis Model\n Pairwise CAR Structure", c(0,3))
dev.off()




