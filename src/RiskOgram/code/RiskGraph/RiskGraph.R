library(RMySQL)
library(gplots)

con <- dbConnect(dbDriver("MySQL"), dbname="proj_dz_risk_indian",username="",  host="", password='');

## Functions

mergeViaRowNames <- function(a, b, aName="", bName="", all.a=F, all.b=F) {
	# This function takes in two matrices and then merges them by their rownames
	# It can take in optional arguments that pre-pend text to the column names of the parts, useful for keeping track of them
	a <- data.frame(a)
	b <- data.frame(b)
	if (aName != "") {colnames(a) <- paste(aName, colnames(a), sep="")}
	if (bName != "") {colnames(b) <- paste(bName, colnames(b), sep="")}
	a <- data.frame(mindex=rownames(a),a )
	b <- data.frame(mindex=rownames(b), b)
	m <- merge(a, b, by.x="mindex", by.y="mindex", all.x=all.a, all.y=all.b)
	rownames(m) <- m[,1]
	m <- m[,2:(dim(m)[2])]
	return(m)
}

preToPost <- function(LR, preTestProb) {
  preTestOdds <- preTestProb/(1-preTestProb);
  postTestOdds <- LR * preTestOdds
  return(postTestOdds/(1+postTestOdds))
}


drawRiskOverallGraph <- function(Disease, Prevalence, Cumulative, SNPCount, ratioLimit=1.5, txtScale=50, prevSize=.1, legcex=.5, showOnlyAggregate=F, segments=F, showNumerics=F, segColor="gray50", ...) {
  #  Takes in four vectors, and then
  Data <- data.frame(Disease=Disease, Prevalence=Prevalence, Cumulative=Cumulative, SNPCount=SNPCount)
  Data$LogPrev <- log10(Data$Prevalence)
  Data$LogCum <- log10(Data$Cumulative)
   #Data$Color <- rep("blue", NROW(Data))
  #Data$Color[(Data$LogCum - Data$LogPrev) > log10(ratioLimit) ] <- "darkred"
  #Data$Color[(Data$LogPrev - Data$LogCum) > log10(ratioLimit) ] <- "darkgreen"
  Data$Color[Data$LogCum > Data$LogPrev ] <- "orange"
  Data$Color[Data$LogPrev > Data$LogCum] <- "blue"
  Data$Color[Data$LogPrev == Data$LogCum] <- "black"
  par(oma=c(0, 0, 0, 0));
  plot.new()
  if (! showNumerics) {
    Data$RightText <- Data$Disease
  } else {
    Data$RightText <- paste(format(Data$Disease), format(Data$SNPCount), sep="  ")

  }
  #print(Data$RightText[1:20])
  TextWidth <- max(strwidth(Data$RightText))*txtScale
  par(mar=c(0.5, 0.2, 0.5, TextWidth))
  CaptionHeight = 1
  minProb <- floor(min(c(Data$LogPrev, Data$LogCum, Data$LogMax)))
  plot.window(xlim=c(minProb, 0), ylim=c(-CaptionHeight*.5, NROW(Data)))

  ticks <- rep(seq(minProb, 0, 1), c(rep(2, abs(minProb)), 1))  + c(rep(c(0, log10(5)), abs(minProb)), 0)
  mticks <-  seq(minProb, 0, 1)
  
  #ticks <- seq(minProb, 0, 1)
  y <- 1:NROW(Data);
  if (segments) {  segments(minProb, y, 0, y, lty="dotted", col=segColor) }

  
  Row <- 1
  while(Row <= NROW(Data)) {
    # plot dotted line
    if (!  showOnlyAggregate) {
      lines(c(Data$LogPrev[Row], Data$LogCum[Row]), c(Row, Row), col=Data$Color[Row], lty="dotted", ...)
      lines(c(Data$LogPrev[Row], Data$LogMax[Row]), c(Row, Row), col=Data$Color[Row], lty="solid", ...)
      lines(c(Data$LogPrev[Row], Data$LogPrev[Row]), c(Row+prevSize, Row-prevSize), lwd=3, col="black")
    #points(Data$LogPrev[Row], Row, cex=1.5, col="black", pch="|", lwd=3)
    #points(Data$LogCum[Row], Row, col="black", pch=20)
      points(Data$LogMax[Row], Row, col="black", pch=4)
    } else {
       lines(c(Data$LogPrev[Row], Data$LogCum[Row]), c(Row, Row), col=Data$Color[Row], lty="solid", ...)
       #lines(c(Data$LogPrev[Row], Data$LogPrev[Row]), c(Row+prevSize, Row-prevSize), lwd=3, col="black")
       ## Hershey font triangle version:
       #if (Data$LogPrev[Row] < Data$LogCum[Row]) {
       #  trig <- "\\#H0853" 
       #} else {trig <- "\\#H0855" }
       #text(Data$LogCum[Row], Row, trig, col="black", vfont=c("sans serif", "plain"), cex=1.5)

       #polygon triangle version
       if (Data$LogPrev[Row] < Data$LogCum[Row]) {
         trigp <-  Data$LogPrev[Row] + (prevSize/2)
       } else {
          trigp <-  Data$LogPrev[Row] -(prevSize/2)
       }
       polygon(c(Data$LogPrev[Row], Data$LogPrev[Row], trigp), c(Row+prevSize, Row-prevSize, Row), col="black")
       
    }
    
    Row <- Row+1;
    }
  # vertical lines
  lines(c(0, 0), c(0, NROW(Data)));  # line for the right side
  for (i in ticks[1:(length(ticks)-1)]) {
    lines(c(i, i),  c(1-.25, NROW(Data)), lty="dotted", col=segColor)
  }

  
  mtext(Data$RightText, at=y, adj=0, side=4, las=2, cex=legcex, family="mono")
  axis(1, at=mticks, paste(labels=10^(mticks)*100, "%", sep=""), pos=0.5, cex.axis=.8, mex=0.5)

}


drawOverallLinear<- function(Disease, Prevalence, Cumulative, SNPCount, ratioLimit=1.5, txtScale=50, prevSize=.1, legcex=.5, showOnlyAggregate=F, segments=F, showNumerics=F, segColor="gray50", ...) {
  Data <- data.frame(Disease=Disease, Prevalence=Prevalence, Cumulative=Cumulative,  SNPCount=SNPCount)
  Data$Color[Data$Cumulative > Data$Prevalence ] <- "orange"
  Data$Color[Data$Prevalence > Data$Cumulative] <- "blue"
  Data$Color[Data$Prevalence == Data$Cumulative] <- "black"
  par(oma=c(0, 0, 0, 0));
  plot.new()
  if (! showNumerics) {
    Data$RightText <- Data$Disease
  } else {
    Data$RightText <- paste(format(Data$Disease), format(Data$SNPCount), sep="  ")

  }
  TextWidth <- max(strwidth(Data$RightText))*txtScale
  par(mar=c(0.5, 0.2, 0.5, TextWidth))
  CaptionHeight = 1
  plot.window(xlim=c(0, 1), ylim=c(-CaptionHeight*.5, NROW(Data)))
  
  ticks <- seq(0, 1, .1)

  y <- 1:NROW(Data);
  if (segments) {  segments(0, y, 1, y, lty="dotted", col=segColor) }

  Row <- 1
  while(Row <= NROW(Data)) {
    #lines(c(Data$Prevalence[Row], Data$Cumulative[Row]), c(Row, Row), col=Data$Color[Row], lty="solid", lwd=3)
    lines(c(Data$Prevalence[Row], Data$Cumulative[Row]), c(Row, Row), col=Data$Color[Row], lty="solid", ...)
    if (Data$Prevalence[Row] < Data$Cumulative[Row]) {
      trigp <-  Data$Prevalence[Row] + (prevSize/10)
    } else {
      trigp <-  Data$Prevalence[Row] -(prevSize/10)
    }
    polygon(c(Data$Prevalence[Row], Data$Prevalence[Row], trigp), c(Row+prevSize, Row-prevSize, Row), col="black")
      
    Row <- Row+1;
  }
  
  # vertical lines
  lines(c(1, 1), c(0, NROW(Data)));  # line for the right side
  for (i in ticks[1:(length(ticks)-1)]) {
    lines(c(i, i),  c(1-.25, NROW(Data)), lty="dotted", col=segColor)
  }

  
  mtext(Data$RightText, at=y, adj=0, side=4, las=2, cex=legcex, family="mono", col="black", font=1)
  axis(1, at=ticks, paste(labels=(ticks)*100, "%", sep=""), pos=0.5, cex.axis=.8, mex=0.5)

}


colFunction <- function(x) {
  if (x > 2) {return("black")}
  if (x == 2) {return("gray30")}
  else {return("gray70")}
}

getMaxAbs <- function(v) {
  maxabs <- max(abs(v))
  v[abs(v)==maxabs][1]
}

getMaxAbsLog <- function(v) {
  maxabs <- max(abs(log(v)))
  v[abs(log(v))==maxabs][1]
}

riskOGramSingleCondition <- function(Prevalence, RiskDF, condition, txtScale=30, segColor="gray50") {
  # RiskDF <- riskAlleles[riskAlleles$broad_phenotype == "Myocardial infarction",]
  PreTestOdds <- Prevalence/(1-Prevalence)
  
  RiskDF <- RiskDF[order(RiskDF$study_cnt_total, RiskDF$sample_size, decreasing=T),]
  RiskDF$LeftText <- paste(RiskDF$symbol, paste("rs", RiskDF$dbSNP, sep=""), RiskDF$genotype)
  RiskSF <- RiskDF[,c("LeftText", "LR", "study_cnt_total", "sample_size")]
  RiskSF <- rbind(list("Prevalence", PreTestOdds, NA, NA), RiskSF)
  #RiskSF$LogLR <- log10(RiskSF$LR)
  RiskSF$CumLR <- cumprod(RiskSF$LR)
  RiskSF$CumProb <- RiskSF$CumLR/(1+RiskSF$CumLR)


  # Only reverse before drawing!
  RiskSF <- RiskSF[NROW(RiskSF):1,]


  # Generate the margin text
  OverallLeftText <- c(RiskSF$LeftText, "Genotype Test")

  # Note that these are reversed!!!
  rtLR <- c(format(RiskSF$LR, digit=2, width=6)[1:(NROW(RiskSF)-1)], "", "LR")
  rtSC <- c(format(RiskSF$study_cnt_total, width=6)[1:(NROW(RiskSF)-1)], "", "Studies")
  rtSamp <- c(format(RiskSF$sample_size, width=6)[1:(NROW(RiskSF)-1)], "", "Samples")
  rtProb <- c(paste(format(RiskSF$CumProb*100, digit=2, width=7), "%", sep=""), "Probability")

  
  OverallRightText <- paste(format(rtLR, justify="right"), format(rtSC, justify="right"), format(rtSamp, justify="right"), format(rtProb, justify="right"))

  
  par(oma=c(0, 0, 0, 0));
  plot.new()
  
  par(mar=c(1, max(strwidth(OverallLeftText))*txtScale, 2, max(strwidth(OverallRightText))*txtScale ))
  CaptionHeight = 1
  minProb <- floor(min(log10(RiskSF$CumProb)))
  plot.window(xlim=c(minProb, 0), ylim=c(-CaptionHeight*.5, NROW(RiskSF)+CaptionHeight*.5))
  #mticks <- seq(minProb, 0, 1)

  ticks <- rep(seq(minProb, 0, 1), c(rep(2, abs(minProb)), 1))  + c(rep(c(0, log10(5)), abs(minProb)), 0)
  print(ticks)
  
  
  y <- 1:NROW(RiskSF);
  segments(minProb, y, 0, y, lty="dotted", col=segColor)

   # Vertical lines
  for (i in ticks) {
    lines(c(i, i),  c(1-.25, NROW(RiskSF)), lty="dotted", col=segColor)
    #lines(c(i, i),  c(1-.25, NROW(RiskSF)), lty="dotted", col=segColor)
  }

  boxScale <- log10(max(RiskSF$sample_size, na.rm=TRUE))


  Row <- 1;    
  Value <- 0;
  ValueBefore <- 0;
  while(Row <= NROW(RiskSF)) {
   # print(Row)
        Value <- log10(RiskSF$CumProb[Row]);
        if(Row == NROW(RiskSF)) {Symbol<-19; Cex<-1;Colour<-"blue"}
        else {Symbol<-15; Cex<-log10(RiskSF$sample_size[Row])/boxScale;Colour<-colFunction(RiskSF$study_cnt_total[Row])}
        points(Value, Row, cex=Cex*4, col=Colour, pch=Symbol);
        if(Row > 1) {
          ## connect with lines   
            lines(c(ValueBefore, Value), c(Row-1, Row));      
        }
        Row <- Row+1;
        ValueBefore <- Value;
    }
      
  #lines(c(1, 1), c(1, NROW(RiskSF));
 
  
  mtext(OverallLeftText, at=c(y, NROW(RiskSF)+.75), adj=1, side=2, las=2, family="mono", font=1)
  mtext(OverallRightText, at=c(y, NROW(RiskSF)+.75), adj=0, side=4, las=2,  family="mono", font=1)
  mtext(condition, side=3, font=1)
  
  #axis(1, at=mticks, paste(labels=10^(mticks)*100, "%", sep=""), pos=0.5, cex.axis=.8, mex=0.5)
  axis(1, at=ticks, paste(labels=10^(ticks)*100, "%", sep=""), pos=0.5, cex.axis=.6, mex=0.5)

}

riskOGramForCondition <- function(condition, riskA=riskAlleles, preTP=preTestProbs, ...) {
  riskDF <- riskA[riskA$broad_phenotype == condition,]
  Prevalence <- preTP$probability[preTP$broad_phenotype == condition]
  #print(Prevalence)
  #print(riskDF)
  if (length(Prevalence) < 1 | NROW(riskDF)==0 ) {print ("Nothing to plot")}
  else {
    riskOGramSingleCondition(Prevalence, riskDF, condition, ...) 
  }

}


########################################################################
# Look at overall risk


# Get the risk data
riskAlleles <- dbGetQuery(con, "select * from indian_snps")


riskSummaryCt <- tapply(riskAlleles$study_cnt, riskAlleles$broad_phenotype, sum)
riskSummarySNP <- tapply(riskAlleles$study_cnt, riskAlleles$broad_phenotype, length)
riskSummaryLRprod <- tapply(riskAlleles$LR, riskAlleles$broad_phenotype, prod)
riskSummaryMaxLR <- tapply(riskAlleles$LR, riskAlleles$broad_phenotype, getMaxAbsLog)
riskSummary <- mergeViaRowNames(riskSummaryCt, riskSummaryLRprod, "StudyCount", "LR")
riskSummary <- mergeViaRowNames(riskSummary, riskSummarySNP, "", "SNPCount")
riskSummary <- mergeViaRowNames(riskSummary, riskSummaryMaxLR, "", "MaxLR")
colnames(riskSummary) <- c("StudyCount", "LR", "SNPCount", "MaxLR")
riskSummary$Disease <- rownames(riskSummary)

# update the table name for the pretest probability
preTestProbs <- dbGetQuery(con, "select * from proj_patient_risk.pretest_indian_female")
riskSummary <- merge(riskSummary, preTestProbs, by.x="Disease", by.y="broad_phenotype", all.x=F, all.y=F)
riskSummary$PostTest <- preToPost(riskSummary$LR, riskSummary$probability)
riskSummary$PostTestWithMax <- preToPost(riskSummary$MaxLR, riskSummary$probability)
write.table(riskSummary, "RiskSummary.txt", row.names=F, sep="\t", quote=F)

# get the risk over all disease
pdf("indian_OverallRisk.pdf", height=7, width=7)
#multRisk <- riskSummary[riskSummary$SNPCount > 1,]
multRisk <- riskSummary[riskSummary$Disease != "Creutzfeldt-Jakob disease" & riskSummary$Disease != "Crohn's disease",]
multRisk <- multRisk[order(multRisk$PostTest),]
drawRiskOverallGraph(multRisk$Disease, multRisk$probability, multRisk$PostTest, multRisk$SNPCount, txtScale=25, prevSize=.2, showOnlyAggregate=T, legcex=.65, lwd=3, segments=T, segColor="gray60", showNumerics=T)
#drawOverallLinear(multRisk$Disease, multRisk$probability, multRisk$PostTest, multRisk$SNPCount, txtScale=25, prevSize=.2, showOnlyAggregate=T, legcex=.65, lwd=3, segments=T, segColor="gray60", showNumerics=T)
dev.off()

# get the risk over all disease
pdf("indian_OverallRisk_LinearScale.pdf", height=7, width=7)
#multRisk <- riskSummary[riskSummary$SNPCount > 1,]
multRisk <- riskSummary[riskSummary$Disease != "Creutzfeldt-Jakob disease" & riskSummary$Disease != "Crohn's disease",]
multRisk <- multRisk[order(multRisk$PostTest),]
#drawRiskOverallGraph(multRisk$Disease, multRisk$probability, multRisk$PostTest, multRisk$SNPCount, txtScale=25, prevSize=.2, showOnlyAggregate=T, legcex=.65, lwd=3, segments=T, segColor="gray60", showNumerics=T)
drawOverallLinear(multRisk$Disease, multRisk$probability, multRisk$PostTest, multRisk$SNPCount, txtScale=25, prevSize=.2, showOnlyAggregate=T, legcex=.65, lwd=3, segments=T, segColor="gray60", showNumerics=T)
dev.off()

pdf("Indian_uterine_leiomyoma_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Uterine leiomyoma", txtScale=70)
dev.off()

pdf("Indian_Bipolar_riskOgram.pdf", height=7, width=12)
riskOGramForCondition("Bipolar disorder", txtScale=70)
dev.off()

pdf("Indian_Esophageal_riskOgram.pdf", height=7, width=12)
riskOGramForCondition("Esophageal cancer", txtScale=70)
dev.off()

pdf("Indian_aneurysm_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Intracranial aneurysm", txtScale=70)
dev.off()

pdf("Indian_gastric_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Gastric cancer", txtScale=70)
dev.off()

pdf("Indian_obesity_riskOgram.pdf", height=7, width=12)
riskOGramForCondition("Obesity", txtScale=70)
dev.off()

pdf("Indian_fibrillation_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Atrial fibrillation", txtScale=70)
dev.off()

#pdf("Indian_Myopia_riskOgram.pdf", height=3, width=12)
#riskOGramForCondition("Myopia", txtScale=70)
#dev.off()

pdf("Indian_Cardiomyopathy_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Dilated cardiomyopathy", txtScale=70)
dev.off()

pdf("Indian_Bladder_riskOgram.pdf", height=5, width=12)
riskOGramForCondition("Bladder cancer", txtScale=70)
dev.off()

pdf("Indian_Stroke_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Stroke", txtScale=70)
dev.off()

pdf("Indian_BCC_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Basal cell carcinoma", txtScale=70)
dev.off()

#pdf("Indian_Crohn_riskOgram.pdf", height=16, width=12)
#riskOGramForCondition("Crohn's disease", txtScale=70)
#dev.off()

pdf("Indian_Myeloproliferative_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Myeloproliferative disorders", txtScale=70)
dev.off()

pdf("Indian_Narcolepsy_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Narcolepsy", txtScale=70)
dev.off()

pdf("Indian_MSA_riskOgram.pdf", height=3, width=12)
riskOGramForCondition("Multiple system atrophy", txtScale=70)
dev.off()


# update the disease names for RiskOGram
#pdf("Snyder_Glaucoma_riskOgram.pdf", height=3, width=12)
#riskOGramForCondition("Open angle glaucoma", txtScale=70)
#dev.off()

#pdf("Snyder_Dyslipidemia_riskOgram.pdf", height=3, width=12)
#riskOGramForCondition("Dyslipidemia", txtScale=70)
#dev.off()

#pdf("Snyder_CAD_riskOgram.pdf", height=10, width=10)
#riskOGramForCondition("Coronary artery disease", txtScale=60)
#dev.off()

#pdf("Snyder_BCC_riskOgram.pdf", height=5, width=10)
#riskOGramForCondition("Basal cell carcinoma", txtScale=60)
#dev.off()

#pdf("Snyder_Hypertriglyceridemia_riskOgram.pdf", height=3, width=10)
#riskOGramForCondition("Hypertriglyceridemia", txtScale=60)
#dev.off()

#pdf("Snyder_T2D_riskOgram.pdf", height=10, width=12)
#riskOGramForCondition("Type 2 diabetes", txtScale=70)
#dev.off()

#pdf("Snyder_Alzheimer_riskOgram.pdf", height=7, width=12)
#riskOGramForCondition("Alzheimer's disease", txtScale=70)
#dev.off()

pdf("Indian_Parkinson_riskOgram.pdf", height=7, width=12)
riskOGramForCondition("Parkinson's disease", txtScale=70)
dev.off()



#pdf("Snyder_Hypertension_riskOgram.pdf", height=3, width=12)
#riskOGramForCondition("Hypertension", txtScale=70)
#dev.off()

#pdf("Snyder_Obesity_riskOgram.pdf", height=7, width=12)
#riskOGramForCondition("Obesity", txtScale=70)
#dev.off()

#pdf("Snyder_Psoriasis_riskOgram.pdf", height=7, width=10)
#riskOGramForCondition("Psoriasis", txtScale=60)
#dev.off()

#pdf("Snyder_AMD_riskOgram.pdf", height=7, width=10)
#riskOGramForCondition("Age related macular degeneration", txtScale=60)
#dev.off()

#pdf("Snyder_Prostate_riskOgram.pdf", height=12, width=10)
#riskOGramForCondition("Prostate cancer", txtScale=60)
#dev.off()

#pdf("Snyder_Stroke_riskOgram.pdf", height=3, width=10)
#riskOGramForCondition("Stroke", txtScale=60)
#dev.off()

# height=3 when SNP_cnt<=5
# height=5 when 5<SNP_cnt<=10
# height=7 when 10<SNP_cnt<=20
# height=10 when 20<SNP_cnt
