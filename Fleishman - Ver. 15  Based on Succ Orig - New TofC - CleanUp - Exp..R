######################################### "Fleishman - Ver. 15  Based on Succ Orig - New TofC - CleanUp - Exp."


############################  11/19/2021.  Start of all work.

############################  12/13/2021 Built directly from "Fleishman Transform for Skewed Data - Simplified 2"
############################  Using the "Constants" from the R script to calculate the lowest level of Stand. Kurtosis

############################  12/18/2021 Built directly from "Fleishman Transform for Skewed Data - External Constants"
############################  I want to generate an approx. null distribution from the measured 4 moments.  This, of course,
############################  would be quite crude in the case of say n=10 or so.

############################  12/23/2021 Built directly from "Fleishman Transform - External Constants - Toy 1.R"
############################   Trying to get a DYNAMIC creation of the "External Constants" for Null Hypothesis gen.

############################   12/26/2021 Built directly from "Fleishman Transform - External Constants - Real 1"
############################    I know want to start the process of actually DYNAMICALLY determining the Fleishman Transforma constants

############################   Uses the "calc_lower_skurt" function to get the Fleishman transforma. "constants"

############################   12/30/2021  "uilt directly from Fleishman Transform - External Constants - Real 2 - Dynamic Constants 1
############################     Now I want to change to function "nonnormvar1" to get the Fleishman transforma. "constants"

############################    1/5/2022  Using first of Bobee's 1975 Skew Corrections

############################    I now want to try Bobee's Third Skew Correction Methods (the one that depends on both SS & Sample Skew

############################   1/8/2022 Built dirctly from "Fleishman Transform - External Constants - Real 4 - Bobee's Third Skew Method"
############################    I want that script to be left as a "re-entry" type of point in case this version runs into problems
############################   1/9/2022 I'm "tieing off" this script to serve as a re-entry point.f

############################   1/9/2022 Built directly on "Fleishman Transform - External Constants - Real 5 - Bobee's Third - Clean Up"
############################   I want to try a new idea on the min. kurtosis for the Calc2

############################   1/15/2022 I'm tying this script off since it work quite well and I want to move on.

############################   1/15/2022 Built directly from:  "Fleishman Transform - External Constants - Real 6 - Bobee's Third "

############################   1/18/2022 Built directly from: "Fleishman Transform - External Constants - Real 7 - Cleanup Time"
############################   As a "re-entry" point since I need to do a major rennovation on the gen. of NData & NullData

############################   1/22/2022 Built directly from:  "Fleishman Transform - External Constants - Real 8 - Cleanup Time"
############################    I just wanted to leave the previous version behind as a re-entry point.

############################   1/28/2022 Built directly from:  "Fleishman Transform - External Constants - Real 9 - Cleanup Time"
############################   I want "..... Real 9" to be a possible re-entry point if needed

############################   1/30/2022 Built directly from:  "Fleishman Transform - External Constants - Real 10 - Cleanup Time"
############################    I want this "..... Real 10" to be a possible re-entry point if needed

############################   2/1/2022 Built directly from:  "Fleishman Transform - External Constants - Real 11 - Cleanup Time"
############################   I want this ".... Real 11" to be a possible re-entry point if needed

###########################    2/5/2022 Built directly from:  "Fleishman Transform - External Constants - Real 12- Specific Kurtosi"
###########################     This is the second of my two-sample scripts.  The first one has some sore of bug in it.

###########################    2/5/2022 Built directly from:  "Fleishman Transform - Two Sample Version Take 2"
###########################     I want a good re-entry point if necessary

###########################    2/9/2022 Built directly from:  "Fleishman Transform - Two Sample Version Take 3"
###########################    I want to really, really move on to the two-sample case but want to leave this one-sample success
###########################    story behind as a re-entry point even though its name " . . .Two Sample Version Take 3" turned out
###########################    to aspirational.  I had to regroup back to the one-sample case.

###########################    2/10/2022  Built directly from "Fleishman Transform - TRUE Two Sample Version Take 1"
###########################     I have already made the t Test portion into a 2-sample test.  I now want to make my tech. 2-sample

###########################    2/12/2022 Built directly from "Fleishman Transform - True Two Sample Version - Now incl. my tech as well"
###########################     I wanted a safe re-entry point using the prev. ver as I try some ideas re. why I see increases in
###########################     FAR for 2-sample case with skew where Ferenci sees a decrease

###########################   2/13/2022 Built directly from "Fleishman Transform - True Two Sample Version - Now incl. my tech as well - exp1"
###########################    I wanted that script to serve as a "re-entry" point in case I screw this one up.
###########################    In this one, I'm goint to do a bunch of cleanup and checking

###########################   2/18/2022 Built directly from "Fleishman Transform - True Two Sample Version - Now incl. my tech as well - exp2"
###########################   I wanted that prev. ver. to be a "re-entry" point if necess.  This new ver. will have the "Constants"
###########################   in Sec. 2 of the script done using a matrix look-up method.

###########################   2/21/2022 Built directly from "Fleishman Transform - True Two Sample Version - Now incl. my tech as well - exp2"
###########################   I wanted that prev. ver. to be a "re-entry" point if necessary.  This new ver. will have NullData2 depend
###########################   on it's own skew & kurtosis

###########################   2/25/2022 Built directly from "Fleishman Transform - Two Sample Version - Incl. my tech  w Active Constants & NullData2 Sep Skew"
###########################   Want that prev. ver. to be a "re-entry" point if necess.  This new ver. will fully implement my Table of Constants
###########################   for the 2nd portion of my script dealing with NullData, NullData2, etc.

###########################   3/1/2022  Built directly from "Fleishman Transform - Two Sample Version - Incl. my tech  w Active Constants & NullData2 Sep Skew -2"
###########################    I want that prev. ver. to be a "re-entry" point if necess.  This new ver. will incorporate an outer
###########################    loop of "cases" (or some better name) to run multiple scenarios.

###########################   3/3/2022 Built directly from "Fleishman Transform - Two Sample Version - Incl.  Active Constants Outer 'Cases' Loop"
###########################    I want that prev. ver. to be a "re-entry" point.  This new ver. will experiment with the harmonic mean

###########################   3/22/2022 Built directly from "Fleishman Transform - Two Sample Version - Active Constants 
###########################    Outer 'Cases' Loop - Harmonic Means".  I want to co-introduce Wilcox test for median and want a re-entry point.

###########################   4/5/2022 Built directly from "Fleishman Transform - Two Sample Version - Active Constants, Harmonic Means & Medians"
###########################    I want a "re-entry point" as I try to deal with negative skew.

###########################   4/23/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew"
###########################   I want a "re-entry point" as I re-extend the program to the two-sample cases & beyond

###########################   4/28/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew II"
###########################   I want a "re-entry point" as I dig deeper and issues includ. how to deal with negative skews

###########################   5/2/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew III"
###########################   I want to be able to input separate standard deviations for the 2-Sample cases

###########################   5/22/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew IV"
###########################     I want a "re-entry point" to try out some new ideas re. calculating "ScoresBoot" and "Score"

###########################   6/1/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 5"
###########################    As usual, I want a "re-entry pont" to try to properly determine the harmonic mean by fist transforming the data
###########################    to log-normal and then back

###########################   11/4/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 6"
###########################    I want a re-entry point so I can do a number of clean up ops on how I do "LookupRow", etc.

###########################   11/19/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 7"
###########################    I want that earlier ver ( ..... 7) to be a possible re-entry point 

###########################   12/3/2022 Built directly on "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 8"
###########################   I want that earlier ver. to be a possible re-entry point.  In this version, I want to try a different
###########################   method to correct one-sample median testing using the delta betwen the mean and the median.

###########################   12/6/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 9"
###########################    Want to be able to re-enter.  This new version will allow some exps with how I deal with (Med-Mean) correc factor

###########################   12/15/2022 Built directly from "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 10"
###########################    "Fleishman Transform - Two Sample Version - Harmonic Means & Medians - Neg Skew 10" does a pretty good job
###########################    but I want to try some radical ideas re. how I calculate "Score" & "Scores.boot" in my Bootstraps
###########################    This was stimulated by the fact that for my One-Sample Mean test I'm seeing a small upward treand(Skew)

###########################   12/28/22 Built directly from "Fleishman Transform - Exper from '.... - Neg Skew 10'"
###########################   I want a re-entry point so that I can re-explore use of Bobee skew correction in this latest version

###########################   1/6/23 Built directly from "Fleishman - 12.28.22 Ver. 1" which will serve as a "Re-Entry Point"
###########################    I want to clean things up.

###########################   1/12/23  Built directly from "Fleishman - 1.6.23 (Ver. 2)"
###########################    I want a re-entry point.  I want to make some possibly important changes in how I hand both neg.
###########################    skew and kurtosis (always positive of course as a quartic) but less than the value for the normal distrib. of 3.

###########################   1/20/23  Built directly from "Fleishman - 1.12.23 (Ver. 3)"
###########################    I want to exp with rounding instead of ceiling for generating the lookuprow and lookupcol

###########################   1/20/23B.  This is basically a return to "Fleishman - 1.12.23 (Ver. 3)" since rounding worked worse than "ceiling".

###########################   1/31/2023 Built directly from "Fleishman - 1.12.23 (Ver. 5)".  I just want a re-entry point.  I want to
###########################    to do a lot of cleanup.  In particular, all those various attempts at "OneSamp.Adj" and "TwoSamp.Adj"

###########################   2/1/2023 Built directly from "Fleishman - 1.12.23 (Ver. 6)".  Just want a re-entry point.  I successfully
###########################    repaired both OneSamp.Adj & TwoSamp.Adj.  Now, I want to start process of removing all of the old
###########################    "commented out" "scafolding".

###########################   3/20/2023 Built directly from "Fleishman - 1.12.23 (Ver. 7".  I want a re-entry poiont to try
###########################    an idea about killing any finding of stat. sig. if the (corrected) Temp.Skew exceeds a large value

###########################   4/4/2023 Built directly from "Fleishman - 1.12.23 (Ver. 8)".  I want to try a couple of ideas re. looking up
###########################    the experimental kurtosis in the table of coefficients I prev. created.  I want a re-entry point.

###########################   4/5/2023 Built directly from "Fleishman - 1.12.23 (Ver. 8)".  I am in the process of generating a new "Table
###########################     of coefficients running from ExcessKurt = -1.2 up to 6 using the R statement: "seq(-1.2,6,by=.1)".  I want the
###########################     "Fleishman - 1.12.23 (Ver. 8)" to be a possible re-entry point.

###########################   4/10/2023 Built directly from "Fleishman - 1.12.23 (Ver. 9) (Experimental wrt ExcessKurt & New Table of Coeffs)"
###########################    I want to do some experiments with TempSkew, etc. and need a re-entry point.

###########################   4/14/2023 Built directly from "Fleishman - 1.12.23 (Ver. 10) (Experimental re. TempSkew, etc.)"
###########################    I want to do some additional experiments with dealing with high skew and low (e.e. 30) Sample Size
###########################    I (obviously) need a "re-entry" point.

###########################   6/1/2023.  Built directly from Fleishman 11.  I want to
###########################   keep working on 2-Samp Median.

###########################   6/16/2023 Copied directly from Posit Cloud Ver. 12

###########################   9/7/23 Copied directly from "Fleishman - 1.12.23 (Ver. 12) (Experimental re. TempSkew, etc.) Modified on 8.25.23"
###########################  Final try to fix 2 Sample Median by adding/subtracting a TempSkew based correction at the Nulldata Level

###########################   9/16/2023 Built directly from "Fleishman - Ver. 13 (Experimental re. TempSkew, etc.) Modified on 9.7.23". I want
###########################   correct several errors I found from a review in Ocean City.  I want the older version to be a re-entry poinnt in
###########################   case of any disaster

###########################   9/19/2023  Built directly from "Fleishman - Ver. 14 Modified on 9.16.23".  I want this "Experimental" version
###########################    to try out restoring the actual calculation of TempKurt and TempKurt2 to try to fix the problem I'm having
###########################    for 2 Sample Median test @ Skew = 0

###########################   4/30/24.  Back to better notes here.  This is directly from "Fleishman - Ver. 15  Based on Succ Orig - New TofC - Experimental"
###########################    I just want to clean up how I switch between mu "Artifiial Data", the LogNormal and the Exponential.

###########################   5/7/2024 Directly from "Fleishman - Ver. 15  Based on Succ Orig - New TofC - Experimental - Cleanup".

###########################   5/14/2024 Directly from "Fleishman - Ver. 15  Based on Succ Orig - New TofC - Experimental - CleanUp"
###########################   I just want to try a little experiment in which I restrict a Skew=4 run to >= 0



rm(list=ls())  #  This command removes any "environment" variables so as not to confuse what's true at the beginning of a script"



library(meta) #  This loads the DerSimonian MA Method from "R"

library(SimMultiCorrData)

library(psych)

library(DescTools)

library(nonpar)

library(moments)

RunMatrix1 <- rep(0,1)                
dim(RunMatrix1) <- c(1,1)
RunMatrix1 <- read.csv(file="C:/Users/Larry2/Documents/FleishmanRead1.txt", header = FALSE, sep = " ", quote = "\"",dec = ".", fill = FALSE)  

Cases <- RunMatrix1[1,1]
#print("Cases")
#print(Cases)

RunMatrix2 <- rep(0,Cases*12)

RunMatrix2 <- read.csv(file="C:/Users/Larry2/Documents/FleishmanRead2.txt", header = FALSE, sep = "\t", quote = "\"",dec = ".", fill = FALSE) 

#RunMatrix3 <- rep(0,NumRuns*40*2)
#dim(RunMatrix3) <- c(2*NumRuns,40)
#RunMatrix3 <- read.csv(file="C:/users/Larry2/Documents/StudySizes.txt", header = FALSE, sep = "\t", quote = "\"",dec = ".", fill = FALSE)

write("",file="C:/Users/Larry2/Documents/FAR.txt",append=FALSE)
write("",file="C:/Users/Larry2/Documents/FAR.t.txt",append=FALSE)
write("",file="C:/Users/Larry2/Documents/FAR.wilcox.txt",append=FALSE)

#write("",file="Harmonic.Real.Reps.txt",ncolumns=1,append=FALSE)
#write("",file="Harmonic.var.Reps.txt",ncolumns=1,append=FALSE)

#write("",file="C:/users/Larry2/documents/MeanofNData.txt",ncolumns=1,append=FALSE)
#write("",file="C:/users/Larry2/documents/AvgTempSkew.txt",append=FALSE)
#write("",file="C:/users/Larry2/documents/AvgDerivedSkew.txt",append=FALSE)
#write("",file="C:/users/Larry2/documents/AvgTempVar.txt",append=FALSE)

for( Combo in 1:Cases) {
  ####################################################################################################################
  #  Start of Actual Fleishman Transformation Work
  ####################################################################################################################
  
  TwoSample <- RunMatrix2[Combo, 1]  #  TwoSample = 1 makes this a two-sample test.  TwoSample = 0 makes this a one-sample test
  
  Method <- RunMatrix2[Combo, 2] #  Mean = 1, Harmonic Mean = 2, Median = 3
  
  Distrib2Var <- 1  #  For the two-sample test, this is the variation of the second distrib.  I.e. DataNorm2
  
  LowKurts <- rep(0,43)

  LowKurts <- read.csv(file="C:/LMP/Statistics & Probability/Normal Distribution - Poss 2021 Paper/Tables of Constants & of LowKurts/LowKurts.txt", header = FALSE)

  LowKurts <- as.numeric(unlist(LowKurts))

  dim(LowKurts) <- c(43,1)
  
  SS <- RunMatrix2[Combo, 4]
  
  SkewEst <- RunMatrix2[Combo, 5]
  
  ExcessKurt <- RunMatrix2[Combo, 6] 
  
  VarEst1 <- RunMatrix2[Combo, 7]
  VarEst2 <- RunMatrix2[Combo, 8]
  
  
  Means1 <- RunMatrix2[Combo, 9]
  
  Means2 <- RunMatrix2[Combo, 10]
  
  SkewLim.Unreal <- 4.2
  
  Distrib <- RunMatrix2[Combo, 11]
  DistParam <- RunMatrix2[Combo, 12]
  
  #TwoSamp.Adj.1 <- RunMatrix2[Combo, 12]
  
  
  #Harmonic.adj <- 30  #  6/7/2022 Based on a series of experiments
  Harmonic.adj <- 500  #  1/17/2023 Based on a NEW series of experiments
  
  
  MeanAdj <- 0   #  To handle the neg. numbers for harmonic mean calcs
  Means1 <- Means1 + MeanAdj
  Means2 <- Means2 + MeanAdj
  
  
  StartReps <- 1
  
  MaxReps <- RunMatrix2[Combo, 3]
  
  
  MaxBoots <- 1000  

	Adjust.Fudge <- 1
	# Exp. on 9/2/2023
	#Adjust.Fudge <- 0

  
  LookUpRow <- 10*SkewEst + 1
  #LookUpRow <- ceiling(LookUpRow)  ###  attic
  LookUpRow <- round(LookUpRow)
  LookUpRow <- max(1,LookUpRow)
  LookUpRow <- min(43,LookUpRow)
  
  
  #LowLookUp <- 10*SkewEst + 1  ###############  Put necess. "guardrails" around LowLookUp  #### Exp on 11/4/2022
  LowLookUp <- LowKurts[LookUpRow,1]
  
  
  #KurtEst <- LowLookUp + ExcessKurt  ##### BetaKurtEst <- ExcessKurt   ##### Beta
  KurtEst <- ExcessKurt   ##### Beta
  
  
  TableofConstantsII <- readRDS("C:/LMP/Statistics & Probability/Normal Distribution - Poss 2021 Paper/Tables of Constants & of LowKurts/TableofConstants")  
  
  dim(TableofConstantsII) <- c(41,104,4)
  
  BobeeA <-  1 + 6.51/SS + 20.2/SS^2
  BobeeB <- 1.48/SS + 6.77/SS^2
  
  #LowKurts <- as.numeric(unlist(LowKurts))  #  2/1/22 Put this here to remind me how to change a "list" var. to a numeric one
  
  
  # Normal distribution with Fleishman's PMT:
  BaseKurtosis <- 3 #  For the normal distribution
  
  BaseKurtosis <- BaseKurtosis + 15
  
  NumTrys <- 1
  
  #write("",file="C:/users/Larry2/documents/Test-Normal.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/TValues.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/CritVal1.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/CritVal2.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/IsItSS.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/IsTTestSS.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/pvalues.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/NormData.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/NData.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/Normal.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/TinStdUnits.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/Score.txt",ncolumns=1,append=FALSE)
  write("",file="C:/users/Larry2/documents/Score.reps.txt",ncolumns=1,append=FALSE)

	write("",file="C:/users/Larry2/documents/AvgNull.median.txt",ncolumns=1,append=FALSE)

  
  
  
  #write("",file="C:/users/Larry2/documents/VarofNData.reps.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/MeanofNData.reps.txt",ncolumns=1,append=FALSE)
  write("",file="C:/users/Larry2/documents/VarScoresBoot.reps.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/VarNullData.reps.txt",ncolumns=1,append=FALSE)
  
  #write("",file="C:/users/Larry2/documents/TempVar.reps.txt",ncolumns=1,append=FALSE)
  
  write("",file="C:/users/Larry2/documents/TempSkew.reps.txt",ncolumns=1,append=FALSE)
  
  write("",file="ScoresBoot.reps.txt",ncolumns=1,append=FALSE)
  #write("",file="C:/users/Larry2/documents/Score.reps.txt",ncolumns=1,append=FALSE)
  
  
  write("",file="C:/users/Larry2/documents/TempKurt.reps.txt",ncolumns=1,append=FALSE)
  
  #write("",file="C:/users/Larry2/documents/NullKurt.reps.txt",ncolumns=1,append=FALSE)
  
  write("",file="C:/users/Larry2/documents/KurtAbove.reps.txt",ncolumns=1,append=FALSE)
  
  #write("",file="C:/users/Larry2/documents/NullSkew.reps.txt",ncolumns=1,append=FALSE)
  
  #write("",file="C:/users/Larry2/documents/BadCalcValues.reps.txt",ncolumns=3,append=FALSE)
  #write("",file="C:/users/Larry2/documents/GiveUp.reps.txt",ncolumns=3,append=FALSE)
  #write("",file="C:/users/Larry2/documents/OldKurt.reps.txt",ncolumns=1,append=FALSE)
  
  
  
  CExt <- rep(0,4)
  dim(CExt) <- c(4,1)
  
  Delta <- rep(0,6)
  dim(Delta) <- c(6,1)
  
  
  SkewEst <- abs(SkewEst)
  SkewEst <- 100*SkewEst
  SkewEst <- ceiling(SkewEst)  ### Attic
  SkewEst <- SkewEst/100
  
  OrigMeanEst <- 0
  OrigVarEst <- 1
  
  
  NDataVarAdj <- 1 #  12/27/21 temp. exper. to fix Skew=3 FARs.
  
  
  DataNorm2 <- rep(0,SS)
  dim(DataNorm2) <- c(SS,1)
  
  FAR <- 0
  FAR.t <- 0
  FAR.wilcox <- 0
  
  #KurtosisEst2 <- KurtosisEst
  
  
  Score.reps <- rep(0,MaxReps)
  dim(Score.reps) <- c(MaxReps,1)
  
  ScoresBoot <- rep(0,MaxBoots)
  dim(ScoresBoot) <- c(MaxBoots,1)
  
  DataNorm <- rep(0,SS)
  dim(DataNorm) <- c(SS,1)
  
  DataNorm2 <- rep(0,SS)
  dim(DataNorm2) <- c(SS,1)
  
  NData <- rep(0,SS)
  dim(NData) <- c(SS,1)
  
  NData2 <- rep(0,SS)
  dim(NData2) <- c(SS,1)
  
  Xform.Data <- rep(0,SS)
  dim(Xform.Data) <- c(SS,1)
  
  temp1 <- rep(0,SS)
  dim(temp1) <- c(SS,1)
  
  ThreshAdj <- 1  

  IgnoreFlag <- 0
  
  
  SkewNeg.Flag <- 0
  SkewNeg2.Flag <- 0
	SkewNegBoth.Flag <- 0
  
  #OneSamp.Adj <- 0
  
  AvgTempSkew <- 0
  AvgDerivedSkew <- 0
  AvgTempVar <- 0
  
  
  
  RealReps <- (MaxReps-StartReps) + 1
  
  	NoCount <- 0
	#Duds <- 0
  
  for(reps in StartReps:MaxReps) { ## Apple
    print("reps")
    print(reps)

	


    Ignore.Flag <- 0
    
    
    Skew.Flag <- 0
    
    Skew.Unreal <- 0
    
    TempSkew <- 0
    Reverse <- 0
    
    #Con2 <- rep(0,4)
    #dim(Con2) <- c(4,1)
    
    Seed <- runif(1,min=0001,max=9999) 
    Seed <- round(Seed,0)
    Seed <- max(1,Seed)
    
    #print("Seed")
    #print(Seed)
    
    
    LookUpRow <- 10*SkewEst + 1
    #LookUpRow <- ceiling(LookUpRow)  ## Attic
    LookUpRow <- round(LookUpRow)
    LookUpRow <- max(1,LookUpRow)
    LookUpRow <- min(41,LookUpRow)
    
    LookUpCol <- KurtEst + 4  ############### Beta
    
    LookUpCol <- round(LookUpCol)
    LookUpCol<- max(1,LookUpCol)
    LookUpCol <- min(104,LookUpCol)
    
    
    CExt <- TableofConstantsII[LookUpRow,LookUpCol,]
    
    
    dim(CExt) <- c(4,1)
    
    LookUpRow2 <- 10*SkewEst + 1
    #LookUpRow2 <- ceiling(LookUpRow2)  #### Attic
    LookUpRow2 <- round(LookUpRow2)
    LookUpRow2 <- max(1,LookUpRow2)
    LookUpRow2 <- min(41,LookUpRow2)
    
    LookUpCol2 <- KurtEst + 4
    
    LookUpCol2 <- round(LookUpCol2)
    LookUpCol2<- max(1,LookUpCol2)
    LookUpCol2 <- min(104,LookUpCol2)
    
    
    CExt.2 <- TableofConstantsII[LookUpRow2,LookUpCol2,]
    
    dim(CExt.2) <- c(4,1)

############################ Rats

	
   
    	DataNorm[1:SS,1] <- rnorm(SS,0,1)

	DataNorm2[1:SS,1] <- rnorm(SS,0,1)

	if(Distrib == 0) {
	    
    	NData[,1] <- CExt[1] + CExt[2]*DataNorm[,1] + CExt[3]*DataNorm[,1]^2 + CExt[4]*DataNorm[,1]^3

	#for(x in 1:SS) {

		#if(NData[x,1] < 0) NData[x,1] <- -NData[x,1]

	#}

	#NData[which(NData[ < 0,1)]) <- -NData[which(NData[ < 0),1])
	
	NData[,1] <- NData[,1] * sqrt(VarEst1)    
    	
   	NData[,1] <- NData[,1] + Means1
	
	NData2[,1] <- CExt.2[1] + CExt.2[2]*DataNorm2[,1] + CExt.2[3]*DataNorm2[,1]^2 + CExt.2[4]*DataNorm2[,1]^3

	#NData2[which(NData2[ < 0,1]) <- -NData2[which(NData2[ < 0),1])

	#for(x in 1:SS) {

		#if(NData2[x,1] < 0) NData2[x,1] <- -NData2[x,1]

	#}


	NData2[,1] <- NData2[,1] * sqrt(VarEst2)    
    	
   	NData2[,1] <- NData2[,1] + Means2


	}	### End of Artificial Data Section	

	############################### Log Normal Section

	if(Distrib == 1) {


	LogParam <- DistParam

	NData[,1] <- exp(LogParam*DataNorm)

    	#NData[,1] <- log(NData[,1])

	NData[,1] <- NData[,1] * sqrt(VarEst1)    
    	
   	NData[,1] <- NData[,1] + Means1 


	NData2[,1] <- exp(LogParam*DataNorm2)	

	#NData2[,1] <- log(NData2[,1])  
	   
	NData2 <- NData2 * sqrt(VarEst2)
    
    	NData2[,1] <- NData2[,1] + Means2

	}  #### End of LogNormal Section

	####################################   Exponential Section

	if(Distrib == 2) {

	RateParam <- DistParam

	NData[,1] <- rexp(SS,RateParam)

    	#NData[,1] <- log(NData[,1])
    
    	NData[,1] <- NData[,1] * sqrt(VarEst1)    
    	
   	NData[,1] <- NData[,1] + Means1
        
    	NData2[,1] <- rexp(SS,rate=RateParam)	

	#NData2[,1] <- log(NData2[,1])    
    
	NData2[,1] <- NData2[,1] * sqrt(VarEst2)
    
    	NData2[,1] <- NData2[,1] + Means2

	}  ## End of Exponential section
	    
    
    temp1 <- var(NData)
    #write(var(NData),file="C:/users/Larry2/documents/VarofNData.reps.txt",ncolumns=1,append=TRUE)
    #write(mean(NData),file="C:/users/Larry2/documents/MeanofNData.reps.txt",ncolumns=1,append=TRUE)
    
    
    StatSig <- 0
    
    z_prime <- 0
    
    
    ValidFlag <- 0
    
    for (i in 1:NumTrys) {
      
      #TempSeed <- TempSeed + i
      
    }
    #if(ValidFlag == 0) stop("I tried NumTrys for a valid.pdf distribution but failed")
    
    N_exp <- rnorm(SS,0,1)
    
    harmonic.mean(NData)
    sd(NData)
    skewness(NData)
    kurtosis(NData)
    
    
    RealNorm <- rnorm(SS,0,1)
    
    temp1 <- ( harmonic.mean(NData) ) / (sd(NData)/sqrt(SS))
    
    #write(temp1,file="C:/users/Larry2/documents/TValues.txt",ncolumns=1,append=TRUE)
    
    temp1 <- rnorm(SS,0,1)
    #write(temp1,file="C:/users/Larry2/documents/NormData.txt",ncolumns=1,append=TRUE)
    #write(NData,file="C:/users/Larry2/documents/NData.txt",ncolumns=1,append=TRUE)
    #write(temp1,file="C:/users/Larry2/documents/Normal.txt",ncolumns=1,append=TRUE)
    
    #TempVar <- mad(NData)^2
    
    TempVar <- var(NData)   ## Banana
    #TempVar <- VarEst1
    TempVar2 <- var(NData2) ## Banana
    #TempVar2 <- VarEst2
    
    
    if( TwoSample == 0 ) {
      
      y_starLMP1 <- qt(.025,(SS-1))
      y_starLMP2 <- qt(.975,(SS-1))
      
    } else if ( TwoSample == 1) {
      
      y_starLMP1 <- qt(.025,(2*SS-2))
      y_starLMP2 <- qt(.975,(2*SS-2))
      
    }
    
    StatSig.t <- 0
    
    if(TwoSample == 0) {
      
      Result.t <- t.test(NData)
      if(Result.t$p.value <= .05) FAR.t <- FAR.t + 1
      
      Test.Wilcox <- wilcox.test(NData,mu=0,exact=TRUE)
      if(Test.Wilcox$p.value <= .05) FAR.wilcox <- FAR.wilcox + 1
      
      
    } else if(TwoSample == 1 ) {
      
      PooledSD <- sqrt(  (SS-1)*TempVar + (SS-1)*TempVar2  )
      PooledSD <- PooledSD / sqrt( 2*SS -2)
      
      #PooledSD <- sqrt(  SS*sd(NData)^2 + SS*sd(NData2)^2  )
      #PooledSD <- PooledSD / sqrt(2*SS)
      
      Result.t <- t.test(NData,NData2)
      if(Result.t$p.value <= .05) FAR.t <- FAR.t + 1
      
      
      Test.Wilcox <- wilcox.test(NData,NData2)
      if(Test.Wilcox$p.value <= .05) FAR.wilcox <- FAR.wilcox + 1
      
      
    }
    
    #write(StatSig.t,file="C:/users/Larry2/documents/IsTTestSS.txt",ncolumns=1,append=TRUE)
    
    
    
    
    #OldSkew <- TempSkew
    
    Flag.Skew <- 0
    
    for(ignore in 1:1) {  #  temp. for loop to not execute a bunch of steps
      
      if(IgnoreFlag == 1) break
      
      #OrigSkew <- TempSkew #### 12/6/22.  "Rescues" signed version of TempSkew for later use and reporting
      
      #if(TempSkew < 0) {
        #TempSkew <- -TempSkew  # Works reasonably OK.  Commenting it out is an exp. on 1/12/23
        #TempSkew <- 0
        #SkewNeg.Flag <- 1
      #}
      
      #TempSkew <-  (  TempSkew * (BobeeA + BobeeB*TempSkew^2) )   ## Banana

	#if(TempSkew < 0) {
        #TempSkew2 <- -TempSkew2  # Commenting it out is an exper. on 1/12/2023
        #TempSkew2 <- 0
        #SkewNeg.Flag <- 1
      #}

      
      #if(TempSkew >= SkewLim.Unreal) Skew.Unreal <- 1 ## 4/10/23 Exp
      
      #if(SS == 20) TempSkew <- 1.1*TempSkew #### Simple exp. on 1/4/23
      
      lambda <- log(2)/3
      
      Var.Corr <- exp(lambda*abs(TempSkew) )
      
      
      q3 <- quantile(NData,.75)
      q1 <- quantile(NData,.25)
      a <- min(NData)
      b <- max(NData)
      
      SD.Star <- .5*( (b-a)/3.735 + (q3-q1)/1.26 )
      #SD.Star <- (b-a)/3.735 
      
      denom1 <- 2 * qnorm( (SS-.375) / (SS+.1) )
      
      #SD.Star <- (b-a) / denom1
      
      #SD.Star <- Var.Star*1.5
      
      #SD.Star <- sqrt(VarEst) + rnorm(1,0,.25)
      #SD.Star <- sqrt(VarEst) 
      
      
      #print( SD.Star )
      #TempVar <- SD.Star^2
      
     
      
      TempSkew <- skewness(NData)   
	TempSkew <- TempSkew * (BobeeA + BobeeB*TempSkew^2)  

      TempSkew2 <- skewness(NData2) 
	TempSkew2 <- TempSkew2 * (BobeeA + BobeeB*TempSkew2^2)	

	if(TempSkew < 0) {
        #TempSkew2 <- -TempSkew2  # Commenting it out is an exper. on 1/12/2023
        #TempSkew2 <- 0
        SkewNeg.Flag <- 1
      }


 	if(TempSkew2 < 0) {
        #TempSkew2 <- -TempSkew2  # Commenting it out is an exper. on 1/12/2023
        #TempSkew2 <- 0
        SkewNeg2.Flag <- 1
      }
    
    
    
    AvgTempSkew <- AvgTempSkew + TempSkew
    
    
    	#Adjust1 <-  -.3*(1-exp(-1* TempSkew ))
		#y = 0.0285x2 - 0.2332x + 0.0049

	
		
	if(TempSkew >= 0) Adjust1 <- .0285*TempSkew^2 - .2332*TempSkew + .0049
	#if(TempSkew < 0) Adjust1 <- .0285*TempSkew^2 - .2332*TempSkew + .0049
	
	if(TempSkew < 0) Adjust1 <- -1 * ( .0285*TempSkew^2 - .2332*abs(TempSkew) + .0049 )


 	if(TempSkew2 >= 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049
	#if(TempSkew2 < 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049

	if(TempSkew2 < 0) Adjust2 <- -1 * (.0285*TempSkew2^2 - .2332*abs(TempSkew2) + .0049)
      
      if( (TwoSample==1) & (TempSkew2 >= SkewLim.Unreal) ) Skew.Unreal <- 1  ## 4/10/23 Exp
      
      
      ############################### Experiment GOOSE 
      
      
      if(Flag.Skew == 1) {
        #TempSkew <- 0
        #print(TempSkew)
        #stop("Temp Skew Was Neg")
      }
      
            
      AvgDerivedSkew <- AvgDerivedSkew + TempSkew
      AvgTempVar <- AvgTempVar + TempVar      
      
      
      TempKurt <- kurtosis(NData)
      
      write(TempKurt,file="C:/users/Larry2/documents/TempKurt.reps.txt",ncolumns=1,append=TRUE)
      
      
      #print("TempKurt")
      #print(TempKurt)
      #print("TempSkew")
      #print(TempSkew)
      
      #  Exp. on 4/10/22
      
      if(Method == 1) TempMeans <- mean(NData)
      if(Method == 2) TempMeans <- harmonic.mean(NData)
      if(Method == 3) TempMeans <- median(NData)
      
      if(Method == 1) TempMeans2 <- mean(NData2)
      if(Method == 2) TempMeans2 <- harmonic.mean(NData2)
      if(Method == 3) TempMeans2 <- median(NData2)
      
      
      #write(TempVar,file="C:/users/Larry2/documents/TempVar.reps.txt",ncolumns=1,append=TRUE)
      
      write(TempSkew,file="C:/users/Larry2/documents/TempSkew.reps.txt",ncolumns=1,append=TRUE)
      
      
      OldKurt <- kurtosis(NData)
      #write(OldKurt,file="C:/users/Larry2/documents/OldKurt.reps.txt",ncolumns=1,append=TRUE)
      
      RandSamp <- rep(0,SS)
      dim(RandSamp) <- c(SS,1)
      
      RandSamp2 <- rep(0,SS)
      dim(RandSamp2) <- c(SS,1)
      
      NullData <- rep(0,SS)
      dim(NullData) <- c(SS,1)
      
      NullData2 <- rep(0,SS)
      dim(NullData2) <- c(SS,1)
      
      
      
      MinDelta <- min(Delta[1:6])
      
      GoodC <- 0
      ExistFlag <- 0
      #remove(Calc2)
      
      TryFlag <- 1
      
      Kurt.Flag <- 0
      
      #print("Skew.Flag")  #### Apple
      #print(Skew.Flag)
      #print("Temp.Skew")
      #print(TempSkew)
      
      for(attempts in 1:1) {
        
        Seed2 <- runif(1,min=0001,max=9999)
        Seed2 <- round(Seed2)
        Seed2 <- max(1,Seed2)
        Seed <- 1234
        
        LookUpRow <- 10* abs(TempSkew)+ 1  ###########  12/20/22 I need to make this statement the same here as it is for LookUpRow2  ### Beta
        #  Exp below	
	#LookUpRow <- ceiling(LookUpRow)
	LookUpRow <- round(LookUpRow)

      LookUpRow <- max(1,LookUpRow)
      LookUpRow <- min(41,LookUpRow)        
       
        
#########################################  Tiger
        

		KurtAbove1 <- kurtosis(NData) - 3 - LowKurts[LookUpRow,1]  #  1/3/23 Exp  ###########  Beta

		#print("KurtAbove1")
		#print(KurtAbove1)

        #KurtAbove1 <- kurtosis(NData) - 3 # 1/3/23 Exp.
        #KurtAbove1 <- 1 * (  kurtosis(NData) - 3 )
        
        
        #KurtAbove1 <- 0  ### Part of an exp. on 4/4/2023 Saratoga
        
        
        
        
        
        LookUpCol <- KurtAbove1 + 4  ##  Temp Exp on 9/19/23
        #LookUpCol <- 13

		#print("KurtAbove1")
		#print(KurtAbove1)

	#####################################################  Puma

	write(KurtAbove1,file="C:/users/Larry2/documents/KurtAbove.reps.txt",ncolumns=1,append=TRUE)
		


        
        LookUpCol <- round(LookUpCol)
        LookUpCol<- max(1,LookUpCol)
        LookUpCol <- min(104,LookUpCol)

		print("Final LookUpCol")
		print(LookUpCol)
        
        #print("Final LookUpCol")
        #print(LookUpCol)
        #print("kurtosis(NData)")
        #print(kurtosis(NData))
        
        
        Con2 <- TableofConstantsII[LookUpRow,LookUpCol,]
        
        Con2 <- unlist(Con2)
        dim(Con2) <- c(4,1)
        
        
      LookUpRow2 <- 10*abs(TempSkew2) + 1
     	#  Exp below
 	#LookUpRow2 <- ceiling(LookUpRow2)
	LookUpRow2 <- round(LookUpRow2)  

      LookUpRow2 <- max(1,LookUpRow2)
      LookUpRow2 <- min(41,LookUpRow2)
      LookUpRow2 <- round(LookUpRow2)

        #if(TempSkew2 > SkewLim.Unreal) {                       ### Apple
        #Skew.Flag <- 1
        #RealReps <- RealReps - 1
        #Duds <- Duds + 1
        #break
        #}
        
        
        KurtAbove2 <- kurtosis(NData2) - 3 - LowKurts[LookUpRow2,1]  #  12/20/22 Exp
	#KurtAbove2 <- kurtosis(NData2) - 3
        
        
        LookUpCol2 <- KurtAbove2 + 4  
        #LookUpCol2 <- 13
        
        LookUpCol2 <- round(LookUpCol2)  
        
        LookUpCol2<- max(1,LookUpCol2)
        LookUpCol2 <- min(104,LookUpCol2)
        
        
        Con2.2 <- TableofConstantsII[LookUpRow2,LookUpCol2,]
        
        dim(Con2.2) <- c(4,1)	
        
        
      }  #  End of for(attempts in 1:1)
      
      
      if(Skew.Flag == 1) {
        print("reps")	  
        print(reps)
        
        Skew.Flag <- 0
        next	
        #print("I'm in the if statement 'if(Skew.Flag ==1' ")
      }
      
      StatSig <- 0
      
      
      #PooledSD <- sqrt(  (SS-1)^sd(NData)^2 + (SS-1)*sd(NData2)^2  )
      #PooledSD <- PooledSD / sqrt( 2*SS -2)
      PooledSD <- mean(sd(NData),sd(NData2))
      
      VarFudge <- 1
      
      #Harmonic.adj <- 100
      
      NData.mod <- NData + Harmonic.adj
      
      AvgNull.sd <- 0
      AvgNull2.sd <- 0
      AvgNull.mean <- 0
      AvgNull.median <- 0
      AvgNull2.median <- 0
      
      MeanCorrec <- .5
      
      
      for(boots in 1:MaxBoots) {  #############################  "Boots" Loop  #####################################
        
        RandSamp[1:SS,1] <- rnorm(SS,0,1 )
        RandSamp2[1:SS,1] <- rnorm(SS,0,1)  

		if( (SkewNeg.Flag == 1) & (SkewNeg2.Flag == 1) ) SkewNegBoth.Flag <- 1
		else SkewNegBoth.Flag <- 0
        
        #if(SkewNegBoth.Flag == 0)   NullData[,1] <-  Con2[1,1] + Con2[2,1]*RandSamp[,1] + Con2[3,1]*RandSamp[,1]^2 + Con2[4,1]*RandSamp[,1]^3 
        if(SkewNeg.Flag == 0) NullData[,1] <-  Con2[1,1] + Con2[2,1]*RandSamp[,1] + Con2[3,1]*RandSamp[,1]^2 + Con2[4,1]*RandSamp[,1]^3  
		else if(SkewNeg.Flag == 1) NullData[,1] <-  Con2[1,1] + -1*Con2[2,1]*RandSamp[,1] + Con2[3,1]*RandSamp[,1]^2 + -1*Con2[4,1]*RandSamp[,1]^3                                                                                                               
        
        #NullData[,1] <- NullData[,1] * sqrt(VarEst)
        
        NullData[,1] <- NullData[,1] * sqrt( TempVar )
	## 9/7/2023 Exp below
	  #NullData[,1] <- NullData[,1] - Adjust1
        
        
        #NullData[,1] <- NullData[,1] * (2*Half.IQR)
        
               
        if(Flag.Skew == 1) {
          #print(TempSkew)
          #NullData[,1] <- NullData[,1] * 10	
        }
        
        #if(SkewNegBoth.Flag == 0)      NullData2[,1] <-  Con2.2[1,1] + Con2.2[2,1]*RandSamp2[,1] + Con2.2[3,1]*RandSamp2[,1]^2 + Con2.2[4,1]*RandSamp2[,1]^3 
       if(SkewNeg2.Flag == 0) NullData2[,1] <-  Con2.2[1,1] + Con2.2[2,1]*RandSamp2[,1] + Con2.2[3,1]*RandSamp2[,1]^2 + Con2.2[4,1]*RandSamp2[,1]^3 	
        else if(SkewNeg2.Flag == 1) NullData2[,1] <-  Con2.2[1,1] + -1*Con2.2[2,1]*RandSamp2[,1] + Con2.2[3,1]*RandSamp2[,1]^2 + -1*Con2.2[4,1]*RandSamp2[,1]^3 


        NullData2[,1] <- NullData2[,1] * sqrt( TempVar2 )
		## 9/7/2023 Exp below
 		#NullData2[,1] <- NullData2[,1] - Adjust2	     
       
        
        
        
        
        AvgNull.sd <- AvgNull.sd + sd(NullData[,1])
        AvgNull2.sd <- AvgNull2.sd + sd(NullData2[,1])
        
        AvgNull.median <- AvgNull.median + median(NullData[,1])
        AvgNull2.median <- AvgNull2.median + median(NullData2[,1])	
        
        
        AvgNull.mean <- AvgNull.mean + mean(NullData)
        
        OneSamp.Adj <- -.3*(1-exp(-1*TempSkew))
        
        
        #TwoSamp.Adj.1 <- -.47*(  1-exp(-.75 * TempSkew) ) 
        
        
        #TwoSamp.Adj.1 <- -.47*(  1-exp(-.75 * mean(TempSkew,TempSkew2) ) )
        
        #TwoSamp.Adj.2 <- +.47*(  1-exp(-.75 * TempSkew2) )
        
        
        #TwoSamp.Adj.1 <- -2*(  1-exp(-.3 * TempSkew) ) 
        
        TwoSamp.Adj.1 <- 0
        
        TwoSamp.Adj.2 <- 0
        
        
        
        if(TwoSample == 0) {
          
          NullData.mod <- NullData + Harmonic.adj
          
          if( Method == 1) ScoresBoot[boots] <- mean(NullData)/ sd(NullData) 
          
          if( Method == 2) {  #  Rho
            Harmonic.s <- var(1/NullData)
            Harmonic.m = 1/SS * sum(1/NullData)
            Harmonic.var <- (1/SS) * Harmonic.s / Harmonic.m^4
            
            
            ScoresBoot[boots] <- ( harmonic.mean(NullData.mod) - Harmonic.adj ) / sd(NullData.mod)
            #ScoresBoot[boots] <- ( harmonic.mean(NullData.mod) - Harmonic.adj ) / sqrt(Harmonic.var)
            
          }
          
          NullDataCompare <- NullData - mean(NullData[,1])	
          
          
          IQR1 <- IQR(NullData.mod)
          
          if( Method == 3) ScoresBoot[boots] <- (  median(NullData) + OneSamp.Adj ) / sd(NullData)  # At least temp. using "TwoSame.Adj" for this one-sample case
          
          
        } else if(TwoSample == 1) {  #####################  Two Sample is TRUE
          
          NullData.mod <- NullData + Harmonic.adj
          NullData2.mod <- NullData2 + Harmonic.adj
          
          if( Method == 1) ScoresBoot[boots] <- (mean(NullData) - mean(NullData2)) / (mean(sd(NullData),sd(NullData2) ) )	
          
          if( Method == 2) {
            
            Xform.Data[,1] <- exp( NullData[,1])
            Harmonic.Mean.Xform <- harmonic.mean(Xform.Data)
            Harmonic.Mean.Real <- log(Harmonic.Mean.Xform)
            Xform.Data[,1] <- exp( NullData2[,1])
            Harmonic.Mean.Xform2 <- harmonic.mean(Xform.Data)
            Harmonic.Mean.Real2 <- log(Harmonic.Mean.Xform2)
            
            
            Harmonic.Mean.Real <-  harmonic.mean(NullData + Harmonic.adj)
            Harmonic.Mean.Real <- Harmonic.Mean.Real - Harmonic.adj
            Harmonic.Mean.Real2 <- harmonic.mean(NullData2 + Harmonic.adj)
            Harmonic.Mean.Real2 <- Harmonic.Mean.Real2 - Harmonic.adj
            
            
            ScoresBoot[boots] <- (Harmonic.Mean.Real - Harmonic.Mean.Real2) / mean(sd(NullData),sd(NullData2) )
            
          }
          
         
          
          if( Method == 3) {	
            
           

		Alpha <- .1 
            
           		
		if(TempSkew >= 0) Adjust1 <- .0285*TempSkew^2 - .2332*TempSkew + .0049
		#if(TempSkew < 0) Adjust1 <- .0285*TempSkew^2 - .2332*TempSkew + .0049
	
		if(TempSkew < 0) Adjust1 <- -1 * ( .0285*TempSkew^2 - .2332*abs(TempSkew) + .0049 )


 		if(TempSkew2 >= 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049
		#if(TempSkew2 < 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049

		if(TempSkew2 < 0) Adjust2 <- -1 * (.0285*TempSkew2^2 - .2332*abs(TempSkew2) + .0049)
		
            ################### 4/22/24 Exp. Below

		#Adjust1 <- 0 * Adjust1
		#Adjust2 <- 0 * Adjust2
		

		Adjust1 <- Adjust1 * Adjust.Fudge
		Adjust2 <- Adjust2 * Adjust.Fudge

		#  Exp on 9/19/23
		#if(  (TempSkew >= 0) & (TempSkew2 < 0) )  next
		#if(  (TempSkew < 0) & (TempSkew2 >= 0) )  next 

################################ Zebra


            #ScoresBoot[boots] <- ( (median(NullData) + Adjust1 ) / sd(NullData) )
            #- ( (median(NullData2) + Adjust2) / sd(NullData2)   ) 

	############# Exp. on 4/13/24

	a <- median(NullData)
	b <- Adjust1
	c <- sd(NullData)
	d <- median(NullData2)
	e <- Adjust2
	f <- sd(NullData2)


 	 #ScoresBoot[boots] <- (a+b)/c  - (d+e)/f
		
	
	ScoresBoot[boots] <- (a)/c  - (d)/f

#################################### Exp on 4/20/24
	#ScoresBoot[boots] <- median( NullData-NullData2 )



          }	
        }
        
        write(ScoresBoot[boots],file="C:/Users/Larry2/Documents/ScoresBoot.reps.txt",ncolumns=1,append=TRUE)
        
        
        VarNullData.reps <- var(NullData)
        #write(VarNullData.reps,file="C:/users/Larry2/documents/VarNullData.reps.txt",ncolumns=1,append=TRUE)
        
        NullData.Kurt <- kurtosis(NullData)
        #write(NullData.Kurt,file="C:/users/Larry2/documents/NullKurt.reps.txt",ncolumns=1,append=TRUE)
        
        NullData.Skew <- skewness(NullData)
        
        #write(NullData.Skew,file="C:/users/Larry2/documents/NullSkew.reps.txt",ncolumns=1,append=TRUE)
        
        
      }  ###################################################  End of "for(boots in 1:MaxBoots" above  #################################################
      
      AvgNull.sd <- AvgNull.sd / MaxBoots
      AvgNull2.sd <- AvgNull2.sd/MaxBoots
      
      AvgNull.median <- AvgNull.median / MaxBoots
      AvgNull2.median <- AvgNull2.median / MaxBoots
      
      
      AvgNull.mean <- AvgNull.mean / MaxBoots
      
      if(TwoSample == 0) {
        
        NData.mod <- NData + Harmonic.adj
        
        ### 4/9/22 Exp below  #######  PEACH	
        
        
       if( Method == 1) Score <- mean(NData)/sd(NData) 	    
        
       if( Method == 2) {  #  Rho
          
          Harmonic.Mean.Real <-  harmonic.mean(NData + Harmonic.adj)
          Harmonic.Mean.Real <- Harmonic.Mean.Real - Harmonic.adj
          Harmonic.Mean.Real2 <- harmonic.mean(NData2 + Harmonic.adj)
          Harmonic.Mean.Real2 <- Harmonic.Mean.Real2 - Harmonic.adj
          
          
          Harmonic.s <- var(1/NData)
          Harmonic.m = 1/SS * sum(1/NData)
          Harmonic.var <- (1/SS) * Harmonic.s / Harmonic.m^4
          
          #write(Harmonic.Mean.Real,file="C:/users/Larry2/documents/Harmonic.Real.Reps.txt",ncolumns=1,append=TRUE)
          #write(Harmonic.var,file="C:/users/Larry2/documents/Harmonic.var.Reps.txt",ncolumns=1,append=TRUE)
          
          Score <- ( harmonic.mean(NData.mod) - Harmonic.adj ) / sd(NData)		
          
        }	
        
        if( Method == 3) Score <-( median(NData) + OneSamp.Adj  ) / (sd(NData) )  
        
        
      } else if(TwoSample == 1) {  ############################  Two Sample is TRUE
        
        NData.mod <- NData + Harmonic.adj
        NData2.mod <- NData2 + Harmonic.adj	
        #if( Method == 1) Score <- (mean(NData) - mean(NData2) ) ### 12/27/22 This method works
        if( Method == 1) Score <- (mean(NData) - mean(NData2) ) / mean( sd(NData),sd(NData2)   )
        
        if( Method == 2) {
          
          Xform.Data[,1] <- exp( NData[,1])
          Harmonic.Mean.Xform <- harmonic.mean(Xform.Data)
          Harmonic.Mean.Real <- log(Harmonic.Mean.Xform)
          Xform.Data[,1] <- exp( NData2[,1])
          Harmonic.Mean.Xform2 <- harmonic.mean(Xform.Data)
          Harmonic.Mean.Real2 <- log(Harmonic.Mean.Xform2)
          
          Harmonic.Mean.Real <-  harmonic.mean(NData + Harmonic.adj)
          Harmonic.Mean.Real <- Harmonic.Mean.Real - Harmonic.adj
          Harmonic.Mean.Real2 <- harmonic.mean(NData2 + Harmonic.adj)
          Harmonic.Mean.Real2 <- Harmonic.Mean.Real2 - Harmonic.adj
          
          #write(Harmonic.Mean.Real,file="C:/users/Larry2/documents/Harmonic.Real.Reps.txt",ncolumns=1,append=TRUE)
          #write(mean(NData),file="C:/users/Larry2/documents/MeanofNData.txt",ncolumns=1,append=TRUE)
          
          
          Score <- (Harmonic.Mean.Real - Harmonic.Mean.Real2) / mean( sd(NData),sd(NData2) )
          
          
        }
        
        
        if( Method == 3) {
          
                      
          
         	if(TempSkew >= 0) Adjust1 <- .0285*TempSkew^2 - .2332*TempSkew + .0049
		if(TempSkew < 0) Adjust1 <- -1 * (.0285*TempSkew^2 - .2332*abs(TempSkew) + .0049)
		#if(TempSkew < 0) Adjust1 <-  (.0285*TempSkew^2 - .2332*TempSkew + .0049)



      	if(TempSkew2 >= 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049
		#if(TempSkew2 < 0) Adjust2 <- .0285*TempSkew2^2 - .2332*TempSkew2 + .0049


		if(TempSkew2 < 0) Adjust2 <- -1 * (.0285*TempSkew2^2 - .2332*abs(TempSkew2) + .0049) 

############################### 4/22/24 Exp. Below

		Adjust1 <- 1 * Adjust1
		Adjust2 <- 1 * Adjust2

#############################################################  Exp on 4/23/24
		Adjust1 <- AvgNull.median
		Adjust2 <- AvgNull2.median
	
	
		
		Adjust1 <- Adjust1 * Adjust.Fudge
		Adjust2 <- Adjust2 * Adjust.Fudge

#####################################   Hyena
	
		###################### Exp on 4/13/24          
      a <- median(NData)
	b <- Adjust1
	c <- sd(NData)
	d <- median(NData2)
	e <- Adjust2
	f <- sd(NData2)

	#Score <- (a)/c  - (d)/f


	Score <- (a+b)/c  - (d+e)/f  
		
	
		
        }
        
      }
      
      write(Score,file="C:/users/Larry2/documents/Score.reps.txt",ncolumns=1,append=TRUE)
      
      
      Score.reps[reps] <- Score
      temp1 <- Score.reps[reps]
      
      #write(temp1,file="C:/users/Larry2/documents/Score.reps.txt",ncolumns=1,append=TRUE)
      
      if(Skew.Unreal ==1) {
        #write("Skew.Unreal")
        #write(Skew.Unreal)
        CritVal1 <- quantile(ScoresBoot,.025)
        CritVal2 <- quantile(ScoresBoot,.975)
        
      } else {
        
        CritVal1 <- quantile(ScoresBoot,.025)
        CritVal2 <- quantile(ScoresBoot,.975)
      }
      
      #write(Score,file="C:/users/Larry2/documents/Score.txt",ncolumns=1,append=TRUE)
      
      if( (Score <=  ThreshAdj * CritVal1) | (Score >=  ThreshAdj * CritVal2) ){
        
        #if( (Score >=  ThreshAdj * CritVal2) ){
        
        FAR <- FAR + 1

	if(Ignore.Flag ==1) FAR <- FAR - 1        
        
      }
      
      #write(StatSig,file="C:/users/Larry2/documents/IsItSS.txt",ncolumns=1,append=TRUE)
      
    }  # End of special "ignore" loop
    
    #if(Flag.Skew == 1) next	
    
    
    
  }  #  End "reps" overall loop
  
  
  
  #Duds <- MaxReps - RealReps
  
  #print("Duds")
  #print(Duds)
  
  Divisor <- RealReps - NoCount
  
  AvgTempSkew <- AvgTempSkew / Divisor
  AvgDerivedSkew <- AvgDerivedSkew / Divisor
  AvgTempVar <- AvgTempVar / Divisor
  
  
  
  #write(AvgTempSkew,file="C:/users/Larry2/documents/AvgTempSkew.txt",ncolumns=2,append=TRUE)
  #write(AvgDerivedSkew,file="C:/users/Larry2/documents/AvgDerivedSkew.txt",ncolumns=2,append=TRUE)
  #write(AvgTempVar,file="C:/users/Larry2/documents/AvgTempVar.txt",ncolumns=2,append=TRUE)

	write(AvgNull.median,file="C:/users/Larry2/documents/AvgNull.median.txt",ncolumns=2,append=TRUE)

  
  temp1 <- FAR/Divisor
  #write("Hi, I'm trying to print out FAR",file="C:/users/Larry2/documents/FAR.txt",ncolumns=2,append=TRUE)
  
  write(temp1,file="C:/Users/Larry2/Documents/FAR.txt",append=TRUE)
  
  #print("SS")
  #print(SS)
  #print("FAR")
  #print(FAR/Divisor)
  
  temp1 <- FAR.t/reps #  FAR.t will always have data from all reps
  
  write(temp1,file="C:/Users/Larry2/Documents/FAR.t.txt",append=TRUE)
  
  #print("FAR.t")
  #print(temp1)
  
  temp1 <- FAR.wilcox/reps  #  FAR.wilcox will always have data from all reps
  
  write(temp1,file="C:/Users/Larry2/Documents/FAR.wilcox.txt",append=TRUE)
  
  # print("FAR.wilcox")
  #print(temp1)
  
  #print("Divisor")
  # print(Divisor)
  
}  ############################  End of overall "Combo" Loop

#print(Score)
#print(CritVal1)
#print(CritVal2)

#print("Skewness of NData")
#print(skewness(NData))

#print("Skewness of NullData")
#print(skewness(NullData))

#print(var(NData))
#print(var(NullData))

#print("AvgNull.mean")
#print(AvgNull.mean)

#print("mean(NData)")
#print(mean(NData))

FAR/reps
FAR.t/reps
FAR.wilcox/reps

#skew(NData)
#kurtosis(NData)
#var(NData)

#skew(NData2)
#kurtosis(NData2)
#var(NData2)

#hist(NData[,1])

#hist(ScoresBoot)
#Score
skewness(NData)
kurtosis(NData)
skewness(NData2)
kurtosis(NData2)

hist(NData,breaks=seq(-30,30,by=.1))
