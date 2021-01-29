library(data.table)
library(caret)
library(rlist)
library(doParallel)
library(ggplot2)
library(gtools)
library(fastDummies)
library(biglm)
library(car)
library(latex2exp)
library(patchwork)
library(stats)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(cobalt)
library(scales)
library(patchwork)
library(stringr)
library(zoo)
library(elasticnet)
library(ranger)
library(nnet)
library(xgboost)

# Helper Functions --------------------------------------------------------

is.even <- function(x) {
  if (x%%2 == 0) {
    return(T)
  } else {
    return(F)
  }
}

# - upper median for lower bound
# - lower median for upper bound
# implemented as 50th element for lower bound and 51st element for upper bound
upperMedian <- function(x) {
  ordered <- x[order(x)]
  if (is.even(length(x)))  {
    med <- ordered[floor(length(x)/2)+1]
  } else {
    med <- median(ordered)
  }
  return(med)
}

lowerMedian <- function(x) {
  ordered <- x[order(x)]
  if (is.even(length(x)))  {
    med <- ordered[length(x)/2]
  } else {
    med <- median(ordered)
  }
  return(med)
}

theme_new <- function(base_size = 12,
                      base_family = "LM Roman 10",
                      base_line_size = base_size / 85,
                      base_rect_size = base_size / 85){
  theme_classic(base_size = base_size, 
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
  theme(
    plot.title = element_text(
      color = rgb(0, 0, 0, maxColorValue = 255), 
      hjust = 0.5,
      size = 10),
    axis.title = element_text(
      color = rgb(0, 0, 0, maxColorValue = 255),
      size = 10),
    axis.text = element_text(
      color = rgb(0,0,0, maxColorValue = 255),
      size = 8),
    
    complete = TRUE
  )
}

plot_gates <- function(gatesObj, blpObj, title=NULL) {
  
  blp <- blpObj[[2]]
  gates <- gatesObj[[1]]
  numGroups <- length(gates)
  
  tgates <- data.table(t(gates))
  tgates[, group := seq(1,numGroups)]
  
  plt <- ggplot(tgates, aes(factor(group), y=ATE)) +
    scale_y_continuous(labels=comma) +
    geom_hline(aes(yintercept = 0), color = "darkslategray4" , size=0.25) +
    geom_hline(data=blp, aes(yintercept=ateLowerBound, linetype = "type1"), size=0.25) +
    geom_hline(data=blp, aes(yintercept=ateUpperBound, linetype = "type1"), size=0.25) +
    geom_hline(data=blp, aes(yintercept=ate, linetype = "type2"), size=0.25) + 
    scale_linetype_manual(name = "",
                          values = c("type1" = "dotted", "type2" = "dashed"),
                          labels = c("90% CI", "ATE"),
                          guide = guide_legend(reverse = TRUE)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), size=0.5, width=.5, col = I("#2A2A2A")) +
    geom_point(size=2, col = I("#2A2A2A")) +
    xlab(TeX("Groups based on CATE proxy \\textit{$S(Z)$}")) +
    ylab("GATES") +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
    theme(legend.position = c(0.8,.5)) + 
    ggtitle(title) +
    theme_new()
  
  return(plt)
  
}

# ACTUAL ESTIMATION -------------------------------------------------------

create_cddf_object <- function(datatable, treatVar, outcomeVar, groupVar=NULL, 
                               method, colsNotToUse=NULL, numSplits=100, 
                               preprocess = c("center", "scale"),
                               sigLevel=0.05, seed=123, numGroups=10, quantiles=NULL,
                               CLANvars=NULL, returnModel=F, propensityScores,
                               neuralNetMaxit=100, verbose=0) {
  
  ###
  # takes data.table as input
  # treatVar needs to be named "treatment"; needs to be numeric, not factor
  # groupVar needs to be factor to ensure sample splitting is working as 
  # expected; default: no stratified splitting
  # works with the following methods: ranger, enet, xgbTree, nnet, pcaNNet
  # numGroups argument partitions observations into equally sized groups
  # quantiles argument partitions observations based on quantiles
  # propScores needs to be numeric vector with propensity scores corresponding
  # to the data
  ###
  
  if (method %in% c("nnet", "pcaNNet")) {
    neuralNet <- TRUE
  } else {
    neuralNet <- FALSE
  }
  
  # get group variable as the first variable; needed later on for within transformation
  setcolorder(datatable, groupVar)
  
  ### STEP -1: CONTROLING TRAINING AND FIXING / SETTING VARIOUS PARAMETERS
  
  # set trainControl
  trControl <- trainControl(method="adaptive_cv", # 5-fold cross-validation
                              number = 10,
                              adaptive = list(min=5, alpha=0.05,
                                              method="gls", complete=FALSE),
                              search = "random",
                              savePredictions = "none",
                              returnData = F,
                              returnResamp = "none",
                              trim = T) 
  
  # fix seed to replicate results
  set.seed(seed)
  
  # save number of rows in dataset
  rowNumber <- nrow(datatable)
  
  # create empty list objects for results
  blpList <- list()
  gatesList <- list()
  clanList <- list()
  lambdaBlpList <- list()
  lambdaGatesList <- list()

  ### STEP 0: FIX NUMBER OF SPLITS AND SIGNIFICANCE LEVEL
  numSplits <- numSplits
  sigLevel <- sigLevel
  
  
  ### STEP 1: PROPENSITY SCORES SKIPPED
  
  ### STEP 2: SPLIT DATA IN HALF
  
  # iterate over number of splits
  for (S in 1:numSplits) {
    
    numOfGroups <- numGroups
    set.seed(seed)
    seed <- seed+1
    
    if (verbose > 0) {
      print("Preprocessing...")
    }
    
    
    datatable[, propScore := propensityScores]
    
    # sample randomly from indices without replacement
    if (!is.null(groupVar)) {
      drawn <- createDataPartition(y = datatable[[groupVar]], times = 1, p = .5, list = F)
    } else {
      drawn <- createDataPartition(y = datatable[[treatVar]], times=1, p=.5, list=F)
    }
    
    # split data into auxiliary and main sample
    auxSample  <- datatable[drawn]
    mainSample <- datatable[-drawn]
    
    # normalize outcome variable to be between 0 and 1 when used with Neural Network
    if (neuralNet == T) {
      # formula: z = (x - min(x)) / (max(x) - min(x))
      auxMin <- min(auxSample[, ..outcomeVar])
      auxMax <- max(auxSample[, ..outcomeVar])
      out <- auxSample[, ..outcomeVar]
      set(auxSample, j=outcomeVar, value=((out - auxMin) / (auxMax - auxMin)))
    }
    
    # TRAIN MACHINE LEARNING MODELS
    # using auxiliary sample
    
    # save and remove propensity scores
    auxSample[, propScore := NULL]
    propScores <- mainSample[["propScore"]]
    mainSample[, propScore := NULL]
    
    # A: split data into treatment and control group  
    auxTreat <- auxSample[get(treatVar) == 1]
    auxControl <- auxSample[get(treatVar) == 0]
    
    # save treatment indicator and remove from data
    treatIndicatorAux <- as.vector(as.matrix(auxSample[, ..treatVar]))
    treatIndicatorMain <- as.vector(as.matrix(mainSample[, ..treatVar]))
    auxSample[, (treatVar) := NULL]
    mainSample[, (treatVar) := NULL] # remove treatment indicator from mainSample
    
    # save outcome variable and predictors separately
    # auxiliary sample
    y_auxTreat <- as.vector(as.matrix(auxTreat[, ..outcomeVar]))
    x_auxTreat <- auxTreat[, !c(treatVar, outcomeVar, groupVar), with=F]
    y_auxControl <- as.vector(as.matrix(auxControl[, ..outcomeVar]))
    x_auxControl <- auxControl[, !c(treatVar, outcomeVar, groupVar), with=F]

    # main sample
    x_main <- mainSample[, !c(outcomeVar, groupVar), with=F]
    y_main <- as.vector(as.matrix(mainSample[, ..outcomeVar]))
    
    # remove colsNotToUse from training data
    if (!is.null(colsNotToUse)) {
      x_auxTreat <- x_auxTreat[, !colsNotToUse, with=F]
      x_auxControl <- x_auxControl[, !colsNotToUse, with=F]
      x_main <- x_main[, !colsNotToUse, with=F]
    }
    
    # B: train and tune separate models for both groups
    # train treatment group model
    if (verbose > 0) {
      print("Training Treatment Group Model...")
    }
    
    if (method == "enet") {
      modelTreat <- train(x = x_auxTreat,
                          y = y_auxTreat,
                          method = method,
                          trControl = trControl,
                          preProcess = preprocess,
                          tuneLength = 15)
      
    } else if (method == "xgbTree") {
      modelTreat <- train(x = x_auxTreat,
                          y = y_auxTreat,
                          method = method,
                          trControl = trControl,
                          preProcess = preprocess, 
                          tuneLength = 15,
                          nthread = detectCores()-2)
      
    } else if (method == "ranger") {
      modelTreat <- train(x = x_auxTreat,
                          y = y_auxTreat,
                          method = method,
                          trControl = trControl,
                          preProcess = preprocess,
                          tuneLength = 15, 
                          num.threads = detectCores()-2)
      
      
    } else if (method == "nnet" | method == "pcaNNet") {
      modelTreat <- train(x = x_auxTreat,
                          y = y_auxTreat,
                          method = method,
                          trControl = trControl,
                          preProcess = preprocess, # scale data to be between 0 and 1
                          linout = neuralNet,
                          tuneLength = 15,
                          trace=FALSE)
    }
    
    
    # train control group model
    if (verbose > 0) {
      print("Training Control Group Model...")
    }
      
    if (method == "enet") {
      modelControl <- train(x = x_auxControl,
                            y = y_auxControl,
                            method = method,
                            trControl = trControl,
                            tuneLength = 15,
                            preProcess = "range") # scale data to be between 0 and 1)
      
    } else if (method == "xgbTree") {
      modelControl <- train(x = x_auxControl,
                            y = y_auxControl,
                            method = method,
                            trControl = trControl,
                            tuneLength = 15,
                            preProcess = "range", # scale data to be between 0 and 1
                            nthread = detectCores()-2)
      
    } else if (method == "ranger") {
      modelControl <- train(x = x_auxControl,
                            y = y_auxControl,
                            method = method,
                            trControl = trControl,
                            tuneLength = 15,
                            preProcess = "range", # scale data to be between 0 and 1
                            num.threads = detectCores()-2)
      
      
    } else if (method == "nnet" | method == "pcaNNet") {
      modelControl <- train(x = x_auxControl,
                            y = y_auxControl,
                            method = method,
                            trControl = trControl,
                            tuneLength = 15, 
                            preProcess = "range", # scale data to be between 0 and 1
                            linout = neuralNet,
                            trace=FALSE)
    }
    
    # C: predict BCA and calculate CATE proxy
    # using main sample
    bca <- predict(modelControl,
                   x_main)
    
    pred_treatModel <- predict(modelTreat,
                               x_main)
    
    # reverse normalization of outcome variable again using min and max from auxiliary sample
    if (neuralNet == T) {
      # formula: x = z*(max(x) - min(x)) + min(x)
      bca <- bca * (auxMax - auxMin) + auxMin
      pred_treatModel <- pred_treatModel * (auxMax - auxMin) + auxMin
    }
    
    # calculate proxy treatment effects S(Z)
    cateProxy <- pred_treatModel - bca
    
    # save results in data.table
    proxy <- data.table(bca = bca,
                        cate = cateProxy)
    setnames(proxy, c("bca", "cate"))
    
    # calculate expectation of S(Z), E[S(Z)]
    expectationS <- mean(cateProxy)
    
    # add other relevant variables back to the data
    proxy[, treatment := treatIndicatorMain] # treatment indicator
    proxy[, cateMinusExp := cate - expectationS]
    
    
    ### ESTIMATE BLP PARAMETERS USING WEIGHTED OLS 
    # using main sample
    
    if (verbose > 0) {
      print("Estimate BLP...")
    }
    
    # prepare data
    proxy[, propScore := propScores] # add propensity scores again
    proxy[, treat := treatment - propScore] # calculate [D-p(Z)]
    proxy[, outcome := y_main]
    proxy[, weights := (1 / (propScore*(1-propScore)))]
    
    # remove observations with infinity weights
    estSample <- proxy[weights < Inf]
    
    # estimate BLP parameters
    blp_weightedLinReg <- biglm(outcome ~ bca + cate + treat + treat:cateMinusExp,
                                weights = ~ weights,
                                data = estSample)
    
    ### ESTIMATE GATES PARAMETERS
    if (verbose > 0) {
      print("Estimate GATES...")
    }
    
    if (!is.null(numOfGroups) & !is.null(quantiles)) {
      
      err <- "Only one argument of numOfGroups or quantiles can be used, not both."
      print(err)
      return(err)
      
    } else if (is.null(numOfGroups)) {
      
      # find quantiles
      numOfGroups <- length(quantiles)-1
      q <- quantcut(proxy$cate, q=quantiles)
      
      # add factor variable for quantiles to proxy and create dummies from it
      proxy[, group := q]
      proxy[, group := factor(group, labels = 1:numOfGroups)]
      proxy <- dummy_cols(proxy, select_columns = "group", remove_selected_columns = F)
      
    } else if (is.null(quantiles)) {
      
      # find quantiles
      q <- quantcut(proxy$cate, q=numOfGroups)
      
      # add factor variable for quantiles to proxy and create dummies from it
      proxy[, group := q]
      proxy[, group := factor(group, labels = 1:numOfGroups)]
      proxy <- dummy_cols(proxy, select_columns = "group", remove_selected_columns = F)
    }
    
    # remove observations with infinity weights
    estSample <- proxy[weights < Inf]
    
    # prepare GATES estimation
    groupCols <- names(proxy)[(ncol(proxy)-(numOfGroups-1)):ncol(proxy)]
    groupInteractions <- paste0(groupCols, ":treat")
    groupColsCollapsed <- paste(groupInteractions, collapse = " + ")
    frmlaAux <- paste0("outcome", " ~ bca + cate")
    frmla <- formula(paste(frmlaAux, groupColsCollapsed, sep = " + "))
    
    # weighted GATES regression
    gates_weightedLinReg <- biglm(formula = frmla, 
                                  weights =  ~ weights,
                                  data = estSample)
    
    
    # ESTIMATE CLAN PARAMETERS
    if (verbose > 0) {
      print("Classification Analysis...")
    }
    
    # convert group column back to numeric
    proxy[, group := as.numeric(group)]
    
    if (!is.null(CLANvars)) {
      
      # add group column to main sample
      mainSample[, group := proxy$group]
      
      # characterize least (min) and most (max) affected group
      minGroupNum <- proxy[, min(group)]
      maxGroupNum <- proxy[, max(group)]
      
      # set group parameters to missing for all groups except most and least affected
      mainSample[((group != minGroupNum) & (group != as.integer(maxGroupNum))), group := NaN]
      
      # calculate t-tests
      testResults <- lapply(mainSample[, (CLANvars), with=F],
                            function(x) t.test(x ~ mainSample$group, var.equal = TRUE))
      
    } else if (is.null(CLANvars)) {
      
      # create empty list object
      testResults <- "No CLAN arguments in function call."
    }
    
    if (verbose > 0) {
      print("Compute performance measures...")
    }
      
    ### COMPUTE PERFORMANCE MEASURES FOR THE ML METHODS
    
    # BLP Criterion
    
    # save coefficient and variance of S(Z) in variables
    beta2 <- as.numeric(coef(blp_weightedLinReg)[5])
    varSZ <- proxy[, stats::var(cate)]
    
    # compute performance measure based on BLP
    lambda_blp <- abs(beta2)^2 * varSZ
    
    
    # GATES Criterion
    
    # save coefficients
    gates_coeffs <- as.numeric(coef(gates_weightedLinReg)[4:(3+numOfGroups)])
    
    
    # calculate relative frequencies
    relGroupFreq <- proxy[, .N / nrow(proxy), by=group]
    setorder(relGroupFreq, by=-group)
    
    # convert to numeric vector
    relGroupFreq <- as.numeric(relGroupFreq$V1)
    
    # calculate actual performance measure
    lambda_gates <- sum((gates_coeffs**2) * relGroupFreq)
    
    
    ### ADD ALL RESULTS TO RESULTS LISTS
    blpList[S] <- list(blp_weightedLinReg)
    gatesList[S] <- list(gates_weightedLinReg)
    clanList[S] <- list(testResults)
    lambdaBlpList[S] <- list(lambda_blp)
    lambdaGatesList[S] <- list(lambda_gates)
    
    # print indicator
    statement <- paste0("Completed iteration: ", S, "/", numSplits)
    print(statement)
    
  }
  
  # combine all results in single list
  resultsList <- list("BLP_Results" = blpList, 
                      "GATES_Results" = gatesList, 
                      "CLAN_Results" = clanList, 
                      "lambda_BLP" = lambdaBlpList, 
                      "lambda_GATES" = lambdaGatesList)
  
  
  
  # return results list object
  return(resultsList)
  
}
extract_blp <- function(cddfObject, sigLevel = 0.05) {
  
  # initialize empty vectors to improve speed
  intercept <- rep(NaN, length(cddfObject["BLP_Results"]))
  bca <- rep(NaN, length(cddfObject["BLP_Results"]))
  cate <- rep(NaN, length(cddfObject["BLP_Results"]))
  beta1 <- rep(NaN, length(cddfObject["BLP_Results"]))
  beta2 <- rep(NaN, length(cddfObject["BLP_Results"]))
  
  beta1_lower <- rep(NaN, length(cddfObject["BLP_Results"])) 
  beta1_upper <- rep(NaN, length(cddfObject["BLP_Results"])) 
  beta1_pvalue <- rep(NaN, length(cddfObject["BLP_Results"]))
  beta2_lower <- rep(NaN, length(cddfObject["BLP_Results"]))
  beta2_upper <- rep(NaN, length(cddfObject["BLP_Results"]))
  beta2_pvalue <- rep(NaN, length(cddfObject["BLP_Results"]))
  
  # loop over every estimation results object and extract coefficients
  for (blp in 1:length(cddfObject[["BLP_Results"]])) {
    
    # fill coefficients vectors iteratively
    intercept[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat[, "Coef"]["(Intercept)"]
    bca[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat[, "Coef"]["bca"]
    cate[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat[, "Coef"]["cate"]
    beta1[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat[, "Coef"]["treat"]
    beta2[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat[, "Coef"]["treat:cateMinusExp"]
    
  }
  
  # loop over every estimation results object and extract confidence intervals and p-values
  for (blp in 1:length(cddfObject[["BLP_Results"]])) {
    
    # fill vectors iteratively
    beta1_lower[blp] <- confint(cddfObject[["BLP_Results"]][[blp]], parm="treat", level = 1-sigLevel)[1]
    beta1_upper[blp] <- confint(cddfObject[["BLP_Results"]][[blp]], parm="treat", level = 1-sigLevel)[2]    
    beta1_pvalue[blp] <- summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat["treat","p"]
    
    beta2_lower[blp] <- confint(cddfObject[["BLP_Results"]][[blp]], parm="treat:cateMinusExp", level = 1-sigLevel)[1]
    beta2_upper[blp] <- confint(cddfObject[["BLP_Results"]][[blp]], parm="treat:cateMinusExp", level = 1-sigLevel)[2]
    beta2_pvalue[blp] <-  summary(cddfObject[["BLP_Results"]][[blp]])[2][1]$mat["treat:cateMinusExp","p"]
    
  }
  
  # combine single vectors to results data.table
  blpResults <- data.table(intercept = intercept,
                           control_bca = bca,
                           control_cate = cate,
                           ate = beta1,
                           ateLower = beta1_lower,
                           ateUpper = beta1_upper,
                           atePValue = beta1_pvalue,
                           hte = beta2,
                           hteLower = beta2_lower,
                           hteUpper = beta2_upper,
                           htePValue = beta2_pvalue)
  
    # calculate median of point estimates
  beta1 <- blpResults[, median(ate)]
  beta2 <- blpResults[, median(hte)]
  
  ### get upper and lower median of confidence interval
  if (is.even(nrow(blpResults))) {
    
    lower <- floor(nrow(blpResults)/2)
    upper <- floor(nrow(blpResults)/2) + 1
    
    beta1Lower <- blpResults[order(ateLower)][upper]$ateLower
    beta1Upper <- blpResults[order(ateUpper)][lower]$ateUpper
    
    beta2Lower <- blpResults[order(hteLower)][upper]$hteLower
    beta2Upper <- blpResults[order(hteUpper)][lower]$hteUpper
    
  } else {
    
    beta1Lower <- blpResults[, median(ateLower)]
    beta1Upper <- blpResults[, median(ateUpper)]
    
    beta2Lower <- blpResults[, median(hteLower)]
    beta2Upper <- blpResults[, median(hteUpper)]
    
  }
  
  
  ### get adjusted p-values
  
  # beta1 - ATE coefficient
  b1p <- blpResults[order(atePValue)]
  if (is.even(nrow(b1p))) {
    beta1PValue <- b1p[floor(nrow(b1p) / 2)]$atePValue
  } else {
    beta1PValue <- b1p[, median(atePValue)]
  }
  
  # beta2 - HTE coefficient
  b2p <- blpResults[order(htePValue)]
  if (is.even(nrow(b2p))) {
    beta2PValue <- b2p[floor(nrow(b2p)/ 2)]$htePValue
  } else {
    beta2PValue <- b2p[, median(htePValue)]
  }
  
  # combine results in final data.table
  finalResults <- data.table(ate = beta1,
                             ateLowerBound = beta1Lower,
                             ateUpperBound = beta1Upper,
                             ateAdjP = beta1PValue * 2,
                             hte = beta2,
                             hteLowerBound = beta2Lower,
                             hteUpperBound = beta2Upper,
                             hteAdjP = beta2PValue * 2)
  
    return(list(blpDT = blpResults, 
              VEINResults = finalResults))
}

extract_gates <- function(cddfObject, sigLevel = 0.05, numOfGroups = 10) {
  
  # get number of coefficients
  numCoefficients <- numOfGroups
  
  # initialize empty data.table
  gatesCoefs <- data.table()
  diffs <- c()
  pValsDiff <- c()
  ftests <- c()
  
  # loop over all splits
  for (S in 1:length(cddfObject[["GATES_Results"]])) {
    
    # initialize empty vector
    coefs <- c()
    pVals <- c()
    i <- 0
    
    # loop over all gates coefficients
    for (group in 4:(3+numCoefficients)) {
      
      i <- i+1
      
      # extract single coefficients and p-values
      coefs[i] <- coefficients(cddfObject[["GATES_Results"]][[S]])[group]
      pVals[i] <- summary(cddfObject[["GATES_Results"]][[S]])[2][1]$mat[(group),"p"]
      
    }
    
    pVals <- t(as.matrix(pVals))
    coefs <- t(coefs)
    
    # extract confidence intervals
    lower <- as.matrix(confint(cddfObject[["GATES_Results"]][[S]], level = 1-sigLevel)[4:(numOfGroups+3)])
    upper <- as.matrix(confint(cddfObject[["GATES_Results"]][[S]], level = 1-sigLevel)[(7+numOfGroups):(6+numOfGroups+numOfGroups)])
    
    # combine coefficients, lower and upper CI and p-values
    coefficients <- c(coefs, lower, upper, pVals)
    
    # rbind vectors
    gatesCoefs <- rbind(gatesCoefs, data.frame(as.list(coefficients)), use.names=F)
    
    # compare coefficients of most and least affected group and perform t-test / Chi-square test
    coefNames <- names(coef(cddfObject[["GATES_Results"]][[S]]))
    toTest <- coefNames[c(4, (3+numOfGroups))]
    toTest <- paste(toTest, collapse = " = ")
    diff <- coefs[1] - coefs[numOfGroups]
    pValDiff <- linearHypothesis(cddfObject[["GATES_Results"]][[S]], toTest)[, "Pr(>Chisq)"][[2]]
    
    # test if all coefficients are equal to zero
    ftest <- linearHypothesis(cddfObject[["GATES_Results"]][[S]], coefNames[4:(numOfGroups+3)])[,"Pr(>Chisq)"][[2]]
    
    # store results
    diffs[S] <- diff
    pValsDiff[S] <- pValDiff
    ftests[S] <- ftest
    
  }
  
  # create data.table and name columns accordingly
  colNames <- paste("group", seq(1:numOfGroups), sep="_")
  lowerNames <- paste("lower_group", seq(1:numOfGroups), sep = "_")
  upperNames <- paste("upper_group", seq(1:numOfGroups), sep = "_")
  pValNames <- paste("pVal_group", seq(1:numOfGroups), sep="_")
  setnames(gatesCoefs, old = names(gatesCoefs), new = c(colNames, lowerNames, upperNames, pValNames))
  
  # extract medians of point estimates and confidence intervals -> VEIN Results
  pointEst <- gatesCoefs[, lapply(.SD, median), .SDcols = names(gatesCoefs)[! (grepl("lower", names(gatesCoefs)) | grepl("upper", names(gatesCoefs)) | grepl("pVal", names(gatesCoefs)))]]
  lower <- gatesCoefs[, lapply(.SD, upperMedian), .SDcols = names(gatesCoefs)[grepl("lower", names(gatesCoefs))]]
  upper <- gatesCoefs[, lapply(.SD, lowerMedian), .SDcols = names(gatesCoefs)[grepl("upper", names(gatesCoefs))]]
  pValues <- gatesCoefs[, lapply(.SD, lowerMedian), .SDcols = names(gatesCoefs)[grepl("pVal", names(gatesCoefs))]]*2
  
  # combine results and return results
  results <- cbind(pointEst, lower, upper, pValues)
  
  results <- data.table(results)
  res <- data.table(group = results[, 1:numOfGroups])
  res <- rbindlist(list(res, results[, (1+numOfGroups):(2*numOfGroups)]), use.names = F)
  res <- rbindlist(list(res, results[, (1+numOfGroups*2):(3*numOfGroups)]), use.names = F)
  res <- rbindlist(list(res, results[, (1+numOfGroups*3):(4*numOfGroups)]), use.names = F)
  resDF <- data.frame(res, row.names = c("ATE", "lower", "upper", "p-Value")) 
  
  diff <- median(diffs)
  pValDiff <- lowerMedian(pValsDiff)*2
  ftest <- lowerMedian(ftests)*2
  
  tests <- data.table(diffATE = diff,
                      pValue = pValDiff,
                      pValueFTest = ftest)
  
  return(list(resDF, tests))
  
}

extract_performance_metrics <- function(cddfObject) {
  
  lambda_blp <- median(unlist(cddfObject[["lambda_BLP"]]))
  gates_blp <- median(unlist(cddfObject[["lambda_GATES"]]))
  
  result <- data.table(lambdaBlp = lambda_blp,
                       gatesBlp = gates_blp)
  
  return(result)
  
}

extract_test_results <- function(cddfObject) {
  
  # clanVars
  clanVars <- names(cddfObject[["CLAN_Results"]][[1]])
  
  # initialize empty results data.table
  results <- data.table()
  results[, (clanVars) := rep(NaN, 6)]
  
  for (var in clanVars) {  
    
    tmpRes <- data.table(groupMin = numeric(),
                         groupMax = numeric(),
                         difference = numeric(),
                         pValue = numeric(),
                         lower = numeric(),
                         upper = numeric())
    
    for (split in 1:length(cddfObject[["CLAN_Results"]])) {
      
      # extract means and p-values / CI
      groupmin_mean <- cddfObject[["CLAN_Results"]][[split]][[var]]$estimate[[1]]
      groupmax_mean <- cddfObject[["CLAN_Results"]][[split]][[var]]$estimate[[2]]
      diff <- groupmin_mean - groupmax_mean
      p_val <- cddfObject[["CLAN_Results"]][[split]][[var]]$p.value
      lower <- cddfObject[["CLAN_Results"]][[split]][[var]]$conf.int[[1]]
      upper <- cddfObject[["CLAN_Results"]][[split]][[var]]$conf.int[[2]]
      
      # save in data.table and rbind results
      tmp <- data.table(groupMin = groupmin_mean,
                        groupMax = groupmax_mean,
                        difference = diff,
                        pValue = p_val,
                        lower = lower,
                        upper = upper)
      
      tmpRes <- rbindlist(list(tmpRes, tmp))
      
    }
    
    # calculate medians for results and save in data.table
    groupMin <- tmpRes[, median(groupMin)]
    groupMax <- tmpRes[, median(groupMax)]
    diff <- tmpRes[, median(diff)]
    pValue <- tmpRes[, lowerMedian(pValue)*2]
    lower <- tmpRes[, upperMedian(lower)]
    upper <- tmpRes[, lowerMedian(upper)]
    
    vec <- c(groupMin, groupMax, diff, pValue, lower, upper)
    
    #return(list(result, clanVars, vec))
    # append results to results vector
    for (row in 1L:6L) {
      set(results, i = row, j = var, value = vec[[row]])
    }
    
  }
  
  # add description of values to data.table
  results[, parameter := c("groupMin", "groupMax", "difference", "pValue", "lower", "upper")]
  setcolorder(results, "parameter")
  
  return(results)
  
}




