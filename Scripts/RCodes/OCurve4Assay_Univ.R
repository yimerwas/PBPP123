
OCurve4Assay_Univ <- function(StanModel, ModelName, SampleSize = seq(3, 20, by = 2), Specs, Confidence, Coverage, SeedNum) {
  
  # Helper function for rounding
  round2 <- function(x, n) {
    posneg = sign(x)
    z = abs(x) * 10^n
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z / 10^n
    z * posneg
  }
  
  # Set tolerance interval parameters
  Betat <- Coverage      # Content
  Gammat <- Confidence    # Confidence
  
  # Specifications
  specs <- Specs
  
  # Extract posterior samples from StanModel
  PosteriorSamples <- extract(StanModel)
  K <- SampleSize
  Kmax <- max(K)
  
  # Posterior distribution of residual SD
  RandomError <- PosteriorSamples$sigma
  X1 <- matrix(NA, ncol = Kmax, nrow = length(RandomError))
  
  set.seed(SeedNum)
  
  # Generate random errors
  for (i in 1:length(RandomError)) {
    X1[i, ] <- rnorm(Kmax, 0, RandomError[i])
  }
  
  # Calculate Residual Tolerance Intervals
  ETM <- NULL
  for (j in 1:length(K)) {
    etm <- matrix(NA, nrow = length(RandomError), ncol = 3)
    etm[, 1] <- K[j]
    kvalue <- K.factor(K[j], alpha = (1 - Gammat), P = Betat, side = 2, method = "EXACT", m = 50)
    
    for (i in 1:length(RandomError)) {
      etm[i, 2:3] <- mean(X1[i, 1:K[j]]) + c(-1, 1) * sd(X1[i, 1:K[j]]) * kvalue
      
      if (i %% 200 == 0) {
        cat("N =", K[j], "(", j, "/", length(K), "),", "chain(", i, "/", length(RandomError), ")", "\n")
      }
    }
    
    ETM <- rbind(ETM, etm)
  }
  
  # Save the tolerance object
  saveRDS(ETM, paste("Results/", paste0(ModelName, "_ETM.rds"), sep = ""))
  
  # Calculate probabilities of success (PoS)
  Betas <- seq(from = specs[1], to = specs[2], by = 0.1)
  POS <- NULL
  
  for (beta in Betas) {
    ETMData <- as.data.frame(ETM)
    ETMData$mean <- beta
    names(ETMData) <- c("k", "lower", "upper", "mean")
    ETMData$lower <- ETMData$lower + beta
    ETMData$upper <- ETMData$upper + beta
    ETMData$success <- ifelse((ETMData$lower >= specs[1]) & (ETMData$upper <= specs[2]), 1, 0)
    
    pos <- ETMData %>%
      group_by(mean, k) %>%
      dplyr::summarize(pos = mean(success))
    
    POS <- rbind(POS, pos)
  }
  
  # Save success/failure object for true batch means (tolerance intervals)
  saveRDS(POS, paste("Results/", paste0(ModelName, "_POS.rds"), sep = ""))
  
  # Recall the tolerance object
  ETM <- readRDS(paste("Results/", paste0(ModelName, "_ETM.rds"), sep = ""))
  
  # Read success/failure of true batch means (tolerance intervals)
  POS <- readRDS(paste("Results/", paste0(ModelName, "_POS.rds"), sep = ""))
  
  # Results graphs
  # Generate the distribution of a future random batch
  RanFun <- function(x, y) {
    rnorm(1, mean = x, sd = y)
  }
  
  # Process mean
  ProcessMean <- PosteriorSamples$b_Intercept
  Mu_ProcessMean <- mean(ProcessMean)
  SigmaBatch <- PosteriorSamples$sd_b
  RandomBatchMean <- mapply(RanFun, ProcessMean, SigmaBatch)
  Process_Mu_CL <- quantile(RandomBatchMean, c(0.025, 0.975))
  
  # Determine POS for the process mean
  POS_Mu_ProcessMean <- POS %>%
    mutate(pos = round2(pos, 2)) %>%
    filter(mean %in% seq(95, 105, by = 0.1) & pos >= 0.95) %>%
    filter(mean == round2(Mu_ProcessMean, 1))
  
  # Create OCurvePlot
  OCurvePlot <- POS %>%
    mutate(n = as.factor(k)) %>%
    ggplot(aes(x = mean, y = pos, group = n, col = n)) +
    geom_line() + 
    labs(x = "True batch mean", y = "PoS", title = ModelName) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(breaks = seq(specs[1], specs[2], 1)) +
    geom_vline(aes(xintercept = Mu_ProcessMean), linetype = "dashed", col = "blue") +
    geom_hline(aes(yintercept = 0.95), col = "red") +
    annotate("text", x = (Mu_ProcessMean + 0.2), y = 0.45, label = "Process mean", color = "blue", angle = 90) +
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = Process_Mu_CL[1], xmax = Process_Mu_CL[2], alpha = 0.1, fill = "brown") +
    annotate("text", x = (Process_Mu_CL[1] - 0.2), y = 0.45, label = "Lower 95% PL", color = "brown", angle = 90) +
    annotate("text", x = (Process_Mu_CL[2] + 0.2), y = 0.45, label = "Upper 95% PL", color = "brown", angle = 90) +
    annotate("text", x = 104, y = 0.9, label = ifelse(is.infinite(min(POS_Mu_ProcessMean$k)), paste("Size>", 40), paste("Size=", min(POS_Mu_ProcessMean$k))), color = "blue", angle = 0, size = 3) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  # Get minimum sample size
  POS_MinimumSample <- POS_Mu_ProcessMean %>% filter(k == min(k))
  
  return(list(OCurvePlot, POS_MinimumSample$k))
}