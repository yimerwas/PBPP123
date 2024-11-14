#==================================================================
# Description: Running simulations for partial borrowing power priors
#==================================================================

## libraries
library(data.table)
library(rstan)
library(cmdstanr)

## cmd arguments
options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
simulation <- args[[1]]   ## S1, S2, S3, S4, S5
nsim <- as.integer(args[[2]])
ndone <- as.integer(args[[3]])

## stan models
# set_cmdstan_path(path = "/home/cmdstan-2.33.1")
blmm <- cmdstan_model(exe_file = "stan/BLMM_allScenarios")
blmm_ca0 <- stan_model(file = "stan/BLMM_normalizing_constant.stan")

## helper functions

#' Simulate observations from a linear mixed model 
#' @param B integer, number of batches
#' @param N integer, total number of observations (must be multiple of B)
#' @param beta0 process mean
#' @param sigmab numeric, inter-batch standard deviation
#' @param sigma0 numeric, intra-batch standard deviation  
#' @param seed (optional) integer random seed
#' @param name character, name to identify dataset
#' @return data.table with columns:
#' * \code{name}: dataset name
#' * \code{i}: integer batch number
#' * \code{j}: integer observation index within batch
#' * \code{y0}: simulated assay value
#' * \code{b0}: simulated random batch effect
#' @example
#' simLMM(10, 10, 100.4, 1.35, 1.07)
#' @author jchau
#' @export
simLMM <- function(N, B, beta0, sigmab, sigma0, seed, name = "current") {
  
  if(!missing(seed)) {
    set.seed(seed)
  }
  
  if(N %% B != 0) {
    stop("N must be a multiple of B")
  }
  
  b0 <- rnorm(B, mean = 0, sd = sigmab)
  y0 <- rnorm(N, beta0 + rep(b0, each = N / B), sd = sigma0)
  
  data.table(
      name = name,
      i = rep(seq_len(B), each = N / B),
      j = rep(seq_len(N / B), times = B),
      y0 = y0,
      b0 = b0
  )
}

#' Calculate normalizing constant C(a0) 
#'
#' @param model rstan object of class \code{stanmodel} 
#' @param a0 numeric vector of a0's at which to evaluate C(a0)
#' @param data data.table with historic data, contains columns \code{i} and \code{y0}
#' @param ... additional parameters passed to \code{rstan::sampling}
#' @return numeric log marginal likelihood
#' @author jchau
#' @importFrom bridgesampling bridge_sampler
#' @export
calc_ca0 <- function(model, a0 = c(0.5, 1), data, ...) {
  vapply(a0, FUN = function(a0_fixed) {
        ## posterior samples
        stanfit <- sampling(
            object = model,
            data = list(
                N = nrow(data),
                B = max(data[["i"]]),
                Y0 = data[["y0"]],
                J0 = data[["i"]],
                a0_fixed = a0_fixed
            ),
            chains = 2,
            iter = 2000,
            warmup = 1000,
            cores = 1,
            refresh = 0,
            show_messages = FALSE,
            init = function() {
              list(
                  sigma0 = 1.07,
                  beta00 = 100.4,
                  sigma_b0 = 1.35
              )
            },
            ...
        )
        ## log marginal likelihood
        invisible(
            capture.output(
                ca0 <- tryCatch({
                      bridgesampling::bridge_sampler(stanfit, method = "normal")$logml
                    }, error = function(e) e)
            )
        )
        if(!inherits(ca0, "error")) {
          return(ca0) 
        } else {
          return(NA_real_)
        }
      },
      FUN.VALUE = numeric(1)
  )
}

#' Extract a0, sigma statistics from cmdstanr-fit
#' @param fit R6-object of class \code{CmdStanMCMC}
#' @param simulation character simulation name, (S1, S3, S3, S4 or S5)
#' @param run integer simulation run
#' @param scenario character simulation scenario identifier
#' @return data.table
#' @author jchau
#' @export
getStats <- function(fit, simulation, run, scenario) {
  
  stats <- capture.output(fit$print(c("a0", "sigma", "beta00", "sigma_b0", "beta0", "sigma_b"), digits = 4), type = "output")
  stats <- as.data.table(read.table(text = paste(stats, collapse = "\n"), header = TRUE, sep = ""))
  
  ## fixed a0
  if(grepl("^PBPP1", scenario)) {
    a0 <- as.numeric(sub("^.+\\=", "", scenario))
    set(stats, i = 1L, j = c("mean", "median", "q5", "q95"), value = a0)
    set(stats, i = 1L, j = c("sd", "mad", "ess_bulk", "ess_tail"), value = 0)
    set(stats, i = 1L, j = "rhat", value = 1)
  }
  
  return(
      cbind(
          data.table(
              name = simulation,
              run = run,
              scenario = scenario
          ),
          stats
      )
  )
}

## initialize parameters
N <- c(100L, 100L)
B <- c(10L, 10L)
beta0 <- c(100.4, 100.4)
sigma_b <- c(1.35, 1.35)

if(simulation == "S1") {
  ## - matching parameter values
  sigma <- c(1.07, 1.07)
} else if(simulation == "S2") {
  ## - sigma0 > sigma
  sigma <- c(1.16, 0.99)
} else if(simulation == "S3") {
  ## - sigma0 < sigma
  sigma <- c(0.99, 1.16)
} else if(simulation == "S4") {
  ## - sigma0 << sigma
  sigma <- c(0.2, 4)
} else if(simulation == "S5") {
  ## - sigma0 >> sigma
  sigma <- c(4, 0.2)
} else {
  stop("Unknown simulation setting")
}

## run simulations
for(i in (ndone + 1):nsim) {
  
  ## simulate data
  historic <- simLMM(N[1], B[1], beta0[1], sigma_b[1], sigma[1], name = "historic")
  current <- simLMM(N[2], B[2], beta0[2], sigma_b[2], sigma[2], name = "current")
  
  ## calculate ln(C(0.5)) and ln(C(1)) used in PBPP3
  ln_c <- calc_ca0(
      model = blmm_ca0,
      a0 = c(0.5, 1),
      data = historic
  )
  
  ## initialize stan data
  data0 <- list(
      N = N,
      B = B,
      Y0 = historic[["y0"]],
      Y1 = current[["y0"]],
      J0 = historic[["i"]],
      J1 = current[["i"]],
      scenario = 0L,
      a0_fixed = 1,
      a0_shape = c(1, 1),
      ln_c05 = 0,
      ln_c1 = 0
  )
  
 for(a0 in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
   
   cat(sprintf("* Scenario 1, run %d, a0 = %g\n", i, a0))
   
   ## PBPP1
   samples <- tryCatch({
         blmm$sample(
             data = c(
                 list(
                     N = N,
                     B = B,
                     Y0 = historic[["y0"]],
                     Y1 = current[["y0"]],
                     J0 = historic[["i"]],
                     J1 = current[["i"]],
                     scenario = 1L,
                     a0_fixed = a0
                 ), 
                 data0[c("a0_shape", "ln_c05", "ln_c1")]
             ),
             show_messages = FALSE,
             show_exceptions = FALSE,
             refresh = 0,
             chains = 2,
             parallel_chains = 1,
             iter_warmup = 1000,
             iter_sampling = 1000
         )
       }, error = function(e) e)
   
   if(!inherits(samples, "error")) {
     stats <- getStats(samples, simulation, i, sprintf("PBPP1: a0=%g", a0))
     fwrite(stats, file = sprintf("results/%s_results.csv", simulation), append = TRUE)
   }
   
 }
  
  for(beta in c(1, 5, 10, 50)) {
    
    a0_shape <- c(beta, ifelse(beta < 50, beta, 20))
   cat(sprintf("* Scenario 2, run %d, beta(%d, %d)\n", i, a0_shape[1], a0_shape[2]))
   
   ## PBPP2
   samples <- tryCatch({
         blmm$sample(
             data = c(
                 list(
                     N = N,
                     B = B,
                     Y0 = historic[["y0"]],
                     Y1 = current[["y0"]],
                     J0 = historic[["i"]],
                     J1 = current[["i"]],
                     scenario = 2L,
                     a0_shape = a0_shape
                 ), 
                 data0[c("a0_fixed", "ln_c05", "ln_c1")]
             ),
             show_messages = FALSE,
             show_exceptions = FALSE,
             chains = 2,
             parallel_chains = 1,
             iter_warmup = 1000,
             iter_sampling = 1000
         )
       }, error = function(e) e)
   
   if(!inherits(samples, "error")) {
     stats <- getStats(samples, simulation, i, sprintf("PBPP2: beta(%d, %d)", a0_shape[1], a0_shape[2]))
     fwrite(stats, file = sprintf("results/%s_results.csv", simulation), append = TRUE)
   }
   
    ## PBPP3
    cat(sprintf("* Scenario 3, run %d, beta(%d, %d)\n", i, a0_shape[1], a0_shape[2]))
    
    if(!any(is.na(ln_c))) {
      
      samples <- tryCatch({
            blmm$sample(
                data = c(
                    list(
                        N = N,
                        B = B,
                        Y0 = historic[["y0"]],
                        Y1 = current[["y0"]],
                        J0 = historic[["i"]],
                        J1 = current[["i"]],
                        scenario = 3L,
                        a0_shape = a0_shape,
                        ln_c05 = ln_c[1],
                        ln_c1 = ln_c[2]
                    ), 
                    data0[c("a0_fixed")]
                ),
                show_messages = FALSE,
                show_exceptions = FALSE,
                chains = 2,
                parallel_chains = 1,
                iter_warmup = 1000,
                iter_sampling = 1000
            )
          }, error = function(e) e)
      
      if(!inherits(samples, "error")) {
        stats <- getStats(samples, simulation, i, sprintf("PBPP3: beta(%d, %d)", a0_shape[1], a0_shape[2]))
        fwrite(stats, file = sprintf("results/%s_results_3.csv", simulation), append = TRUE)
      }
      
    }
    
  }
  
}





















