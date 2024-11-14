Approximate_Ca0 = function(StanData_Historical,
                   StanModel,
                   NumberOfIterations = 5000,
                   Burnin = 1000,
                   Chains = 3,
                   Max_treedepth = 15,
                   Thin = 3,
                   Adapt_delta = 0.999,
                   a0_Increment = 0.05, 
                   seed = 1234567890){
  
  d <- data.frame(a0 = seq(0,1,by = a0_Increment), C = NA)
  for(i in seq(0,1,by = a0_Increment)){
    n = 0
    print(paste0(i * 100,"%: a0 = ", i))
    seed = seed
    success = F
    while(!success){
      seed = seed + 1
      n = n + 1
      if(n > 30){stop("Approximation of the normalising constant not possible without divergent transitions. Try increasing Adapt_delta above ", Adapt_delta, " or choosing a more informative prior distribution for the between-cluster SD")}
      StanData_Historical$a0=i

        result = suppressWarnings(rstan::sampling(StanModel, data = StanData_Historical, refresh = 0, control = list(adapt_delta = Adapt_delta, max_treedepth = Max_treedepth), cores = 3, iter = NumberOfIterations, thin = Thin, seed = seed, warmup = Burnin))
        
      t <- rstan::get_sampler_params(result, inc_warmup = F)
      divergent <- sum(t[[1]][,"divergent__"],t[[2]][,"divergent__"],t[[3]][,"divergent__"],t[[4]][,"divergent__"])
      success = ifelse(divergent == 0,T,F)
    }
    set.seed(seed)
    ddpcr::quiet(d$C[d$a0 == i] <- bridgesampling::bridge_sampler(result, method = "normal")$logml)
    d$divergent[d$a0 == i] <- divergent
  }
  
  g <- mgcv::gam(C ~ s(a0), data = d)
  a0_grid = seq(0,1,length = 10000)
  C_grid = predict(g, data.frame(a0 = a0_grid))
  return(cbind.data.frame(a0_grid=a0_grid, C_grid=C_grid))
}