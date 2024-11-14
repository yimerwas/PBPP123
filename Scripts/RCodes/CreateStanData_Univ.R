CreateStanData_Univ <- function(BrmsModel_Historical, BrmsModel_Current, a0, random=FALSE, normalized=FALSE){
  
  BrmsModel_Historical_StanData <- standata(BrmsModel_Historical)
  BrmsModel_Current_StanData <- standata(BrmsModel_Current)
  
  
  if(random==FALSE){
    
    StanData <- list(
      # Historical data
      N0 = BrmsModel_Historical_StanData$N,
      Y0 = BrmsModel_Historical_StanData$Y,
      N0_1 = BrmsModel_Historical_StanData$N_1,
      M0_1 = BrmsModel_Historical_StanData$M_1,
      J0_1 = BrmsModel_Historical_StanData$J_1,
      Z0_1_1 = BrmsModel_Historical_StanData$Z_1_1,
      prior_only0 = BrmsModel_Historical_StanData$prior_only,
      a0 = a0,
      #Current data
      N = BrmsModel_Current_StanData$N,
      Y = BrmsModel_Current_StanData$Y,
      N_1 = BrmsModel_Current_StanData$N_1,
      M_1 = BrmsModel_Current_StanData$M_1,
      J_1 = BrmsModel_Current_StanData$J_1,
      Z_1_1 = BrmsModel_Current_StanData$Z_1_1,
      prior_only = BrmsModel_Current_StanData$prior_only
    )}
  
    else {
        StanData <- list(
          # Historical data
          N0 = BrmsModel_Historical_StanData$N,
          Y0 = BrmsModel_Historical_StanData$Y,
          N0_1 = BrmsModel_Historical_StanData$N_1,
          M0_1 = BrmsModel_Historical_StanData$M_1,
          J0_1 = BrmsModel_Historical_StanData$J_1,
          Z0_1_1 = BrmsModel_Historical_StanData$Z_1_1,
          prior_only0 = BrmsModel_Historical_StanData$prior_only,
          a01 = a0[1],
          a02 = a0[2],
          #Current data
          N = BrmsModel_Current_StanData$N,
          Y = BrmsModel_Current_StanData$Y,
          N_1 = BrmsModel_Current_StanData$N_1,
          M_1 = BrmsModel_Current_StanData$M_1,
          J_1 = BrmsModel_Current_StanData$J_1,
          Z_1_1 = BrmsModel_Current_StanData$Z_1_1,
          prior_only = BrmsModel_Current_StanData$prior_only)  
  }
  
  return(StanData)
}
