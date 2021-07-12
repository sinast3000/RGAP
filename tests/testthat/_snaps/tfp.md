# TFP model

    Call:
    TFPmodel(tsl = tsList, cycle = "RAR2")
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    Object is a valid object of class TFPmodel.

# TFP MLE fit

    Call:
    fitTFP(model = model)
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic  p-value
      cA        5.34e-01       1.53e-01        3.49 0.000478
      cSigma    8.62e-05       2.26e-05        3.82 0.000133
      cTau      8.70e+00       3.26e+00        2.67 0.007610
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdOmega    5.02e-03       3.85e-03        1.30   0.193
      tdPhi      8.94e-01       1.02e-01        8.77   0.000
      tdSigma    7.77e-06       6.68e-06        1.16   0.245
    
    cubs
              Coefficient Standard Error t-statistic p-value
      cuC0       2.717588       2.96e-01       9.168 0.00000
      cuConst    0.002091       9.14e-03       0.229 0.81915
      cuSigma    0.000153       5.62e-05       2.731 0.00632
      RMSE: 0.0289
      R2: 0.367
      Box-Ljung test: X-squared = 1.77, df = 4, p-value = 0.777
    
             loglik             AIC             BIC             HQC signal-to-noise 
          192.86842      -367.73683      -353.48516      -362.76262         0.09014 

# TFP bayesian fit

    Call:
    fitTFP(model = model, method = "bayesian", R = 1000, thin = 2, 
        MLEfit = fit)
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    
    	MCMC estimation results
    
    cycle
                 Mean   Median 85% HPDI-LB 85% HPDI-UB
      cA     0.425372 0.427395    1.85e-01    6.55e-01
      cSigma 0.000155 0.000121    6.72e-05    2.28e-04
      cTau   7.930563 7.538559    2.84e+00    1.14e+01
    
    trend
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      tdOmega 6.60e-03 6.32e-03    7.60e-04    1.04e-02
      tdPhi   9.82e-01 9.84e-01    9.73e-01    9.90e-01
      tdSigma 6.37e-07 6.02e-07    3.65e-07    8.41e-07
    
    cubs
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      cuC0    1.40e+00 1.40e+00    1.378517    1.419370
      cuConst 8.76e-05 5.68e-05   -0.000871    0.000891
      cuSigma 4.37e-04 4.24e-04    0.000266    0.000551
    signal-to-noise 
           0.004112 

