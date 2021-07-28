# NAWRU model

    Call:
    NAWRUmodel(tsl = tsList, exoType = exoType)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    Object is a valid object of class NAWRUmodel.

# NAWRU MLE fit

    Call:
    fitNAWRU(model = model)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic  p-value
      cPhi1        1.270         0.1303       9.744 0.00e+00
      cPhi2       -0.403         0.5044      -0.799 4.24e-01
      cSigma       0.211         0.0411       5.145 2.68e-07
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdSigma     0.00152        0.00189       0.807    0.42
    
    phillips curve
              Coefficient Standard Error t-statistic  p-value
      pcC0      -3.48e-03       1.47e-03     -2.3645 1.81e-02
      pcConst    6.62e-05       2.08e-03      0.0319 9.75e-01
      pcSigma    1.21e-04       2.26e-05      5.3623 8.22e-08
      pcddws     9.83e-01       9.68e-02     10.1583 0.00e+00
      RMSE: 0.0116
      R2: 0.627
      Box-Ljung test: X-squared = 21, df = 10, p-value = 0.0209
    
             loglik             AIC             BIC             HQC signal-to-noise 
          1.389e+02      -2.618e+02      -2.452e+02      -2.554e+02       7.207e-03 

# NAWRU bayesian fit

    Call:
    fitNAWRU(model = model, method = "bayesian", R = 1000, thin = 2, 
        MLEfit = fit)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    
    	MCMC estimation results
    
    cycle
                Mean  Median 85% HPDI-LB 85% HPDI-UB
      cPhi1   0.7906  0.7872       0.705      0.8958
      cPhi2  -0.0726 -0.0698      -0.156      0.0242
      cSigma  0.2391  0.2338       0.175      0.2993
    
    trend
                 Mean   Median 85% HPDI-LB 85% HPDI-UB
      tdSigma 0.00048 0.000475    0.000369    0.000599
    
    phillips curve
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      pcC0    -0.07524 -0.07483     -0.1257     -0.0336
      pcConst  0.00155  0.00189     -0.0253      0.0322
      pcSigma  0.05929  0.05869      0.0464      0.0743
      pcddws   0.05198  0.05482     -0.6225      0.6758
              MRMSE signal-to-noise 
           0.046210        0.002006 

