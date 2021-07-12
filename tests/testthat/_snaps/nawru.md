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
      cPhi1        1.269         0.1303       9.737 0.00e+00
      cPhi2       -0.403         0.5024      -0.801 4.23e-01
      cSigma       0.211         0.0411       5.144 2.69e-07
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdSigma      0.0015        0.00186       0.807    0.42
    
    phillips curve
              Coefficient Standard Error t-statistic  p-value
      pcC0      -3.47e-03       1.47e-03     -2.3564 1.85e-02
      pcConst    6.61e-05       2.07e-03      0.0319 9.75e-01
      pcSigma    1.22e-04       2.27e-05      5.3588 8.38e-08
      pcddws     9.86e-01       9.69e-02     10.1823 0.00e+00
      RMSE: 0.0116
      R2: 0.627
      Box-Ljung test: X-squared = 7.05, df = 4, p-value = 0.133
    
             loglik             AIC             BIC             HQC signal-to-noise 
          1.389e+02      -2.618e+02      -2.452e+02      -2.554e+02       7.108e-03 

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
      cPhi1   0.7904  0.7870       0.705      0.8958
      cPhi2  -0.0727 -0.0698      -0.156      0.0242
      cSigma  0.2391  0.2338       0.175      0.2993
    
    trend
                 Mean   Median 85% HPDI-LB 85% HPDI-UB
      tdSigma 0.00048 0.000475    0.000373    0.000603
    
    phillips curve
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      pcC0    -0.07597 -0.07513     -0.1257     -0.0345
      pcConst  0.00167  0.00211     -0.0262      0.0306
      pcSigma  0.05928  0.05868      0.0464      0.0743
      pcddws   0.04965  0.05394     -0.6225      0.6758
    signal-to-noise 
           0.002007 

