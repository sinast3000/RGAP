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
    fit.TFPmodel(model = model)
    
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
      cA        5.34e-01       1.53e-01        3.50 0.000472
      cSigma    8.64e-05       2.26e-05        3.82 0.000136
      cTau      8.70e+00       3.26e+00        2.67 0.007621
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdOmega    5.01e-03       3.85e-03        1.30   0.193
      tdPhi      8.95e-01       1.01e-01        8.86   0.000
      tdSigma    7.69e-06       6.58e-06        1.17   0.243
    
    cubs
              Coefficient Standard Error t-statistic p-value
      cuC0       2.717101       0.295610        9.19 0.00000
      cuConst    0.002106       0.009158        0.23 0.81810
      cuSigma    0.000153       0.000056        2.74 0.00614
      RMSE: 0.0289
      R2: 0.367
      Box-Ljung test: X-squared = 2.17, df = 7.2, p-value = 0.956
    
             loglik             AIC             BIC             HQC signal-to-noise 
          192.86826      -367.73652      -353.48485      -362.76230         0.08896 

# TFP bayesian fit

    Call:
    fit.TFPmodel(model = model, method = "bayesian", R = 1000, thin = 2, 
        MLEfit = f)
    
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
      cA     0.416761 0.409588    1.64e-01    6.56e-01
      cSigma 0.000127 0.000101    4.93e-05    1.79e-04
      cTau   8.078838 7.588844    3.09e+00    1.24e+01
    
    trend
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      tdOmega 7.81e-03 7.33e-03    3.61e-04    1.42e-02
      tdPhi   9.15e-01 9.37e-01    8.44e-01    9.90e-01
      tdSigma 9.28e-06 7.16e-06    2.02e-06    1.56e-05
    
    cubs
                   Mean    Median 85% HPDI-LB 85% HPDI-UB
      cuC0     1.40e+00  1.40e+00    1.376785    1.422306
      cuConst -1.75e-05 -8.19e-06   -0.000946    0.001007
      cuSigma  5.44e-04  5.05e-04    0.000316    0.000787
              MRMSE signal-to-noise 
            0.02290         0.07336 

