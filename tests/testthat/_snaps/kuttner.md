# Kuttner model

    Call:
    KuttnerModel(tsl = tsList)
    
    	State space model object of class KuttnerModel
    
    cycle 				AR2
    trend 				RW1
    inflation equation
      cycle lags 			1
      error term			ARMA(0,3)
      exogenous variables		gdpGL1
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	30
      period 			1991 - 2020
      frequency 			annual
    
    Object is a valid object of class KuttnerModel.

# Kuttner MLE fit

    Call:
    fitKuttner(model = model, parRestr = parRestr)
    
    	State space model object of class KuttnerModel
    
    cycle 				AR2
    trend 				RW1
    inflation equation
      cycle lags 			1
      error term			ARMA(0,3)
      exogenous variables		gdpGL1
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	30
      period 			1991 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic  p-value
      cPhi1     1.511900       1.56e-01       9.686 0.00e+00
      cPhi2    -0.571461       1.92e+00      -0.297 7.66e-01
      cSigma    0.000139       1.93e-05       7.195 6.23e-13
    
    trend
              Coefficient Standard Error t-statistic  p-value
      tSigma     5.75e-05       5.11e-05        1.12 2.61e-01
      tdConst    1.99e-02       2.92e-03        6.83 8.55e-12
    
    inflation equation
                    Coefficient Standard Error t-statistic  p-value
      inflC1             -0.107         0.7716      -0.138 0.890049
      inflConst          -0.134         0.0993      -1.351 0.176731
      inflErrGamma1      -0.219         0.2478      -0.884 0.376952
      inflErrGamma2      -0.292         0.1438      -2.033 0.042017
      inflErrGamma3      -0.490         0.2366      -2.071 0.038324
      inflGdpGL1          6.707         4.6484       1.443 0.149090
      inflSigma           0.378         0.1026       3.681 0.000232
      RMSE: 0.638
      R2: 0.385
      Box-Ljung test: X-squared = 2.53, df = 4, p-value = 0.64
    
             loglik             AIC             BIC             HQC signal-to-noise 
            50.8208        -77.6415        -60.8272        -72.2625          0.4146 

