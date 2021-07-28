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
      cPhi1     1.512436       1.56e-01       9.671 0.00e+00
      cPhi2    -0.571865       1.93e+00      -0.296 7.67e-01
      cSigma    0.000139       1.71e-05       8.160 4.44e-16
    
    trend
              Coefficient Standard Error t-statistic  p-value
      tSigma     5.63e-05       5.08e-05        1.11 2.67e-01
      tdConst    1.99e-02       2.93e-03        6.80 1.03e-11
    
    inflation equation
                    Coefficient Standard Error t-statistic  p-value
      inflC1             -0.105         0.7723      -0.136 0.891552
      inflConst          -0.144         0.0986      -1.459 0.144622
      inflErrGamma1      -0.220         0.2431      -0.907 0.364592
      inflErrGamma2      -0.292         0.1444      -2.023 0.043054
      inflErrGamma3      -0.488         0.2316      -2.107 0.035143
      inflGdpGL1          7.179         4.6153       1.556 0.119825
      inflSigma           0.378         0.1028       3.677 0.000236
      RMSE: 0.637
      R2: 0.385
      Box-Ljung test: X-squared = 2.55, df = 6, p-value = 0.862
    
             loglik             AIC             BIC             HQC signal-to-noise 
            50.8338        -77.6676        -60.8532        -72.2885          0.4045 

