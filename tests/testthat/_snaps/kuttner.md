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
             Coefficient Standard Error t-statistic p-value
      cPhi1     1.443973       1.41e-01      10.219 0.00000
      cPhi2    -0.521264       1.23e+00      -0.424 0.67172
      cSigma    0.000198       5.81e-05       3.397 0.00068
    
    trend
              Coefficient Standard Error t-statistic  p-value
      tSigma     1.81e-05       9.55e-06        1.90 5.79e-02
      tdConst    1.98e-02       2.63e-03        7.55 4.29e-14
    
    inflation equation
                    Coefficient Standard Error t-statistic  p-value
      inflC1            -0.0973          0.803      -0.121 0.903533
      inflConst         -0.1234          0.100      -1.228 0.219464
      inflErrGamma1     -0.2143          0.253      -0.848 0.396673
      inflErrGamma2     -0.2959          0.143      -2.074 0.038065
      inflErrGamma3     -0.4938          0.244      -2.027 0.042707
      inflGdpGL1         6.1713          4.693       1.315 0.188498
      inflSigma          0.3749          0.102       3.687 0.000227
      RMSE: 0.638
      R2: 0.384
      Box-Ljung test: X-squared = 2.71, df = 6, p-value = 0.844
    
             loglik             AIC             BIC             HQC signal-to-noise 
           51.06918       -78.13836       -61.32400       -72.75930         0.09163 

