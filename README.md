# business-school-econometrics

Multi-linear regressions in Javascript

Returns :
* coefficients estimates
* estimates covariance matrix
* t-test
* F-test
* chi2-test

Diagnostic functions for :
* Durbin-Watson and Breusch-Godfrey autocorrelation test
* Jarque Bera normality test

```javascript
// wrap data series into a matrix
var Xs = [];
for (var t=0; t<Ys.length; t++) {
	Xs.push( [F1s[t], F2s[t], F3s[t], F3s[t]] );
}

// perform regression
var intercept = true;
MR = MultiRegression(Xs, Ys, intercept);

// main regression output
MR.beta;
MR.stdErr;
MR.betaCov;
MR.R2, MR.SSRes, MR.SSTot, MR.SSExp;
MR.F, MR.F_df1, MR.F_df2, MR.F_pval;
MR.DW;

// run Breusch-Godfrey test
var lags = 3;
var BG = BreuschGodfrey(MR, lags);
BG.LM;
BG.pval;
```