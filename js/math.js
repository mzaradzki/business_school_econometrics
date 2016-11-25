
var gaussian = function() {
    // Box Muller method
    var u = 1 - Math.random(); // Subtraction to flip [0, 1) to (0, 1].
    var v = 1 - Math.random();
    var g = Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
    return g;
}

var DurbinWatson = function(vector) {
    var dwNum = 0;
    for (var i=1; i<vector.length; i++) {
        var d = (vector[i]-vector[i-1]);
        dwNum+= d*d;
    }
    return dwNum / numeric.norm2Squared(vector);
}

/*
// DEPRECATED : based on simple_statistics module
function multilinearRegressionLine(mb) {
    // Return a function that computes a `y` value for each
    // x value it is given, based on the values of `b` and `a`
    // that we just computed.
    return function(x) {
        return mb.b + (mb.m * x);
    };
}
*/
function multilinearRegressionLine(beta, addConstantToX) {
    return function(x) {
        if (addConstantToX) {
            return numeric.dot(beta, [1].concat(x));
        }
        else
        {
            return numeric.dot(beta, x);
        }
    };
}
var MultiRegression = function(Xs, Ys, addConstantToX) {
    // A * beta = V
    var fullXs = numeric.clone(Xs);
    if (addConstantToX) {
        fullXs = [];
        for (var r=0; r<Xs.length; r++) {
            var newrow = [1];
            fullXs.push( newrow.concat(Xs[r]) );
        }
    }
    var A = numeric.dot(numeric.transpose(fullXs), fullXs);
    var V = numeric.dot(numeric.transpose(fullXs), Ys);
    var Beta = numeric.solve(A, V);
    var predictions = numeric.dot(fullXs, Beta);
    var predictionsMean = numeric.sum(predictions)/  predictions.length;
    var SSExp = numeric.norm2Squared(numeric.sub(predictions, predictionsMean));
    var YsMean = numeric.sum(Ys) / predictions.length;
    var SSTot = numeric.norm2Squared(numeric.sub(Ys, YsMean));
    var residuals = numeric.sub(Ys, predictions);
    var SSRes = numeric.norm2Squared(residuals);
    var dof = (Ys.length-Beta.length);
    var s = Math.sqrt(SSRes / dof);
    var BetaCovar = numeric.mul(s*s, numeric.inv(A));
    var R2 = 1-SSRes/SSTot;
    var DW = DurbinWatson(residuals);
    var F = ((SSTot-SSRes)/(Beta.length-1)) / (SSRes/(Ys.length-Beta.length));
    //var regFun = multilinearRegressionLine(Beta, addConstantToX);
    return {
            addConstantToX:addConstantToX,
            beta:Beta,
            dof:dof,
            regressors: A,
            residuals:residuals,
            SSRes:SSRes, SSTot:SSTot, SSExp:SSExp,
            R2:R2,
            F:F, F_df1:Beta.length-1, F_df2:Ys.length-Beta.length, F_pval:1-jStat.centralF.cdf(F, Beta.length-1, Ys.length-Beta.length),
            stdErr:s, betaCov:BetaCovar,
            DW:DW,
        };
}
var BreuschGodfrey = function(MR, lags) {
    console.log('WARNING : needs code review and testing')
    var Xs = [];
    var Ys = [];
    for (var o=lags; o<MR.residuals.length; o++) {
        var row = [];
        for (var l=1; l<=lags; l++) {
            row.push( residuals[o-l] );
        }
        Ys.push( residuals[o] );
        Xs.push( row.concat(MR.regressors[o]) );
    }
    var auxMR = MultiRegression(Xs, Ys, false);
    var LM = (Xs.length-lags) * auxMR.R2;
    var pval = 1-jStat.chisquare.cdf(LM, lags);
    return {LM:LM, pval:pval, dof:lags};
}
var testMultiRegression = function(addConstantToX) {
    var Xs = [[1,0], [1,0.5],[2,0],[-1,0],[1,-1],[1,2]];
    var Ys = numeric.dot(Xs, [2,3]);
    for (var yi=0; yi<Ys.length; yi++) {
        Ys[yi] = Ys[yi] + gaussian()*0.5;
    }
    return MultiRegression(Xs, Ys, addConstantToX);
}

/*
// DEPRECATED : based on simple_statistics module
var simAndReg = function(_Xlist, _intercept, _beta, _errstd) {
    var M = [];
    for (var i=0; i<_Xlist.length; i++) {
        var y = 1*_intercept + _beta*_Xlist[i] + _errstd*gaussian(); // WARNNIG : *1 for casting purpose
        M.push([_Xlist[i], y]);
    }
    //console.log(M[0]);
    var regcoeff = ss.linearRegression(M);
    //console.log(regcoeff);
    var regfun = ss.linearRegressionLine(regcoeff);
    return {XY:M, regcoeff:regcoeff, regfun:regfun};
}
*/
var simAndReg = function(_Xlist, _intercept, _beta, _errstd) {
    var M = [];
    var Ys = [];
    for (var i=0; i<_Xlist.length; i++) {
        var y = 1*_intercept + _beta*_Xlist[i] + _errstd*gaussian(); // WARNNIG : *1 for casting purpose
        Ys.push(y);
        M.push([_Xlist[i], y]);
    }
    //
    var mr = MultiRegression(numeric.transpose([_Xlist]), Ys, true);
    var regfun = multilinearRegressionLine(mr.beta, mr.addConstantToX);
    return {XY:M, Y:Ys, MR:mr, regfun:regfun};
}

