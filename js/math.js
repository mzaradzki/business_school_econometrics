var gaussian = function() {
    // Box Muller method
    var u = 1 - Math.random(); // Subtraction to flip [0, 1) to (0, 1].
    var v = 1 - Math.random();
    var g = Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
    return g;
}
var MultiReg = function(Xs, Ys, addConstant) {
    // A * beta = V
    var fullXs = numeric.clone(Xs);
    if (addConstant) {
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
    var residuals = numeric.sub(Ys, predictions);
    var dof = (Ys.length-Beta.length);
    var SSR = numeric.norm2(residuals);
    var s = Math.sqrt(SSR / dof);
    var BetaCovar = numeric.mul(s*s, numeric.inv(A));
    return {beta:Beta, dof:dof, SSR:SSR, stdErr:s, betaCov:BetaCovar};
}
var testMultiReg = function(addConstant) {
    var Xs = [[1,0], [1,0.5],[2,0],[-1,0],[1,-1],[1,2]];
    var Ys = numeric.dot(Xs, [2,3]);
    for (var yi=0; yi<Ys.length; yi++) {
        Ys[yi] = Ys[yi] + gaussian()*0.5;
    }
    return MultiReg(Xs, Ys, addConstant);
}
var simAndReg = function(_Xlist, _intercept, _beta, _errstd) {
    var M = [];
    for (var i=0; i<_Xlist.length; i++) {
        var y = 1*_intercept + _beta*_Xlist[i] + _errstd*gaussian(); // WARNNIG : *1 for casting purpose
        M.push([_Xlist[i], y]);
    }
    //console.log(M[0]);
    var regcoeff = ss.linearRegression(M);
    var regfun = ss.linearRegressionLine(regcoeff);
    return {XY:M, regcoeff:regcoeff, regfun:regfun};
}