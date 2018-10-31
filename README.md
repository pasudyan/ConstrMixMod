# ConstrMixMod
R package for [Flexible Mixture Modeling on Constrained Spaces](https://arxiv.org/abs/1809.09238)

The functions here are used to generate figures in the paper. The main function to draw inference with Gibbs sampling is: 
```r
infConstrMixMod(mix, Xdata, constr, K, thr, constrType, niwpr, ...)
```
The package also provides estimation of posterior predictive distribution for a given test set through:
```r
postPredDist(mix, constrType, res, test, pts, realDens, ...)
postPredDistMv(mix, constrType, res, test, constr, thinsc, ...)
```

