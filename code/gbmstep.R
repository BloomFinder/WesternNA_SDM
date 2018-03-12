# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  July 2017
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('gbmstep3'),
                   packages='dismo',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(family='bernoulli',
                                      gbm.y=1,
                                      gbm.x=2:18,
                                      tree.complexity=5,
                                      max.trees=5000,
                                      n.trees=50,
                                      n.folds = 10,
                                      learning.rate=0.005,
                                      bag.fraction = 0.75),
                   fitFunction = 'dismo::gbm.step',
                   settingRules = function(x='sdmVariables',f='fitSettings') {
                     #if (!is.null(userSetting)) fitSettings <- .assign(fitSettings,userSettings)
                     #else if (x@distribution == 'ab') fitSettings[['distribution']] <- "poisson"
                     #else if (x@distribution == 'n') fitSettings[['distribution']] <- "multinomial"
                     if (x@distribution %in% c('poisson','multinomial')) {
                       f[['family']] <- x@distribution
                     }
                     #if (x@number.of.records[1] < 55) {
                     #  f[['n.minobsinnode']] <- max(floor(x@number.of.records[1] / 5),3)
                     #  if (x@number.of.records[1] < 30) f[['bag.fraction']] <- 0.9
                     #}
                     list(fitSettings=f)
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(n.trees=6000,type='response'),
                   predictFunction='predict.gbm',
                   #------ metadata (optional):
                   title='Boosted Regression Trees',
                   creator='Babak Naimi',
                   authors=c('Greg Ridgeway'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "A working guide to boosted regression trees",
                                          author = as.person("J. Elith [aut], J. R. Leathwick [aut], T. Hastie [aut]"),
                                          year = "2008",
                                          journal = "Journal of Animal Ecology",
                                          number="77",
                                          pages="802-813",
                                          publisher="Wiley Online Library")
                   ),
                   description='Fits Boosting regression trees (BRT), called also generalized boosting regression model (GBM). Boosting is the process of iteratively adding basis functions in a greedy fashion so that each additional basis function further reduces the selected loss function [see the help for gbm function in gbm package]'
)
#------------