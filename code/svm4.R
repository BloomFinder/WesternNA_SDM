# Author: Ian Breckheimer
# Date (last update):  March 2018
# Version 1
# Licence GPL v3

predict2 <- function(...) {predict(...)[,2]}

#-------------
methodInfo <- list(name=c('svm4','SVM4','ksvm4'),
                   packages='kernlab',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(x='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(type='C-svc',kernel='rbfdot',epsilon=0.1,
                                      prob.model=TRUE,shrinking=TRUE,C=50),
                   fitFunction = 'ksvm',
                   settingRules = function(x='sdmVariables',f='fitSettings') {
                     if (x@distribution == 'multinomial') f[['type']] <- 'C-svc'
                     list(fitSettings=f)
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='probabilities'),
                   predictFunction='predict2',
                   #------ metadata (optional):
                   title='Two-class Support Vector Machines with probability output',
                   creator='Babak Naimi, modified by Ian Breckheimer',
                   authors=c('Alexandros Karatzoglou'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "LIBSVM: a library for Support Vector Machines",
                                          author = as.person("C. Chih-Chung [aut], L. Chih-Jen [aut]"),
                                          year='2015',
                                          journal = "http://www.csie.ntu.edu.tw/~cjlin/libsvm"
                   )
                   ),
                   description='Support Vector Machines are an excellent tool for classification, novelty detection, and regression'
)