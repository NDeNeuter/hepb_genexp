setwd('/Users/pietermeysman/Data/Immunology/GOA/HepB gen exp/reads/nicolas_try/day3response')

res.alls <- t(read.csv('day3R-logfolds.csv', row.names = 1))
response <- factor(read.csv('response.csv',row.names=1)[,1])

pcares<-prcomp(res.alls,center=TRUE,scale=TRUE)
require('ggplot2')
ggplot(data.frame(pcares$x),aes(x=PC1,y=PC2,label=response)) + geom_text(aes(colour=response))

data <- data.frame(pcares$x[,1:2],response)
levels(data$response) <- make.names(levels(factor(data$response)))

data <- data.frame(res.alls,response)
levels(data$response) <- make.names(levels(factor(data$response)))


require('caret')
fitControl<-trainControl(method="cv",number=5,repeats=10,savePredictions="final",classProbs=TRUE,summaryFunction=twoClassSummary)
fit<-train(response ~ ., data=data, method="glmnet",trControl=fitControl, metric='ROC')
require("ROCR")
rocpred <- prediction(predictions = fit$pred$X1, labels=fit$pred$obs)
perf <- performance(rocpred, measure = "tpr", x.measure = "fpr")
plot(perf) 
abline(0,1)
