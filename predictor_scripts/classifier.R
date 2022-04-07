
# This is the count table, normalized and only for the day 0 samples
table <- data.frame(read.table('day0_normcounts.txt',header=T,row.names=1))
table <- table[,]

# Compress features with PCA
pcares<-prcomp(table[,-13596],center=TRUE,scale=TRUE)

# Write PCA loadings to file
#for (i in (1:10)){
#  write.table(colnames(table)[pcares$rotation[,i] < -0.01], file = paste0("predictor_output/PC",i,"_-001.txt",collapse = NULL))
#}


# Data set for Naive Bayes classifier

data <- data.frame(pcares$x[,1:5])
data["Response"] = as.factor(table$Response)

#Naive Bayes

values <- matrix(0,34,1)

require('klaR')

# Leave-one-out cross validation

for(i in 1:34){
  fit<-NaiveBayes(Response ~ ., data=data[-i,], usekernel=TRUE, prior=c(20/34,14/34))
  #klaR::plot(fit)
  pred<-predict(fit, newdata=data[i,])
  values[i]<-pred$posterior[,2]
}
# ROC curve

require("ROCR")
rocpred <- prediction(predictions = values, labels=data$Response)
perf <- performance(rocpred, measure = "tpr", x.measure = "fpr")
plot(perf) 
abline(0,1)
performance(rocpred,measure="auc")@y.values[[1]]
# AUCROC = 0.6857143 (in R 4.1.1)

# PR curve
perfPR <- performance(rocpred, measure = "prec", x.measure = "rec")
plot(perfPR) 

# Confusion matrix
predfit<-values>0.5
table(predfit,data$Response)

# Accuracy
performance(rocpred, measure = "acc")@y.values[[1]]

fit<-NaiveBayes(Response ~ ., data=data, usekernel=TRUE)
#plot(fit)

# With the granulocytes

cellcount <- read.table('../data/cellcounts.txt',header=T,row.names=1,sep="\t")
total <- merge(cellcount,table[,13596, drop=FALSE],by="row.names",all.x=FALSE)

data.gr <- data.frame(pcares$x[,1:5],total$GRA0)
data.gr["Response"] = as.factor(table$Response)

#Naive Bayes
values.gr <- matrix(0,34,1)

require('klaR')
for(i in 1:34){
  fit<-NaiveBayes(Response ~ ., data=data.gr[-i,], usekernel=TRUE, prior=c(20/34,14/34))
  #klaR::plot(fit)
  pred<-predict(fit, newdata=data.gr[i,])
  values.gr[i]<-pred$posterior[,2]
}
require("ROCR")
rocpred.gr <- prediction(predictions = values.gr, labels=data.gr$Response)
perf.gr <- performance(rocpred.gr, measure = "tpr", x.measure = "fpr")
plot(perf.gr) 
abline(0,1)
performance(rocpred.gr,measure="auc")@y.values[[1]]
#AUCROC = 0.7178571 (in R 4.1.1)

perfPR.gr <- performance(rocpred.gr, measure = "prec", x.measure = "rec")
plot(perfPR.gr) 

predfit.gr<-values.gr>0.5
table(predfit.gr,data.gr$Response)

performance(rocpred.gr, measure = "acc")@y.values[[1]]

