
# This is the count table, normalized and only for the day 0 samples
table <- data.frame(read.table('day0_normcounts.txt',header=T,row.names=1))

# Compress features with PCA
pcares<-prcomp(table[,-13596],center=TRUE,scale=TRUE)

#for (i in (1:10)){
#  write.table(colnames(table)[pcares$rotation[,i] < -0.01], file = paste0("predictor_output/PC",i,"_-001.txt",collapse = NULL))
#}


#Naive Bayes

data <- data.frame(pcares$x[,1:10],table$Response)
colnames(data)[11] <- "Response"
data$Response <- factor(data$Response)
levels(data$Response) <- make.names(levels(factor(data$Response)))


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
# AUCROC = 0.661654

perfPR <- performance(rocpred, measure = "prec", x.measure = "rec")
plot(perfPR) 

predfit<-values>0.5
table(predfit,data$Response)

performance(rocpred, measure = "acc")@y.values[[1]]

fit<-NaiveBayes(Response ~ ., data=data, usekernel=TRUE)
plot(fit)

# With the granulocytes

cellcount <- read.table('../data/cellcounts.txt',header=T,row.names=1,sep="\t")
total <- merge(cellcount,table[,13596, drop=FALSE],by="row.names",all.x=FALSE)

data <- data.frame(pcares$x[,1:10],total$GRA0,table$Response)
colnames(data)[12] <- "Response"
data$Response <- factor(data$Response)
levels(data$Response) <- make.names(levels(factor(data$Response)))

#Naive Bayes
values <- matrix(0,34,1)

require('klaR')
for(i in 1:34){
  fit<-NaiveBayes(Response ~ total.GRA0, data=data[-i,], usekernel=FALSE, prior=c(20/34,14/34))
  #klaR::plot(fit)
  pred<-predict(fit, newdata=data[i,])
  values[i]<-pred$posterior[,2]
}
require("ROCR")
rocpred <- prediction(predictions = values, labels=data$Response)
perf <- performance(rocpred, measure = "tpr", x.measure = "fpr")
plot(perf) 
abline(0,1)
performance(rocpred,measure="auc")@y.values[[1]]
#AUCROC = 0.6954887

perfPR <- performance(rocpred, measure = "prec", x.measure = "rec")
plot(perfPR) 

predfit<-values>0.5
table(predfit,data$Response)

performance(rocpred, measure = "acc")@y.values[[1]]

