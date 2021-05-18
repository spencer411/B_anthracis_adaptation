library(randomForest)
library(spm)

Data = read.csv("Data.csv")
set.seed(12345)
Train = createDataPartition(Data$Clade1, p=0.75, times=1)
Train_Data = Data[Train$Resample1,]
Valid_Data = Data[-Train$Resample1,]

a=c()
a=matrix(NA,ncol=4,nrow=50*20)
for (n in seq(100,5000,by = 100)) {
  for (m in 1:20) {
    set.seed(12345)
    rf.cv = RFcv(Train_Data[,-1], Train_Data[,1], mtry=m, ntree=n , cv.fold=5)
    a[m+((n/100)-1)*20,1] = n
    a[m+((n/100)-1)*20,2] = m
    a[m+((n/100)-1)*20,3] = rf.cv$kappa
    a[m+((n/100)-1)*20,4] = rf.cv$ccr
    print(n);print(m)
  }
}

Results = as.data.frame(a)
colnames(Results)=c("ntree","mtry","Kappa","ccr")
view(Results)

set.seed(12345)
RF1 = randomForest(Clade1 ~ ., data=Train_Data, importance=T, ntree=400, mtry=9)
RF1

predV = predict(RF1, Valid_Data, type="class")
table(predV, Valid_Data$Clade1)
mean(predV == Valid_Data$Clade1)

Imp = importance(RF1)