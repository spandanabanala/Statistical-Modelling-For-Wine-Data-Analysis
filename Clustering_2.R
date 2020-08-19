data1 = read.csv("WineV1.csv")
data2 = read.csv("WineV2.csv")

DataWine1 <- data1[ -c(1,3:4) ]
DataWine1 <- DataWine1[ -c(4:8) ]
head(DataWine1)
DataWine2 <- data2[ -c(1,3:4) ]
DataWine2 <- DataWine2[ -c(4:11) ]
head(DataWine2)

DataWine1$price <- ifelse(is.na(DataWine1$price), mean(DataWine1$price, na.rm=TRUE), 
                          DataWine1$price)
head(DataWine1)
DataWine2$price <- ifelse(is.na(DataWine2$price), mean(DataWine2$price, na.rm=TRUE), 
                          DataWine2$price)
head(DataWine2)
CombineData <- rbind(DataWine1, DataWine2)
head(CombineData)

US_Country <- CombineData[CombineData$country == "US",]

US_Country <- US_Country[ -c(1) ]
write.csv(US_Country, "C:/Users/LENOVO/Desktop/test/USfinal.csv", row.names = TRUE)
plot(US_Country)
#install.packages("mclust")
library(mclust)

US_Country[2] = log(US_Country[2],2)

write.csv(US_Country, "C:/Users/LENOVO/Desktop/test/ussss.csv", row.names = TRUE)


toy_data_3 <- apply(US_Country, 2, scale)

toy_data_3
fit <- Mclust(toy_data_3)

print(fit)

summary(fit)

fit$parameters$pro

fit$parameters$mean

plot(fit, what = "classification")

plot(fit, what = "uncertainty")

head(fit$classification)

head(fit$uncertainty)

plot(fit, what = "BIC")

fit$BIC

#fit2 <- Mclust(US_Country, G = 21, modelNames = "VEE")
#plot(fit2, what = "classification")


fit3 <- Mclust(toy_data_3, G = 7, modelNames = "VVV")
plot(fit3, what = "classification")
plot(fit3, what = "uncertainty")

table(fit$classification, fit3$classification)



