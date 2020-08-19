library(dplyr)
library(ggplot2)
library(tidyverse)
library(imputeTS)
library(GGally)
library(Boruta)

wine_full <- read.csv("F:/Trinity/asm/winemag-data-130k-v2.csv/winemag-data-130k-v2.csv")
head(wine_full)


us_wines <- subset(df, country =='US' )
us_wines <- droplevels(us_wines)
summary(us_wines)
dim(us_wines)
##us_wines <-  wine_full%>%
##  filter(country == 'US')

## missing data imputation
cat('Rows and Columns -> ')
dim(us_wines)
cat('Complete cases -> ', sum(complete.cases(us_wines)))
cat('Total null values -> ' , sum(is.na(us_wines)))
list_na <- colnames(us_wines)[apply(us_wines, 2, anyNA)]
cat('Columns with null values -> ')
list_na

nulls <- sapply(us_wines, function(x) sum(is.na(x)))
cat('null count for each column -> ')
nulls
cat('Using imputeTS to fill NA values with mean ... ')

mean_fill <- us_wines%>%
                na_mean()
cat('Missing values after filling with mean -> ', sum(is.na(mean_fill)))
us_wines <- mean_fill 

## visualizing the distribution of points
ggplot(data = us_wines, aes(x= points, colour = I('black'), fill = I('#099DD9')))+
  geom_histogram(binwidth = 1)+
  labs(x = "Points", y= "Frequency", title = "Distribution of points")


## sampling dataset
ntotal <- nrow(us_wines)
nsamp <- 500 ## subset to sample
ind_samp <- sample(1:ntotal, nsamp, replace = FALSE) ## which datapoints to sample
wine_samp <- us_wines[ind_samp, ] ## which datapoints to sample
wine_samp <- droplevels(wine_samp)

wine_samp <- within(wine_samp, rm("X","country","taster_name","taster_twitter_handle"))

# Perform Boruta search
boruta_output <- Boruta(points ~ ., data=wine_samp, doTrace=1)  

boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)

# Variable Importance Scores
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
head(imps2[order(-imps2$meanImp), ],10)  # descending sort

# Plot variable importance
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")


## model 1 based on the features selected
model <- lm(points~price+province+region_2+winery+designation, data = wine_samp)
summary(model)

plot(model, which = 2)

yhat <- predict(model)
plot(yhat, wine_samp$points)