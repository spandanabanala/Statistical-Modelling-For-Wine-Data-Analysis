#Read the Wine data from two files by pasting in the project folder
DataWine = read.csv("WineV2.csv")
Data1 = read.csv("WineV1.csv")

#Removing the unwanted columns from the data and keeping only 4 columns Country,Points,price and Variety
DataWine1 <- DataWine[ -c(1,3:4) ]
DataWine2 <- DataWine1[ -c(4:9) ]
DataWine3 <- DataWine2[-c(5)]
Data2 <- Data1[ -c(1,3:4) ]
Data3 <- Data2[ -c(4:6) ]
Data4 <- Data3[-c(5)]

#Appending the data from two files and storing in Datawine3
DataWine3 <- rbind(Data4, DataWine3)


#Removing the null rows
DataWine3=DataWine3[rowSums(is.na(DataWine3)) == 0,]

#Filtering only the data with price equals to 15
DataWine3 <- DataWine3[DataWine3$price == 15,]


#Filtering only South Africa acountry
SA_Country <- DataWine3[DataWine3$country == "South Africa",]
#Filtering only Chile country
Ch_Country <- DataWine3[DataWine3$country == "Chile",]
#Filtering only Sauvigon Blanc from South Africa
SouthAfricaSB <- SA_Country[SA_Country$variety == "Sauvignon Blanc",]
#Filtering only Chardonnay from Chile
ChileCh <- Ch_Country[Ch_Country$variety == "Chardonnay",]


#Appending the filtered wine data..Southafrica is placed in top of Chile
FinalData <- rbind(SouthAfricaSB, ChileCh)


#Copying the country in separate column so as to get number for country factor value
FinalData$Country1 <- c(FinalData$country)



#1.a.i

#got this code from Hierarchical models case study 
library(ggplot2)
ggplot(FinalData) + geom_boxplot(aes(x = reorder(Country1, points, median), points,
                                fill = reorder(Country1, points, median)), show.legend=FALSE)


ggplot(FinalData, aes(x = reorder(Country1, Country1, length))) + stat_count()
ggplot(FinalData, aes(points)) + stat_bin()

ggplot(data.frame(size = tapply(FinalData$points, FinalData$Country1, length),
                  mean_score = tapply(FinalData$points, FinalData$Country1, mean)), aes(size, mean_score)) +
  geom_point() + xlab("Wine sample size") + ylab("Mean Score") +
  ggtitle("Effect size versus sample size")


sapply(FinalData, class)
is.na(FinalData[1, 1])

class(FinalData$Country1)
Test1 <- sort(table(FinalData$points), decreasing = TRUE)
head(FinalData)


#Copying the compare_m_gibbs function from Hierarchical models case study 
compare_m_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400,
                            a0 = 1, b0 = 50, alpha0 =1, beta0 = 50, maxiter = 5000)
{
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  ###
  ### starting values
  #m <- nlevels(ind)
  ### Changed this line to 2 as only ywo countries are selected
  m <- 2
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  alphan <- alpha0 + sum(n_m)/2
  ###
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  ### MCMC algorithm
  for(s in 1:maxiter)
  {
    # sample new values of the thetas
    for(j in 1:m)
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    betan <- beta0 + ss/2
    tau_w <- rgamma(1, alphan, betan)
    #sample a new value of mu
    taum <- m * tau_b + tau0
    mum <- (mean(theta) * m * tau_b + mu0 * tau0) / taum
    mu <- rnorm(1, mum, 1/ sqrt(taum))
    # sample a new value of tau_b
    am <- a0 + m/2
    bm <- b0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, am, bm)
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

#Running the gibbs model for our filtered data
Result <- compare_m_gibbs(FinalData$points, FinalData$Country1)
head(Result)
#write.csv(Result, "C:/Users/LENOVO/Desktop/test/FinalData3.csv", row.names = TRUE)

#Applying the function as given in case study
apply(Result$params, 2, mean) #mean of columns fit params

apply(Result$params, 2, sd) #SD of columns 


## within school standard variation
mean(1/sqrt(Result$params[, 2]))

sd(1/sqrt(Result$params[, 2]))

## between school standard variation
mean(1/sqrt(Result$params[, 3]))

sd(1/sqrt(Result$params[, 3]))



theta_hat <- apply(Result$theta, 2, mean) ## get basic posterior summary
theta_hat
names(theta_hat) <- 1:2 ## keep track of different country wines
# 1 refers to Souverign blanc and 2 refers to Chardonny
tc <- sort(theta_hat, decreasing = TRUE) ## which wine did best and worst?
#this gives result for which wine performs better
tc




#1.a.ii

#Copying the compare_2_gibbs from the Hierarchical methods case study
compare_2_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400, del0 = 0, gamma0 = 1/400,
                            a0 = 1, b0 = 50, maxiter = 5000)
{
  y1 <- y[ind == 1]
  y2 <- y[ind == 2]
  n1 <- length(y1)
  n2 <- length(y2)
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  for(s in 1 : maxiter)
  {
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    ##update mu
    taun <- tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    ##update del
    gamman <- tau0 + tau*(n1 + n2)
    deln <- ( del0 * tau0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

#need this packckage..runn it only once
install.packages("MCMCpack")
library(MCMCpack)

#Running the gibbs_2_compare function to get the probabilty result
Result1 <- compare_2_gibbs(FinalData$points, FinalData$Country1)
plot(as.mcmc(Result1))

raftery.diag(as.mcmc(Result1))

apply(Result1, 2, mean)

apply(Result1, 2, sd)

mean(1/sqrt(Result1[, 3]))
sd(1/sqrt(Result1[, 3]))




y1_sim <- rnorm(5000, Result1[, 1] + Result1[, 2], sd = 1/sqrt(Result1[, 3]))
y2_sim <- rnorm(5000, Result1[, 1] - Result1[, 2], sd = 1/sqrt(Result1[, 3]))
ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))

#this gives us the probabilty that sauvorign blanc wine has more probabilty of being better
mean(y1_sim > y2_sim)
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0)




