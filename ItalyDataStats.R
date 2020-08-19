DataWineItaly = read.csv("WineData.csv")


DataWineItaly1 <- DataWineItaly[ -c(1,3:4,7,9:14) ]

Italy_Country <- DataWineItaly1[DataWineItaly1$country == "Italy",]
Italy_Country$price <- ifelse(is.na(Italy_Country$price), mean(Italy_Country$price, na.rm=TRUE), 
                              Italy_Country$price)
write.csv(Italy_Country, "C:/Users/LENOVO/Desktop/test/ItalyData.csv", row.names = TRUE)
Italy_Country <- Italy_Country[Italy_Country$price < 20,]
colnames(Italy_Country)
names(Italy_Country)[names(Italy_Country) == "region_1"] <- "Region"

Italy_Country$Region[Italy_Country$Region == ""] <- NA
write.csv(Italy_Country, "C:/Users/LENOVO/Desktop/test/ItalyData2.csv", row.names = TRUE)
ItalyFinal <- Italy_Country %>% group_by(Region) %>% filter(n() > 4)

write.csv(ItalyFinal, "C:/Users/LENOVO/Desktop/test/Italyfinal.csv", row.names = TRUE)
AverageWinePoints <- mean(ItalyFinal$points)

ItalyFinal$region_numeric <- c(ItalyFinal$Region)
head(ItalyFinal)
ItalyFinal[!is.na(ItalyFinal$Region), ]
write.csv(ItalyFinal, "C:/Users/LENOVO/Desktop/test/Italyfinal123.csv", row.names = TRUE)
Deleted_NA = read.csv("Italyfinal123.csv")



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
  m <- 1333
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
Result <- compare_m_gibbs(Deleted_NA$points, Deleted_NA$region_numeric)
write.csv(Result, "C:/Users/LENOVO/Desktop/test/Result.csv", row.names = TRUE)


head(Result)


