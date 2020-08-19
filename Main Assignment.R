#################### Loading packages and preparing data ##################

## Set Working Directory
setwd('~/R/ASM main assessment')

## Load Libraries
library(data.table)
library(ggplot2)
library(MCMCpack)
library(mclust)

## Load Data
# df1 <- fread("winemag-data_first150k.csv")
# df1$tab = 1
df2 <- fread("winemag-data-130k-v2.csv")
# df2$tab = 2

## Add missing columns for binding
# df1$taster_name = NA
# df1$taster_twitter_handle = NA 
# df1$title = NA

## Bind 2 data sets into one
# df = rbind(df1,df2)
df <- df2

## Free up memory by removing initial data sets
# rm(df1,df2)
rm(df2)

## Summarise Data
summary(df)
head(df)

##################### EDA and Data Cleaning ####################

## Looking at the distinct countries
## Required country values (South Africa,Chile)
sort(unique(df$country))


## Varieties in South Africa
## Picking "Sauvignon Blanc". Ignoring "Sauvignon Blanc-Chenin Blanc", "Sauvignon Blanc-Semillon"
sort(unique(df[df$country == "South Africa"]$variety))

## Varieties in Chile
## Picking "Chardonnay" and ignoring "Chardonnay-Viognier"
sort(unique(df[df$country == "Chile"]$variety))


## Filter records for Sauvignon Blanc in Sout Africa and Chardonnay in Chile
df2 <- df[(df$country == "South Africa" & df$variety == "Sauvignon Blanc") | (df$country == "Chile" & df$variety == "Chardonnay"),c("country","points","price","variety")]
df2 <- df2[complete.cases(df2), ]

## Looking at distribution of point values by variety
ggplot(df2) + geom_jitter(aes(variety, points, shape = df2$variety)) + geom_boxplot(aes(variety, points, fill = variety, alpha = 0.1))


################### Difference in 2 Groups ######################

## t test to see if there is a significant difference in the means of points of the 2 varieties of wine
## The test has a p-value of less than 0.05 which means the hypothesis that there isn't a difference in the means is rejected and there exists a significant evidence of difference in the means of two wine varieties
t.test(points ~ variety, data=df2, var.equal = TRUE)

## Function definition for mean comparison
compare_2_gibbs <- function(y, ind, mu0 = 85, tau0 = 1/10, del0 = 10, gamma0 = 1/400, 
                            a0 = 1, b0 = 50, maxiter = 5000)
{
  y1 <- y[ind == sort(unique(ind))[1]] # Filtering for 1st category
  y2 <- y[ind == sort(unique(ind))[2]] # Filtering for 2nd category
  
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
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    
    ##update del
    gamman <-  tau0 + tau*(n1 + n2)
    deln <- ( del0 * tau0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}


fit <- compare_2_gibbs(df2$points, df2$variety)
acf(fit)

plot(as.mcmc(fit))
raftery.diag(as.mcmc(fit))

y1_sim <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

y_sim_diff = y1_sim - y2_sim
ggplot(data.frame(y_sim_diff)) + stat_bin(aes(y_sim_diff))

## Difference in ratings
print(paste0("The difference in the rating of Chardonnay and Sauvignon Blanc",toString(round(mean(y_sim_diff),3))))

## The confidence Intervals for the difference values are
quantile(y_sim_diff,c(0.025,0.5,0.975))

print(paste0("Probability of Chardonnay being better than Sauvignon Blanc is ",toString(mean(y1_sim > y2_sim))))
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0)


########## For wines more than 15$ price
df3 <- df2[df2$price == 15,]

t.test(points ~ variety, data=df3, var.equal = TRUE)

fit2 <- compare_2_gibbs(df3$points, df3$variety)

plot(as.mcmc(fit2))
raftery.diag(as.mcmc(fit2))

y1_sim <- rnorm(5000, fit2[, 1] + fit2[, 2], sd = 1/sqrt(fit2[, 3]))
y2_sim <- rnorm(5000, fit2[, 1] - fit2[, 2], sd = 1/sqrt(fit2[, 3]))

y_sim_diff = y1_sim - y2_sim
ggplot(data.frame(y_sim_diff)) + stat_bin(aes(y_sim_diff))

## Difference in ratings
print(paste0("The difference in the rating of Chardonnay and Sauvignon Blanc",toString(round(mean(y_sim_diff),3))))

## The confidence Intervals for the difference values are
quantile(y_sim_diff,c(0.025,0.5,0.975))

print(paste0("Probability of Chardonnay being better than Sauvignon Blanc is ",toString(mean(y1_sim > y2_sim))))
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0)


################## Wines by region ##################

df4 <- df[(df$country == "Italy") & df$price < 20 & df$region_1 != "",c("country","points","price","variety","region_1")]
df4_r <- data.table(table(df4$region_1))
df4 <- df4[(df4$region_1 %in% df4_r[(df4_r$N > 3),V1]),]
df4 <- df4[complete.cases(df4),]
# df4_m <- df4[,list(mean_points=mean(points)), by = region_1]

ggplot(df4) + geom_boxplot(aes(x = reorder(region_1, points, median), points, 
                               fill = reorder(region_1, points, median)), show.legend=FALSE)

ggplot(df4, aes(x = reorder(region_1, region_1, length))) + stat_count()
ggplot(df4, aes(points)) + stat_bin(bins = 30)

ggplot(data.frame(size = tapply(df4$points, df4$region_1, length), 
                  mean_score = tapply(df4$points, df4$region_1, mean)), aes(size, mean_score)) + 
  geom_point() + xlab("Region sample size") + ylab("Mean Points") + 
  ggtitle("Effect size versus sample size")

## Adding small random noise to avoid 0 variance
df4$points <- df4$points + rnorm(nrow(df4),1,1)/10000

compare_m_gibbs <- function(y, ind, mu0 = 85, tau0 = 1/10, 
                            a0 = 1, b0 = 50, alpha0 =1, beta0 = 50, maxiter = 500)
{
  
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  ###
  
  ### starting values
  m <- nlevels(ind)
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

df4$region_1 <- as.factor(df4$region_1)
df4 <- droplevels(df4)

regions <- levels(df4$region_1)

fit3 <- compare_m_gibbs(df4$points, df4$region_1)
apply(fit3$params, 2, mean)
apply(fit3$params, 2, sd)
mean(1/sqrt(fit3$params[, 3])) 


theta_hat <- apply(fit3$theta, 2, mean) ## get basic posterior summary
names(theta_hat) <- regions # 1:nlevels(df4$region_1) ## keep track of different schools
sort(theta_hat, decreasing = TRUE) ## which schools did best and worst?

theta_ci <- apply(fit3$theta, 2, quantile, prob = c(0.025, .975)) ## upper/lower bounds for thetas
df_error <- data.frame(lower = theta_ci[1, ], upper = theta_ci[2, ], mean = theta_hat, 
                       region = regions)
ggplot(df_error, aes(x = reorder(region, mean), mean)) + geom_errorbar(aes(ymin = lower, ymax = upper))

## reformat samples for ggplot
theta_df <- data.frame(samples = as.numeric(fit3$theta), 
                       region = rep(1:ncol(fit3$theta), each = nrow(fit3$theta))) 

ggplot(theta_df) + geom_boxplot(aes(x = reorder(region, samples, median), samples, 
                                    fill = reorder(region, samples, median)), show.legend=FALSE)

ggplot(data.frame(size = tapply(df4$points, df4$region_1, length), theta_hat = theta_hat), 
       aes(size, theta_hat)) + geom_point()

ggplot(data.frame(ybar = tapply(df4$points, df4$region_1, mean), theta_hat = theta_hat), 
       aes(ybar, theta_hat)) + geom_point()

df4_t <- data.table(theta_hat)
df4_t$region <- regions

## Regions in Italy with ratings more than the average rating
unique(df4_t[df4_t$theta_hat > mean(fit3$params[,1]),])