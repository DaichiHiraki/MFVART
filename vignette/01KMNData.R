library(fredr)
library(MFVART)
setwd("/home/hoanguc3m/Dropbox/WP16/Code/MFVART/Data/")
DatM_name <- "FredMD2024-10.csv"
DatQ_name <- "FredQD2024-10.csv"
headers <- read.csv(DatM_name, header = F, nrows = 1, as.is = T)
Monthly <- read.csv(DatM_name, skip = 2, header = F)
colnames(Monthly) <- headers

headers <- read.csv(DatQ_name, header = F, nrows = 1, as.is = T)
Quarterly <- read.csv(DatQ_name, skip = 3, header = F)
colnames(Quarterly) <- headers
# Remove Transform:,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2,2,2,5,5,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1,2,1,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,2,6,6,5,6,6,7,6,6,6,2,5,2,5,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,2,6,6,6,1

lag <- 12
t_max <- nrow(Monthly)
Time <- as.Date(Monthly$sasdate, format = "%m/%d/%Y")[(lag+1):t_max]

Var1 <- diff(log(Monthly$INDPRO), lag = lag) * 100
Var2 <- diff(log(Monthly$CPIAUCSL), lag = lag) * 100
#plot(Var2, type ="l")
Var3 <- (Monthly$UNRATE)[(lag+1):t_max]
plot(Var3, type ="l")

Var4 <- log(Monthly$VIXCLSx)[(lag+1):t_max]
plot(Var3, Var4)
Monthly <- data.frame(Time = Time, 
                      IP = Var1, 
                      Inflation = Var2, 
                      UNRATE = Var3,
                      VIX = Var4)

GDP <- Quarterly$GDPC1
GDP_growth <- 100*diff(log(GDP),lag = 1) # Using 1 lag instead of 4 lags.
Time_q <- as.Date(Quarterly$sasdate, format = "%m/%d/%Y")[2:nrow(Quarterly)]

Quarterly <- data.frame(Time = Time_q, GDP_growth = GDP_growth)
Data <- merge(x = Monthly, y = Quarterly, by = "Time", all.x = T)

adj_lag <- 43-lag
y_all <- as.matrix(Data[adj_lag:nrow(Data),2:6])
y_all[1,5] <- Data[adj_lag-4,6]
y_all[2,5] <- Data[adj_lag-1,6]

y0 <- y_all[1:2,]
y <- y_all[3:nrow(y_all),]
p = 2
K = 5
t_max <- nrow(y)
vars::VARselect( na.omit(y))

# load("/home/hoanguc3m/MEGA/WP16/RData/03KMN.RData")
# head(y)


prior <- get_prior(y, p = p, dist="Gaussian", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain1 <- BMFVAR.novol(y, K = K, p = p, dist = "Gaussian", 
                       y0 = y0, prior = prior, inits = inits)
B <- matrix(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:30], nrow = 5)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[31:40]
A <- t(A)
A_inv <- solve(A)
A_inv %*% t(A_inv)

prior <- get_prior(y, p = p, dist="Student", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain2 <- BMFVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)

plot(apply(Chain1$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:45])
abline(a = 0, b = 1)
apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[46]

prior <- get_prior(y, p = p, dist="OT", SV = F, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain3 <- BMFVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
plot(apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain2$mcmc$param, MARGIN = 2, FUN = median)[1:45])
abline(a = 0, b = 1)
apply(Chain3$mcmc$param, MARGIN = 2, FUN = median)[46:50]

#############################################
prior <- get_prior(y, p = p, dist="cG", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain4 <- BMFVAR.SV(y, K = K, p = p, dist = "cG", y0 = y0, prior = prior, inits = inits)
# 
# prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(5))
# inits <- get_init(prior)
# # prior$ymiss_prior <- 1
# inits <- get_init(prior, samples = 3000, burnin = 00, thin = 1)
# # prior$ymiss_prior <- 4

# Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
# plot(get_post(Chain5$mcmc$param, "B"))
# 
# mcmc_tmp <- matrix(NA, nrow = nrow(Chain5$mcmc$y_miss), 
#                    ncol = ncol(Chain5$mcmc$y_miss))
# 
# for (i in c(1:nrow(Chain5$mcmc$y_miss))){
#   # Sys.sleep(0.05)
#   # plot(as.numeric(Make_bandSparse(t_max, aggregation="triangular") %*%
#   #              Chain5$mcmc$y_miss[i,])[600:t_max], type = "l", main = i)
#   mcmc_tmp[i,] <- as.numeric(Make_bandSparse(t_max, aggregation="triangular") %*% 
#                                Chain5$mcmc$y_miss[i,])
# }
# 
# plot(Chain5$mcmc$y_miss[,690])
# plot(Chain5$mcmc$y_miss[,691])
# plot(Chain5$mcmc$y_miss[,692])
# plot(Chain5$mcmc$y_miss[,693])
# plot(Chain5$mcmc$y_miss[,694])
# 
# 
# plot(mcmc_tmp[,690])
# plot(mcmc_tmp[,691]) # Obs
# plot(mcmc_tmp[,692])
# plot(mcmc_tmp[,693])
# plot(mcmc_tmp[,694]) # Obs
# plot(mcmc_tmp[,695])
# plot(mcmc_tmp[,696])
# plot(mcmc_tmp[,697]) # Obs

############################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
# prior$ymiss_prior <- 4
Chain5 <- BMFVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)

prior <- get_prior(y, p = p, dist="Student", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain6 <- BMFVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
median(Chain6$mcmc$param[,41])

apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[1:46]
plot(apply(Chain5$mcmc$param, MARGIN = 2, FUN = median)[1:45],
     apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,42:46)])
abline(a = 0, b = 1)

prior <- get_prior(y, p = p, dist="OT", SV = T, aggregation = "triangular", idq = c(5))
inits <- get_init(prior)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain7 <- BMFVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[41:45]
# nu1       nu2       nu3       nu4       nu5 
# 14.367142 19.200400  8.738173  5.907335 19.278837 

B <- matrix(apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[1:30], nrow = 5)
A = diag(1, K, K)
A[upper.tri(A, diag=FALSE)] <- apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[31:40]
A <- t(A)
A_inv <- solve(A)
A_inv %*% t(A_inv)

plot(apply(Chain6$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,42:46)],
     apply(Chain7$mcmc$param, MARGIN = 2, FUN = median)[c(1:40,46:50)])
abline(a = 0, b = 1)

# save.image("/home/hoanguc3m/MEGA/WP16/RData/03KMN.RData")
############################################

prior <- get_prior(y, p = p, dist="Gaussian", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain8 <- BMFVAR.FSV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)



prior <- get_prior(y, p = p, dist="Student", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain9 <- BMFVAR.FSV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
plot(apply(get_post(Chain6$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain9$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain6$mcmc$y_miss), 2, mean),
     apply(Chain9$mcmc$y_miss, 2, mean))

prior <- get_prior(y, p = p, dist="OT", SV = T, 
                   aggregation = "triangular", idq = c(5), r = 2)
inits <- get_init(prior, samples = 60000, burnin = 10000, thin = 10)
Chain10 <- BMFVAR.FSV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)

plot(apply(get_post(Chain7$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain10$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain7$mcmc$y_miss), 2, mean),
     apply(Chain10$mcmc$y_miss, 2, mean))


plot(apply(get_post(Chain5$mcmc$param, element = "B"), 2, mean),
     apply(get_post(Chain8$mcmc$param, element = "B"), 2, mean))
abline(a = 0, b = 1)
plot(apply((Chain5$mcmc$y_miss), 2, mean),
     apply(Chain8$mcmc$y_miss, 2, mean))

plot(apply((Chain8$mcmc$y_miss), 2, mean),
     apply(Chain6$mcmc$y_miss, 2, mean))


plot(Chain1)
plot(Chain5)
plot(Chain8)

plot(Chain2)
plot(Chain6)
plot(Chain9)

plot(Chain3)
plot(Chain7)
plot(Chain10)

c(Chain5$esttime, Chain6$esttime, Chain7$esttime)
c(Chain8$esttime, Chain9$esttime, Chain10$esttime)

colMeans(get_post(Chain7$mcmc$param, element = "nu"))
colMeans(get_post(Chain10$mcmc$param, element = "nu"))

colMeans(get_post(Chain10$mcmc$param, element = "nu"))

save.image("/home/hoanguc3m/MEGA/WP16/RData/03KMN.RData")
