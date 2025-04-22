# eBH vs eBH+ simulations

library(mvtnorm)


# e-value z test

ztestE <- function(x, a) {
	evalue <- exp(sum(a*x-a^2/2))
	return(evalue)
}

# eBH

eBH <- function(evalues, alpha=0.05){
	
	l <- length(evalues)
	Sevalues <- sort(evalues, decreasing=T)
	sorted_indices <- order(evalues, decreasing=T)
	
	temp <- rep(NA, times=l)
	
	for (k in 1:l) {
		left <- Sevalues[k]
		right <- 1 / alpha * l / k
		temp[k] <- (left>=right)
	}
	
	if (sum(temp)==0) {
		R <- "emptyset"
		rejections <- 0
		indices <- 0
	} else if (sum(temp)>0) {
		R <- 1:max(which(temp==T))
		rejections <- length(R)
		indices <- sorted_indices[R]
	} else {print("faliure: something wrong with the rejected set")}
	
	
	out <- list()
	out$indices <- indices
	out$evalues <- evalues[indices]
	out$rejections <- rejections
	
	return(out)
	
}

# our method taking the max over subsets

eBHplus <- function(evalues, alpha=0.05) {
	
library(sets)
	
Sevalues <- sort(evalues, decreasing=T)
sorted_indices <- order(evalues, decreasing=T)


l <- length(Sevalues)
temp <- as.numeric(1:l)
powerset <- set_power(temp)

sumconsistence <- rep(NA, times=length(powerset))

for (i in 1:length(powerset)) {
	
	U <- unlist(as.vector(powerset)[[i]])
	
	left <- rep(NA, times=length(powerset))
	right <- rep(NA, times= length(powerset))
	
	for (j in 1:length(powerset)) {
	
		S <- unlist(as.vector(powerset)[[j]])
		
		if (length(S)>0) {
			left[j] <- 1 / length(S) * sum(Sevalues[S])
		} else if (length(S)==0) {
			left[j] <- 0
		} else {print("failure in length of S"); break}
		
		if (length(U)>0) {
			right[j] <- length(intersect(U, S)) / length(U) / alpha
		} else if (length(U)==0) {
			right[j] <- 0
		} else {print("failure in length of U"); break}
		
		}
	sumconsistence[i] <- sum(left>=right)==length(powerset)	
	
}

ind <- max(which(sumconsistence==T))
R <- unlist(as.vector(powerset)[[ind]])

if (is.null(R)) {
	R <- "emptyset" 
	rejections <- 0
	indices <- 0
	} else {
		rejections <- length(R)
		indices <- sorted_indices[R]
		}

out <- list()
out$indices <- indices
out$evalues <- evalues[indices]
out$rejections <- rejections

return(out)

}


# Positive dependence

pi0 <- 0.75
m <- 8
n0 <- round(m*pi0)

A <- 0.125

# means
mu <- c(rep(A, times=m-n0), rep(0, times=n0))

rho <- 0.8

# covariance matrix
# Note that BH does not control the FDR if rho < 0
Sigma <- matrix(nrow=length(mu), ncol=length(mu))
diag(Sigma) <- 1
for (i in 1:length(mu)) {
	for (j in 1:length(mu)) {
		if (j != i) {
		Sigma[i, j] <- rho^{abs(i-j)}
		}
	}
}

alpha <- 0.05

a <- A

nsim <- 1000

rej.maxSC <- rej.eBH <- FDP.maxSC <- FDP.eBH <- rep(NA, times=nsim)


for (k in 1:nsim) {
	
	Z <- rmvnorm(100, mean=mu, sigma=Sigma)
	
	evaluesH <- rep(NA, times=length(mu))
	for (j in 1:length(mu)) {
		evaluesH[j] <- ztestE(Z[,j], a = a)
		}
	
	out.maxSC <- eBHplus(evaluesH, alpha=alpha)
	rej.maxSC[k] <- out.maxSC$rejections

	fp.maxSC <- sum(out.maxSC$indices > m-n0)
	FDP.maxSC[k] <- fp.maxSC / max(1, out.maxSC$rejections)

	
	out.eBH <- eBH(evaluesH, alpha=alpha)
	rej.eBH[k] <- out.eBH$rejections

 	fp.eBH <- sum(out.eBH$indices > m-n0)
 	FDP.eBH[k] <- fp.eBH / max(1, out.eBH$rejections)
	

	
	print(k)
	
}

TotalS.FDP <- mean(FDP.maxSC)
TotaleBH.FDP <- mean(FDP.eBH)
PowerS <- mean(rej.maxSC) / (m-n0)
PowereBH <- mean(rej.eBH) / (m-n0)

x <- c(0.125, 0.25, 0.375, 0.5)


plot(x, PowerS, type="l", col="blue")
points(x, PowerS, pch=3, col="blue")
lines(x, PowereBH, col="red")
points(x, PowereBH, pch=2, col="red")



# Independence


pi0 <- 0.75
m <- 8
n0 <- round(m*pi0)

A <- 0.125

# means
mu <- c(rep(A, times=m-n0), rep(0, times=n0))

Sigma <- diag(1, 8, 8)

alpha <- 0.05

a <- A

nsim <- 1000

rej.maxSC <- rej.eBH <- FDP.maxSC <- FDP.eBH <- rep(NA, times=nsim)


for (k in 1:nsim) {
	
	Z <- rmvnorm(100, mean=mu, sigma=Sigma)
	
	evaluesH <- rep(NA, times=length(mu))
	for (j in 1:length(mu)) {
		evaluesH[j] <- ztestE(Z[,j], a = a)
		}
	
	out.maxSC <- eBHplus(evaluesH, alpha=alpha)
	rej.maxSC[k] <- out.maxSC$rejections

	fp.maxSC <- sum(out.maxSC$indices > m-n0)
	FDP.maxSC[k] <- fp.maxSC / max(1, out.maxSC$rejections)

	
	out.eBH <- eBH(evaluesH, alpha=alpha)
	rej.eBH[k] <- out.eBH$rejections

 	fp.eBH <- sum(out.eBH$indices > m-n0)
 	FDP.eBH[k] <- fp.eBH / max(1, out.eBH$rejections)
	

	
	print(k)
	
}

TotalS.FDP <- mean(FDP.maxSC)
TotaleBH.FDP <- mean(FDP.eBH)
PowerS <- mean(rej.maxSC) / (m-n0)
PowereBH <- mean(rej.eBH) / (m-n0)

TotalS.FDP <- c(mean(FDP.maxSC), TotalS.FDP)
TotaleBH.FDP <- c(mean(FDP.eBH), TotaleBH.FDP)
PowerS <- c(mean(rej.maxSC) / (m-n0), PowerS)
PowereBH <- c(mean(rej.eBH) / (m-n0), PowereBH)

x <- c(0.125, 0.25, 0.375, 0.5)

# Negative dependence

pi0 <- 0.75
m <- 8
n0 <- round(m*pi0)

A <- 0.25

# means
mu <- c(rep(A, times=m-n0), rep(0, times=n0))

rho <- -0.8

# covariance matrix
# Note that BH does not control the FDR if rho < 0
Sigma <- matrix(nrow=length(mu), ncol=length(mu))
diag(Sigma) <- 1
for (i in 1:length(mu)) {
	for (j in 1:length(mu)) {
		if (j != i) {
		Sigma[i, j] <- rho / (m-1)
		}
	}
}

alpha <- 0.05

a <- A

nsim <- 1000

rej.maxSC <- rej.eBH <- FDP.maxSC <- FDP.eBH <- rep(NA, times=nsim)


for (k in 1:nsim) {
	
	Z <- rmvnorm(100, mean=mu, sigma=Sigma)
	
	evaluesH <- rep(NA, times=length(mu))
	for (j in 1:length(mu)) {
		evaluesH[j] <- ztestE(Z[,j], a = a)
		}
	
	out.maxSC <- eBHplus(evaluesH, alpha=alpha)
	rej.maxSC[k] <- out.maxSC$rejections

	fp.maxSC <- sum(out.maxSC$indices > m-n0)
	FDP.maxSC[k] <- fp.maxSC / max(1, out.maxSC$rejections)

	
	out.eBH <- eBH(evaluesH, alpha=alpha)
	rej.eBH[k] <- out.eBH$rejections

 	fp.eBH <- sum(out.eBH$indices > m-n0)
 	FDP.eBH[k] <- fp.eBH / max(1, out.eBH$rejections)
	

	
	print(k)
	
}

TotalS.FDP <- mean(FDP.maxSC)
TotaleBH.FDP <- mean(FDP.eBH)
PowerS <- mean(rej.maxSC) / (m-n0)
PowereBH <- mean(rej.eBH) / (m-n0)

TotalS.FDP <- c(mean(FDP.maxSC), TotalS.FDP)
TotaleBH.FDP <- c(mean(FDP.eBH), TotaleBH.FDP)
PowerS <- c(mean(rej.maxSC) / (m-n0), PowerS)
PowereBH <- c(mean(rej.eBH) / (m-n0), PowereBH)

x <- c(0.125, 0.25, 0.375, 0.5)


