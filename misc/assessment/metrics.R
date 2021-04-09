#!/bin/env Rscript

require(lpSolve)

PongDist <- function(x,y,K){
  thresh <- 0.1/K
  arr <- x-y
  arr2 <- which(x+y > thresh)
  arr <- arr[arr2]
  return(1 - sqrt(sum(arr^2) / (2*length(arr))))
}

PongMatch <- function(Theta_true,Theta){
  # Matches Theta to Theta_true
  # Returns reordered version of Theta
  dist_mat <- apply(Theta_true,2,function(x){apply(Theta,2,function(y){PongDist(x,y,K=ncol(Theta))})})
  col_ord <- apply(lp.assign(dist_mat, direction = "max")$solution,2,which.max)
  return(Theta[,col_ord])
}

MAE <- function(Theta_true,Theta){
  # Function to calculate MAE
  Theta_permute <- PongMatch(Theta_true,Theta)
  err <- median(apply(abs((as.matrix(Theta_true)-as.matrix(Theta_permute))),1,median))
  return(err)
}

ErrMinPerm <- function(Theta_true,Theta){
	# Function to calculate RMSE
	Theta_permute <- PongMatch(Theta_true,Theta)
	err <- norm(as.matrix(Theta_true) - as.matrix(Theta_permute),'F') / sqrt(nrow(Theta)*ncol(Theta))
	return(err)
}

JSD <- function(Theta_true,Theta_unmatched,eps=1e-9){
  # Calculates Jensen-Shannon Divergence
  Theta <- PongMatch(Theta_true,Theta_unmatched)
  Theta_true2 <- Theta_true
  Theta2 <- Theta
  if (any(Theta_true2 == 0)){
    Theta_true2[Theta_true2 == 0] <- eps
    Theta_true2 <- t(apply(Theta_true2,1,function(x){x/sum(x)}))
  }
  if (any(Theta2 == 0)){
    Theta2[Theta2 == 0] <- eps
    Theta2 <- t(apply(Theta2,1,function(x){x/sum(x)}))
  }
  M <- 0.5 * (Theta_true2 + Theta2)
  KL_mat1 <- t(sapply(1:nrow(Theta_true2),function(i){as.numeric(Theta_true2[i,]*log(Theta_true2[i,]/M[i,]))}))
  KL_mat2 <- t(sapply(1:nrow(Theta2),function(i){as.numeric(Theta2[i,]*log(Theta2[i,]/M[i,]))}))
  KL_mat <- (KL_mat1 + KL_mat2) * 0.5
  return(mean(rowSums(KL_mat*is.finite(KL_mat),na.rm = T)))
}

KL <- function(Theta_true,Theta,eps=1e-9){
  # Calculate KL Divergence
  Theta_true2 <- Theta_true
  Theta2 <- PongMatch(Theta_true,Theta)
  if (any(Theta_true2 == 0)){
    Theta_true2[Theta_true2 == 0] <- eps
    Theta_true2 <- t(apply(Theta_true2,1,function(x){x/sum(x)}))
  }
  if (any(Theta2 == 0)){
    Theta2[Theta2 == 0] <- eps
    Theta2 <- t(apply(Theta2,1,function(x){x/sum(x)}))
  }
  KL_mat <- t(sapply(1:nrow(Theta_true2),function(i){as.numeric(Theta_true2[i,]*log(Theta_true2[i,]/Theta2[i,]))}))
  return(mean(rowSums(KL_mat*is.finite(KL_mat),na.rm = T)))
}

CalculateMetrics <- function(Theta_true,Theta){
	# Function to calculate metrics
	theta_permute <- Theta
	kld <- KL(Theta_true,theta_permute)
	jsd <- JSD(Theta_true,theta_permute)
	rmse <- ErrMinPerm(Theta_true,Theta)
	metrics <- c(kld,jsd,rmse)
	names(metrics) <- c("KL","JSD","RMSE")
	return(metrics)
}

