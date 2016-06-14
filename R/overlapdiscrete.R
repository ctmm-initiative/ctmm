overlapdiscrete <- function(P, Q, method="bc", dA = 1.0) {
	switch(method,
	       bc  = BhattacharyyaCoefficientDiscrete(P, Q, dA),
	       bd  = BhattacharyyaDistanceDiscrete(P, Q, dA),
	       hd  = HellingerDistanceDiscrete(P, Q, dA),
	       tvd = TotalVariationDistanceDiscrete(P, Q, dA),
	       kld = KullbackLeiblerDivergence(P, Q, dA),
	       jsd = JensenShannonDistanceDiscrete(P, Q, dA),
	       BhattacharyyaCoefficientDiscrete(P, Q, dA))
}

BhattacharyyaCoefficientDiscrete <- function(P, Q, dA = 1.0) {
	sum(sqrt(P * Q)) * dA
}

BhattacharyyaDistanceDiscrete    <- function(P, Q, dA = 1.0) {
	-log(BhattacharyyaCoefficientDiscrete(P, Q, dA))
}

HellingerDistanceDiscrete        <- function(P, Q, dA = 1.0) {
	sqrt(1 - BhattacharyyaCoefficientDiscrete(P, Q, dA))
}

TotalVariationDistanceDiscrete   <- function(P, Q, dA = 1.0) {
	0.5 * sum(abs(P - Q)) * dA
}

KullbackLeiblerDivergence        <- function(P, Q, dA = 1.0) {
	#KLD is only defined if P == 0 implies Q == 0
	#if we assume that the support of P, Q are identical, then
	#the misalignment is likely due to simple numerical instability
	#which we correct for here
	if (any(xor(P == 0, Q == 0))) {
		P <- P + .Machine$double.xmin
		Q <- Q + .Machine$double.xmin
	}
	sum(P * log2(P / Q)) * dA
}

JensenShannonDistanceDiscrete    <- function(P, Q, dA = 1.0) {
	if (any(xor(P == 0, Q == 0))) {
		P <- P + .Machine$double.xmin
		Q <- Q + .Machine$double.xmin
	}
	M <- 0.5 * (P + Q)
	0.5 * KullbackLeiblerDivergence(P, M, dA = dA) + 0.5 * KullbackLeiblerDivergence(Q, M, dA = dA)
}