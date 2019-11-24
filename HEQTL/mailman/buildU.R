buildU <- function(m, sigma) {
	S = length(sigma)
	if (m <= 1) {
		return (t(matrix(sigma)))
	} else {
		Usub = buildU(m - 1, sigma)
		up = t(matrix(rep(sigma, each = ncol(Usub))))
		down = do.call(cbind, replicate(S, Usub, simplify = FALSE))
		return (rbind(up, down))
	}
}
