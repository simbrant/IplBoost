
test_that("Compute S0",
{
library(survival)
data("pbc")
t = pbc$time[1:100]/365
d = pbc$status[1:100]
mat = pbc[1:100, c(4, 7)]

d <- d[order(t)]
mat <- as.matrix(mat[order(t), ])
t <- t[order(t)] + rnorm(length(t), sd = 0.01)
lms <- seq(0, round(0.75*max(t)), round(0.75*max(t))/3)
betas <- matrix(seq(1, length(lms)*dim(mat)[2], 1), ncol=dim(mat)[2])
risk.s <- exp(mat %*% t(betas))
tmp = risk.s
expect_equal(as.matrix(compute_S0(tmp, t, length(t), length(lms))),
             as.matrix(apply(risk.s[order(-t), ], 2, cumsum)[order(-(1:length(t))), ]),
             check.attributes=FALSE)})

betas <- matrix(rnorm(length(lms)*dim(mat)[2]), ncol=dim(mat)[2])
risk.s <- exp(mat %*% t(betas))
diff <- ((as.matrix(compute_S0(as.matrix(risk.s), t, length(t), length(lms))) - as.matrix(apply(risk.s[order(-t), ], 2, cumsum)[order(-(1:length(t))), ])))
mean(apply(diff, 1, function(x)mean(abs(x))))
