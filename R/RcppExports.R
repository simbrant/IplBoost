# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.compute_u_j <- function(j, status, mat, times, S0, S1j, n, S, lms, w) {
    .Call(`_IplBoost_compute_u_j`, j, status, mat, times, S0, S1j, n, S, lms, w)
}

.compute_negI_j <- function(j, status, times, S0, S1j, S2j, n, S, lms, w, lambda) {
    .Call(`_IplBoost_compute_negI_j`, j, status, times, S0, S1j, S2j, n, S, lms, w, lambda)
}

.compute_ipl <- function(times, status, mat, betas, lms, w, S, n, p) {
    .Call(`_IplBoost_compute_ipl`, times, status, mat, betas, lms, w, S, n, p)
}

.scale_columns <- function(mat, vec, m, n) {
    invisible(.Call(`_IplBoost_scale_columns`, mat, vec, m, n))
}

.compute_S0 <- function(risk, n, S) {
    .Call(`_IplBoost_compute_S0`, risk, n, S)
}

.compute_S1_j <- function(j, risk, mat, n, S) {
    .Call(`_IplBoost_compute_S1_j`, j, risk, mat, n, S)
}

.compute_S2_j <- function(j, risk, mat, n, S) {
    .Call(`_IplBoost_compute_S2_j`, j, risk, mat, n, S)
}

