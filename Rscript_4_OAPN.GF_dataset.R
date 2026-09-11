========================Code _for _real_life_data==========================
  # Distribution Comparison for Given Data
  # OAPN: Odd Alpha Power Normal distribution
  # Normal baseline with generating family:
  #     G(x) = alpha^F0(x) * F0(x) /(1 - F0(x) + alpha^F0(x) * F0(x))
  #
  # Closely related Normal-baseline competitors included:
  # Normal, Exponentiated Normal, Alpha-Power Transformed Normal,
  # Proportional-Odds / Marshall-Olkin Normal, Odd Log-Logistic Normal,
  # Generalized Odd Log-Logistic Normal, Kumaraswamy Normal, Beta Normal,
  # Transmuted Normal, and Topp-Leone Normal.
  #
  # Outputs:
  #   OAPN_model_comparison.csv       includes estimated parameters column
  #   OAPN_parameter_estimates.csv    long-format parameter estimates
  #   OAPN_all_plots.pdf              PDF, CDF, survival, hazard, QQ, PP, PIT
  # ========================================================================

rm(list = ls())
set.seed(123)

# =========================================================
# Packages
# =========================================================
pkgs <- c("goftest")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}

# =========================================================
# Data
# Replace x below with your own data vector if needed.
# Since the Normal baseline has support on the full real line, the code does
# NOT restrict x to positive values.
# =========================================================
x <- c(0.0094, 0.0500, 0.4064, 4.6307, 7.1645, 7.2316, 8.2616, 9.2662, 9.3812, 9.5223, 9.8783,
       10.4791, 11.0760, 11.3250, 11.5284, 11.9226, 12.0294, 12.5381, 12.8049, 13.4615, 13.8530,
       5.1741, 5.8808, 6.3348, 10.4077, 10.0192, 9.9346, 12.1835, 12.0740, 12.354)
x <- na.omit(as.numeric(x))
x <- x[is.finite(x)]
n <- length(x)

if (n < 3) stop("At least three finite observations are required.")

eps <- 1e-12
x_ks <- jitter(x, factor = 1e-8)

# =========================================================
# Numerical helpers
# =========================================================
clip01 <- function(u, eps = 1e-12) {
  pmin(pmax(u, eps), 1 - eps)
}

safe_density <- function(d, eps = 1e-12) {
  d[!is.finite(d) | d <= eps] <- eps
  d
}

safe_cdf <- function(F, eps = 1e-12) {
  F[!is.finite(F)] <- NA
  clip01(F, eps)
}

fmt_par <- function(p, digits = 6) {
  paste0(names(p), "=", formatC(as.numeric(p), digits = digits, format = "fg"), collapse = "; ")
}

# =========================================================
# Baseline Normal functions
# =========================================================
dN0 <- function(t, mu, sigma) {
  dnorm(t, mean = mu, sd = sigma)
}

pN0 <- function(t, mu, sigma) {
  pnorm(t, mean = mu, sd = sigma)
}

# =========================================================
# Distribution functions: Normal-baseline generated models
# Each function takes t and named parameter vector p.
# =========================================================

# --------------------- Normal ----------------------------
pNormal <- function(t, p) {
  safe_cdf(pN0(t, p["mu"], p["sigma"]), eps)
}

dNormal <- function(t, p) {
  safe_density(dN0(t, p["mu"], p["sigma"]), eps)
}

# --------------------- OAPN -----------------------------
# User-specified generating family:
# G = alpha^F * F /(1 - F + alpha^F * F)
# Density multiplier:
# alpha^F * {1 + F(1-F)log(alpha)} / {1 - F + alpha^F F}^2
# For global monotonicity, require alpha >= exp(-4).
pOAPN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  alpha <- p["alpha"]
  A <- exp(log(alpha) * F)
  G <- A * F / (1 - F + A * F)
  safe_cdf(G, eps)
}

dOAPN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  alpha <- p["alpha"]
  L <- log(alpha)
  A <- exp(L * F)
  den <- (1 - F + A * F)^2
  mult <- A * (1 + F * (1 - F) * L) / den
  safe_density(f * mult, eps)
}

# --------------------- Exponentiated Normal / Power Normal ----------------
pEN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  safe_cdf(F^p["a"], eps)
}

dEN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  safe_density(p["a"] * f * F^(p["a"] - 1), eps)
}

# --------------------- Alpha-Power Transformed Normal ----------------------
pAPTN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  alpha <- p["alpha"]
  L <- log(alpha)
  if (abs(L) < 1e-7) return(safe_cdf(F, eps))
  G <- expm1(L * F) / expm1(L)
  safe_cdf(G, eps)
}

dAPTN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  alpha <- p["alpha"]
  L <- log(alpha)
  if (abs(L) < 1e-7) return(safe_density(f, eps))
  safe_density((L * exp(L * F) / expm1(L)) * f, eps)
}

# --------------------- Proportional-Odds / Marshall-Olkin Normal -----------
pMON <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  alpha <- p["alpha"]
  G <- alpha * F / (1 - (1 - alpha) * F)
  safe_cdf(G, eps)
}

dMON <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  alpha <- p["alpha"]
  safe_density(alpha * f / (1 - (1 - alpha) * F)^2, eps)
}

# --------------------- Odd Log-Logistic Normal -----------------------------
pOLLN_lit <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  S <- 1 - F
  a <- p["a"]
  G <- F^a / (F^a + S^a)
  safe_cdf(G, eps)
}

dOLLN_lit <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  S <- 1 - F
  f <- dN0(t, p["mu"], p["sigma"])
  a <- p["a"]
  den <- (F^a + S^a)^2
  safe_density(a * f * F^(a - 1) * S^(a - 1) / den, eps)
}

# --------------------- Generalized Odd Log-Logistic Normal -----------------
pGOLLN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  a <- p["a"]
  theta <- p["theta"]
  U <- F^(a * theta)
  V <- (1 - F^theta)^a
  G <- U / (U + V)
  safe_cdf(G, eps)
}

dGOLLN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  a <- p["a"]
  theta <- p["theta"]
  U <- F^(a * theta)
  V <- (1 - F^theta)^a
  den <- (U + V)^2
  safe_density(a * theta * f * F^(a * theta - 1) * (1 - F^theta)^(a - 1) / den, eps)
}

# --------------------- Kumaraswamy Normal ---------------------------------
pKwN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  a <- p["a"]
  b <- p["b"]
  G <- 1 - (1 - F^a)^b
  safe_cdf(G, eps)
}

dKwN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  a <- p["a"]
  b <- p["b"]
  safe_density(a * b * f * F^(a - 1) * (1 - F^a)^(b - 1), eps)
}

# --------------------- Beta Normal ----------------------------------------
pBetaN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  safe_cdf(pbeta(F, shape1 = p["a"], shape2 = p["b"]), eps)
}

dBetaN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  safe_density(f * dbeta(F, shape1 = p["a"], shape2 = p["b"]), eps)
}

# --------------------- Transmuted Normal ----------------------------------
pTN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  lambda <- p["lambda"]
  G <- (1 + lambda) * F - lambda * F^2
  safe_cdf(G, eps)
}

dTN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  lambda <- p["lambda"]
  safe_density(f * (1 + lambda - 2 * lambda * F), eps)
}

# --------------------- Topp-Leone Normal ----------------------------------
pTLN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  a <- p["a"]
  G <- (F * (2 - F))^a
  safe_cdf(G, eps)
}

dTLN <- function(t, p) {
  F <- clip01(pN0(t, p["mu"], p["sigma"]), eps)
  f <- dN0(t, p["mu"], p["sigma"])
  a <- p["a"]
  safe_density(2 * a * f * (1 - F) * (F * (2 - F))^(a - 1), eps)
}

# =========================================================
# Fitting helpers
# =========================================================
Models <- c()
LL <- c()
NPAR <- c()
PARAMS <- list()
PDF_FUN <- list()
CDF_FUN <- list()
CONV <- c()

add_model <- function(res) {
  Models <<- c(Models, res$model)
  LL <<- c(LL, res$LL)
  NPAR <<- c(NPAR, res$k)
  PARAMS[[length(Models)]] <<- res$par
  PDF_FUN[[length(Models)]] <<- res$pdf
  CDF_FUN[[length(Models)]] <<- res$cdf
  CONV <<- c(CONV, res$convergence)
}

safe_run <- function(expr) {
  tryCatch(eval(expr), error = function(e) {
    cat("Skipped model due to error:", conditionMessage(e), "\n")
    NULL
  })
}

fit_model <- function(model, trans, start, lower, upper, pdf, cdf,
                      nstarts = 35, maxit = 20000) {
  nll <- function(eta) {
    p <- trans(eta)
    dens <- suppressWarnings(pdf(x, p))
    if (any(!is.finite(dens)) || any(dens <= 0)) return(1e100)
    val <- -sum(log(pmax(dens, eps)))
    if (!is.finite(val)) 1e100 else val
  }
  
  starts <- list(start)
  for (i in seq_len(nstarts - 1)) {
    starts[[length(starts) + 1]] <- runif(length(start), lower, upper)
  }
  
  best <- NULL
  for (st in starts) {
    fit <- tryCatch(
      optim(st, nll, method = "L-BFGS-B", lower = lower, upper = upper,
            control = list(maxit = maxit)),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (is.finite(fit$value) && (is.null(best) || fit$value < best$value)) {
      best <- fit
    }
  }
  
  if (is.null(best)) stop(paste("No successful optimization for", model))
  
  # One final polish from the best solution.
  fit2 <- tryCatch(
    optim(best$par, nll, method = "L-BFGS-B", lower = lower, upper = upper,
          control = list(maxit = maxit), hessian = TRUE),
    error = function(e) NULL
  )
  if (!is.null(fit2) && is.finite(fit2$value) && fit2$value <= best$value + 1e-7) {
    best <- fit2
  }
  
  p_hat <- trans(best$par)
  list(model = model, par = p_hat, LL = -best$value,
       pdf = pdf, cdf = cdf, k = length(p_hat), convergence = best$convergence)
}

# =========================================================
# Starting values and bounds
# =========================================================
mu0 <- mean(x)
sigma0 <- sd(x)
if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- max(1, abs(mu0) * 0.25)

mu_lower <- min(x) - 10 * sigma0
mu_upper <- max(x) + 10 * sigma0
sigma_lower <- max(sigma0 / 100, 1e-5)
sigma_upper <- max(sigma0 * 100, 10)

log_shape_lower <- log(0.03)
log_shape_upper <- log(50)

# OAPN monotonicity condition: alpha >= exp(-4)
log_oapn_alpha_lower <- log(exp(-4) + 1e-6)
log_oapn_alpha_upper <- log(50)

base_lower <- c(mu_lower, log(sigma_lower))
base_upper <- c(mu_upper, log(sigma_upper))
base_start <- c(mu0, log(sigma0))

# =========================================================
# Fit models
# =========================================================

# Normal baseline
res <- safe_run(quote({
  p_hat <- c(mu = mean(x), sigma = sqrt(mean((x - mean(x))^2)))
  LL_hat <- sum(log(pmax(dNormal(x, p_hat), eps)))
  list(model = "Normal", par = p_hat, LL = LL_hat,
       pdf = dNormal, cdf = pNormal, k = length(p_hat), convergence = 0)
}))
if (!is.null(res)) add_model(res)

# OAPN: user-specified distribution
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), alpha = exp(eta[3]))
  fit_model(
    model = "OAPN",
    trans = trans,
    start = c(base_start, log(1.5)),
    lower = c(base_lower, log_oapn_alpha_lower),
    upper = c(base_upper, log_oapn_alpha_upper),
    pdf = dOAPN,
    cdf = pOAPN
  )
}))
if (!is.null(res)) add_model(res)

# Exponentiated Normal / Power Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), a = exp(eta[3]))
  fit_model(
    model = "EN",
    trans = trans,
    start = c(base_start, log(1.2)),
    lower = c(base_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper),
    pdf = dEN,
    cdf = pEN
  )
}))
if (!is.null(res)) add_model(res)

# Alpha-Power Transformed Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), alpha = exp(eta[3]))
  fit_model(
    model = "APTN",
    trans = trans,
    start = c(base_start, log(1.5)),
    lower = c(base_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper),
    pdf = dAPTN,
    cdf = pAPTN
  )
}))
if (!is.null(res)) add_model(res)

# Proportional-Odds / Marshall-Olkin Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), alpha = exp(eta[3]))
  fit_model(
    model = "MON",
    trans = trans,
    start = c(base_start, log(1.5)),
    lower = c(base_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper),
    pdf = dMON,
    cdf = pMON
  )
}))
if (!is.null(res)) add_model(res)

# Odd Log-Logistic Normal from the literature
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), a = exp(eta[3]))
  fit_model(
    model = "OLLN",
    trans = trans,
    start = c(base_start, log(1.2)),
    lower = c(base_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper),
    pdf = dOLLN_lit,
    cdf = pOLLN_lit
  )
}))
if (!is.null(res)) add_model(res)

# Generalized Odd Log-Logistic Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]),
                           a = exp(eta[3]), theta = exp(eta[4]))
  fit_model(
    model = "GOLLN",
    trans = trans,
    start = c(base_start, log(1.2), log(1.2)),
    lower = c(base_lower, log_shape_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper, log_shape_upper),
    pdf = dGOLLN,
    cdf = pGOLLN,
    nstarts = 45
  )
}))
if (!is.null(res)) add_model(res)

# Kumaraswamy Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]),
                           a = exp(eta[3]), b = exp(eta[4]))
  fit_model(
    model = "KwN",
    trans = trans,
    start = c(base_start, log(1.2), log(1.2)),
    lower = c(base_lower, log_shape_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper, log_shape_upper),
    pdf = dKwN,
    cdf = pKwN,
    nstarts = 45
  )
}))
if (!is.null(res)) add_model(res)

# Beta Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]),
                           a = exp(eta[3]), b = exp(eta[4]))
  fit_model(
    model = "BetaN",
    trans = trans,
    start = c(base_start, log(1.2), log(1.2)),
    lower = c(base_lower, log_shape_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper, log_shape_upper),
    pdf = dBetaN,
    cdf = pBetaN,
    nstarts = 45
  )
}))
if (!is.null(res)) add_model(res)

# Transmuted Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), lambda = eta[3])
  fit_model(
    model = "TN",
    trans = trans,
    start = c(base_start, 0),
    lower = c(base_lower, -0.999),
    upper = c(base_upper,  0.999),
    pdf = dTN,
    cdf = pTN
  )
}))
if (!is.null(res)) add_model(res)

# Topp-Leone Normal
res <- safe_run(quote({
  trans <- function(eta) c(mu = eta[1], sigma = exp(eta[2]), a = exp(eta[3]))
  fit_model(
    model = "TLN",
    trans = trans,
    start = c(base_start, log(1.2)),
    lower = c(base_lower, log_shape_lower),
    upper = c(base_upper, log_shape_upper),
    pdf = dTLN,
    cdf = pTLN
  )
}))
if (!is.null(res)) add_model(res)

if (length(Models) == 0) stop("No model was successfully fitted.")

# =========================================================
# Goodness-of-fit statistics
# =========================================================
gof_stats <- function(cdf_fun, pars) {
  Fhat <- function(t) {
    F <- cdf_fun(t, pars)
    clip01(F, eps)
  }
  ks <- tryCatch(unname(ks.test(x_ks, Fhat)$statistic), error = function(e) NA)
  cvm <- tryCatch(unname(goftest::cvm.test(x, Fhat)$statistic), error = function(e) NA)
  ad <- tryCatch(unname(goftest::ad.test(x, Fhat)$statistic), error = function(e) NA)
  c(KS = ks, CVM = cvm, AD = ad)
}

GOF <- do.call(rbind, lapply(seq_along(Models), function(i) {
  gof_stats(CDF_FUN[[i]], PARAMS[[i]])
}))

AIC <- -2 * LL + 2 * NPAR
BIC <- -2 * LL + log(n) * NPAR
CAIC <- -2 * LL + NPAR * (log(n) + 1)
HQIC <- -2 * LL + 2 * NPAR * log(log(n))

comparison <- data.frame(
  Model = Models,
  Estimates = sapply(PARAMS, fmt_par),
  NumPar = NPAR,
  LogLik = LL,
  AIC = AIC,
  BIC = BIC,
  CAIC = CAIC,
  HQIC = HQIC,
  KS = GOF[, "KS"],
  CVM = GOF[, "CVM"],
  AD = GOF[, "AD"],
  Convergence = CONV,
  stringsAsFactors = FALSE
)

comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

cat("\n================ MODEL COMPARISON ================\n")
print(comparison)

#write.csv(comparison, "OAPN_model_comparison.csv", row.names = FALSE)

param_table <- do.call(rbind, lapply(seq_along(Models), function(i) {
  data.frame(Model = Models[i], Parameter = names(PARAMS[[i]]),
             Estimate = as.numeric(PARAMS[[i]]), stringsAsFactors = FALSE)
}))
#write.csv(param_table, "OAPN_parameter_estimates.csv", row.names = FALSE)

# =========================================================
# Robust quantile function from fitted CDF
# =========================================================
q_from_cdf <- function(p, cdf_fun, pars, data_x) {
  sapply(p, function(prob) {
    prob <- min(max(prob, eps), 1 - eps)
    spread <- sd(data_x)
    if (!is.finite(spread) || spread <= 0) spread <- 1
    lower <- min(data_x) - 8 * spread
    upper <- max(data_x) + 8 * spread
    
    f_lower <- cdf_fun(lower, pars) - prob
    f_upper <- cdf_fun(upper, pars) - prob
    
    iter <- 0
    while (is.finite(f_lower) && is.finite(f_upper) && f_lower * f_upper > 0 && iter < 30) {
      lower <- lower - 2 * spread * (iter + 1)
      upper <- upper + 2 * spread * (iter + 1)
      f_lower <- cdf_fun(lower, pars) - prob
      f_upper <- cdf_fun(upper, pars) - prob
      iter <- iter + 1
    }
    
    if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower * f_upper > 0) {
      return(NA_real_)
    }
    
    tryCatch(
      uniroot(function(z) cdf_fun(z, pars) - prob,
              lower = lower, upper = upper, tol = 1e-10)$root,
      error = function(e) NA_real_
    )
  })
}

# =========================================================
# Plot helpers
# =========================================================
S_FUN <- function(t, cdf_fun, pars) {
  pmax(1 - cdf_fun(t, pars), eps)
}

H_FUN <- function(t, pdf_fun, cdf_fun, pars) {
  h <- pdf_fun(t, pars) / S_FUN(t, cdf_fun, pars)
  h[!is.finite(h) | h < 0] <- NA
  h
}

panel_layout <- function(m) {
  nr <- ceiling(sqrt(m))
  nc <- ceiling(m / nr)
  c(nr, nc)
}

plot_overlay <- function(xx, yfun, ylab, main, ylim_quantile = 0.995) {
  cols <- rainbow(length(Models))
  y_list <- lapply(seq_along(Models), function(i) yfun(xx, i))
  y_all <- unlist(y_list)
  y_all <- y_all[is.finite(y_all)]
  y_top <- if (length(y_all) > 0) as.numeric(quantile(y_all, ylim_quantile, na.rm = TRUE)) else 1
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  
  plot(xx, y_list[[1]], type = "l", lwd = 2, col = cols[1],
       ylim = c(0, y_top * 1.05), xlab = "x", ylab = ylab, main = main)
  if (length(Models) > 1) {
    for (i in 2:length(Models)) lines(xx, y_list[[i]], col = cols[i], lwd = 2)
  }
  legend("topright", legend = Models, col = cols, lwd = 2,
         cex = 0.68, bty = "n")
}



############################################################
# COMPLETE OAPN DIAGNOSTIC PLOTS USING dev.new()
############################################################

xx <- seq(min(x), max(x), length.out = 700)

# Generate rainbow colours
cols <- rainbow(length(Models))

# Make OAPN black
cols[grep("OAPN", Models, ignore.case = TRUE)] <- "black"

# Plot order: OAPN first, all others, Normal last
plot_order <- c(
  grep("OAPN", Models, ignore.case = TRUE),
  setdiff(
    seq_along(Models),
    c(
      grep("OAPN", Models, ignore.case = TRUE),
      grep("Normal", Models, ignore.case = TRUE)
    )
  ),
  grep("Normal", Models, ignore.case = TRUE)
)

############################################################
# 1. FITTED PDFs OVER HISTOGRAM
############################################################
windowsFonts(Times = windowsFont("Times New Roman"))

tiff(
  filename = "Figure5_PDF.tiff",
  width = 10,
  height = 7,
  units = "in",
  res = 600,          # Publication quality (>=300 dpi)
  compression = "lzw"
)
par(
  family = "Times",
  cex = 0.9,
  cex.axis = 0.9,
  cex.lab = 0.9,
  cex.main = 0.9,
  cex.sub = 0.9,
  las = 1,
  bty = "l"
)

dev.new(width = 10, height = 7)

hist(
  x,
  probability = TRUE,
  col = "grey85",
  border = "white",
  main = "",
  xlab = "x",
  ylim = c(
    0,
    max(
      density(x)$y,
      sapply(
        seq_along(Models),
        function(i)
          max(
            PDF_FUN[[i]](
              xx,
              PARAMS[[i]]
            ),
            na.rm = TRUE
          )
      )
    ) * 1.35
  )
)

# Draw fitted PDFs
for(i in plot_order)
{
  lines(
    xx,
    PDF_FUN[[i]](
      xx,
      PARAMS[[i]]
    ),
    col = cols[i],
    lwd = 3
  )
}

# Legend
legend(
  "topright",
  inset = c(0, -0.08),
  xpd = TRUE,
  legend = Models[plot_order],
  col = cols[plot_order],
  lwd = 3,
  cex = 0.9,
  bty = "n"
)
############################################################
# Locate the OAPN model
############################################################

OAPN_index <- grep(
  "OAPN",
  Models,
  ignore.case = TRUE
)

############################################################
# 2. FITTED CDF OVER EMPIRICAL CDF
############################################################
windowsFonts(Times = windowsFont("Times New Roman"))
tiff(
  filename = "Figure6_CDF.tiff",
  width = 10,
  height = 7,
  units = "in",
  res = 600,          # Publication quality (>=300 dpi)
  compression = "lzw"
)
par(
  family = "Times",
  cex = 0.9,
  cex.axis = 0.9,
  cex.lab = 0.9,
  cex.main = 0.9,
  cex.sub = 0.9,
  las = 1,
  bty = "l"
)
dev.new(width = 10, height = 7)
plot(
  ecdf(x),
  main = "",
  xlab = "x",
  ylab = "Cumulative Distribution Function",
  verticals = TRUE,
  do.points = FALSE
)
lines(
  xx,
  CDF_FUN[[OAPN_index]](
    xx,
    PARAMS[[OAPN_index]]
  ),
  col = "black",
  lwd = 3
)

legend(
  "bottomright",
  legend = Models[OAPN_index],
  col = "black",
  lwd = 3,
  cex = 0.9,
  bty = "n"
)
############################################################
# Arrange models: OAPN first, Normal last
############################################################

plot_order <- c(
  grep("OAPN", Models, ignore.case = TRUE),
  setdiff(
    seq_along(Models),
    c(
      grep("OAPN", Models, ignore.case = TRUE),
      grep("Normal", Models, ignore.case = TRUE)
    )
  ),
  grep("Normal", Models, ignore.case = TRUE)
)

############################################################
# 6. PP PLOTS
############################################################
windowsFonts(Times = windowsFont("Times New Roman"))

tiff(
  filename = "Figure7_PP.tiff",
  width = 10,
  height = 7,
  units = "in",
  res = 600,          # Publication quality (>=300 dpi)
  compression = "lzw"
)
par(
  family = "Times",
  cex = 0.9,
  cex.axis = 0.9,
  cex.lab = 0.9,
  cex.main = 0.9,
  cex.sub = 0.9,
  las = 1,
  bty = "l"
)
lay <- panel_layout(length(Models))

dev.new(width = 12, height = 8)

par(
  mfrow = lay,
  mar = c(4,4,2.5,1)
)

emp_p <- (seq_len(n)-0.5)/n

x_sorted <- sort(x)
for(j in seq_along(plot_order))
{
  i <- plot_order[j]
  
  Ffit <- clip01(
    CDF_FUN[[i]](
      x_sorted,
      PARAMS[[i]]
    ),
    eps
  )
  
  plot(
    Ffit,
    emp_p,
    pch = 16,
    cex = 0.9,
    main = Models[i],
    xlab = "Fitted CDF",
    ylab = "Empirical CDF"
  )
  
  abline(
    0,
    1,
    col = "red",
    lwd = 2
  )
}

par(mfrow = c(1,1))

############################################################
# 7. PIT HISTOGRAMS
############################################################
windowsFonts(Times = windowsFont("Times New Roman"))
tiff(
  filename = "Figure8_PIT.tiff",
  width = 10,
  height = 7,
  units = "in",
  res = 600,          # Publication quality (>=300 dpi)
  compression = "lzw"
)
par(
  family = "Times",
  cex = 0.9,
  cex.axis = 0.9,
  cex.lab = 0.9,
  cex.main = 0.9,
  cex.sub = 0.9,
  las = 1,
  bty = "l"
)
lay <- panel_layout(length(Models))

dev.new(width = 12, height = 8)

par(
  mfrow = lay,
  mar = c(4,4,2.5,1)
)

for(j in seq_along(plot_order))
{
  i <- plot_order[j]
  
  pit <- clip01(
    CDF_FUN[[i]](
      x,
      PARAMS[[i]]
    ),
    eps
  )
  
  hist(
    pit,
    breaks = seq(0, 1, length.out = 11),
    col = "lightblue",
    border = "white",
    main = Models[i],
    xlab = "PIT Values",
    ylab = "Frequency"
  )
  
  abline(
    h = length(pit)/10,
    col = "red",
    lwd = 2
  )
}

par(mfrow = c(1,1)) 
