# Mixed-model analysis of high-activity bouts (HABs): Supplemental Fig. 8 and Supplemental Table 2
#
# Fits negative-binomial generalized linear mixed models with glmmTMB (Brooks et al., 2017, R Journal 9:378-400):
#   count ~ 0 + session + <injury / drug terms per time point> + offset(log(bins/120)) + (1 | larva)
# - random intercept for each larva (repeated recordings of the same individuals)
# - fixed effect for each recording session (plate x time point), because the HAB threshold is set within each recording
# - offset for the number of one-minute bins in a recording (3-h recordings = 180 bins; 103-min recording = 103 bins)
# Segments without treatment test injury x time; post-treatment segments test injury x drug across time points.
# Per-time-point rate ratios are Holm-corrected across the time points of each comparison within a segment.
#
# Input: the *_long.csv files in data/ (one row per larva x time point).
# Usage:  Rscript supp_fig8_models.R            (writes results_R.csv and results_R_tests.csv)
# Needs:  install.packages("glmmTMB")

suppressPackageStartupMessages(library(glmmTMB))
data_dir <- "data"

prep <- function(d) {
  tps <- unique(d$tp[order(d$tp_order)])
  d$tp      <- factor(d$tp, levels = tps)
  d$session <- factor(d$session)
  d$larva   <- factor(d$larva)
  d$off     <- log(d$bins / 120)
  has_drug  <- !all(is.na(d$drug))
  for (i in seq_along(tps)) {
    ind <- as.numeric(d$tp == tps[i])
    d[[paste0("T", i)]] <- ind * d$tbi
    if (has_drug) {
      d[[paste0("D", i)]]  <- ind * d$drug
      d[[paste0("TD", i)]] <- ind * d$tbi * d$drug
    }
  }
  list(d = d, tps = as.character(tps))
}

fit_nb <- function(rhs, d) {
  f <- as.formula(paste("count ~ 0 + session +", rhs, "+ offset(off) + (1 | larva)"))
  glmmTMB(f, family = nbinom2, data = d)
}

terms_for <- function(prefix, n) paste0(prefix, seq_len(n))

# rate ratio at each time point for the sum of the coefficients named by `prefixes`
contrast <- function(m, prefixes, tps, label) {
  b <- fixef(m)$cond
  V <- vcov(m)$cond
  rows <- lapply(seq_along(tps), function(i) {
    w <- setNames(numeric(length(b)), names(b))
    w[paste0(prefixes, i)] <- 1
    est <- sum(w * b)
    se  <- sqrt(as.numeric(t(w) %*% V[names(b), names(b)] %*% w))
    data.frame(tp = tps[i], est = est, se = se)
  })
  r <- do.call(rbind, rows)
  r$p      <- 2 * pnorm(-abs(r$est / r$se))
  r$p_holm <- p.adjust(r$p, "holm")
  data.frame(contrast = label, tp = r$tp, RR = exp(r$est),
             lo = exp(r$est - 1.96 * r$se), hi = exp(r$est + 1.96 * r$se),
             p = r$p, p_holm = r$p_holm)
}

lrt <- function(m_small, m_big, test, d) {
  chi2 <- as.numeric(2 * (logLik(m_big) - logLik(m_small)))
  df   <- attr(logLik(m_big), "df") - attr(logLik(m_small), "df")
  data.frame(test = test, chi2 = chi2, df = df, p = pchisq(chi2, df, lower.tail = FALSE),
             larvae = nlevels(d$larva), n = nrow(d))
}

injury_by_time <- function(d) {
  P <- prep(d); d <- P$d; tps <- P$tps; n <- length(tps)
  m0 <- fit_nb("tbi", d)
  m1 <- fit_nb(paste(terms_for("T", n), collapse = " + "), d)
  list(tab = contrast(m1, "T", tps, "TBI vs sham"),
       test = lrt(m0, m1, "injury x time", d))
}

injury_by_drug <- function(d) {
  P <- prep(d); d <- P$d; tps <- P$tps; n <- length(tps)
  base <- paste(c(terms_for("T", n), terms_for("D", n)), collapse = " + ")
  m_no   <- fit_nb(base, d)
  m_full <- fit_nb(paste(base, "+", paste(terms_for("TD", n), collapse = " + ")), d)
  tab <- rbind(contrast(m_full, "T",           tps, "TBI vs sham (vehicle)"),
               contrast(m_full, c("D", "TD"),  tps, "Drug vs vehicle (injured)"),
               contrast(m_full, "D",           tps, "Drug vs vehicle (sham)"))
  list(tab = tab, test = lrt(m_no, m_full, "injury x drug (time-varying)", d))
}

jobs <- data.frame(
  panel = c("A", "B", "B", "C", "C", "D", "D"),
  file  = c("fig7", "fig8_pre", "fig8_post", "fig9_pre", "fig9_post", "fig10_pre", "fig10_post"),
  kind  = c("time", "time", "drug", "time", "drug", "time", "drug"),
  stringsAsFactors = FALSE)

tabs <- list(); tests <- list()
for (k in seq_len(nrow(jobs))) {
  d <- read.csv(file.path(data_dir, paste0(jobs$file[k], "_long.csv")),
                colClasses = c(plate = "character", session = "character", tp = "character"))
  res <- if (jobs$kind[k] == "time") injury_by_time(d) else injury_by_drug(d)
  res$tab  <- cbind(panel = jobs$panel[k], experiment = d$experiment[1], segment = d$segment[1], res$tab)
  res$test <- cbind(panel = jobs$panel[k], experiment = d$experiment[1], segment = d$segment[1], res$test)
  tabs[[k]] <- res$tab; tests[[k]] <- res$test
  cat(sprintf("%s %-5s %-18s %-28s chi2(%d) = %.2f, p = %.3g\n", jobs$panel[k], d$experiment[1], d$segment[1],
              res$test$test, res$test$df, res$test$chi2, res$test$p))
}
write.csv(do.call(rbind, tabs),  "results_R.csv", row.names = FALSE)
write.csv(do.call(rbind, tests), "results_R_tests.csv", row.names = FALSE)

# Values from the original analysis (Python, adaptive Gauss-Hermite quadrature, 15 nodes).
# glmmTMB uses the Laplace approximation, so small differences in the last decimal are expected.
expected <- data.frame(
  segment = c("Fig7 whole", "Fig8 pre", "Fig8 post", "Fig9 pre", "Fig9 post", "Fig10 pre", "Fig10 post"),
  chi2    = c(41.60, 14.36, 29.12, 16.87, 27.08, 27.93, 28.52),
  df      = c(13, 6, 7, 6, 8, 5, 9))
cat("\nExpected (original analysis):\n"); print(expected)
