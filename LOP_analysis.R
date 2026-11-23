# =====================================================================
# PSY310 Solo Project - Levels-of-processing experiment: full analysis
# Put this script in the SAME folder as your 8 CSV files, then in RStudio:
#   Session > Set Working Directory > To Source File Location
# then click "Source" (or run line by line with Ctrl+Enter).
# =====================================================================

rm(list = ls())
if (!requireNamespace("geepack", quietly = TRUE)) install.packages("geepack")
library(geepack)

# ---------------------------------------------------------------------
# 1. READ DATA AND SCORE RECALL
# ---------------------------------------------------------------------
files <- sort(list.files(pattern = "^[0-9]+_LOP_experiment.*\\.csv$"))
cat("Files found:", length(files), "\n")
stopifnot(length(files) == 8)

PRACTICE <- c("spoon", "snake")     # practice words = intrusions, not recall
trials <- list()
intrusions <- data.frame()

for (f in files) {
  pid <- as.integer(sub("^([0-9]+)_.*", "\\1", f))
  d <- read.csv(f, na.strings = c("", "NA"), stringsAsFactors = FALSE)
  m <- d[d$condition %in% c("deep", "shallow"), ]      # the 16 main trials
  items <- m$word

  recalled_words <- NULL                                # NULL = no recall data
  if ("recall_textbox.text" %in% names(d)) {
    txt <- d$recall_textbox.text[!is.na(d$recall_textbox.text)]
    if (length(txt) > 0) {
      toks <- unlist(strsplit(gsub("\\\\n", "\n", txt[1]), "[\r\n]+"))
      toks <- unique(trimws(tolower(toks)))
      toks <- toks[toks != ""]
      hit <- character(0); n_prac <- 0; n_other <- 0
      for (t in toks) {
        if (t %in% items) {
          hit <- c(hit, t)
        } else if (t %in% PRACTICE) {
          n_prac <- n_prac + 1
        } else {
          dist <- adist(t, items)                       # accept 1-letter typos (boel -> bowl)
          if (min(dist) <= 1) hit <- c(hit, items[which.min(dist)]) else n_other <- n_other + 1
        }
      }
      recalled_words <- unique(hit)
      intrusions <- rbind(intrusions, data.frame(P = pid, practice = n_prac, other = n_other))
    }
  }

  corr <- m$key_resp.corr; corr[is.na(corr)] <- 0       # timeouts count as incorrect
  trials[[length(trials) + 1]] <- data.frame(
    P = pid, word = m$word, cond = m$condition, corr = corr,
    answered = !is.na(m$key_resp.keys), rt = m$key_resp.rt,
    recalled = if (is.null(recalled_words)) NA_integer_ else as.integer(m$word %in% recalled_words)
  )
}
T <- do.call(rbind, trials)
Tr <- T[!is.na(T$recalled), ]                           # trials with recall data

cat("\nParticipants with recall data:", paste(sort(unique(Tr$P)), collapse = ", "), "\n")
cat("Timeouts:", sum(!T$answered), "of", nrow(T), "trials\n")
print(table(T$cond[!T$answered]))
cat("\nIntrusions (practice words / other):\n"); print(intrusions)

# ---------------------------------------------------------------------
# 2. HELPERS
# ---------------------------------------------------------------------
wide <- function(df, y, id = "P") {                     # long -> one row per participant (or per word)
  w <- reshape(df[, c(id, "cond", y)], idvar = id, timevar = "cond", direction = "wide")
  names(w) <- sub(paste0(y, "\\."), "", names(w)); w
}
describe <- function(x, label) {
  n <- length(x); se <- sd(x) / sqrt(n); ci <- qt(.975, n - 1) * se
  cat(sprintf("%-22s n=%d  M=%.3f  SD=%.3f  SE=%.3f  95%%CI [%.3f, %.3f]  Mdn=%.2f\n",
              label, n, mean(x), sd(x), se, mean(x) - ci, mean(x) + ci, median(x)))
}
paired_report <- function(a, b, label) {
  cat("\n---", label, "---\n")
  describe(a, "Condition 1 (deep)"); describe(b, "Condition 2 (shallow)")
  tt <- t.test(a, b, paired = TRUE); print(tt)
  dd <- a - b; n <- length(dd); dz <- mean(dd) / sd(dd); df <- n - 1
  tobs <- unname(tt$statistic)
  lo <- uniroot(function(nc) pt(tobs, df, ncp = nc) - .975, c(-30, 30))$root
  hi <- uniroot(function(nc) pt(tobs, df, ncp = nc) - .025, c(-30, 30))$root
  cat(sprintf("Cohen's dz = %.3f, 95%% CI [%.3f, %.3f]; Hedges gz = %.3f\n",
              dz, lo / sqrt(n), hi / sqrt(n), dz * (1 - 3 / (4 * df - 1))))
  cat("Differences >0 / =0 / <0:", sum(dd > 0), "/", sum(dd == 0), "/", sum(dd < 0), "\n")
  print(shapiro.test(dd))
  invisible(tt)
}

# ---------------------------------------------------------------------
# 3. PRIMARY ANALYSIS: RECALL (by participant, paired t-test)
# ---------------------------------------------------------------------
Rw <- wide(aggregate(recalled ~ P + cond, data = Tr, FUN = sum), "recalled")
print(Rw)
paired_report(Rw$deep, Rw$shallow, "PRIMARY: recall (words out of 8), paired t-test")
cat("\nMean proportion recalled: deep =", round(mean(Rw$deep) / 8, 3),
    " shallow =", round(mean(Rw$shallow) / 8, 3), "\n")
cat("Total words recalled: deep =", sum(Rw$deep), " shallow =", sum(Rw$shallow), "\n")

# non-parametric checks
cat("\nWilcoxon signed-rank (R uses a normal approximation when there are ties/zeros,\n",
    "so this p may differ slightly from the exact p in the report):\n")
print(suppressWarnings(wilcox.test(Rw$deep, Rw$shallow, paired = TRUE, exact = FALSE)))
dd <- Rw$deep - Rw$shallow; flips <- expand.grid(rep(list(c(1, -1)), length(dd)))
perm <- apply(flips, 1, function(s) abs(mean(s * dd)))
cat("Exact sign-flip permutation p =", round(mean(perm >= abs(mean(dd)) - 1e-12), 4), "\n")

# ---------------------------------------------------------------------
# 4. SECONDARY / EXPLORATORY RECALL ANALYSES
# ---------------------------------------------------------------------
# (a) answered trials only
Ra <- wide(aggregate(recalled ~ P + cond, data = Tr[Tr$answered, ], FUN = mean), "recalled")
paired_report(Ra$deep, Ra$shallow, "Recall (proportion), answered trials only")

# (b) by items (16 words; proportion of participants recalling each word)
It <- wide(aggregate(recalled ~ word + cond, data = Tr, FUN = mean), "recalled", id = "word")
print(It)
cat("\n--- By-item paired t-test ---\n"); print(t.test(It$deep, It$shallow, paired = TRUE))
di <- It$deep - It$shallow
cat(sprintf("dz = %.3f; items deep>shallow / tied / deep<shallow: %d / %d / %d\n",
            mean(di) / sd(di), sum(di > 0), sum(di == 0), sum(di < 0)))
print(suppressWarnings(wilcox.test(It$deep, It$shallow, paired = TRUE, exact = FALSE)))

# (c) trial-level logistic regression (GEE, clustered on participant)
Tr$deep <- as.integer(Tr$cond == "deep"); Tr <- Tr[order(Tr$P), ]
g <- geeglm(recalled ~ deep, id = P, data = Tr, family = binomial, corstr = "exchangeable")
print(summary(g))
b <- coef(summary(g))["deep", ]
cat(sprintf("Odds ratio = %.2f, 95%% CI [%.2f, %.2f]\n", exp(b["Estimate"]),
            exp(b["Estimate"] - 1.96 * b["Std.err"]), exp(b["Estimate"] + 1.96 * b["Std.err"])))

# (d) by list (participants 1-4 = List A, 5-8 = List B)
Rw$list <- ifelse(Rw$P <= 4, "A", "B"); Rw$diff <- Rw$deep - Rw$shallow
cat("\nMean deep-minus-shallow difference by list:\n")
print(aggregate(diff ~ list, data = Rw, FUN = function(x) c(mean = mean(x), n = length(x))))

# ---------------------------------------------------------------------
# 5. ENCODING TASK (manipulation check): accuracy and response time, n = 8
# ---------------------------------------------------------------------
Aw <- wide(aggregate(corr ~ P + cond, data = T, FUN = mean), "corr")
Tw <- wide(aggregate(rt ~ P + cond, data = T[T$answered, ], FUN = mean), "rt")
paired_report(Aw$deep, Aw$shallow, "Encoding accuracy (proportion correct)")
paired_report(Tw$deep, Tw$shallow, "Encoding response time (s, answered trials)")

# ---------------------------------------------------------------------
# 6. SENSITIVITY / POWER
# ---------------------------------------------------------------------
cat("\nPower at n = 7, dz = 0.75:", round(power.t.test(n = 7, delta = .75, sd = 1, type = "paired")$power, 3), "\n")
cat("Smallest dz detectable with 80% power at n = 7:",
    round(power.t.test(n = 7, power = .8, sd = 1, type = "paired")$delta, 2), "\n")
cat("n needed for 80% power at dz = 0.75:",
    ceiling(power.t.test(delta = .75, power = .8, sd = 1, type = "paired")$n), "\n")

# ---------------------------------------------------------------------
# 7. FIGURES (saved as PNG in the same folder)
# ---------------------------------------------------------------------
pair_plot <- function(a, b, ylab, ylim, file, scale = 1) {
  png(file, width = 1400, height = 1200, res = 300); par(mar = c(3, 4, 1, 1))
  plot(NA, xlim = c(.6, 2.4), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab, bty = "l")
  axis(1, at = 1:2, labels = c("Deep (living?)", "Shallow (letter E?)"))
  segments(1, a, 2, b, col = "grey60"); points(rep(1, length(a)), a, pch = 16, col = "grey50")
  points(rep(2, length(b)), b, pch = 16, col = "grey50")
  for (i in 1:2) {
    x <- list(a, b)[[i]]; m <- mean(x); ci <- qt(.975, length(x) - 1) * sd(x) / sqrt(length(x))
    arrows(i, m - ci, i, m + ci, angle = 90, code = 3, length = .06, lwd = 2, col = c("#1f4e79", "#b5651d")[i])
    points(i, m, pch = 19, cex = 1.6, col = c("#1f4e79", "#b5651d")[i])
  }
  dev.off()
}
pair_plot(Rw$deep / 8, Rw$shallow / 8, "Proportion recalled (of 8)", c(0, 1), "fig1_recall.png")
pair_plot(Aw$deep, Aw$shallow, "Encoding accuracy", c(.3, 1.05), "fig2a_accuracy.png")
pair_plot(Tw$deep, Tw$shallow, "Mean response time (s)", c(1, 2.5), "fig2b_rt.png")

# ---------------------------------------------------------------------
# 8. EXPORT A WIDE FILE FOR JASP (one row per participant)
# ---------------------------------------------------------------------
names(Rw)[names(Rw) == "deep"] <- "Recall_Deep"; names(Rw)[names(Rw) == "shallow"] <- "Recall_Shallow"
names(Aw)[2:3] <- c("Acc_Deep", "Acc_Shallow"); names(Tw)[2:3] <- c("RT_Deep", "RT_Shallow")
J <- Reduce(function(x, y) merge(x, y, by = "P", all = TRUE), list(Rw[, c("P", "Recall_Deep", "Recall_Shallow")], Aw, Tw))
J <- merge(data.frame(P = 1:8), J, by = "P", all.x = TRUE)
write.csv(J, "LOP_wide_for_JASP.csv", row.names = FALSE, na = "")
print(J)
cat("\nDone. Figures and LOP_wide_for_JASP.csv saved in:", getwd(), "\n")
