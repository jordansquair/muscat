plotMeanDisp <- function(x, filter = TRUE, q = 0.01) {
    rd <- rowData(x)
    df <- data.frame(
        mean = rd$AveLogCPM,
        disp = rd$tagwise.dispersion,
        trend = rd$trended.dispersion)
    if (filter)
        df <- filter(df,
            disp > quantile(disp,   q),
            disp < quantile(disp, 1-q),
            mean < quantile(mean, 1-q))
    p <- ggplot(df, aes(mean, disp)) +
        geom_point(aes(col = "royalblue"), size = 2, fill = "white") +
        geom_point(size = 1, col = "white", fill = "white") +
        geom_line(size = 1, show.legend = FALSE,
            aes(col = "tomato", mean, trend)) +
        geom_hline(size = 0.5, lty = 2, show.legend = FALSE,
            aes(col = "black", yintercept = metadata(x)$common.dispersion)) +
        scale_y_log10(labels = function(u) format(u, scientific = TRUE)) + 
        theme_bw() + theme(
            panel.border = element_rect(size = 1),
            plot.background = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(rep(0,4), "mm"),
            axis.ticks = element_line(size = 0.4),
            axis.text = element_text(color = "black"),
            axis.title = element_blank(),
            legend.position = c(1,1),
            legend.justification = c(1,1),
            legend.key.size = unit(0, "mm"),
            legend.background = element_blank(),
            legend.margin = margin(-1,1,0,0,"mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "lightgrey", size = 0.4)) +
        scale_color_identity(NULL, guide = "legend",
            breaks = c("royalblue", "tomato", "black"),
            labels = c("tagwise", "trended", "common")) +
        guides(color = guide_legend(override.aes = list(size = 3)))
    ggMarginal(p, type = "histogram", bins = 50, fill = "grey", col = "white")
}
library(MASS)

y <- rowData(x)
# q <- quantile(y$tagwise_disp, 0.01)
# y <- y[y$tagwise_disp > q, ]

ltd <- log(td <- y$tagwise_disp)
hist(td, n = 50, prob = TRUE)

par(mfrow = c(1, 2))
hist(y$mean_logCPM, n = 50)
fit <- fitdistr(ltd, "normal")
ests <- as.list(fit$estimate)
xs <- seq(min(ltd), max(ltd), l = 1e3)
hist(ltd, n = 50, prob = TRUE)
lines(xs, dnorm(xs, ests$mean, ests$sd), col = "red")

fit <- fitdistr(td, "t")
ests <- as.list(fit$estimate)
xs <- seq(min(ltd), max(ltd), l = 1e3)
hist(ltd, n = 50, prob = TRUE)
#lines(xs, rnorm(xs, ests$shape, ests$rate), col = "red")

fun <- approxfun(density(ltd))
lines(xs, fun(xs), col = "red")

plotMeanDisp(x)
