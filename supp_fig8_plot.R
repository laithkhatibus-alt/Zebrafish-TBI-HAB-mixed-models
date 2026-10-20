# Supplemental Fig. 8: model-estimated rate ratios from results_R.csv
# Run after supp_fig8_models.R, from the same folder (Session > Set Working Directory > To Source File Location).
# Needs ggplot2:  install.packages("ggplot2")
# Writes Supplemental_Fig8.pdf (vector) and Supplemental_Fig8.tiff (600 dpi) into this folder.

library(ggplot2)

res <- read.csv("results_R.csv", stringsAsFactors = FALSE)

# x position: midpoint of each recording, in hours after injury
res$time <- sapply(strsplit(res$tp, "-"), function(z) mean(as.numeric(z)))

# three series per panel, offset slightly so their error bars don't overlap
series_labels <- c("Injury effect: TBI vs sham (vehicle groups after treatment)",
                   "Drug effect, injured larvae: TBI + drug vs TBI + vehicle",
                   "Drug effect, sham larvae: sham + drug vs sham + vehicle")
res$series <- ifelse(grepl("^TBI vs sham", res$contrast), series_labels[1],
               ifelse(res$contrast == "Drug vs vehicle (injured)", series_labels[2], series_labels[3]))
res$series <- factor(res$series, levels = series_labels)
res$x <- res$time + c(0, -0.55, 0.55)[as.integer(res$series)]

# filled symbol = Holm-adjusted p < 0.05; open symbol = not significant
res$fill_key <- ifelse(res$p_holm < 0.05, as.character(res$series), "ns")

titles <- c(A = "A   No treatment (Fig. 7)",
            B = "B   Delayed indomethacin (Fig. 8)",
            C = "C   Delayed retigabine (Fig. 9)",
            D = "D   Delayed MCC950 (Fig. 10)")
res$panel_lab <- factor(titles[res$panel], levels = titles)

# dark phases (grey) and time of drug addition (dashed line) for each panel
dark <- data.frame(panel = c("A", "A", "B", "B", "C", "C", "D", "D"),
                   xmin  = c(12, 36, 12, 36, 13, 36, 13, 36),
                   xmax  = c(22, 45, 22, 45, 22, 46, 22, 46))
dark$panel_lab <- factor(titles[dark$panel], levels = titles)
drug <- data.frame(panel = c("B", "C", "D"), x = c(26.5, 26.5, 22.5),
                   label = c("indomethacin added", "RTG added", "MCC950 added"))
drug$panel_lab <- factor(titles[drug$panel], levels = titles)

cols   <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), series_labels)
shapes <- setNames(c(21, 22, 24), series_labels)
fills  <- c(cols, ns = "white")

p <- ggplot(res, aes(x = x, y = RR)) +
  geom_rect(data = dark, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),  # ymin = 0, not -Inf: on a log axis 0 maps to the panel floor
            inherit.aes = FALSE, fill = "grey92") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_vline(data = drug, aes(xintercept = x), linetype = "dashed", linewidth = 0.4) +
  geom_text(data = drug, aes(x = x + 0.5, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 2.6) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = series), linewidth = 0.6) +
  geom_point(aes(colour = series, shape = series, fill = fill_key), size = 2, stroke = 0.8) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_fill_manual(values = fills, guide = "none") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(fill = unname(cols))),
         shape  = guide_legend(ncol = 1, override.aes = list(fill = unname(cols)))) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30),
                labels = c("0.01", "0.03", "0.1", "0.3", "1", "3", "10", "30")) +
  scale_x_continuous(breaks = seq(0, 54, 6), limits = c(-1, 59), expand = c(0, 0)) +
  facet_wrap(~ panel_lab, ncol = 1, scales = "free_y") +
  labs(x = "Time after injury (h; midpoint of each recording)", y = "Rate ratio",
       caption = "Filled symbols: Holm-adjusted p < 0.05. Open symbols: not significant.") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", size = 9),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 7.5))

ggsave("Supplemental_Fig8.pdf",  p, width = 7, height = 9.5)
ggsave("Supplemental_Fig8.tiff", p, width = 7, height = 9.5, dpi = 600, compression = "lzw")
print(p)
