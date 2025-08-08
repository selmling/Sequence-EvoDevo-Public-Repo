library("lme4")
library("here")
library("broom")
library("scales")
library("emmeans")
library("cowplot")
library("tidytext")
library("tidyverse")
library("kableExtra")
library("ggbeeswarm")

# ---- helper functions

theme_Elm <- function(base_size = 10){ 
  \
  theme(
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = base_size * 1.1),
    legend.text = element_text(size = base_size),
    legend.title = element_text(size = base_size * 1.1),
    plot.title = element_text(size = base_size * 1.3, hjust = 0.5),
    plot.subtitle = element_text(size = base_size * 1.2, hjust = 0.5),
    plot.caption = element_text(size = base_size * 0.9),
    plot.margin = unit(c(.25, .25, .25, .25), "cm")
  )
}

custom_format <- function(x) {
  x %>%
    as.numeric() %>% # Ensure numeric input
    formatC(format = "f", drop0trailing = TRUE) %>% # Drop trailing zeros
    gsub(pattern = "^0\\.", replacement = ".", .) %>% # Remove leading zero for positive numbers
    gsub(pattern = "^-0\\.", replacement = "-.", .) # Remove leading zero for negative numbers
}

# ---- data import

inf_syl_dat <- read_csv(here("Study 2", "data", "infant_sequence_data.csv")) %>%
  mutate(syllables_per_second = syllable_total/duration,
         age = factor(age, levels = c("5 months", "10 months")))

# ---- descriptive statistics

# number of subjects
inf_syl_dat %>%
  ungroup() %>%
  summarize(subject_count = n_distinct(sub)) %>%
  kbl("pipe")

desc_stats <- inf_syl_dat %>%
  filter(seq_cat == "sequence") %>%  # Only sequences
  group_by(age) %>%                  # Group by age
  summarize(
    Mean_Syllable_Total = mean(syllable_total, na.rm = TRUE),
    SD_Syllable_Total = sd(syllable_total, na.rm = TRUE),
    Mean_Duration = mean(duration, na.rm = TRUE),
    SD_Duration = sd(duration, na.rm = TRUE)
  )

desc_stats %>%
  kbl("pipe")

# difference in sequence duration from 5 to 10 months

m0 <- inf_syl_dat %>%
  filter(seq_cat == "sequence") %>% 
  lmer(duration ~ age + (1 | sub), data = .)

summary(emmeans(m0, "age", contr = "pairwise"), infer = TRUE)$contrast %>%
  mutate(p.value = format(p.value, scientific = TRUE, digits = 3)) %>%
  kbl("pipe")

# ---- visualize

# Syllables per second

seq_syl_x_sec_sumstats <- inf_syl_dat %>%
    filter(seq_cat == "sequence") %>%
    group_by(sub, age) %>%
    summarize(mean_syll_per_sec = mean(syllables_per_second),
              std_syll_per_sec = sd(syllables_per_second))

seq_syl_x_sec_diffs <- pivot_wider(seq_syl_x_sec_sumstats,
                                id_cols = c("sub"),
                                names_from = "age",
                                values_from = c("mean_syll_per_sec", "std_syll_per_sec")) %>%
                    mutate(mean_syll_per_sec_dif = `mean_syll_per_sec_10 months` - `mean_syll_per_sec_5 months`,
                           std_syll_per_sec_dif = `std_syll_per_sec_10 months` - `std_syll_per_sec_5 months`) %>% 
                    select(sub, mean_syll_per_sec_dif, std_syll_per_sec_dif)

seq_comp_stat_label <- data.frame(
  age = 1.5,
  mean_syll_per_sec = 4.5
)

font_sz <- 10

# ---- Figure 4A: Syllables per second over development

p4A <- seq_syl_x_sec_sumstats %>%
    mutate(age = factor(age, levels = c("5 months", "10 months"),
           labels = c("5 mo's", "10 mo's"))) %>%
    ggplot(aes(x = age, y = mean_syll_per_sec)) +
    stat_boxplot(geom = "errorbar", size = 1, width = 0.35) +
    geom_boxplot(size = 1, outlier.shape = 18, width = 0.5,
                 outlier.size = 2, outlier.alpha = .75) +
    geom_line(aes(group = sub), position = "identity",
              alpha = .15, size = 1) +
    geom_point(size = 3.5, alpha = .15,
               stroke = 0, shape = 19) +
    geom_text(data = seq_comp_stat_label, label = "***", size = 6, hjust = .5) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(title = "Sequences",
         y = "Average syllables per second",
         x = "") +
    theme_classic() +
    theme_Elm() +
    theme(
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)))

# Figure 4A Statistics
# ---------------------

inf_syl_dat <- inf_syl_dat %>%
  group_by(age, sub) %>%
  mutate(syll_per_sec = ifelse(seq_cat != "sequence", NA, syllable_total/duration))

inf_syl_dat_dev <- inf_syl_dat %>%
  group_by(age, sub) %>%
  summarize(Mean = mean(syll_per_sec, na.rm = TRUE),
            SD = sd(syll_per_sec, na.rm = TRUE),
            n = n())

inf_syl_dat_dev %>%
  group_by(age) %>%
  summarize(Average = mean(Mean, na.rm = TRUE))

m1 <- inf_syl_dat_dev %>%
  lmer(Mean ~ age + (1 | sub), data = .)

kbl(summary(emmeans(m1, "age", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# ---- Figure 4B: Developmental change in syllables per second

width <- .55

y_lab <- expression(atop("Syllables per second " ~ Delta, "(10 months - 5 months)"))

# Calculate the position for the annotation
p4Blab_dat_c <- seq_syl_x_sec_diffs %>% ungroup() %>% 
  summarize(x = "sequence", mean_syll_per_sec_dif = 3.2)

p4B <- seq_syl_x_sec_diffs %>%
    ggplot(aes("", mean_syll_per_sec_dif)) +
        geom_boxplot(size = 1, color = "grey", outlier.shape = NA,
                     width = 0.5, fill = "grey") +
        stat_summary(fun = median, aes(ymin=..y.., ymax=..y..),
                     geom='errorbar', color = "white", size=2,
                     width = width, linetype = "solid") +
        geom_beeswarm(size = 2, cex = 3) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "red", size = 1) +
        coord_cartesian(ylim = c(-1, 3.25)) +
        geom_text(data = p4Blab_dat_c, label = "*",size=6,hjust = .5) +
        labs(x = "",
             y = y_lab) +
        theme_classic() +
        theme_Elm() +
        theme(axis.title.y = element_text(size = font_sz),
              axis.text.x = element_text(vjust = 0.5),
              legend.position = "none",
              plot.margin = margin(l = 0.5, r = 0.5, t = 0.25, b = 0.25, unit = "cm"))

# Figure 4B Statistics
# ---------------------

signs_dat <- seq_syl_x_sec_diffs %>%
  filter(!is.na(mean_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(mean_syll_per_sec_dif > 0, "positive", "negative"))

signs_dat_stats <- signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    num_positive = positive,
    total = positive + negative)

sign_test <- binom.test(signs_dat_stats$num_positive, signs_dat_stats$total, p = 0.5, alternative = "two.sided")

tidy(sign_test) %>%
  mutate(
    `Positive Changes` = signs_dat_stats$num_positive,
    `Negative Changes` = signs_dat_stats$num_negative
  ) %>%
  select(`Positive Changes`, `Negative Changes`, p.value) %>%
  kbl("pipe")


# ---- Figure 4C: Developmental change in syllables per second stability

y_lab <- expression(atop("Syllables per second SD" ~ Delta, "(10 months - 5 months)"))

p4C <- seq_syl_x_sec_diffs %>%
    ggplot(aes("", std_syll_per_sec_dif)) +
        geom_boxplot(size = 1, color = "grey", outlier.shape = NA,
                     fill = "grey", width = 0.5) +
        stat_summary(fun = median, aes(ymin=..y.., ymax=..y..),
                     geom='errorbar', color = "white", size=2,
                     width = width, linetype = "solid") +
        geom_beeswarm(size = 2, cex = 3) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "red", size = 1) +
        coord_cartesian(ylim = c(-1, 2.25)) +
        labs(x = "",
             y = y_lab) +
        theme_classic() +
        theme_Elm() +
        theme(axis.title.y = element_text(size = font_sz),
              axis.text.x = element_text(vjust = 0.5),
              legend.position = "none",
              plot.margin = margin(l = 0.5, r = 0.5, t = 0.25, b = 0.25, unit = "cm"))

# Figure 4C Statistics
# ---------------------

SD_signs_dat <- seq_syl_x_sec_diffs %>%
  filter(!is.na(std_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(std_syll_per_sec_dif > 0, "positive", "negative"))

SD_signs_dat_stats <- SD_signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    num_positive = positive,
    total = positive + negative)

SD_sign_test <- binom.test(SD_signs_dat_stats$num_positive, SD_signs_dat_stats$total,
                        p = 0.5, alternative = "two.sided")

tidy(SD_sign_test) %>%
  mutate(
    `Positive Changes` = SD_signs_dat_stats$num_positive,
    `Negative Changes` = SD_signs_dat_stats$num_negative
  ) %>%
  select(`Positive Changes`, `Negative Changes`, p.value) %>%
  kbl("pipe")

# ---- caregiver responses predict sequence compression

cg_resp_dat <- read_csv(here("Study 2", "data", "infant_cg_response_data.csv"))

inf_seq_count <- inf_syl_dat %>%
  mutate(seq_cat = as.factor(seq_cat)) %>%
  group_by(age, sub, seq_cat, .drop = FALSE) %>%
  summarize(n = n())

# developmental change in syllables per second

inf_syl_dat_dev <- inf_syl_dat %>%
  group_by(age, sub) %>%
  mutate(syll_per_sec = ifelse(seq_cat != "sequence", NA, syllable_total/duration)) %>% 
  filter(seq_cat == "sequence") %>% # only sequences
  summarize(Mean = mean(syll_per_sec, na.rm = TRUE),
            SD = sd(syll_per_sec, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(age), values_from = c(Mean, SD)) %>%
  drop_na(`Mean_5 months`) %>% # drop infants without sequences at 5 mo's
  mutate(syll_per_sec_dev_diff = `Mean_10 months` - `Mean_5 months`,
         seq_comp_SD_dev_diff =  `SD_10 months` - `SD_5 months`)


# zero sequence infants (28, 48, 72, 76, 85)
inf_seq_count %>%
  filter(n == 0) %>%
  kbl("pipe")

c0dat <- cg_resp_dat %>%
  filter(age == "5 months",
         seq_cat == "sequence") %>%
  group_by(sub) %>%
  summarize(prop_response = mean(elicited_response)) %>%
  left_join(inf_syl_dat_dev)

p4D <- c0dat %>%
  ggplot(aes(prop_response, syll_per_sec_dev_diff)) +
  geom_smooth(color = "red", fill = "red", method = lm,
              alpha = .2, size = 1.5) +
  geom_point(shape = 20, stroke = 2.5, size = .01) +
  coord_cartesian(ylim = c(-1, 3)) +
  scale_x_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.1)) +
  # geom_text(aes(x,y), label = "*",data = stats_lab_1,
  #           size=6, hjust = .75, vjust = .75) +
  theme_classic() +
  labs(#title="Feedback 40dph",
    y=str_wrap("Sequence compression from 5 to 10 months", 24),
    x= "Prop 5 month old sequences\nthat elicited caregiver response") + 
  theme_Elm()

# Figure 4D Statistics
# ---------------------

seq_cg_cor <- tidy(cor.test(c0dat$prop_response, c0dat$syll_per_sec_dev_diff, method = "spearman"))

seq_cg_cor %>% kbl("pipe")

# ---- Sequence compression predicts MCDI

mcdi_summary <- read_csv(here("Study 2", "data", "mcdi_summary.csv"))

seq_mcdi <- inf_syl_dat_dev %>%
  left_join(mcdi_summary, by = c("sub" = "subject_id"))

p4E <- seq_mcdi %>%
  ggplot(aes(x = syll_per_sec_dev_diff, y = `Words Understood Percentile-sex`)) +
  geom_smooth(color = "red", fill = "red", method = lm,
              alpha = .2, size = 1.5) +
  geom_point(shape = 21, stroke = 2, size = 1) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme_classic(12) +
  labs(
    x = "Sequence compression\nfrom 5 to 10 months",
    y = "Vocabulary %tile at 18 months\n(MCDI Comprehension)"
  ) +
  theme_Elm()

# Figure 4E Statistics

seq_mcdi_cor <- tidy(cor.test(seq_mcdi$syll_per_sec_dev_diff,
                              seq_mcdi$`Words Understood Percentile-sex`,
                              method = "spearman"))

seq_mcdi_cor %>% kbl("pipe")

# ---- combine panels

# Adjust margins for consistent spacing
p4A <- p4A + theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5))
p4B <- p4B + theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
                   axis.text.x = element_text(margin = margin(t = 5)))
p4C <- p4C + theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
                   axis.text.x = element_text(margin = margin(t = 5)))
p4D <- p4D + theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5))
p4E <- p4E + theme(plot.margin = margin(t = 5, r = 5, b = 20, l = 5))

# Align A (p4A) with B (p4B) vertically and their x-axes
aligned_p4A_p4B <- align_plots(p4A, p4B, align = "v", axis = "tblr")

# Align B (p4B) with E (p4E) vertically
aligned_p4B_p4E <- align_plots(aligned_p4A_p4B[[2]], p4E2, align = "v", axis = "lr")

# Combine p4B and p4C into a single plot
combined_p4B_p4C <- plot_grid(
  aligned_p4B_p4E[[1]],  p4C, 
  ncol = 2, labels = c("B", "C"), label_size = 11
)

# Combine the aligned p4B_p4C with p4E (E)
col2 <- plot_grid(
  combined_p4B_p4C, aligned_p4B_p4E[[2]],
  ncol = 1, labels = c("", "E"), label_size = 11
)

# Align p4A (A) with p4D (D)
aligned_p4A_p4D <- align_plots(
    aligned_p4A_p4B[[1]], p4D,
    align = "v", axis = "tblr"
)

# Combine aligned p4A with p4D into col1
col1 <- plot_grid(
  aligned_p4A_p4D[[1]], aligned_p4A_p4D[[2]],
  ncol = 1, labels = c("A", "D"), label_size = 11
)

# Combine col1 and col2 into the final layout
final_plot <- plot_grid(
  col1, col2, ncol = 2
)

final_plot
