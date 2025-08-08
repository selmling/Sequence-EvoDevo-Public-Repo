library("here")
library("lme4")
library("broom")
library("readxl")
library("scales")
library("cowplot")
library("emmeans")
library("tidytext")
library("tidyverse")
library("kableExtra")
library("ggbeeswarm")

# ---- helper functions

font_sz <- 12

theme_Elm <- function(){
  theme(
    text = element_text(size=font_sz),
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.ticks.length = unit(-3.5, "pt"),
    legend.position = "none"
  )
}

custom_format <- function(x) {
  x %>%
    as.numeric() %>% # Ensure numeric input
    formatC(format = "f", drop0trailing = TRUE) %>% # Drop trailing zeros
    gsub(pattern = "^0\\.", replacement = ".", .) %>% # Remove leading zero for positive numbers
    gsub(pattern = "^-0\\.", replacement = "-.", .) # Remove leading zero for negative numbers
}

RmType <- function(string) {
  sub("._", "", string)
}

# ---- wrangle

zb_syl_dat <- read_csv(here("Study 3", "data", "zb_seg_dat_counts.csv"))

wrong_names <- c("Y_O", "Y", "R_DBL")

zb_syl_dat <- zb_syl_dat %>%
  mutate(bird_ID = ifelse(bird_ID %in% wrong_names,
                                       paste0(bird_ID, "-"),bird_ID))

zb_syl_dat_dev <- zb_syl_dat %>%
  group_by(age, condition, bird_ID, Brother_ID) %>%
  summarize(Mean = mean(syll_per_sec),
            SD = sd(syll_per_sec),
            n = n()) %>%
  pivot_wider(names_from = c(age), values_from = c(Mean, SD)) %>%
  mutate(syll_per_sec_dev_diff = Mean_ninety - `Mean_forty-five`,
         seq_comp_SD_dev_diff = SD_ninety- `SD_forty-five`) %>%
  select(condition, bird_ID, Brother_ID, n, Mean_ninety, `Mean_forty-five`,
         syll_per_sec_dev_diff, SD_ninety, `SD_forty-five`, seq_comp_SD_dev_diff)

# ---- visualize

font_sz <- 12

zb_syl_dat_sum <- zb_syl_dat %>%
    mutate(age = recode(age, `forty-five` = "45 dph", `ninety` = "90 dph")) %>%
    group_by(condition, age) %>%
    summarise(Mean = mean(motif_duration))

zb_syl_dat <- zb_syl_dat %>%
  mutate(age = recode(age, `forty-five` = "45 dph", `ninety` = "90 dph"))

zb_stats_lab <-  data.frame(species = "Zebra finch", mean_syll_per_sec = 20, age = 1.5)

seq_syl_x_sec_sumstats <- zb_syl_dat %>%
    group_by(bird_ID, condition, Brother_ID, age) %>%
    summarize(mean_syll_per_sec = mean(syll_per_sec),
              std_syll_per_sec = sd(syll_per_sec))

seq_syl_x_sec_diffs <- pivot_wider(seq_syl_x_sec_sumstats,
                                id_cols = c("bird_ID", "condition", "Brother_ID"),
                                names_from = "age",
                                values_from = c("mean_syll_per_sec", "std_syll_per_sec")) %>%
                    mutate(mean_syll_per_sec_dif = `mean_syll_per_sec_90 dph` - `mean_syll_per_sec_45 dph`,
                           std_syll_per_sec_dif = `std_syll_per_sec_90 dph` - `std_syll_per_sec_45 dph`) %>% 
                    select(bird_ID, condition, mean_syll_per_sec_dif, std_syll_per_sec_dif)

zb_stats_lab <- data.frame(
  condition = "Contingent",
  mean_syll_per_sec = 20,
  age = 1.5
)

# ---- Figure 5A: Zebra finch developmental change in syllables per second 

fig5a <- seq_syl_x_sec_sumstats %>%
  ggplot(aes(x = age, y = mean_syll_per_sec)) +
  stat_boxplot(geom = "errorbar", size = 1, width = 0.35) +
  geom_boxplot(size = 1, outlier.shape = 18, width = 0.75,
                 outlier.size = 2, outlier.alpha = .75) +
  geom_line(aes(group = bird_ID), position = "identity",
            alpha = .15, size = 1) +
  geom_point(size = 3, alpha = .15, stroke = 0, shape = 19) +
  geom_text(data = zb_stats_lab, aes(label = "**", y = mean_syll_per_sec), size = 6, hjust = 0.5) +
  coord_cartesian(ylim = c(4, 20)) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  labs(y = "Average syllables per second", x = "") +
  facet_wrap(~condition) + # Hide facet labels
  theme_classic() +
  theme_Elm()

# Figure 5A Statistics
# ---------------------

# Contingent
lm1AL <- seq_syl_x_sec_sumstats %>%
       mutate(age = factor(age, levels = c("90 dph", "45 dph"))) %>%
       filter(condition == "Contingent") %>%
       lmer(mean_syll_per_sec ~ age + (1 | bird_ID), data = ., REML= FALSE)
       
kbl(summary(emmeans(lm1AL, "age", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# Yoked
lm1AR <- seq_syl_x_sec_sumstats %>%
       mutate(age = factor(age, levels = c("90 dph", "45 dph"))) %>%
       filter(condition == "Yoked") %>%
       lmer(mean_syll_per_sec ~ age + (1 | bird_ID), data = ., REML= FALSE)
       
kbl(summary(emmeans(lm1AR, "age", contr = "pairwise"), infer = TRUE)$contrast, "pipe")


# ---- Figure 5B: Zebra finch developmental change in syllables per second

y_lab <- expression(atop("Syllables per second" ~ Delta, "(90 - 45 dph)"))

stat_label <- data.frame(mean_syll_per_sec_dif=c(8),condition = "Contingent")

fig5b <- seq_syl_x_sec_diffs %>%
    ggplot(aes(condition, mean_syll_per_sec_dif, color = condition)) +
        geom_boxplot(size = 1, color = "grey", outlier.shape = NA,
                     fill = "grey", width = 0.5) +
        stat_summary(fun = median, aes(ymin=..y.., ymax=..y..),
                     geom='errorbar', color = "white", size=2,
                     width = .5, linetype = "solid") +
        geom_beeswarm(size = 2, cex = 3)  +
        scale_color_manual(values = c(Contingent = "#45ABF5", Yoked = "#F52F2C")) +
        scale_x_discrete(labels = c("Contingent" = "C", "Yoked" = "YC")) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "black", size = 1) +
        geom_text(data = stat_label,label = "*",size=6,color="black") +
        labs(x = "",
             y = y_lab) +
        theme_classic() +
        theme_Elm() +
        theme(axis.title.y = element_text(size = font_sz),
              axis.text.x = element_text(vjust = 0.5),
              legend.position = "none")

# Figure 5B Statistics
# ---------------------

# Contingent
contingent_signs_dat <- seq_syl_x_sec_diffs %>%
  filter(condition == "Contingent") %>%
  filter(!is.na(mean_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(mean_syll_per_sec_dif > 0, "positive", "negative"))

c_signs_dat_stats <- contingent_signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    num_positive = positive,
    total = positive + negative) %>%
  mutate(condition = "Contingent") %>%
  select(condition, num_positive, total)

contingent_wilcox <- seq_syl_x_sec_sumstats %>%
  filter(condition == "Contingent") %>%
  group_by(condition) %>%
  do(wilcox.test(.$mean_syll_per_sec[.$age == "45 dph"], 
                 .$mean_syll_per_sec[.$age == "90 dph"], 
                 paired = TRUE) %>% 
       tidy()) %>%
  select(condition, p.value)

c_signs_dat_stats %>%
  left_join(contingent_wilcox, by = "condition") %>%
  kbl("pipe")

# Yoked
yoked_signs_dat <- seq_syl_x_sec_diffs %>%
  filter(condition == "Yoked") %>%
  filter(!is.na(mean_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(mean_syll_per_sec_dif > 0, "positive", "negative"))

yc_signs_dat_stats <- yoked_signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    num_positive = positive,
    total = positive + negative) %>%
  mutate(condition = "Yoked") %>%
  select(condition, num_positive, total)

# drop birds that don't have both ages
yoked_birds <- seq_syl_x_sec_sumstats %>%
  group_by(bird_ID) %>%
  summarise(has_both_ages = all(c("45 dph", "90 dph") %in% age), .groups = "drop") %>%
  filter(has_both_ages) %>% 
  pull(bird_ID)

yoked_wilcox <- seq_syl_x_sec_sumstats %>%
  filter(condition == "Yoked", bird_ID %in% yoked_birds) %>%
  group_by(condition) %>%
  do(wilcox.test(.$mean_syll_per_sec[.$age == "45 dph"], 
                 .$mean_syll_per_sec[.$age == "90 dph"], 
                 paired = TRUE) %>% 
       tidy()) %>%
  select(condition, p.value)

yc_signs_dat_stats %>%
  left_join(yoked_wilcox, by = "condition") %>%
  kbl("pipe")

# ---- Figure 5C: Zebra finch developmental change in stability in syllables per second

y_lab <- expression(atop("Syllables per second SD" ~ Delta, "(90 - 45 dph)"))

stat_label <- data.frame(std_syll_per_sec_dif=c(-1.5), condition="Contingent")

fig5c <- seq_syl_x_sec_diffs %>%
    ggplot(aes(condition, std_syll_per_sec_dif, color = condition)) +
        geom_boxplot(size = 1, color = "grey", outlier.shape = NA,
                     fill = "grey", width = 0.5) +
        stat_summary(fun = median, aes(ymin=..y.., ymax=..y..),
                     geom='errorbar', color = "white", size=2,
                     width = 0.5, linetype = "solid") +
        geom_beeswarm(size = 2, cex = 3)  +
        scale_color_manual(values = c(Contingent = "#45ABF5", Yoked = "#F52F2C")) +
        scale_x_discrete(labels = c("Contingent" = "C", "Yoked" = "YC")) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "black", size = 1) +
        geom_text(data = stat_label,label = "**",size=6,color="black") +
        labs(x = "",
             y = y_lab) +
        theme_classic() +
        theme_Elm() +
        theme(axis.title.y = element_text(size = font_sz),
              axis.text.x = element_text(vjust = 0.5),
              legend.position = "none")

# Figure 5C Statistics
# ---------------------

# Contingent
contingent_sd_signs_dat <- seq_syl_x_sec_diffs %>%
  filter(condition == "Contingent") %>%
  filter(!is.na(std_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(std_syll_per_sec_dif > 0, "positive", "negative"))

c_sd_signs_dat_stats <- contingent_sd_signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    # num_positive = positive,
    total = negative) %>%
  mutate(condition = "Contingent") %>%
  select(condition, num_negative, total)

contingent_sd_wilcox <- seq_syl_x_sec_sumstats %>%
  filter(condition == "Contingent") %>%
  group_by(condition) %>%
  do(wilcox.test(.$std_syll_per_sec[.$age == "45 dph"], 
                 .$std_syll_per_sec[.$age == "90 dph"], 
                 paired = TRUE) %>% 
       tidy()) %>%
  select(condition, p.value)

c_sd_signs_dat_stats %>%
  left_join(contingent_sd_wilcox, by = "condition") %>%
  kbl("pipe")

# Yoked
yoked_sd_signs_dat <- seq_syl_x_sec_diffs %>%
  filter(condition == "Yoked") %>%
  filter(!is.na(std_syll_per_sec_dif)) %>%
  mutate(direction = ifelse(std_syll_per_sec_dif > 0, "positive", "negative"))

yc_sd_signs_dat_stats <- yoked_sd_signs_dat %>%
  group_by(direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  summarize(num_negative = negative,
    num_positive = positive,
    total = positive + negative) %>%
  mutate(condition = "Yoked") %>%
  select(condition, num_positive, total)

# drop birds that don't have both ages
yoked_birds <- seq_syl_x_sec_sumstats %>%
  group_by(bird_ID) %>%
  summarise(has_both_ages = all(c("45 dph", "90 dph") %in% age), .groups = "drop") %>%
  filter(has_both_ages) %>% 
  pull(bird_ID)

yoked_sd_wilcox <- seq_syl_x_sec_sumstats %>%
  filter(condition == "Yoked", bird_ID %in% yoked_birds) %>%
  group_by(condition) %>%
  do(wilcox.test(.$std_syll_per_sec[.$age == "45 dph"], 
                 .$std_syll_per_sec[.$age == "90 dph"], 
                 paired = TRUE) %>% 
       tidy()) %>%
  select(condition, p.value)

yc_sd_signs_dat_stats %>%
  left_join(yoked_sd_wilcox, by = "condition") %>%
  kbl("pipe")

# ---- Predicting sequence compression

pback_dat <- read_xlsx(here("Study 3", "data", "BySubject SPSS Behavioral Data FINAL excel.xlsx"))

sub_map <- read_xlsx(here("Study 3", "data","subID_map.xlsx"))

sub_map <- sub_map %>%
  rename(MaleNum = blind_ID)

# bring in subjectID into pback data
pback_dat <- pback_dat %>% left_join(sub_map)

# summary
seg_dat_sum_stats <- zb_syl_dat %>%
  group_by(condition, bird_ID, Brother_ID, age) %>%
  summarize(mean_syll_per_sec = mean(syll_per_sec),
            median_syll_per_sec = median(syll_per_sec),
            SD_syll_per_sec = sd(syll_per_sec),
            SE_syll_per_sec = sd(syll_per_sec)/sqrt(n())) %>%
  arrange(condition, Brother_ID)

# join average response
Avg_ALLContProp <- pback_dat %>%
  select(bird_ID, condition, Avg_ALLContProp)

Avg_Overall_Similarity <- pback_dat %>%
  select(bird_ID, condition, Overall_Similarity)

seg_dat_sum_stats <- seg_dat_sum_stats %>% left_join(Avg_ALLContProp)

seg_dat_sum_stats <- seg_dat_sum_stats %>% left_join(Avg_Overall_Similarity)

zb_syl_dat_dev <- zb_syl_dat_dev %>% left_join(Avg_ALLContProp)

zb_syl_dat_dev <- zb_syl_dat_dev %>% left_join(Avg_Overall_Similarity)

# ---- Figure 5E: Predicting overall similarity to tutor from sequence rate

fig5e <- zb_syl_dat_dev %>%
  ggplot(aes(syll_per_sec_dev_diff, Overall_Similarity)) +
  geom_point(aes(color=condition), shape = 21, stroke = 2, size = 1) + 
  geom_smooth(aes(color=condition,fill=condition),method=lm) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  coord_cartesian(xlim = c(-5,10)) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  facet_grid(.~condition) +
  theme_classic() +
  labs(title="Vocal maturity",
       y="Overall song similarity",
       x=str_wrap("Sequence compression development (90 - 45 dph)",32)) +
  theme(legend.position = "none")


# Figure 5E Statistics
# ---------------------

cor_meth <- "spearman"

zb_syl_dat_dev %>%
  filter(condition == "Contingent") %>%
  ungroup() %>%
  do(tidy(cor.test(.$syll_per_sec_dev_diff,
                   .$Overall_Similarity,
                   method = cor_meth))) %>%
  mutate(condition = "Contingent") %>%
  kbl("pipe")

zb_syl_dat_dev %>%
  filter(condition == "Yoked") %>%
  ungroup() %>%
  do(tidy(cor.test(.$syll_per_sec_dev_diff,
                   .$Overall_Similarity,
                   method = cor_meth))) %>%
  mutate(condition = "Yoked") %>%
  kbl("pipe")

# ---- Figure 5D: predicting sequence rate from responses at developmental timepoints

pback_dat <- pback_dat %>%
  mutate(T5_ContProp = coalesce(T5_ContProp, T5_AccContVidProp),
         T10_ContProp = coalesce(T10_ContProp, T10_AccContVidProp),
         T15_ContProp = coalesce(T15_ContProp, T15_AccContVidProp),
         T20_ContProp = coalesce(T20_ContProp, T20_AccContVidProp),
         T25_ContProp = coalesce(T25_ContProp, T25_AccContVidProp))

ContProp <- pback_dat %>%
  select(bird_ID, condition, T5_ContProp, T10_ContProp, T15_ContProp, T20_ContProp, T25_ContProp)

zb_syl_dat_dev <- zb_syl_dat_dev %>% left_join(ContProp)

fig5d <- zb_syl_dat_dev %>%
  # filter(bird_ID != "Mg") %>% # based on boxplot, `Mg` is outlier in avg responses
  ggplot(aes(T5_ContProp, syll_per_sec_dev_diff)) +
  geom_point(aes(color=condition), shape = 20, stroke = 3) + 
  geom_smooth(aes(color=condition,fill=condition),method=lm) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  coord_cartesian(ylim = c(-5,10)) +
  scale_x_continuous(limits=c(0,1),labels = label_number(accuracy = 0.1)) +
  facet_grid(.~condition) +
  theme_classic() +
  labs(title="Feedback 40dph",
       y=str_wrap("Sequence compression development (90 - 45 dph)",32),
       x=str_wrap("Prop. of songs received response at 40 dph",28)) +
  theme(legend.title = element_blank(),
        legend.position = "none")


# Figure 5D Statistics
# ---------------------

cor_meth <- "spearman"

zb_syl_dat_dev %>%
  filter(condition == "Contingent") %>%
  ungroup() %>%
  do(tidy(cor.test(.$T5_ContProp,
                   .$syll_per_sec_dev_diff,
                   method = cor_meth))) %>%
  mutate(condition = "Contingent") %>%
  kbl("pipe")

zb_syl_dat_dev %>%
  filter(condition == "Yoked") %>%
  ungroup() %>%
  do(tidy(cor.test(.$T5_ContProp,
                   .$syll_per_sec_dev_diff,
                   method = cor_meth))) %>%
  mutate(condition = "Yoked") %>%
  kbl("pipe")

# ---- Figure 5: combined plot

# align fig5a (A) and (D) vertically
aligned_fig5a_fig5d <- align_plots(fig5a, fig5d, align = "hv", axis = "tblr")

# align fig5b (C) and fig5e (E) vertically
aligned_fig5b_fig5e <- align_plots(fig5b, fig5e, align = "v", axis = "tblr")

# Combine fig5b and fig5c into a single plot
aligned_fig5b_fig5c <- plot_grid(
  aligned_fig5b_fig5e[[1]],  fig5c, 
  ncol = 2, labels = c("B", "C"), label_size = 11
)

# Combine the aligned fig5c_fig5e with fig5d (E)
col2 <- plot_grid(
  aligned_fig5b_fig5c, aligned_fig5b_fig5e[[2]],
  ncol = 1, labels = c("", "E"), label_size = 11
)

col1 <- plot_grid(
  fig5a, aligned_fig5a_fig5d[[2]], align = "v", axis = "tblr",
  ncol = 1, labels = c("A", "D"), label_size = 11
)

plot_grid(col1, col2, ncol = 2)