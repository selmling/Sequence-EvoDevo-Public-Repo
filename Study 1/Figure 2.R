library("lme4")
library("here")
library("rsvg")
library("broom")
library("purrr")
library("magick")
library("emmeans")
library("cowplot")
library("tidyverse")
library("kableExtra")
library("reticulate")

# ---- helper functions

source_python(here("Study 1", "data_proc", "baserate_proc.py"))
source_python(here("Study 1", "data_proc", "exclude_proc.py"))

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

# ---- data import

alb_data <- read_csv(here("Study 1", "data", "prop_elicited_response_resp_type_points.csv"))

alb_session_durations <- read_csv(here("Study 1", "data", "session_durations.csv"))

# global threshold parameter

threshold <- 2

# wrangle

alb_data <- alb_data %>%
  inner_join(alb_session_durations, by = "sub")

alb_voc_data <- alb_data %>%
  mutate(dur = offset - onset) %>%
  filter(tier == "Infraphonology") %>%
  filter(cat != "Other")

# sample size

subjects <- alb_data %>%
  pull(sub) %>%
  unique()

length(subjects)

# ---- response rate analyses

# frequency per minute of the number of responses across types

freq_per_min_cg_response <- alb_voc_data %>%
    filter(tier == "Infraphonology" & cat != "Other") %>%
    group_by(sub) %>%
    summarise(
        # Total counts for each response type
        verbal_count = sum(response_modality == "Vocal", na.rm = TRUE),
        nonverbal_count = sum(response_modality == "NonVocal", na.rm = TRUE),
        multimodal_count = sum(response_modality == "Both", na.rm = TRUE),
        sess_dur = unique(total_session_duration)  # Session duration in seconds
    ) %>%
    mutate(
        # Calculate frequencies per minute
        verbal_per_minute = verbal_count / (sess_dur / 60),
        nonverbal_per_minute = nonverbal_count / (sess_dur / 60),
        multimodal_per_minute = multimodal_count / (sess_dur / 60)
    )

freq_per_min_response_sumstat <- freq_per_min_cg_response %>%
    summarise(
        mean_verbal_per_minute = mean(verbal_per_minute, na.rm = TRUE),
        sd_verbal_per_minute = sd(verbal_per_minute, na.rm = TRUE),
        mean_nonverbal_per_minute = mean(nonverbal_per_minute, na.rm = TRUE),
        sd_nonverbal_per_minute = sd(nonverbal_per_minute, na.rm = TRUE),
        mean_multimodal_per_minute = mean(multimodal_per_minute, na.rm = TRUE),
        sd_multimodal_per_minute = sd(multimodal_per_minute, na.rm = TRUE)
    )

# display response frequency per minute statistics

freq_per_min_response_sumstat %>% round(2) %>% kbl("pipe")

# verbal more frequent than other response types?

freq_per_min_model_data <- freq_per_min_cg_response %>%
    select(sub, verbal_per_minute, nonverbal_per_minute, multimodal_per_minute) %>%
    pivot_longer(
        cols = c(verbal_per_minute, nonverbal_per_minute, multimodal_per_minute),
        names_to = "response_type",
        values_to = "frequency_per_minute"
    ) %>%
    mutate(
        response_type = factor(response_type, levels = c("verbal_per_minute", "nonverbal_per_minute", "multimodal_per_minute"))
    )

# verbal vs nonverbal freq per min difference
lm01 <- freq_per_min_model_data %>%
    filter(response_type %in% c("verbal_per_minute", "nonverbal_per_minute")) %>%
    lmer(frequency_per_minute ~ response_type + (1|sub), data=., REML= FALSE)

summary(emmeans(lm01, "response_type", contr = "pairwise"), infer = TRUE)$contrast %>%
    mutate(`p.value` = ifelse(`p.value` < 0.0001, "< .0001")) %>%
    pull(`p.value`) %>%
    kbl("pipe")

# verbal vs multimodal freq per min difference
lm02 <- freq_per_min_model_data %>%
    filter(response_type %in% c("verbal_per_minute", "multimodal_per_minute")) %>%
    lmer(frequency_per_minute ~ response_type + (1|sub), data=., REML= FALSE)

summary(emmeans(lm02, "response_type", contr = "pairwise"), infer = TRUE)$contrast %>%
    mutate(`p.value` = ifelse(`p.value` < 0.0001, "< .0001")) %>%
    pull(`p.value`) %>%
    kbl("pipe")

# verbal responses

alb_verbal_sumstats <- alb_voc_data %>%
  group_by(sub, sequence) %>%
  summarize(total_vocalizations = n(),
            verbal_response_elicited = sum(response_modality == "Vocal", na.rm = TRUE),
            proportion_verbally_responded = verbal_response_elicited / total_vocalizations,
            .groups = "drop")

alb_verbal_avg_stats <- alb_verbal_sumstats %>%
  group_by(sequence) %>%
  summarize(mean = mean(proportion_verbally_responded, na.rm = TRUE),
            std = sd(proportion_verbally_responded, na.rm = TRUE),
            .groups = "drop")

alb_verbal_avg_stats

# difference between sequence and individual syllable CG verbal response rates
lm1 <- alb_verbal_sumstats %>%
    lmer(proportion_verbally_responded ~ sequence + (1|sub), data=., REML= FALSE)

kbl(summary(emmeans(lm1, "sequence", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# non-verbal responses

alb_nonverb_sumstats <- alb_voc_data %>%
  group_by(sub, sequence) %>%
  summarize(total_vocalizations = n(),
            nonverbal_response_elicited = sum(response_modality == "NonVocal", na.rm = TRUE),
            proportion_nonverbal_responded = nonverbal_response_elicited / total_vocalizations,
            .groups = "drop")

alb_nonverb_avg_stats <- alb_nonverb_sumstats %>%
  group_by(sequence) %>%
  summarize(mean = mean(proportion_nonverbal_responded, na.rm = TRUE),
            std = sd(proportion_nonverbal_responded, na.rm = TRUE),
            .groups = "drop")

alb_nonverb_avg_stats

# difference between sequence and individual syllable CG non-verbal response rates
lm2 <- alb_nonverb_sumstats %>%
    lmer(proportion_nonverbal_responded ~ sequence + (1|sub), data=., REML= FALSE)

kbl(summary(emmeans(lm2, "sequence", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# multimodal responses

alb_multimodal_sumstats <- alb_voc_data %>%
    group_by(sub, sequence) %>%
    summarize(total_vocalizations = n(),
              multimodal_response_elicited = sum(response_modality == "Both", na.rm = TRUE),
              proportion_multimodal_responded = multimodal_response_elicited / total_vocalizations,
              .groups = "drop")

alb_multimodal_avg_stats <- alb_multimodal_sumstats %>%
    group_by(sequence) %>%
    summarize(mean = mean(proportion_multimodal_responded, na.rm = TRUE),
                std = sd(proportion_multimodal_responded, na.rm = TRUE),
                .groups = "drop")

alb_multimodal_avg_stats

# difference between sequence and individual syllable CG multimodal response rates
lm3 <- alb_multimodal_sumstats %>%
    lmer(proportion_multimodal_responded ~ sequence + (1|sub), data=., REML= FALSE)

kbl(summary(emmeans(lm3, "sequence", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# vocal duration followup analysis

alb_voc_dur_stats <- alb_voc_data %>% 
  filter(tier == "Infraphonology", cat != "Other") %>%
  group_by(sub, elicited_response) %>%
  summarize(avg_duration = mean(dur, na.rm = TRUE),
            sd_duration = sd(dur, na.rm = TRUE),
            n = n(),
            .groups = "drop")

# duration summary
alb_voc_dur_stats %>% 
  group_by(elicited_response) %>%
  summarize(mean_duration = mean(avg_duration, na.rm = TRUE),
            sd_duration = sd(avg_duration, na.rm = TRUE),
            n = sum(n))

# duration by response elicitation
alb_voc_dur_stats %>%
  group_by(sub, elicited_response) %>%
  summarise(avg = mean(avg_duration, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = elicited_response, values_from = avg) %>%
  with(t.test(`0`, `1`, paired = TRUE)) %>%
  tidy()

# logistic regression no dur

lm_dur1 <- alb_voc_data %>% 
  filter(tier == "Infraphonology", cat != "Other",
         response_modality == "Vocal" | is.na(response_modality)) %>%
  mutate(elicited_response = as.numeric(elicited_response)) %>%
  glmer(elicited_response ~ sequence + (1 | sub),
        data = ., 
        family = binomial)

kbl(summary(emmeans(lm_dur1, "sequence", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# logistic regression w/ dur

lm_dur2 <- alb_voc_data %>% 
  filter(tier == "Infraphonology", cat != "Other",
         response_modality == "Vocal" | is.na(response_modality)) %>%
  mutate(elicited_response = as.numeric(elicited_response)) %>%
  glmer(elicited_response ~ sequence + dur + (1 | sub),
        data = ., 
        family = binomial)

kbl(summary(emmeans(lm_dur2, "sequence", contr = "pairwise"), infer = TRUE)$contrast, "pipe")

# ---- comparison to chance analyses

# wrangle

# is this the first time running the simultion routines?
first_run <- 0  # put '0' here if already run the sims, and it will load results RDS files

if (first_run == 1) {
    # number of sequences
    alb_seq_data <- alb_voc_data %>%
    filter(cat == "seq")

    seq_n <- alb_seq_data %>%
    count(cat) %>%
    filter(cat == "seq") %>%
    pull(n)

    seq_n

    # number of individual syllables
    alb_isosyl_data <- alb_voc_data %>%
    filter(sequence == 0)

    indsyll_n <- alb_isosyl_data %>%
        count(cat) %>%
        summarize(total = sum(n)) %>%
        pull(total)

    indsyll_n
} else {}

# caregiver verbal response dataset for deriving chance
alb_cg_verbal_data <- alb_data %>%
  filter(tier == "CGR Modality", cat == "Vocal")

# sequence verbal response elicitation rate expected by chance
if (first_run == 1) {
    alb_v_see <- create_expected_elic_rate(alb_cg_verbal_data, alb_seq_data, subjects, threshold, seq_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_v_see_mean <- alb_v_see %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_v_see_mean, here("Study 1", "data", "alb_v_see_mean.rds"))
} else {
    alb_v_see_mean <- readRDS(here("Study 1", "data", "alb_v_see_mean.rds"))
}

# individual syllable verbal response elicitation rate expected by chance
if (first_run == 1) {
    alb_v_isee <- create_expected_elic_rate(alb_cg_verbal_data, alb_isosyl_data, subjects, threshold, indsyll_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_v_isee_mean <- alb_v_isee %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_v_isee_mean, here("Study 1", "data", "alb_v_isee_mean.rds"))
} else {
    alb_v_isee_mean <- readRDS(here("Study 1", "data", "alb_v_isee_mean.rds"))
}

# caregiver non-verbal response dataset for deriving chance
alb_cg_nonverbal_data <- alb_data %>%
  filter(tier == "CGR Modality", cat == "NonVocal")

# sequence nonverbal response elicitation rate expected by chance
if (first_run == 1) {
    alb_nv_see <- create_expected_elic_rate(alb_cg_nonverbal_data, alb_seq_data, subjects, threshold, seq_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_nv_see_mean <- alb_nv_see %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_nv_see_mean, here("Study 1", "data", "alb_nv_see_mean.rds"))
} else {
    alb_nv_see_mean <- readRDS(here("Study 1", "data", "alb_nv_see_mean.rds"))
}

# individual syllable nonverbal response elicitation rate expected by chance
if (first_run == 1) {
    alb_nv_isee <- create_expected_elic_rate(alb_cg_nonverbal_data, alb_isosyl_data, subjects, threshold, indsyll_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_nv_isee_mean <- alb_nv_isee %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_nv_isee_mean, here("Study 1", "data", "alb_nv_isee_mean.rds"))
} else {
    alb_nv_isee_mean <- readRDS(here("Study 1", "data", "alb_nv_isee_mean.rds"))
}

# caregiver multimodal response dataset for deriving chance
alb_cg_multimodal_data <- alb_data %>%
  filter(tier == "CGR Modality", cat == "Both")

# sequence multimodal response elicitation rate expected by chance
if (first_run == 1) {
    alb_mm_see <- create_expected_elic_rate(alb_cg_multimodal_data, alb_seq_data, subjects, threshold, seq_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_mm_see_mean <- alb_mm_see %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_mm_see_mean, here("Study 1", "data", "alb_mm_see_mean.rds"))
} else {
    alb_mm_see_mean <- readRDS(here("Study 1", "data", "alb_mm_see_mean.rds"))
}

# individual syllable multimodal response elicitation rate expected by chance
if (first_run == 1) {
    alb_mm_isee <- create_expected_elic_rate(alb_cg_multimodal_data, alb_isosyl_data, subjects, threshold, indsyll_n) %>%
        mutate(expected_elicitation_rate = as.numeric(expected_elicitation_rate))
    alb_mm_isee_mean <- alb_mm_isee %>%
        summarize(mean_expected_elicitation_rate = mean(expected_elicitation_rate)) %>%
        pull()
    saveRDS(alb_mm_isee_mean, here("Study 1", "data", "alb_mm_isee_mean.rds"))
} else {
    alb_mm_isee_mean <- readRDS(here("Study 1", "data", "alb_mm_isee_mean.rds"))
}

# difference above chance?

alb_verbal_sumstats_diffs <- alb_verbal_sumstats %>% 
    select(sub,sequence,proportion_verbally_responded) %>%
    pivot_wider(names_from = sequence,
                values_from = c(proportion_verbally_responded)) %>%
    rename(sequence_response_rate = `1`,
           isolated_syllable_response_rate = `0`) %>%
    mutate(diff = sequence_response_rate - isolated_syllable_response_rate)

expected_v_diff <- alb_v_see_mean - alb_v_isee_mean

# individual syllable response rate above chance
alb_verbal_sumstats %>%
  filter(sequence == 0) %>%
  summarize(t_test = list(tidy(t.test(proportion_verbally_responded, mu = alb_v_isee_mean)))) %>%
  unnest(t_test)

# sequences response rate above chance
alb_verbal_sumstats %>%
  filter(sequence == 1) %>%
  summarize(t_test = list(tidy(t.test(proportion_verbally_responded, mu = alb_v_see_mean)))) %>%
  unnest(t_test)

# difference between observed and expected difference
tidy(t.test(alb_verbal_sumstats_diffs$diff, mu = expected_v_diff))

# nonverbal response difference from expected chance difference?

alb_nonverb_sumstats_diffs <- alb_nonverb_sumstats %>% 
    select(sub,sequence,proportion_nonverbal_responded) %>%
    pivot_wider(names_from = sequence,
                values_from = c(proportion_nonverbal_responded)) %>%
    rename(sequence_response_rate = `1`,
           isolated_syllable_response_rate = `0`) %>%
    mutate(diff = sequence_response_rate - isolated_syllable_response_rate)

expected_nv_diff <- alb_nv_see_mean - alb_nv_isee_mean

# individual syllable response rate above chance
alb_nonverb_sumstats %>%
  filter(sequence == 0) %>%
  summarize(t_test = list(tidy(t.test(proportion_nonverbal_responded, mu = alb_nv_isee_mean))) ) %>%
  unnest(t_test)

# sequences response rate above chance
alb_nonverb_sumstats %>%
  filter(sequence == 1) %>%
  summarize(t_test = list(tidy(t.test(proportion_nonverbal_responded, mu = alb_nv_see_mean))) ) %>%
  unnest(t_test)

# difference between observed and expected difference
tidy(t.test(alb_nonverb_sumstats_diffs$diff, mu = expected_nv_diff))

# multimodal response difference from expected chance difference?

alb_multimodal_sumstats_diffs <- alb_multimodal_sumstats %>% 
    select(sub,sequence,proportion_multimodal_responded) %>%
    pivot_wider(names_from = sequence,
                values_from = c(proportion_multimodal_responded)) %>%
    rename(sequence_response_rate = `1`,
           isolated_syllable_response_rate = `0`) %>%
    mutate(diff = sequence_response_rate - isolated_syllable_response_rate)

expected_mm_diff <- alb_mm_see_mean - alb_mm_isee_mean

# individual syllable response rate above chance
alb_multimodal_sumstats %>%
  filter(sequence == 0) %>%
  summarize(t_test = list(tidy(t.test(proportion_multimodal_responded, mu = alb_mm_isee_mean))) ) %>%
  unnest(t_test)

# sequences response rate above chance
alb_multimodal_sumstats %>%
  filter(sequence == 1) %>%
  summarize(t_test = list(tidy(t.test(proportion_multimodal_responded, mu = alb_mm_see_mean))) ) %>%
  unnest(t_test)

# difference between observed and expected difference
tidy(t.test(alb_multimodal_sumstats_diffs$diff, mu = expected_mm_diff))

# ---- visualize

# proportion inf voc type elicited verbal responses

alb_verbal_sumstats <- alb_verbal_sumstats %>%
  mutate(sequence = as.factor(sequence))

xlabs <- c("individual\nsyllable", "sequence")

vocal_label <- data.frame(proportion_verbally_responded=c(.01),sequence = c(1.5))

stat_label <- data.frame(proportion_verbally_responded=c(.9),sequence = c(1.5))

# ---- Figure 2A: Caregiver verbal response rate across infant voc types

p1 <- ggplot(alb_verbal_sumstats,aes(x = sequence, y = proportion_verbally_responded)) +
          stat_summary(fun.y=mean, geom="bar", colour="black",
                       fill="black", alpha = .15,
                       size=1, width = 0.35) + 
          geom_jitter(color="red", width = 0.01, alpha = .75) +
          # geom_jitter(color="black", width = 0.01, alpha = .75) +
          stat_summary(fun.data = mean_se, geom = "errorbar", size=1,
                       width=0.075) +
          geom_segment(aes(x=.65, xend=1.4, y=alb_v_isee_mean, yend=alb_v_isee_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
          geom_segment(aes(x=1.65, xend=2.4, y=alb_v_see_mean, yend=alb_v_see_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
          geom_text(data = stat_label,label = "***",size=7,color="black") +
          coord_cartesian(ylim=c(0, 1)) +
          labs(y = "proportion of vocalizations \n that elicited caregiver response",
               x = "",
               title = "Caregiver\nverbal response") +
          theme_classic() +
          scale_x_discrete(labels= xlabs) +
          theme_Elm() +
          theme(plot.title = element_text(hjust = 0.5))

# ---- Figure 2B: Caregiver non-verbal response rate across infant voc types

alb_nonverb_sumstats <- alb_nonverb_sumstats %>%
  mutate(sequence = as.factor(sequence))

xlabs <- c("individual\nsyllable", "sequence")

nonvoc_st_label <- data.frame(proportion_nonverbal_responded=c(.9),sequence = c(1.5),response_modality="NonVocal")

p2 <- ggplot(alb_nonverb_sumstats,aes(x = sequence, y = proportion_nonverbal_responded)) +
          stat_summary(fun.y=mean, geom="bar", colour="black",
                       fill="#8965FC", alpha = .15,
                       size=1, width = 0.35) + 
          geom_jitter(color="#8965FC", width = 0.01, alpha = .75) +
          stat_summary(fun.data = mean_se, geom = "errorbar", size=1,
                       width=0.075) +
          geom_segment(aes(x=.65, xend=1.4, y=alb_nv_isee_mean, yend=alb_nv_isee_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
          geom_segment(aes(x=1.65, xend=2.4, y=alb_nv_see_mean, yend=alb_nv_see_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
        #   geom_text(data = nonvoc_st_label,label = "*",size=7,color="black") +
          coord_cartesian(ylim=c(0, 1)) +
          labs(y = "",# str_wrap("Proportion of vocalizations which elicited caregiver response",36),
               x = "",
               title = "Caregiver\nnon-verbal response") +
          theme_classic() +
          scale_x_discrete(labels= xlabs) +
          theme_Elm() +
          theme(plot.title = element_text(hjust = 0.5))

# ---- Figure 2C: Caregiver multimodal response rate across infant voc types

alb_multimodal_sumstats <- alb_multimodal_sumstats %>%
  mutate(sequence = as.factor(sequence))

xlabs <- c("individual\nsyllable", "sequence")

multimodal_st_label <- data.frame(proportion_multimodal_responded=c(.9),sequence = c(1.5),response_modality="Both")

p3 <- ggplot(alb_multimodal_sumstats,aes(x = sequence, y = proportion_multimodal_responded)) +
          stat_summary(fun.y=mean, geom="bar", colour="black",
                       fill="#DB00B7", alpha = .15,
                       size=1, width = 0.35) + 
          geom_jitter(color="#DB00B7", width = 0.01, alpha = .75) +
          stat_summary(fun.data = mean_se, geom = "errorbar", size=1,
                       width=0.075) +
          geom_segment(aes(x=.65, xend=1.4, y=alb_mm_isee_mean, yend=alb_mm_isee_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
          geom_segment(aes(x=1.65, xend=2.4, y=alb_mm_see_mean, yend=alb_mm_see_mean),
                         col="black",
                         linetype="dotted",
                         size=1) +
        #   geom_text(data = multimodal_st_label,label = "*",size=7,color="black") +
          coord_cartesian(ylim=c(0, 1)) +
          labs(y = "",# str_wrap("Proportion of vocalizations which elicited caregiver response",36),
               x = "",
               title = "Caregiver\nmultimodal response") +
          theme_classic() +
          scale_x_discrete(labels= xlabs) +
          theme_Elm() +
          theme(plot.title = element_text(hjust = 0.5))

# ---- Figure 2D: Verbal response expected difference compared to chance difference

stat_label <- data.frame(diff=c(.6))

p4 <- ggplot(alb_verbal_sumstats_diffs, aes(x = "", y = diff)) +
        geom_jitter(color="red", width = 0.01, alpha = .75) +
     #    geom_jitter(color="black", width = 0.01, alpha = .75) +
        stat_summary(fun.y=mean,geom="point",shape=19,size=3) + 
        stat_summary(fun.data = mean_se, geom = "errorbar", size=1.5, width=0.00) +
          geom_segment(aes(x=.65, xend=1.38, y=expected_v_diff, yend=expected_v_diff),
                         col="black",
                         linetype="dotted",
                         size=1) +
          coord_cartesian(ylim=c(-.75, .75)) +
          geom_text(data = stat_label, label = "*", size=7, color="black") +
          labs(y = str_wrap("elicitation rate difference (sequence - individual syllable)",36),
               x = str_wrap("verbal response",12)) +
          theme_classic() +
          theme_Elm()

# ---- Figure 2E: NonVerbal response expected difference compared to chance difference

p5 <- ggplot(alb_nonverb_sumstats_diffs, aes(x = "", y = diff)) +
        geom_jitter(color="#8965FC", width = 0.01, alpha = .75) +
        stat_summary(fun.y=mean,geom="point",shape=19,size=3) +
        stat_summary(fun.data = mean_se, geom = "errorbar", size=1.5, width=0.00) +
          geom_segment(aes(x=.65, xend=1.38, y=expected_nv_diff, yend=expected_nv_diff),
                         col="black",
                         linetype="dotted",
                         size=1) +
          coord_cartesian(ylim=c(-.75, .75)) +
        #   geom_text(data = stat_label,label = "^",size=4,color="black") +
          labs(y = "",#str_wrap("elicitation rate difference (sequence - isolated syllable)",36),
               x = str_wrap("non-verbal response",12)) +
          theme_classic() +
          theme_Elm()

# ---- Figure 2F: Multimodal response expected difference compared to chance difference

stat_label <- data.frame(diff=c(-.6))

p6 <- ggplot(alb_multimodal_sumstats_diffs, aes(x = "", y = diff)) +
        geom_jitter(color="#DB00B7", width = 0.01, alpha = .75) +
        stat_summary(fun.y=mean,geom="point",shape=19,size=3) +
        stat_summary(fun.data = mean_se, geom = "errorbar", size=1.5, width=0.00) +
          geom_segment(aes(x=.65, xend=1.38, y=expected_mm_diff, yend=expected_mm_diff),
                         col="black",
                         linetype="dotted",
                         size=1) +
          coord_cartesian(ylim=c(-.75, .75)) +
          geom_text(data = stat_label,label = "**",size=7,color="black") +
          labs(y = "",#str_wrap("elicitation rate difference (sequence - isolated syllable)",36),
               x = str_wrap("multimodal response",12)) +
          theme_classic() +
          theme_Elm()

# ---- Figure 2

col1 <- plot_grid(p1, p4, labels = c('A', 'D'), align = "v", nrow=2)

col2 <- plot_grid(p2, p5, labels = c('B', 'E'), align = "v", nrow=2)

col3 <- plot_grid(p3, p6, labels = c('C', 'F'), align = "v", nrow=2)

legend_path <- here("Sequence-Response","figures","legend_1.svg")

legend <- ggdraw() + draw_image(magick::image_read_svg(legend_path))

plot_grid(col1, col2, col3, legend, align = "hv", ncol=4, rel_widths = c(1, 1, 1, .6))
