library("here")
library("Rcpp")
library("lme4")
library("broom")
library("emmeans")
library("cowplot")
library("tidyverse")
library("kableExtra")

# ---- helper functions

sourceCpp(here("Study 1", "data_proc", "temporal_profile.cpp"))

# ---- wrangle

alb_dat <- read_csv(here("Study 1", "data","prop_elicited_response_resp_type_points.csv")) %>%
    select(-c(...1, `Unnamed: 0`))

# event counts

alb_dat <- alb_dat %>% 
    mutate(sub = as.factor(sub),
           category = factor(case_when(
                cat %in% c("seq") ~ "seq",
                cat %in% c("FRV", "CS", "MSFR", "QRV", "CRV", "MSQR") ~ "nonseq",
                cat %in% c("Vocal") ~ "Vocal",
                cat %in% c("NonVocal") ~ "NonVocal",
                cat %in% c("Both") ~ "Both",
                TRUE ~ NA_character_)))

event_counts <- alb_dat %>%
    group_by(category, sub, .drop = FALSE) %>%
    summarise(n = n())

# excluder dyads

Both_excluders <- event_counts %>%
    filter(n == 0 & category == "Both") %>%
    pull(sub) %>%
    unlist()

NonVocal_excluders <- event_counts %>%
    filter(n == 0 & category == "NonVocal") %>%
    pull(sub) %>%
    unlist()

inf_excluders <- event_counts %>% 
    filter(n < 5 & category == "seq" |
           n < 5 & category == "nonseq") %>%
    pull(sub) %>%
    unlist()

# drop "Other" vocalizations

alb_dat <- alb_dat %>%
    filter(cat != "Other") # %>% 
    # filter(!sub %in% inf_excluders)

# drop caregiver behaviors that are not in response to infant vocalizations

alb_dat <- alb_dat %>%
  filter(!is.na(category)) %>%
  arrange(sub, onset) %>%
  mutate(response_eliciting_vocalization_type = if_else(
    tier == "CGR Modality" & lag(elicited_response, default = 0) == 1,
    lag(category),
    NA_character_
  ))

# number of subjects
alb_dat %>%
  ungroup() %>%
  summarize(subject_count = n_distinct(sub)) %>%
  kbl("pipe")

alb_dat <- alb_dat %>%
    filter(!(tier == "CGR Modality" & is.na(response_eliciting_vocalization_type)))

# proportion of caregiver verbal response to vocal types
alb_dat %>% filter(tier == "CGR Modality") %>%
  filter(cat == "Vocal") %>%
  mutate(value = ifelse(is.na(response_eliciting_vocalization_type), "NA", as.character(response_eliciting_vocalization_type))) %>%
  count(value) %>%
  mutate(proportion = n / sum(n))

# average duration of infant vocalizations in frames
alb_dat %>% 
    filter(tier == "Infraphonology") %>% 
    group_by(sub) %>% 
    summarize(mean_dur = mean(dur, na.rm = TRUE)) %>% 
    summarize(grand_mean = mean(mean_dur)) %>% 
    mutate(avg_inf_voc_dur_frames = grand_mean * 29.97)

# average duration of caregiver verbal response in frames
alb_dat %>% 
    filter(category == "Vocal") %>% 
    group_by(sub) %>% 
    summarize(mean_dur = mean(dur, na.rm = TRUE)) %>% 
    summarize(grand_mean = mean(mean_dur)) %>% 
    mutate(avg_cg_verb_response_frames = grand_mean * 29.97)

# ---- temporal profiles

# sequences after caregiver verbal responses

seq_df <- alb_dat %>%
    filter(category == "seq" | category == "Vocal") %>%
    group_by(sub) %>% 
    nest()

seq_temp_profs_ons <- seq_df %>% 
  mutate(temporal_profile = map(data, ~generateTemporalProfile(.x, "onset", "Vocal", -150, 30, "seq")))

seq_temp_profs_offs <- seq_df %>% 
  mutate(temporal_profile = map(data, ~generateTemporalProfile(.x, "offset", "Vocal", -30, 150, "seq")))

# helper
temp_prof_sumstat <- function(df) {
  df %>%
    select(sub, temporal_profile) %>%
    unnest(cols = c(temporal_profile)) %>% 
    group_by(sub, frame_index) %>% 
    summarize(prop_response = mean(target_presence))
}

seq_temp_profs_ons_sumstat <- seq_temp_profs_ons %>%
    temp_prof_sumstat() %>% 
    mutate(sequence = "vocal sequence")

seq_temp_profs_offs_sumstat <- seq_temp_profs_offs %>%
    temp_prof_sumstat() %>% 
    mutate(sequence = "vocal sequence")

# individual syllables after caregiver verbal responses

nonseq_df <- alb_dat %>%
    filter(category == "nonseq" | category == "Vocal") %>%
    group_by(sub) %>% 
    nest()

nonseq_temp_profs_ons <- nonseq_df %>% 
  mutate(temporal_profile = map(data, ~generateTemporalProfile(.x, "onset", "Vocal", -150, 30, "nonseq")))

nonseq_temp_profs_offs <- nonseq_df %>% 
  mutate(temporal_profile = map(data, ~generateTemporalProfile(.x, "offset", "Vocal", -30, 150, "nonseq")))

# unnest and take average over time

nonseq_temp_profs_ons_sumstat <- nonseq_temp_profs_ons %>%
    temp_prof_sumstat() %>% 
    mutate(sequence = "individual syllable")

nonseq_temp_profs_offs_sumstat <- nonseq_temp_profs_offs %>%
    temp_prof_sumstat() %>% 
    mutate(sequence = "individual syllable")

# combine

ons_list <- list(seq_temp_profs_ons_sumstat, nonseq_temp_profs_ons_sumstat)

offs_list <- list(seq_temp_profs_offs_sumstat, nonseq_temp_profs_offs_sumstat)

temp_profs_ons_vocal <- ons_list %>%
    reduce(rbind) %>% 
    mutate(seconds = round(frame_index * (1/30), 2))

temp_profs_offs_vocal <- offs_list %>%
    reduce(rbind) %>% 
    mutate(seconds = round(frame_index * (1/30), 2))

# grand mean

mean_line_vocal <- list(temp_profs_ons_vocal, temp_profs_offs_vocal) %>% 
    reduce(rbind) %>%
    ungroup() %>% 
    summarize(grand_mean = mean(prop_response)) %>% 
    pull()

# visualize

# grand mean post caregiver verbal response

mean_post_line_vocal <- temp_profs_offs_vocal %>%
    ungroup() %>% 
    summarize(grand_mean = mean(prop_response)) %>% 
    pull()

p1 <- ggplot(temp_profs_offs_vocal, aes(x=seconds, y=prop_response, group=sequence)) +
        geom_smooth(size=1.25, aes(linetype = sequence), color = "black", se = TRUE, level = 0.95) +
        scale_linetype_manual(values=c("solid", "dashed"),
                              labels = c("individual\nsyllable", "sequence"),
                              guide = guide_legend(reverse = TRUE, keywidth = 3)) +
        coord_cartesian(ylim = c(0.0, 0.3),
                        xlim = c(-1, 5), expand = TRUE) +
        scale_x_continuous(breaks = 0:5) +
        geom_vline(xintercept = 0, size=1) +
        # geom_hline(yintercept = mean_post_line_vocal, size=.75, linetype= 'dotted') +
        theme_classic() +
        labs(y = "% frames with\ninfant vocalization",
             x = "Time (s)\nfrom caregiver verbal response offset") +
        annotate("rect", xmin = -1, xmax = 0,
                 ymin = -Inf, ymax = Inf,
                 fill = "gray", alpha=.25) +
        theme(text = element_text(size=font_sz),
               axis.text.x = element_text(vjust = 0.5, hjust=.5, colour="black"),
               axis.text.y = element_text(colour="black"),
               axis.ticks.length = unit(-3.5, "pt"),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               legend.position = c(0.75, 0.75),
               plot.title = element_text(hjust = 0.5),
               legend.title=element_blank(),
               legend.key.size = unit(.45, 'cm'),
               legend.text = element_text(size=font_sz),
               legend.key.height = unit(.75, "cm"))

p1

binned_data <- temp_profs_offs_vocal %>%
  mutate(second_bin = floor(seconds),
         second_bin = ifelse(second_bin == 5, 4, second_bin)) %>%  # Bin into 1-second intervals
  group_by(second_bin, sequence, sub) %>%
  summarize(mean_prop_response = mean(prop_response), .groups = 'drop')

model_results <- binned_data %>%
  group_by(second_bin) %>%
  do({
    # Fit model for each bin
    model <- lm(mean_prop_response ~ sequence, data = .)
    tidy(model)  # Extract results in a tidy format
  }) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

model_results <- model_results %>%
  mutate(p_adj = p.adjust(p.value, method = "bonferroni"))

model_results %>%
    select(second_bin,estimate,std.error,statistic,p_adj) %>%
    mutate(estimate = round(estimate, 4),
           std.error = round(std.error, 4),
           statistic = round(statistic, 4),
           p_adj = ifelse(p_adj < 0.0001, "<.0001", round(p_adj, 4))) %>%
    kbl("pipe")

# ---- proportion of caregiver verbal responses that were followed by infant vocal sequences vs individual syllables

library("reticulate")

source_python(here("Study 1", "data_proc", "resp_proc.py"))

threshold <- 2
subjects <- unique(alb_dat$sub)

alb_resp_dat <- alb_dat %>% 
  filter(category != "Both" & category != "NonVocal") %>%
  as.data.frame() %>%  # Convert to data frame for compatibility with Python
  mutate(caregiver_response_elicited_vocalization = prop_inf_response(., threshold, subjects)$caregiver_response_elicited_vocalization) %>%
  as_tibble()

alb_resp_dat <- alb_resp_dat %>%
  group_by(sub) %>%
  mutate(
    CG_resp_elic_voc_type = lead(category, default = NA),
    CG_resp_elic_voc_type = if_else(lead(caregiver_response_elicited_vocalization, default = 0) == 1, CG_resp_elic_voc_type, NA_character_)
  ) %>%
  ungroup()

alb_resp_dat %>% 
    select(sub,onset,offset,cat,tier,category,sequence,
           elicited_response,caregiver_response_elicited_vocalization,
           CG_resp_elic_voc_type) %>%
    view()

CG_resp_voc_type_elicit_stats <- alb_resp_dat %>%
    filter(tier == "CGR Modality") %>%
    group_by(sub) %>%
    summarize(total = n(),
              ind_syl_elicited = sum(CG_resp_elic_voc_type == "nonseq", na.rm = TRUE),
              seqs_elicited = sum(CG_resp_elic_voc_type == "seq", na.rm = TRUE),
              prop_seqs = (seqs_elicited / total),
              prop_ind_syll = (ind_syl_elicited / total)) %>%
    ungroup()

CG_resp_voc_type_elicit_stats %>% 
    summarize(prop_seq_mean = mean(prop_seqs),
              prop_ind_syll_mean = mean(prop_ind_syll))

# total vocalizations elicited after caregiver response
alb_resp_dat %>%
    summarize(sum = sum(caregiver_response_elicited_vocalization, na.rm = TRUE))

# total sequences elicited after caregiver response
alb_resp_dat %>%
    filter(sequence == 1) %>%
    summarize(sum = sum(caregiver_response_elicited_vocalization, na.rm = TRUE))

# total individual syllables elicited after caregiver response
alb_resp_dat %>%
    filter(sequence == 0) %>%
    summarize(sum = sum(caregiver_response_elicited_vocalization, na.rm = TRUE))

# ---- establish chance level infant responding to caregiver response with an individual syllable vs vocal sequence
#
# Since there are more individual syllables than vocal sequences, we would expect a higher proportion of individual syllables to be produced after a caregiver response by chance alone.
# To account for this, we will establish this random chance value, by resampling, as many times as there are responses to infant vocalizations and calculating the proportion of vocal types after those.
# We will repeate this process 1000 times and take the mean as our chance level estimate.
# 

alb_session_durations <- read_csv(here("Study 1", "data", "session_durations.csv")) %>%
  mutate(sub = as.factor(sub))

source_python(here("Study 1", "data_proc", "baserate_proc.py"))

alb_resp_dat <- alb_resp_dat %>%
    left_join(alb_session_durations, by = "sub") %>%
    mutate(sub = as.character(sub),
           dur = offset - onset)

# make sure input df only has caregiver verbal responses to infant vocalizations and infant non-cry vocalizations (sequence or individual syllables)

# sample as many times as there are caregiver verbal responses to infant vocalizations

n <- alb_resp_dat %>%
    filter(category == "Vocal") %>%
    summarize(total = n()) %>%
    pull()

n

alb_cg_data <- alb_resp_dat %>%
    filter(category == "Vocal")

subIDs <- unique(alb_resp_dat$sub)

# chance elicitation of sequences

alb_seq_dat <- alb_resp_dat %>%
    filter(category == "seq")

alb_CG_seq_chance <- derive_expected_elicitation_rate(
    trigger_df = alb_cg_data,
    response_df = alb_seq_dat,
    subjects = subIDs,
    threshold = 2,
    n = n
)

alb_CG_seq_chance <- alb_CG_seq_chance %>%
  mutate(expected_elicitation_rate = unlist(expected_elicitation_rate))

alb_CG_seq_chance_mean <- alb_CG_seq_chance %>%
  summarize(mean = mean(expected_elicitation_rate, na.rm = TRUE)) %>%
  pull()

# chance elicitation of individual syllables

alb_nonseq_dat <- alb_resp_dat %>%
    filter(category == "nonseq")

alb_CG_nonseq_chance <- derive_expected_elicitation_rate(
    trigger_df = alb_cg_data,
    response_df = alb_nonseq_dat,
    subjects = subjects,
    threshold = 2,
    n = n
)

alb_CG_nonseq_chance <- alb_CG_nonseq_chance %>%
  mutate(expected_elicitation_rate = unlist(expected_elicitation_rate))

alb_CG_nonseq_chance_mean <- alb_CG_nonseq_chance %>%
  summarize(mean = mean(expected_elicitation_rate, na.rm = TRUE)) %>%
  pull()

# visualized

stat_label <- data.frame(prop=c(.45), vocal_type="prop_seqs")

p2 <- CG_resp_voc_type_elicit_stats %>%
    pivot_longer(cols = c(prop_seqs, prop_ind_syll), names_to = "vocal_type", values_to = "prop") %>%
    ggplot(.,aes(x = vocal_type, y = prop)) +
      stat_summary(fun.y=mean, geom="bar", colour="black",
                   fill="white", alpha = .15,
                   size=1, width = 0.35) + 
      geom_jitter(color="grey", width = 0.01, alpha = .75) +
      stat_summary(fun.data = mean_se, geom = "errorbar", size=1,
                   width=0.075) +      
      geom_segment(aes(x=.65, xend=1.4, y=alb_CG_nonseq_chance_mean, yend=alb_CG_nonseq_chance_mean),
                      col="black",
                      linetype="dotted",
                      size=1) +
      geom_segment(aes(x=1.65, xend=2.4, y=alb_CG_seq_chance_mean, yend=alb_CG_seq_chance_mean),
                      col="black",
                      linetype="dotted",
                      size=1) +
      geom_text(data=stat_label, label="**", size=7, color="black") +
      coord_cartesian(ylim=c(0, 1)) +
      labs(y = "prop verbal responses\nfollowed by infant vocalization",
            x = "",
            title = "")+
      scale_x_discrete(labels = c("prop_ind_syll" = "individual\nsyllable",
                                  "prop_seqs" = "sequence")) +
      coord_cartesian(ylim = c(0.0, 0.5)) +
      theme_classic() +
        theme_classic() +
        theme(text = element_text(size = font_sz),
              axis.ticks.length = unit(-3.5, "pt"),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(hjust = 0.5, size = font_sz),
              legend.text = element_text(size = font_sz),
              legend.key.size = unit(.45, 'cm'))

# group differences
CG_resp_voc_type_elicit_stats_long <- CG_resp_voc_type_elicit_stats %>%
    pivot_longer(cols = c(prop_seqs, prop_ind_syll), names_to = "vocal_type", values_to = "prop")

lm1 <- lmer(prop ~ vocal_type + (1|sub),data=CG_resp_voc_type_elicit_stats_long, REML= FALSE)
emm1 <- emmeans(lm1,pairwise~vocal_type)
summary(emm1$contrasts)

# difference from chance

tidy(t.test(CG_resp_voc_type_elicit_stats$prop_ind_syll, mu = alb_CG_nonseq_chance_mean))

tidy(t.test(CG_resp_voc_type_elicit_stats$prop_seqs, mu = alb_CG_seq_chance_mean))

plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), align = "h")
