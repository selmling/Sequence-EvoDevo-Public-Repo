from random import sample
import pandas as pd
import numpy as np
import random
import math

# Set a seed for replicability
SEED = 46
random.seed(SEED)
np.random.seed(SEED)

# helper round function

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

# baserate random sample function

def br_rand_samp(df,subjects,threshold,n):
    indat = df.copy()
    outdat = pd.DataFrame(index = np.arange(n), columns = np.arange(2))
    outdat.columns = ['expected sequence response','expected non-sequence response']
    for i, row in outdat.iterrows():
        sub = sample(subjects,1) # sample rand sub
        subject_frame = indat[indat["sub"] == sub[0]]
        time_range = pd.DataFrame({
            "sess_min_time": [subject_frame["onset"].min()],
            "sess_max_time": [subject_frame["session_offset"].max()]})
        rand_time = round_up(random.uniform(
            time_range["sess_min_time"][0],
            time_range["sess_max_time"][0])
                            ,3)
        within_seq_window = subject_frame[subject_frame["sequence"] == 1]
        within_seq_window = within_seq_window[(rand_time > within_seq_window["offset"]) &
                                              (rand_time < within_seq_window["offset"]+threshold)]
        if len(within_seq_window) > 0:
            outdat.loc[i,"expected sequence response"] = 1
        else:
            outdat.loc[i,"expected sequence response"] = 0
        within_nonseq_window = subject_frame[subject_frame["sequence"]==0]
        within_nonseq_window = within_nonseq_window[(rand_time > within_nonseq_window["offset"]) &
                                                    (rand_time < within_nonseq_window["offset"]+threshold)]
        if len(within_nonseq_window) > 0:
            outdat.loc[i,"expected non-sequence response"] = 1
        else:
            outdat.loc[i,"expected non-sequence response"] = 0
    return outdat

def create_base_rate(df,subjects,threshold,n):
    result = outdat = pd.DataFrame(index = np.arange(1000), columns = np.arange(1))
    result.columns = ["expected_sequence_response_ratio"]
    
    for i, row in result.iterrows():
        resp_br = br_rand_samp(df,subjects,threshold,n)
        ns = resp_br["expected non-sequence response"].sum()
        s = resp_br["expected sequence response"].sum()
        ex_cg_seq_rsp_ratio = s/(ns+s)
        result.loc[i,"expected_sequence_response_ratio"] = ex_cg_seq_rsp_ratio
    return result

# expected elicitation rate random sample function

def ee_rand_samp(cg_df, inf_df, subjects, threshold, n):
    infdat = inf_df.copy()
    cgdat = cg_df.copy()
    outdat = pd.DataFrame(index = np.arange(n), columns = np.arange(1))
    outdat.columns = ['expected elicitation']
    for i, row in outdat.iterrows():
        dur = infdat["dur"].sample(n = 1) # randomly sampled vocalization duration
        sub = sample(subjects,1) # sample random caregiver
        subject_frame = cgdat[cgdat["sub"] == sub[0]]
        time_range = pd.DataFrame({
            "sess_min_time": [subject_frame["onset"].min()], 
            "sess_max_time": [subject_frame["session_offset"].max()]})
        if pd.isna(time_range["sess_min_time"][0]) or pd.isna(time_range["sess_max_time"][0]):
            outdat.loc[i, "expected elicitation"] = 0
            continue
        rand_time = round_up(random.uniform(
            time_range["sess_min_time"][0],
            time_range["sess_max_time"][0])
                            ,3)
        rand_onset = rand_time
        rand_offset = rand_time + dur
        within_resp_window = subject_frame[(rand_onset < subject_frame["onset"]) &
                                           ((rand_offset.values[0] + threshold) > subject_frame["onset"])]
        if len(within_resp_window) > 0:
            outdat.loc[i,"expected elicitation"] = 1
        else:
            outdat.loc[i,"expected elicitation"] = 0
    return outdat
    
def create_expected_elic_rate(cg_df, inf_df, subjects, threshold, n):
    result = outdat = pd.DataFrame(index = np.arange(1000), columns = np.arange(1))
    result.columns = ["expected_elicitation_rate"]
    
    for i, row in result.iterrows():
        eliciation_sample = ee_rand_samp(cg_df, inf_df, subjects, threshold, n)
        elicitation_rate = eliciation_sample["expected elicitation"].mean()
        result.loc[i,"expected_elicitation_rate"] = elicitation_rate
    return result

# Generalized expected elicitation rate random sample function

def random_elicitation_sample(trigger_df, response_df, subjects, threshold, n):
    response_data = response_df.copy()
    trigger_data = trigger_df.copy()
    outdat = pd.DataFrame(index = np.arange(n), columns = np.arange(1))
    outdat.columns = ['expected_elicitation']
    for i, row in outdat.iterrows():
        dur = trigger_data["dur"].sample(n=1).values[0]
        sub = sample(subjects,1)
        subject_frame = response_data[response_data["sub"] == sub[0]]
        time_range = pd.DataFrame({
            "sess_min_time": [subject_frame["onset"].min()], 
            "sess_max_time": [subject_frame["session_offset"].max()]})
        rand_time = round_up(random.uniform(
            time_range["sess_min_time"][0],
            time_range["sess_max_time"][0]
        ), 3)
        rand_onset = rand_time
        rand_offset = rand_time + dur
        within_resp_window = subject_frame[(rand_onset < subject_frame["onset"]) &
                                           ((rand_offset + threshold) > subject_frame["onset"])]
        if len(within_resp_window) > 0:
            outdat.loc[i, "expected_elicitation"] = 1
        else:
            outdat.loc[i, "expected_elicitation"] = 0
    return outdat

def derive_expected_elicitation_rate(trigger_df, response_df, subjects, threshold, n):
    result = outdat = pd.DataFrame(index = np.arange(1000), columns = np.arange(1))
    result.columns = ["expected_elicitation_rate"]

    for i, row in result.iterrows():
        print(f"Processing sample {i + 1} of 1000")
        elicitation_sample = random_elicitation_sample(trigger_df, response_df, subjects, threshold, n)
        elicitation_rate = elicitation_sample["expected_elicitation"].mean()
        result.loc[i, "expected_elicitation_rate"] = elicitation_rate
        
    return result