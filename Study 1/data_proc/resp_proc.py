import pandas as pd
import numpy as np

# Set a seed for replicability
SEED = 46
# random.seed(SEED)
np.random.seed(SEED)


def prop_response(df,threshold,subjects):
    result = df.copy()
    result["elicited_response"] = 0
    for sub in subjects:
        subject_frame = result[result["sub"] == sub]
        subject_frame = subject_frame[subject_frame["cat"] != "Other"]
        vocs = subject_frame[subject_frame["tier"] == "VCV"].index.values
        for v_idx in vocs:
            voc_onset = result["onset"][v_idx]
            voc_offset = result["offset"][v_idx]
            within_window = subject_frame[(subject_frame["onset"] >= voc_onset) &
                                          (subject_frame["onset"] <= voc_offset + threshold)]
            cont_parent_row = within_window.index[within_window["tier"] == "Parent"].tolist()
            if len(cont_parent_row) > 0:
                result.loc[cont_parent_row[0]-1, "elicited_response"] = 1
    return result


def prop_inf_response(df, threshold, subjects):
    result = df.copy()
    result["caregiver_response_elicited_vocalization"] = 0
    for sub in subjects:
        print(sub)
        subject_frame = result[(result["sub"] == sub) & (result["cat"] != "Other")]
        caregiver_indices = subject_frame[subject_frame["tier"] == "CGR Modality"].index.values
        for cg_idx in caregiver_indices:
            cg_offset = result["offset"][cg_idx]  # Use caregiver offset as the start of the response window
            
            # Identify rows within the threshold window following the caregiver vocalization
            within_window = subject_frame[(subject_frame["onset"] >= cg_offset) &
                                          (subject_frame["onset"] <= cg_offset + threshold)]
            # Check if an infant response occurs within this window
            infant_rows = within_window.index[within_window["tier"] == "Infraphonology"].tolist()
            if infant_rows:
                result.loc[infant_rows[0], "caregiver_response_elicited_vocalization"] = 1
    return result

def create_result(path,subjects,outname):
    # read in data
    table = pd.read_csv(path)
    result = prop_response(table, 2,subjects)
    result.to_csv("../data/{}".format(outname))

def prop_interrupt(df):
    result = df.copy()
    result["interrupted"] = 0
    for sub in subjects:
        subject_frame = result[result["sub"] == sub]
        subject_frame = subject_frame[subject_frame["cat"] != "Other"]
        vocs = subject_frame[subject_frame["tier"] == "Infraphonology"].index.values
        for v_idx in vocs:
            voc_onset = result["onset"][v_idx]
            voc_offset = result["offset"][v_idx]
            within_window = subject_frame[(subject_frame["onset"] >= voc_onset) & (subject_frame["onset"] <= voc_offset)]
            cont_parent_row = within_window.index[within_window["tier"] == "CGR Vocal Type"].tolist()
            if len(cont_parent_row) > 0:
                result.loc[cont_parent_row[0]-2, "interrupted"] = 1
    return result

def create_interrupts(path):
    # read in data
    table = pd.read_csv(path)
    result = prop_interrupt(table)
    result.to_csv("../data/prop_interrupt_resp_type.csv")

def prop_turns(df,threshold):
    result = df.copy()
    result["conversational_turn"] = 0
    for sub in subjects:
        subject_frame = result[result["sub"] == sub]
        subject_frame = subject_frame[subject_frame["cat"] != "Other"]
        vocs = subject_frame[subject_frame["tier"] == "Infraphonology"].index.values
        for v_idx in vocs:
            voc_onset = result["onset"][v_idx]
            voc_offset = result["offset"][v_idx]
            within_window = subject_frame[(subject_frame["onset"] >= voc_onset - threshold) & (subject_frame["onset"] <= voc_offset + threshold)]
            cont_parent_row = within_window.index[within_window["tier"] == "CGR Vocal Type"].tolist()
            if len(cont_parent_row) > 0:
                result.loc[cont_parent_row[0]-2, "conversational_turn"] = 1
    return result

def create_turns(path):
    # read in data
    table = pd.read_csv(path)
    result = prop_turns(table,5)
    result.to_csv("../data/prop_turns_resp_type.csv")
