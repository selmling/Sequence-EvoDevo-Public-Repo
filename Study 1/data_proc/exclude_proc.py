import pandas as pd
import numpy as np

def less_five_voc(df):
    count_sumstats = df.pivot_table(index=["sub"],
                                    columns='sequence',
                                    aggfunc='size',
                                    fill_value=0)
    
    # who produced less than 5 individual syllables
    under_ind_syll = count_sumstats[count_sumstats[0.0] < 5].reset_index()
    under_ind_syll = under_ind_syll[["sub"]]

    # who produced less than 5 sequences
    under_voc_seq = count_sumstats[count_sumstats[1.0] < 5].reset_index()
    under_voc_seq = under_voc_seq[["sub"]]

    # Combine both under_ind_syll and under_voc_seq using concat
    exclude_list = pd.concat([under_ind_syll, under_voc_seq], ignore_index=True)
    exclude_list["eligibility"] = "exclude"
    
    # Merge to get all exclude indices
    exclude = df.reset_index().merge(exclude_list, on=["sub"]).set_index('index')
    
    indices = exclude.index
    return indices