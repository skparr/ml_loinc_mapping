
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
from collections import defaultdict
import config


# In[2]:

def group_func(group):
    if group.LOINC_NUM.count() >= 2:
        return group


# #### Creates LOINC synonyms for terms with same Component, Property, Time Aspect, System, and Scale Type. 'Key' is the entry with null Method Type and Status == 'ACTIVE'

# In[3]:

def get_loinc_groups():
    if config.print_status == 'Y':
        print('Processing LOINC synonyms')
    raw_loinc_dat = pd.read_csv(config.loinc_file_path, 
        usecols=['LOINC_NUM', 'COMPONENT', 'PROPERTY', 'TIME_ASPCT', 'SYSTEM', 'METHOD_TYP',
           'SCALE_TYP', 'STATUS', 'CLASSTYPE', 'VersionLastChanged'], low_memory=False)
    raw_loinc_dat = raw_loinc_dat[(raw_loinc_dat.CLASSTYPE == 1) | (raw_loinc_dat.CLASSTYPE == 2)].reset_index(drop=True)
    
    raw_loinc_mults = raw_loinc_dat[['LOINC_NUM', 'COMPONENT', 'PROPERTY', 'TIME_ASPCT', 'SCALE_TYP',
    'SYSTEM']].groupby(['COMPONENT', 'PROPERTY', 'TIME_ASPCT', 'SYSTEM', 'SCALE_TYP']) \
    .apply(group_func) \
    .drop('LOINC_NUM', axis=1) \
    .dropna()
    
    raw_loinc_mults = raw_loinc_dat.merge(raw_loinc_mults, how='inner', left_on=['COMPONENT', 'PROPERTY', 'TIME_ASPCT', 
            'SYSTEM', 'SCALE_TYP'], right_on=['COMPONENT', 'PROPERTY', 'TIME_ASPCT', 'SYSTEM', 'SCALE_TYP']) \
    .drop_duplicates() \
    .reset_index(drop=True)
    
    raw_loinc_mults['CANDIDATE_KEYS'] = np.nan
    raw_loinc_mults['LOINC_KEY'] = np.nan
    
    return raw_loinc_mults


# In[4]:

def quant_nans(group):
    if group[group.METHOD_TYP.isnull()].shape[0] == 1:
        _ = group['CANDIDATE_KEYS'] = 1
        _ = group['LOINC_KEY'] = group[group.METHOD_TYP.isnull()].LOINC_NUM.values[0]
    elif not group['METHOD_TYP'].isnull().any():
        if len(group.METHOD_TYP.unique() == 1) and group[group.STATUS == 'ACTIVE'].shape[0] == 1:
            _ = group['CANDIDATE_KEYS'] = 1
            _ = group['LOINC_KEY'] = group[group.STATUS == 'ACTIVE'].LOINC_NUM.values[0]
        else:
            _ = group['CANDIDATE_KEYS'] = 0
    elif group['METHOD_TYP'].isnull().any() and group[(group.METHOD_TYP.isnull()) & (group.STATUS=='ACTIVE')]         .shape[0]==1:
            _ = group['CANDIDATE_KEYS'] = 1
            _ = group['LOINC_KEY'] = group[(group.METHOD_TYP.isnull()) & (group.STATUS=='ACTIVE')]                 .LOINC_NUM.values[0]
    elif group['METHOD_TYP'].isnull().any() and group[(group.METHOD_TYP.isnull()) & (group.STATUS=='DISCOURAGED')]         .shape[0]==1:
            _ = group['CANDIDATE_KEYS'] = 1
            _ = group['LOINC_KEY'] = group[(group.METHOD_TYP.isnull()) & (group.STATUS=='DISCOURAGED')]                 .LOINC_NUM.values[0]
    elif group[group.STATUS == 'DEPRECATED'].shape[0] == group.shape[0]:
        _ = group['CANDIDATE_KEYS'] = group.shape[0]
        max_vers = group['VersionLastChanged'].max()
        if group[(group.VersionLastChanged == max_vers) & (group.METHOD_TYP.isnull())].shape[0] > 1:
            _ = group['LOINC_KEY'] = group[(group.VersionLastChanged == max_vers) & 
                                           (group.METHOD_TYP.isnull())].LOINC_NUM.values[-1]
        else:
            _ = group['LOINC_KEY'] = group[(group.VersionLastChanged == max_vers) & 
                                           (group.METHOD_TYP.isnull())].LOINC_NUM.values[0]
    else:
        _ = group['CANDIDATE_KEYS'] = group[group.METHOD_TYP.isnull()].shape[0]
    return group


# In[3]:

def get_loinc_synonyms():
    raw_loinc_mults = get_loinc_groups()

    raw_loinc_mults = raw_loinc_mults.groupby(['COMPONENT', 'PROPERTY', 'TIME_ASPCT', 'SYSTEM', 'SCALE_TYP'])         .apply(quant_nans)

    raw_loinc_mults = raw_loinc_mults[~raw_loinc_mults.LOINC_KEY.isnull()]
    
    final_loinc_keys = raw_loinc_mults[['LOINC_NUM', 'LOINC_KEY']].reset_index(drop = True)

    final_loinc_keys.to_csv(config.out_dir + 'loinc_synonymns.csv', index=False)
    
    return final_loinc_keys


# In[ ]:



