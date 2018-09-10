
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import os
import time
from collections import defaultdict
import csv
import seaborn as sb
import config
import rpy2.robjects as robjects
from MapLOINCFields import *
from CleanTestsAndSpecimens import *
from APISearchRequests import *
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr


# ### Read in aggregate data & join with parsed/cleaned test name and specimen

# In[2]:

def build_cube():
    if config.print_status == 'Y':
        print('Building analytic data cube')
    agg_source_data = pd.read_csv(config.in_file, sep=config.delim, quoting=csv.QUOTE_NONE, 
              encoding = "ISO-8859-1", keep_default_na=False, na_values=config.missing)
    
    if os.path.exists(config.out_dir + "Cleaned_Lab_Names.csv") and os.path.exists(config.out_dir + "Cleaned_Specimen_Names.csv"):
        cleaned_tests = pd.read_csv(config.out_dir + "Cleaned_Lab_Names.csv", sep="|", quoting=csv.QUOTE_NONE,
                                encoding = "ISO-8859-1", keep_default_na=False,
                                na_values=config.missing)
        cleaned_specimen = pd.read_csv(config.out_dir + "Cleaned_Specimen_Names.csv", sep="|", quoting=csv.QUOTE_NONE,
                                   encoding = "ISO-8859-1", keep_default_na=False,
                                  na_values=config.missing)
    else:
        cleaned_tests, cleaned_specimen = import_source_data()
    
    agg_source_data[config.test_col] = agg_source_data[config.test_col].str.strip().str.upper()
    agg_source_data[config.spec_col] = agg_source_data[config.spec_col].str.strip().str.upper()
    agg_source_data[config.units] = agg_source_data[config.units].str.strip().str.upper()

    if agg_source_data[config.site].dtypes != cleaned_specimen['Site'].dtypes:
        agg_source_data[config.site] = agg_source_data[config.site].astype(str)
        cleaned_specimen.Site = cleaned_specimen.Site.astype(str)
        cleaned_tests.Site = cleaned_tests.Site.astype(str)
    
    joined_dat = agg_source_data.merge(cleaned_tests, how='left', left_on=[config.site, config.test_col], 
                                   right_on=['Site', 'OriginalTestName'])
    joined_dat = joined_dat.merge(cleaned_specimen, how='left', left_on=[config.site, config.spec_col], 
                                 right_on=['Site', 'OriginalSpecimen'])
    joined_dat = joined_dat.drop(['Site_x', 'Site_y', 'OriginalTestName', 'OriginalSpecimen'], axis=1)
    
    ## Get total # of lab results per site, create normalized 'FreqPercent' variable
    joined_dat = joined_dat.merge(pd.Series.to_frame(joined_dat.groupby(config.site)[config.count].sum(), name='TotalCount').reset_index(),
                                 how='inner', left_on=config.site, right_on=config.site)
    joined_dat['FreqPercent'] = joined_dat[config.count]/joined_dat.TotalCount * 100.0
    
    return joined_dat


# ### Read in UMLS mapped names and put in data frame

# In[3]:

def compile_cuis(data):
    master_list = defaultdict(list)
    for i in range(data.shape[0]):
        if data.loc[i, 'SourceTerm'] not in master_list:
            master_list[data.loc[i, 'SourceTerm']].append(data.loc[i, 'CUI'])
        if data.loc[i, 'CUI'] not in master_list[data.loc[i, 'SourceTerm']]:
            master_list[data.loc[i, 'SourceTerm']].append(data.loc[i, 'CUI'])
    return master_list


# In[4]:

def add_cuis_to_cube(dat):
    if config.print_status == 'Y':
        print('Adding UMLS CUIs')
    feature_col_number = config.num_cuis
    if (os.path.exists(config.out_dir + "UMLS_Mapped_Specimen_Names.csv") and 
        os.path.exists(config.out_dir + "UMLS_Mapped_Test_Names.csv")):
        master_spec_UMLS = compile_cuis(pd.read_csv(config.out_dir + "UMLS_Mapped_Specimen_Names.csv", sep="|"))
        master_test_UMLS = compile_cuis(pd.read_csv(config.out_dir + "UMLS_Mapped_Test_Names.csv", sep="|"))
    else:
        test_input, specimen_input = data_setup()
        master_spec_UMLS = compile_cuis(parse_dat(specimen_input, "Specimen"))
        master_test_UMLS = compile_cuis(parse_dat(test_input, "Test"))
        
    for i in range(feature_col_number):
        dat['SpecCUI{0}'.format(i + 1)] = 'NONE'
        dat['TestCUI{0}'.format(i + 1)] = 'NONE'
        
    for j in range(dat.shape[0]):
        for k in range(len(master_spec_UMLS[dat.at[j, 'CleanedSpecimen']])):
            if k < feature_col_number:
                dat.at[j, 'SpecCUI{0}'.format(k + 1)] = master_spec_UMLS[dat.at[j, 'CleanedSpecimen']][k]
        for l in range(len(master_test_UMLS[dat.at[j, 'CleanedTestName']])):
            if l < feature_col_number:
                dat.at[j, 'TestCUI{0}'.format(l + 1)] = master_test_UMLS[dat.at[j, 'CleanedTestName']][l]
    return dat


# ### Map LOINC System tokens to LOINC Long Name Tokens

# In[5]:

short_to_long, parsed_loinc_fields = parse_loinc()


# In[6]:

def map_loinc_system():
    if config.print_status == 'Y':
        print('Mapping LOINC System')
    if os.path.exists(config.out_dir + "LOINC_System_to_Long.csv"):
        system_map = pd.read_csv(config.out_dir + "LOINC_System_to_Long.csv", sep="|")
    else:
        robjects.numpy2ri.activate()
        stringdist = importr('stringdist', lib_loc=config.lib_loc)
        loinc_syst = parsed_loinc_fields[['System', 'LongName']]
        loinc_syst = loinc_syst[(~pd.isnull(loinc_syst.System)) & (loinc_syst.System != '')].reset_index(drop=True)
        loinc_syst.System = loinc_syst.System.str.split(" ")
        loinc_syst.LongName = loinc_syst.LongName.str.split(" ")
        system_tokens = pd.Series([y for x in loinc_syst.System for y in x]).unique()
        longname_tokens = pd.Series([y for x in loinc_syst.LongName for y in x]).unique()
        system_df = pd.DataFrame(0, index=system_tokens, columns=longname_tokens)
        n_rows = loinc_syst.shape[0]
        for i in range(n_rows):
            for j in range(len(loinc_syst.System[i])):
                dists = stringdist.stringdist(loinc_syst.System[i][j], loinc_syst.LongName[i], method = 'jw', p=0)
                bestMatch = loinc_syst.LongName[i][np.argmin(dists)]
                system_df.loc[loinc_syst.System[i][j], bestMatch] = system_df.loc[loinc_syst.System[i][j], bestMatch] + 1
        high_count = system_df.idxmax(axis=1).values
        system_map = pd.DataFrame({'SystemToken': system_tokens, 'SystemMap': high_count})
        system_map.to_csv(config.out_dir + "LOINC_System_to_Long.csv", sep="|", index=False)
    return system_map


# ### Get highest count, longest string for each short name token mapped by counts to long name tokens

# In[7]:

def map_loinc_token_counts():
    ## Get highest count for each short name token mapped by counts to long name tokens (from the MapLOINCFields script)
    idx = short_to_long.groupby(['Token'])['Count'].transform(max) == short_to_long['Count']
    loinc_terms_max = short_to_long[idx].drop('Count', 1).reset_index(drop = True)
    loinc_terms_max['FinalTokenMap'] = np.nan
    ## If TokenMap contains elongation of TokenAbbreviation, populate AcronymMap column
    loinc_terms_max['AcronymnMap'] = np.nan
    for i in range(loinc_terms_max.shape[0]):
        if not pd.isnull(loinc_terms_max.TokenMap[i]):
            token_map_tokens = loinc_terms_max.TokenMap[i].split(" ")
            if len(loinc_terms_max.Token[i]) > 1 and len(token_map_tokens) >= len(loinc_terms_max.Token[i]):
                counter = 0
                string = ""
                for j in range(len(loinc_terms_max.Token[i])):
                    if loinc_terms_max.Token[i][j:j+1] == token_map_tokens[j][0:1]:
                        counter = counter + 1
                        if len(string) < 1:
                            string = token_map_tokens[j]
                        else:
                            string = string + " " + token_map_tokens[j]
                if counter == len(loinc_terms_max.Token[i]):
                    loinc_terms_max.loc[i, 'AcronymnMap'] = string
    return loinc_terms_max


# In[8]:

def group_func(group):
    ## If Token == TokenMap, make this the FinalTokenMap, otherwise use the shortest TokenMap as the key
    if not group['AcronymnMap'].isnull().all():
        _ = group['FinalTokenMap'].fillna(str(group.AcronymnMap.dropna().unique()[0]), inplace=True)
    elif group['Token'].astype(str).any() == group['TokenMap'].astype(str).any():
        _ = group['FinalTokenMap'] = group['Token']
    elif not group['SystemMap'].isnull().all():
        _ = group['FinalTokenMap'] = group['SystemMap']
    else:
        _ = group['FinalTokenMap'] = min(group.TokenMap, key=len)
    return group


# In[9]:

def combine_loinc_mapping():
    system_map_final = map_loinc_system()
    loinc_terms_max = map_loinc_token_counts()
    loincmap = system_map_final.merge(loinc_terms_max, how='outer', left_on='SystemToken', right_on='Token')
    loincmap = loincmap.set_value(loincmap[loincmap['Token'].isnull()].index, 'Token', 
                              loincmap[loincmap['Token'].isnull()]['SystemToken'])
    loincmap = loincmap.set_value(loincmap[loincmap['TokenMap'].isnull()].index, 'TokenMap', 
                              loincmap[loincmap['TokenMap'].isnull()]['SystemMap'])
    
    if config.print_status == 'Y':
        print('Generating LOINC Groups')
    loincmap = loincmap.groupby('Token').apply(group_func)
    loincmap = loincmap[['Token', 'FinalTokenMap']].drop_duplicates().reset_index(drop=True)
    return loincmap


# In[10]:

# Find best string matches between source data test or specimen tokens and mappings from LOINC short name to
# LOINC long name words
def get_matches(data_col, loincmap):
    if config.print_status == 'Y':
        print('String Distance Matching Source Data Terms to LOINC')
    robjects.numpy2ri.activate()
    stringdist = importr('stringdist', lib_loc=config.lib_loc)
    tokenized_list = [data_col[k].split() for k in range(len(data_col))]
    longest_phrase = len(max(tokenized_list, key=len))
    match_matrix_LV = pd.DataFrame(np.nan, index=data_col, columns=range(longest_phrase))
    match_matrix_JW = pd.DataFrame(np.nan, index=data_col, columns=range(longest_phrase))
    rows = len(tokenized_list)
    for i in range(rows):
        if config.print_status == 'Y' and i % 500 == 0:
            print('Matching Term', i, '/', rows)
        for j in range(len(tokenized_list[i])):
            if not pd.isnull(tokenized_list[i][j]):
                dists_LV = stringdist.stringdist(tokenized_list[i][j], loincmap.Token.values, method='lv')
                dists_JW = stringdist.stringdist(tokenized_list[i][j], loincmap.Token.values, method='jw', p=0)
                match_matrix_LV.iloc[i, j] = loincmap.FinalTokenMap[np.argmin(dists_LV)]
                match_matrix_JW.iloc[i, j] = loincmap.FinalTokenMap[np.argmin(dists_JW)]
    return match_matrix_LV, match_matrix_JW


# In[11]:

def concatenate_match_results(input_matrix, data_type):
    n_rows = input_matrix.shape[0]
    n_cols = input_matrix.shape[1]
    for i in range(n_rows):
        for j in range(1, n_cols):
            if not pd.isnull(input_matrix.iloc[i, j]):
                input_matrix.iloc[i, 0] = input_matrix.iloc[i, 0] + " " + input_matrix.iloc[i, j]
    if data_type == 1:
        return pd.DataFrame(input_matrix.iloc[:, 0], columns=['TestNameMap'])
    else:
        return pd.DataFrame(input_matrix.iloc[:, 0], columns=['SpecimenMap'])


# In[12]:

def add_string_distance_features():
    joined_data = build_cube()
    data = add_cuis_to_cube(joined_data)
    
    loincmap = combine_loinc_mapping()
    
    unique_tests = data[~data.CleanedTestName.isnull()].CleanedTestName.unique()
    unique_specimen_types = data[~data.CleanedSpecimen.isnull()].CleanedSpecimen.unique()

    test_match_matrix_LV, test_match_matrix_JW = get_matches(unique_tests, loincmap)
    spec_match_matrix_LV, spec_match_matrix_JW = get_matches(unique_specimen_types, loincmap)

    if config.print_status == 'Y':
        print('Concatenating String Match Results')
    concat_lv_test_match_result = concatenate_match_results(test_match_matrix_LV, 1)
    concat_jw_test_match_result = concatenate_match_results(test_match_matrix_JW, 1)
    concat_lv_spec_match_result = concatenate_match_results(spec_match_matrix_LV, 2)
    concat_jw_spec_match_result = concatenate_match_results(spec_match_matrix_JW, 2)

    concat_lv_test_match_result.columns.values[0] = 'TestNameMapLV'
    concat_jw_test_match_result.columns.values[0] = 'TestNameMapJW'
    concat_lv_spec_match_result.columns.values[0] = 'SpecimenMapLV'
    concat_jw_spec_match_result.columns.values[0] = 'SpecimenMapJW'

    concat_test_match_result = pd.concat([concat_lv_test_match_result, concat_jw_test_match_result], axis=1)
    concat_spec_match_result = pd.concat([concat_lv_spec_match_result, concat_jw_spec_match_result], 
        axis=1)
    
    dat = data.merge(concat_test_match_result, how='left', left_on='CleanedTestName', right_index=True)
    dat = dat.merge(concat_spec_match_result, how='left', left_on='CleanedSpecimen', right_index=True)
    
    loinc_comp_syst = parsed_loinc_fields[['LOINC', 'Component', 'System']]
    loinc_comp_syst = loinc_comp_syst[(~pd.isnull(loinc_comp_syst.System)) & 
        (loinc_comp_syst.System != '')].reset_index(drop=True)
    loinc_comp_syst['ExpandedSystem'] = np.nan
    loinc_comp_syst.ExpandedSystem = loinc_comp_syst.ExpandedSystem.astype(object)
    
    loinc_num_set = loinc_comp_syst.LOINC.unique()

    if config.print_status == 'Y':
        print('Generating LOINC System Field Expansion')
    for i in range(loinc_comp_syst.shape[0]):
        if not pd.isnull(loinc_comp_syst.System[i]):
            loinc_comp_syst.set_value(i, 'System', loinc_comp_syst.System[i].split(" "))
            for j in range(len(loinc_comp_syst.System[i])):
                mapped_term = loincmap.loc[loincmap.Token == loinc_comp_syst.System[i][j], 'FinalTokenMap'].values[0]
                if pd.isnull(loinc_comp_syst.ExpandedSystem[i]):
                    loinc_comp_syst.set_value(i, 'ExpandedSystem', mapped_term)
                else:
                    loinc_comp_syst.loc[i, 'ExpandedSystem'] = loinc_comp_syst.ExpandedSystem[i] + " " + mapped_term
                
    unique_combos = dat[['TestNameMapJW', 'SpecimenMapJW', 'TestNameMapLV', 'SpecimenMapLV']].drop_duplicates().reset_index(drop=True)
    unique_components = loinc_comp_syst.Component.unique()
    unique_system = loinc_comp_syst[~pd.isnull(loinc_comp_syst.ExpandedSystem)].ExpandedSystem.unique()
    
    unique_combos = pd.concat([unique_combos, pd.DataFrame(columns=[['PredictedComponentJW', 'ComponentMatchDistJW', 'PredictedComponentLV', 'ComponentMatchDistLV', 
               'PredictedSystemJW', 'SystemMatchDistJW', 'PredictedSystemLV', 'SystemMatchDistLV']])])
    
    robjects.numpy2ri.activate()
    stringdist = importr('stringdist', lib_loc=config.lib_loc)

    if config.print_status == 'Y':
        print('String Distance Matching to LOINC Component and System')

    nrows = unique_combos.shape[0]

    for i in range(nrows):
        if i % 500 == 0 and config.print_status == 'Y':
            print('Matching', i, '/', nrows)
        matches = stringdist.stringdist(unique_combos.loc[i, 'TestNameMapJW'], unique_components,
            method='jw', p=0)
        bestmatch = np.argmin(matches)
        unique_combos.loc[i, 'PredictedComponentJW'] = unique_components[bestmatch]
        unique_combos.loc[i, 'ComponentMatchDistJW'] = matches[bestmatch]

        matches = stringdist.stringdist(unique_combos.loc[i, 'TestNameMapLV'], unique_components,
            method='lv')
        bestmatch = np.argmin(matches)
        unique_combos.loc[i, 'PredictedComponentLV'] = unique_components[bestmatch]
        unique_combos.loc[i, 'ComponentMatchDistLV'] = matches[bestmatch]

        matches = stringdist.stringdist(unique_combos.loc[i, 'SpecimenMapJW'], unique_system,
            method='jw', p=0)
        bestmatch = np.argmin(matches)
        unique_combos.loc[i, 'PredictedSystemJW'] = unique_system[bestmatch]
        unique_combos.loc[i, 'SystemMatchDistJW'] = matches[bestmatch]

        matches = stringdist.stringdist(unique_combos.loc[i, 'SpecimenMapLV'], unique_system,
            method='lv')
        bestmatch = np.argmin(matches)
        unique_combos.loc[i, 'PredictedSystemLV'] = unique_system[bestmatch]
        unique_combos.loc[i, 'SystemMatchDistLV'] = matches[bestmatch]
        
    dat = dat.merge(unique_combos, how='left', left_on=['TestNameMapLV', 'TestNameMapJW', 'SpecimenMapLV',
       'SpecimenMapJW'], right_on=['TestNameMapLV', 'TestNameMapJW', 'SpecimenMapLV',
       'SpecimenMapJW'])
    
    dat.to_csv(config.out_dir + 'datCube.csv', index=False)
    
    return dat


# In[ ]:



