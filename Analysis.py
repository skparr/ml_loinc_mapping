
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib inline')
import pandas as pd
import numpy as np
from scipy import stats
from hyperopt import fmin, hp, tpe, STATUS_OK, Trials, space_eval
import pickle
from sklearn import preprocessing
from sklearn.utils import shuffle
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import accuracy_score, f1_score, log_loss, classification
import matplotlib.pyplot as plt
import time
import os
from collections import defaultdict
import seaborn as sb
import config
from LOINCSynonyms import *
from DatasetCreationConfigd import *


# In[3]:

import warnings
warnings.filterwarnings(action='ignore', category= classification.UndefinedMetricWarning)


# In[4]:

seed = 12341
N_SPLITS = config.n_splits
TUNING_EVALS = config.tuning_evals


# In[5]:

filepath = config.out_dir


# ### Load LOINC Synonymn dictionary

# In[6]:

def get_loinc_dict(loinc_synonyms):
    loinc_dict = defaultdict()
    for i in range(loinc_synonyms.shape[0]):
        loinc_dict[loinc_synonyms.loc[i, 'LOINC_NUM']] = loinc_synonyms.loc[i, 'LOINC_KEY']
    return loinc_dict


# ### Load Data Cube

# In[7]:

def get_data():
    if os.path.exists(filepath + 'datCube.csv'):
        dat = pd.read_csv(filepath + 'datCube.csv', sep=',', encoding = "ISO-8859-1",
            keep_default_na=False, na_values=['', 'NULL', 'N/A', 'N\A'])
    else:
        dat = add_string_distance_features()
    return dat


# In[8]:

def transform_and_filter_data():
    dat = get_data()
    loinc_synonyms = get_loinc_synonyms()
    
    feature_col_number = config.num_cuis
    X_unfiltered = dat.copy()  
    X_unfiltered = X_unfiltered.merge(loinc_synonyms, how='left', left_on=config.loinc_col, right_on='LOINC_NUM')     .drop('LOINC_NUM', axis=1)
    X_unfiltered = X_unfiltered.set_value(X_unfiltered[X_unfiltered.LOINC_KEY.isnull()].index, 
        'LOINC_KEY', X_unfiltered.LOINC)
    
    cumulct_orig_loinc = pd.Series.to_frame(X_unfiltered.groupby(config.loinc_col)[config.count].sum(), 
    name='CumulCountByOrigLOINC').reset_index()
    rowct_orig_loinc = pd.Series.to_frame(X_unfiltered[config.loinc_col].value_counts(), name='RowCtByOrigLOINC')
    cumulct_loinc_grp = pd.Series.to_frame(X_unfiltered.groupby('LOINC_KEY')[config.count].sum(),
        name='CumulCountByLOINCGrp').reset_index()
    rowct_loinc_grps = pd.Series.to_frame(X_unfiltered['LOINC_KEY'].value_counts(), name='RowCtByLOINCGrp')

    X_unfiltered = X_unfiltered.merge(cumulct_orig_loinc, how='left', left_on=config.loinc_col, right_on=config.loinc_col)
    X_unfiltered = X_unfiltered.merge(cumulct_loinc_grp, how='left', left_on='LOINC_KEY', right_on='LOINC_KEY')
    X_unfiltered = X_unfiltered.merge(rowct_orig_loinc, how='left', left_on=config.loinc_col, right_index=True)
    X_unfiltered = X_unfiltered.merge(rowct_loinc_grps, how='left', left_on='LOINC_KEY', right_index=True)
    
    # Recode NaNs as 'MISSING' so encoder can be used
    X_unfiltered[[config.loinc_col, config.units, 'TestNameMapJW', 'TestNameMapLV', 'SpecimenMapLV',
             'SpecimenMapJW', 'PredictedComponentJW', 'PredictedComponentLV', 
              'PredictedSystemJW', 'PredictedSystemLV', 'LOINC_KEY']] = X_unfiltered[[config.loinc_col, config.units,
              'TestNameMapJW', 'TestNameMapLV', 'SpecimenMapLV',
             'SpecimenMapJW', 'PredictedComponentJW', 'PredictedComponentLV', 
              'PredictedSystemJW', 'PredictedSystemLV', 'LOINC_KEY']].replace(np.nan, 'MISSING')
    
    cols_to_tranform = [config.units, 'TestNameMapJW', 'TestNameMapLV', 'TestNameMapJW', 
                'SpecimenMapLV', 'SpecimenMapJW', 'PredictedComponentJW',
               'PredictedComponentLV', 'PredictedSystemJW', 'PredictedSystemLV']
    
    for i in range(feature_col_number):
        cols_to_tranform.append('TestCUI{0}'.format(i + 1))
        cols_to_tranform.append('SpecCUI{0}'.format(i + 1))
        
    label_encoder_dict = defaultdict(preprocessing.LabelEncoder)

    fitted = X_unfiltered[cols_to_tranform].apply(lambda x: label_encoder_dict[x.name].fit_transform(x))

    ### Create separate LOINC encoder than includes codes from both LOINC and LOINC_KEY columns
    loinc_coder = preprocessing.LabelEncoder().fit(pd.concat([X_unfiltered[config.loinc_col], X_unfiltered.LOINC_KEY]).unique())

    X_unfiltered[cols_to_tranform] = X_unfiltered[cols_to_tranform].apply(lambda x: label_encoder_dict[x.name].transform(x))

    X_unfiltered[[config.loinc_col, 'LOINC_KEY']] = X_unfiltered[[config.loinc_col, 'LOINC_KEY']]         .apply(lambda x: loinc_coder.transform(x))

    ### Filter out rows with missing test name, specimen type, or LOINC code
    unknowns = X_unfiltered[(X_unfiltered.CleanedTestName.isnull()) |
                           (X_unfiltered.CleanedSpecimen.isnull()) |
                           (loinc_coder.inverse_transform(X_unfiltered.LOINC_KEY)== 'MISSING')]

    X_labeled = X_unfiltered.iloc[~X_unfiltered.index.isin(unknowns.index)]
    
    ### Calculate Number of Sites per LOINC code and per LOINC_KEY
    site_cts = X_labeled.copy()
    loinc_site_cts = site_cts.groupby(config.loinc_col)[config.site].nunique()
    joiner = pd.DataFrame(loinc_site_cts.values, index=loinc_site_cts.index, columns=['SiteCountByOrigLOINC'])
    site_cts = site_cts.merge(joiner, how='left', left_on=config.loinc_col, right_index=True)
    loinc_grp_site_cts = site_cts.groupby('LOINC_KEY')[config.site].nunique()
    joiner2 = pd.DataFrame(loinc_grp_site_cts.values, index=loinc_grp_site_cts.index, columns=['SiteCountByLOINCGrp'])
    site_cts = site_cts.merge(joiner2, how='left', left_on='LOINC_KEY', right_index=True)
    
    ### Filter out rows where LOINC group occurs at only one site
    unknowns = pd.concat([unknowns, site_cts[site_cts.SiteCountByLOINCGrp <= config.min_sites_per_loinc_key]         .drop(['SiteCountByOrigLOINC', 'SiteCountByLOINCGrp'], axis=1)])
    
    X_labeled = X_labeled.iloc[~X_labeled.index.isin(unknowns.index)]
    
    ## Filter out rows where cumulative number of test instances is <= 9 for LOINC group or
    ## where row count per LOINC group is < 2
    unknowns = pd.concat([unknowns, X_labeled[(X_labeled.CumulCountByLOINCGrp <= config.min_tests_per_loinc_group) |
        (X_labeled.RowCtByLOINCGrp < config.min_row_count_per_loinc_group)]])
    
    X_labeled = X_labeled.iloc[~X_labeled.index.isin(unknowns.index)].reset_index(drop=True)
    
    return label_encoder_dict, loinc_coder, X_unfiltered, X_labeled, unknowns


# ### Agnostic site-splitting (does not ensure that test labels are present in training data)

# In[8]:

def get_site_splits():
    n_splits = N_SPLITS
    N_rows = X0.shape[0]
    site_list = X0[config.site].unique()
    site_list.sort()

    site_list = shuffle(site_list, random_state=seed)
    working_list = site_list.copy()

    site_splits = []
    for i in range(n_splits):
        site_ct = 0
        test_sites = []
        while site_ct < len(working_list) and (len(test_sites) == 0 or
            X0[X0[config.site].isin(test_sites)].shape[0] / N_rows < ((1 / n_splits) - 0.005)):
            test_sites.append(working_list[site_ct])
            site_ct = site_ct + 1
        if (i < n_splits - 1):
            site_splits.append(test_sites)
            working_list = working_list[site_ct :]
    site_splits.append(working_list)
    
    return site_splits


# In[9]:

def get_indices():
    test_ind = []
    tune_train_ind = []
    tune_test_ind = []
    
    site_splits = get_site_splits()

    for j in range(len(site_splits)):
        np.random.seed(seed=seed)
        test_ind.append(X0[X0[config.site].isin(site_splits[j])].index)
        tune_test_ind.append(np.random.choice(test_ind[j], replace=False, size=int(len(test_ind[j]) * 0.1666666)))
        tune_train_ind.append(test_ind[j][~test_ind[j].isin(tune_test_ind[j])])
    return test_ind, tune_train_ind, tune_test_ind


# In[ ]:




# In[10]:

label_encoder_dict, loinc_coder, X_unfiltered, X_labeled, unknowns = transform_and_filter_data()

## Select columns for labeled and unlabeled/unknown datasets
data_cols = [config.site, config.units, config.mean_col, config.min_col, config.max_col, 
   config.perc_5, config.perc_25, config.median_col, config.perc_75,
   config.perc_95, config.loinc_col, 'FreqPercent', 'TestNameMapLV',
   'TestNameMapJW', 'SpecimenMapLV', 'SpecimenMapJW',
   'ComponentMatchDistJW', 'ComponentMatchDistLV', 'PredictedComponentJW',
   'PredictedComponentLV', 'PredictedSystemJW', 'PredictedSystemLV',
   'SystemMatchDistJW', 'SystemMatchDistLV', 'LOINC_KEY']
for i in range(config.num_cuis):
    data_cols.append('SpecCUI{0}'.format(i + 1))
    data_cols.append('TestCUI{0}'.format(i + 1))

X0 = X_labeled[data_cols]

data_cols.remove('LOINC_KEY')
    
unknowns_analysis = unknowns[data_cols]


# In[ ]:




# In[11]:

## Get indices for hyperparameter tuning so that the test set data is not used during evaluating hyperparameters
test_ind, tune_train_ind, tune_test_ind = get_indices()


# In[ ]:




# ### Create dictionary for hyperparameters for hyperopt package

# In[30]:

spacedict = {'criterion': ['gini', 'entropy'],
           'max_features': np.arange(2, (X0.shape[1] - 3), 2),
           'max_depth': np.arange(5, 35, 5),
           'min_samples_split': np.arange(2, 20, 2),
           'n_estimators': np.array([10, 20, 50, 75, 100, 125, 150, 175, 200])}

space4rf = {key: hp.choice(key, spacedict[key]) for key in spacedict.keys()}


# ## User hyperopt package to tune RF hyperparameters

# In[ ]:

def rf_hyperopt_train_test(rf_params):
    score_rf = []
    if rf_trials.trials[-1]['tid'] % 5 == 0:
        print('Trial: ', rf_trials.trials[-1]['tid'])
    for i in range(N_SPLITS):
        clf = RandomForestClassifier(random_state=seed, n_jobs=-1, **rf_params)
        X_train = X0.iloc[np.concatenate(tune_train_ind[:i] + tune_train_ind[i + 1:])].drop([config.site, config.loinc_col], axis=1)
        y_train = X_train.pop('LOINC_KEY')
        X_test = X0.iloc[np.concatenate(tune_test_ind[:i] + tune_test_ind[i + 1:])].drop([config.site, config.loinc_col], axis=1)
        y_test = X_test.pop('LOINC_KEY')
        clf.fit(X_train, y_train)
        y_preds = clf.predict(X_test)
        score_rf.append(f1_score(y_test, y_preds, labels=clf.classes_, average='weighted'))
        del clf
    return np.mean(score_rf)


# In[31]:

def rf_f(rf_params):
    global rf_best
    f1 = rf_hyperopt_train_test(rf_params)
    if f1 > rf_best:
        rf_best = f1
        print('new best: ', rf_best, rf_params)
    if rf_trials.trials[-1]['tid'] % 5 == 0:
        pickle.dump(rf_trials, open(filepath + 'rf_tuning_trials_final', 'wb'))
    return {'loss': -f1, 'status': STATUS_OK}


# In[32]:

def get_rf_trials():
    rf_best = 0
    try:
        rf_trials = pickle.load(open(filepath + 'rf_tuning_trials_final', 'rb'))
        if rf_trials.trials[len(rf_trials) - 1]['result']['status'] == 'new':
            rf_trials.trials.pop()
        for i in range(len(rf_trials.trials)):
            if (rf_trials.trials[i]['result']['status'] == 'ok' and
                -rf_trials.trials[i]['result']['loss'] > rf_best):
                    rf_best = -rf_trials.trials[i]['result']['loss']
    except FileNotFoundError:
        rf_trials = Trials()
        
    while len(rf_trials) < TUNING_EVALS:
        rf_best = fmin(rf_f, space4rf, algo=tpe.suggest, max_evals=TUNING_EVALS, 
               trials=rf_trials, rstate=np.random.RandomState(seed))
        
    return rf_trials


# In[33]:

rf_trials = get_rf_trials()

rf_final_parms = dict()

for key in rf_trials.best_trial['misc']['vals'].keys():
    rf_final_parms[key] = spacedict[key][rf_trials.best_trial['misc']['vals'][key][0]]


# ## User hyperopt package to tune RF hyperparameters for OVR 

# In[34]:

def ovr_hyperopt_train_test(ovr_params):
    score_ovr = []
    if ovr_trials.trials[-1]['tid'] % 5 == 0:
        print('Trial: ', ovr_trials.trials[-1]['tid'])
    for i in range(N_SPLITS):
        ovr = OneVsRestClassifier(RandomForestClassifier(random_state=seed, n_jobs=-1, **ovr_params))
        X_train = X0.iloc[np.concatenate(tune_train_ind[:i] + tune_train_ind[i + 1:])].drop([config.site, config.loinc_col], axis=1)
        y_train = X_train.pop('LOINC_KEY')
        X_test = X0.iloc[np.concatenate(tune_test_ind[:i] + tune_test_ind[i + 1:])].drop([config.site, config.loinc_col], axis=1)
        y_test = X_test.pop('LOINC_KEY')
        ovr.fit(X_train, y_train)
        y_preds = ovr.predict(X_test)
        score_ovr.append(f1_score(y_test, y_preds, labels=ovr.classes_, average='weighted'))
        del ovr
    return np.mean(score_ovr)


# In[35]:

def ovr_f(ovr_params):
    global ovr_best
    f1 = ovr_hyperopt_train_test(ovr_params)
    if f1 > ovr_best:
        ovr_best = f1
        print('new best: ', ovr_best, ovr_params)
    if ovr_trials.trials[-1]['tid'] % 5 == 0:
        pickle.dump(ovr_trials, open(filepath + 'ovr_tuning_trials_final', 'wb'))
    return {'loss': -f1, 'status': STATUS_OK}


# In[36]:

def get_ovr_trials():
    ovr_best = 0
    try:
        ovr_trials = pickle.load(open(filepath + 'ovr_tuning_trials_final', 'rb'))
        if ovr_trials.trials[len(ovr_trials) - 1]['result']['status'] == 'new':
            ovr_trials.trials.pop()
        for i in range(len(ovr_trials.trials)):
            if (ovr_trials.trials[i]['result']['status'] == 'ok' and
                -ovr_trials.trials[i]['result']['loss'] > ovr_best):
                    ovr_best = -ovr_trials.trials[i]['result']['loss']
    except FileNotFoundError:
        ovr_trials = Trials()
        
    while len(ovr_trials) < TUNING_EVALS:
        ovr_best = fmin(ovr_f, space4rf, algo=tpe.suggest, max_evals=TUNING_EVALS, 
            trials=ovr_trials, rstate=np.random.RandomState(seed))
    
    return ovr_trials


# In[37]:

ovr_trials = get_ovr_trials()

ovr_final_parms = dict()

for key in ovr_trials.best_trial['misc']['vals'].keys():
    ovr_final_parms[key] = spacedict[key][ovr_trials.best_trial['misc']['vals'][key][0]]


# ## Get performance estimates for models after hyperparameters tuned

# In[38]:

def run_cv(X0, unknowns_analysis):
    metric_names = ['Accuracy', 'F1 weighted', 'F1 macro', 'F1 micro']
    output_names = ['accuracy', 'f1_weighted', 'f1_macro', 'f1_micro']

    ovr_accuracy = np.zeros(N_SPLITS)
    ovr_f1_weighted = np.zeros(N_SPLITS)
    ovr_f1_macro = np.zeros(N_SPLITS)
    ovr_f1_micro = np.zeros(N_SPLITS)
    ovr_preds = []
    ovr_unk_preds = []

    rf_accuracy = np.zeros(N_SPLITS)
    rf_f1_weighted = np.zeros(N_SPLITS)
    rf_f1_macro = np.zeros(N_SPLITS)
    rf_f1_micro = np.zeros(N_SPLITS)
    rf_preds = []
    rf_unk_preds = []

    for i in range(N_SPLITS):
        print("RF CV: ", i + 1)
        X_train = X0[~X0[config.site].isin(site_splits[i])].drop([config.site, config.loinc_col], axis=1)
        y_train = X_train.pop('LOINC_KEY')
        X_test = X0[X0[config.site].isin(site_splits[i])].drop([config.site, config.loinc_col], axis=1)
        y_test = X_test.pop('LOINC_KEY')
        X_unk_test = unknowns_analysis[unknowns_analysis[config.site].isin(site_splits[i])].drop([config.site, 
            config.loinc_col], axis=1)
        rf = RandomForestClassifier(criterion=rf_final_parms['criterion'],
            max_features=rf_final_parms['max_features'],
            max_depth=rf_final_parms['max_depth'],
            min_samples_split=rf_final_parms['min_samples_split'],
            n_estimators=rf_final_parms['n_estimators'],
            n_jobs=-1, random_state=seed)

        rf.fit(X_train, y_train)      
        y_pred = rf.predict(X_test)
        rf_preds.append(y_pred)
        y_unk_pred = rf.predict(X_unk_test)
        rf_unk_preds.append(y_unk_pred)

        rf_accuracy[i] = accuracy_score(y_test, y_pred)
        rf_f1_weighted[i] = f1_score(y_test, y_pred, labels=rf.classes_, average='weighted')
        rf_f1_macro[i] = f1_score(y_test, y_pred, labels=rf.classes_, average='macro')
        rf_f1_micro[i] = f1_score(y_test, y_pred, labels=rf.classes_, average='micro')

        del rf

        print("OVR CV: ", i + 1)
        ovr = OneVsRestClassifier(RandomForestClassifier(criterion=ovr_final_parms['criterion'],
            max_features=ovr_final_parms['max_features'],
            max_depth=ovr_final_parms['max_depth'],
            min_samples_split=ovr_final_parms['min_samples_split'],
            n_estimators=ovr_final_parms['n_estimators'],
            n_jobs=-1, random_state=seed), n_jobs= -1)    

        ovr.fit(X_train, y_train)      
        y_pred = ovr.predict(X_test)
        ovr_preds.append(y_pred)
        y_unk_pred = ovr.predict(X_unk_test)
        ovr_unk_preds.append(y_unk_pred)

        ovr_accuracy[i] = accuracy_score(y_test, y_pred)
        ovr_f1_weighted[i] = f1_score(y_test, y_pred, labels=ovr.classes_, average='weighted')
        ovr_f1_macro[i] = f1_score(y_test, y_pred, labels=ovr.classes_, average='macro')
        ovr_f1_micro[i] = f1_score(y_test, y_pred, labels=ovr.classes_, average='micro')

        del ovr
    
    for i in range(len(ovr_preds)):
        X_intermed = X0[X0[config.site].isin(site_splits[i])].copy()
        X_intermed['RFCVPredLOINC'] = rf_preds[i]
        X_intermed['RFCVLOINC'] = np.where(~(X_intermed[config.loinc_col] == X_intermed['LOINC_KEY']) &
             (X_intermed['LOINC_KEY'] == X_intermed['RFCVPredLOINC']), 
                X_intermed[config.loinc_col], X_intermed['RFCVPredLOINC'])
        X_intermed['OVRCVPredLOINC'] = ovr_preds[i]
        X_intermed['OVRCVLOINC'] = np.where(~(X_intermed[config.loinc_col] == X_intermed['LOINC_KEY']) &
             (X_intermed['LOINC_KEY'] == X_intermed['OVRCVPredLOINC']), 
                X_intermed[config.loinc_col], X_intermed['OVRCVPredLOINC'])

        X_unk_intermed = unknowns_analysis[unknowns_analysis[config.site].isin(site_splits[i])].copy()
        X_unk_intermed['RFCVPredLOINC'] = rf_unk_preds[i]
        X_unk_intermed['OVRCVPredLOINC'] = ovr_unk_preds[i]
        if i == 0:
            X_cvs = X_intermed
            X_unk_cvs = X_unk_intermed
        else:
            X_cvs = pd.concat([X_cvs, X_intermed])
            X_unk_cvs = pd.concat([X_unk_cvs, X_unk_intermed])
        
    X_cvs = X_cvs.sort_index()
    X_unk_cvs = X_unk_cvs.sort_index()

    cvs = np.zeros((len(output_names), 2))
    for i in range(len(output_names)):
        cvs[i][0] = vars()['rf_' + output_names[i]].mean()
        cvs[i][1] = vars()['ovr_' + output_names[i]].mean()

    metrics = pd.DataFrame(data=cvs, index=metric_names, columns=['RF', 'OVR'])
    
    return X_cvs, X_unk_cvs, metrics


# In[39]:

if config.run_cv == 'Y':
    X_cvs, X_unk_cvs, cv_metrics = run_cv(X0, unknowns_analysis)


# ## Fit model to full labeled dataset, make predictions, then make predictions on unlabeled dataset

# In[40]:

X_overall = X0.drop([config.site, config.loinc_col, 'LOINC_KEY'], axis=1)
y_overall = X0['LOINC_KEY']

rf_final = RandomForestClassifier(criterion=rf_final_parms['criterion'],
    max_features=rf_final_parms['max_features'],
    max_depth=rf_final_parms['max_depth'],
    min_samples_split=rf_final_parms['min_samples_split'],
    n_estimators=rf_final_parms['n_estimators'],
    n_jobs=-1, random_state=seed)

rf_final.fit(X_overall, y_overall)
rf_preds = rf_final.predict(X_overall)
rf_preds_frame = pd.DataFrame(rf_preds, index=X0.index, columns=['RFFullModelPredLOINCKey'])

ovr_final = OneVsRestClassifier(RandomForestClassifier(criterion=ovr_final_parms['criterion'],
    max_features=ovr_final_parms['max_features'],
    max_depth=ovr_final_parms['max_depth'],
    min_samples_split=ovr_final_parms['min_samples_split'],
    n_estimators=ovr_final_parms['n_estimators'],
    n_jobs=-1, random_state=seed), n_jobs= -1)

ovr_final.fit(X_overall, y_overall)
ovr_preds = ovr_final.predict(X_overall)
ovr_preds_frame = pd.DataFrame(ovr_preds, index=X0.index, columns=['OVRFullModelPredLOINCKey'])

X_unk_overall = unknowns_analysis.drop([config.site, config.loinc_col], axis=1)

rf_unknown_preds = rf_final.predict(X_unk_overall)
rf_unk_preds_frame = pd.DataFrame(rf_unknown_preds, index=X_unk_overall.index, columns=['RFFullModelPredLOINCKey'])

ovr_unknown_preds = ovr_final.predict(X_unk_overall)
ovr_unk_preds_frame = pd.DataFrame(ovr_unknown_preds, index=X_unk_overall.index, columns=['OVRFullModelPredLOINCKey'])


# In[41]:

if config.run_cv == 'Y':
    X0 = X_cvs.copy()
    unknowns_analysis = X_unk_cvs.copy()


# In[64]:

X_final = X0.merge(rf_preds_frame, how='inner', left_index=True, right_index=True)     .merge(ovr_preds_frame, how='inner', left_index=True, right_index=True)
unknowns_final = unknowns_analysis.merge(rf_unk_preds_frame, how='inner', left_index=True, right_index=True)     .merge(ovr_unk_preds_frame, how='inner', left_index=True, right_index=True)


# In[65]:

## If original LOINC code is a member of the predicted LOINC_KEY, retain the original LOINC code as the 
## final label
X_final['RFFullModelFinalLOINCKey'] = np.where(~(X_final[config.loinc_col] == X_final['LOINC_KEY']) &
         (X_final['LOINC_KEY'] == X_final['RFFullModelPredLOINCKey']), 
            X_final[config.loinc_col], X_final['RFFullModelPredLOINCKey'])
X_final['OVRFullModelFinalLOINCKey'] = np.where(~(X_final[config.loinc_col] == X_final['LOINC_KEY']) &
         (X_final['LOINC_KEY'] == X_final['OVRFullModelPredLOINCKey']), 
            X_final[config.loinc_col], X_final['OVRFullModelPredLOINCKey'])


# In[66]:

cols_to_transform = X_final.columns[X_final.columns.isin(label_encoder_dict.keys())]


# In[67]:

## Inverse transforms factor variables back to original data
X_final[cols_to_transform] = X_final[cols_to_transform]     .apply(lambda x: label_encoder_dict[x.name].inverse_transform(x))
unknowns_final[cols_to_transform] = unknowns_final[cols_to_transform]     .apply(lambda x: label_encoder_dict[x.name].inverse_transform(x))

## Inverse transforms factor-converted LOINC codes back to actual LOINC codes
if config.run_cv == 'Y':
    X_final[['RFCVPredLOINC', 'RFCVLOINC', 'OVRCVPredLOINC', 'OVRCVLOINC']] = loinc_coder         .inverse_transform(X_final[['RFCVPredLOINC', 'RFCVLOINC', 'OVRCVPredLOINC', 'OVRCVLOINC']])
    unknowns_final[['RFCVPredLOINC', 'OVRCVPredLOINC']] = loinc_coder         .inverse_transform(unknowns_final[['RFCVPredLOINC', 'OVRCVPredLOINC']])
X_final[[config.loinc_col, 'LOINC_KEY', 'RFFullModelPredLOINCKey', 'OVRFullModelPredLOINCKey',
        'RFFullModelFinalLOINCKey', 'OVRFullModelFinalLOINCKey']] = loinc_coder \
    .inverse_transform(X_final[[config.loinc_col, 'LOINC_KEY', 'RFFullModelPredLOINCKey', 'OVRFullModelPredLOINCKey',
        'RFFullModelFinalLOINCKey', 'OVRFullModelFinalLOINCKey']])
unknowns_final[[config.loinc_col, 'RFFullModelPredLOINCKey', 'OVRFullModelPredLOINCKey']] = loinc_coder     .inverse_transform(unknowns_final[[config.loinc_col, 'RFFullModelPredLOINCKey', 'OVRFullModelPredLOINCKey']])


# In[68]:

X_final = X_labeled[[config.site, config.test_col, 'CleanedTestName', config.spec_col, 'CleanedSpecimen', 
    config.count, config.loinc_col, 'LOINC_KEY']] \
    .merge(X_final.drop([config.site, config.loinc_col, 'LOINC_KEY'], axis=1), 
           how='inner', left_index=True, right_index=True)
    
unknowns_final = unknowns[[config.site, config.test_col, 'CleanedTestName', config.spec_col, 'CleanedSpecimen',
    config.count, config.loinc_col, 'LOINC_KEY']] \
    .merge(unknowns_final.drop([config.site, config.loinc_col], axis=1),
          how='inner', left_index=True, right_index=True)

X_final.to_csv(filepath + 'labeled_data_predictions.csv')
unknowns_final.to_csv(filepath + 'unlabeled_data_predictions.csv')

