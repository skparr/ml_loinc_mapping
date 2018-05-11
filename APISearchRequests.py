
# coding: utf-8

# In[1]:

from Authentication import *
from CleanTestsAndSpecimens import *
import requests
import json
import argparse
import time
from collections import defaultdict
import pandas as pd
import numpy as np
import os
import urllib
import config


# In[3]:

def data_setup():
    if os.path.exists(config.out_dir + "Cleaned_Lab_Names.csv") and os.path.exists(config.out_dir + "Cleaned_Specimen_Names.csv"):
        test_input = pd.read_csv(config.out_dir + "Cleaned_Lab_Names.csv", sep='|')
        specimen_input = pd.read_csv(config.out_dir + "Cleaned_Specimen_Names.csv", sep='|')
    else:
        test_input, specimen_input = import_source_data()    
    return [test_input, specimen_input]


# In[ ]:

def uri_setup():
    apikey = config.api_key
    AuthClient = Authentication(apikey)
    return [uri, AuthClient]


# In[4]:

def get_match(term, searchType, tgt, AuthClient):
    maxCount = 20
    count = 0
    cont = True
    term = urllib.parse.quote(term)
    uri = "https://uts-ws.nlm.nih.gov"
    while cont and count < maxCount:
        try:
            if searchType == "exact":
                content_endpoint = "/rest/search/current?string=" + term + "&searchType=exact&ticket=" + AuthClient.getst(tgt)
            if searchType == "word":
                 content_endpoint = "/rest/search/current?string=" + term + "&searchType=words&ticket=" + AuthClient.getst(tgt)
            r = requests.get(uri+content_endpoint)
            cont = False
        except:
            print("failed in match loop")
            cont = True
            count = count + 1
    items = json.loads(r.text)
    return items["result"]   


# In[5]:

def parse_dat(data, filetype):
    uri, AuthClient = uri_setup()
    
    start = time.time()
    tgt = AuthClient.gettgt()
    ref_col = data.columns[2]
    
    master = defaultdict(list)
        
    for i in range(data.shape[0]):
        if i > 0:
            if ((time.time() - start) > 270):
                start = time.time()
                tgt = AuthClient.gettgt()
            token = data.loc[i, ref_col]
            if token not in master:
                jsonData = get_match(token, "exact", tgt, AuthClient)    
                if jsonData["results"][0]["ui"] == "NONE":
                    jsonData = get_match(token, "word", tgt, AuthClient) 
                    if jsonData["results"][0]["ui"] == "NONE":
                        for term in token.split(sep=" "):
                            jsonData = get_match(term, "exact", tgt, AuthClient)
                            if jsonData["results"][0]["ui"] != "NONE" and len(master[token]) < 3:
                                master[token].append([jsonData["results"][0]["ui"],
                                                                    jsonData["results"][0]["name"],
                                                                    jsonData["results"][0]["rootSource"]])
            if token not in master:
                for j, sets in enumerate(jsonData["results"]):
                    if len(sets) > 2:
                        if j < 3 and len(master[token]) < 3:
                            master[token].append([sets["ui"], sets["name"], sets["rootSource"]])
                    else:
                        master[token].append([sets["ui"], sets["name"], "NONE"])
    master = pd.DataFrame(([keys, items[0], items[1]] for keys in master.keys() for items in master[keys]),
        columns = ['SourceTerm', 'CUI', 'MappedName'])
    if config.write_file_umls_cuis:
        if filetype == 'Specimen':
            master.to_csv(config.out_dir + "UMLS_Mapped_Specimen_Names.csv", sep='|', index=False)
        else:
            master.to_csv(config.out_dir + "UMLS_Mapped_Test_Names.csv", sep='|', index=False)
    return master


# In[ ]:



