
# coding: utf-8

# In[1]:

import csv
import re
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import config


# In[2]:

loincFilePath = config.loinc_file_path


# In[3]:

def clean_terms(dataElement):
    insigWords = ["IN", "FROM", "ON", "OR", "OF", "BY", "AND", "&", "", " "]
    
    data = (dataElement.replace("'", "").replace(",", " ").replace(".", " ")         .replace(":", " ").replace('\t', " ").replace("^", " ").replace("+", " ")         .replace("*", " ").replace("~", " ").replace("(", " ").replace(")", " ")         .replace("!",  " ").replace("[", " ").replace("]", " ")         .replace("_", " ").replace("|", " ").replace('"', " ")         .replace("-", " ").replace("/", " ").replace("\\", " ")         .replace("#", " ").replace("?", " ").replace("%", " ")         .replace("<", " ").replace(">", " ").replace("@", " ")         .replace("=", " ").split(" "))
    
    i = 0
    while i < len(data):
        if data[i] in insigWords:
            data.remove(data[i])
        else:
            data[i] = data[i].strip()
            i = i + 1
    return data


# In[4]:

def add_match(data, shortName, matchName):
    if matchName != "":
        data[shortName][matchName] += 1


# In[5]:

def expand_words(data, shortWords, longWords):
    stringDelta = len(longWords) - len(shortWords)
    for i in range(len(shortWords)):
        match1 = ""
        match2 = ""
        match3 = ""
        j = 0
        for j in range(len(longWords)):
            if stringDelta >= 0:
                if ((i + 0) > j): continue
                if ((i + 2) < j): break
                if ((i + 0) == j):
                    match1 = longWords[j]
                if ((i + 1) == j):
                    match2 = match1 + " " + longWords[j]; 
                if ((i + 2) == j):
                    match3 = match2 + " " + longWords[j]; 
            else:
                if ((i - 2) > j): continue;
                if ((i + 0) < j): break;
                if ((i - 2) == j):
                    match1 = longWords[j]; 
                if ((i - 1) == j):
                    match2 = match1 + " " + longWords[j]; 
                if ((i + 0) == j):
                    match3 = match2 + " " + longWords[j];
        add_match(data, shortWords[i], match1);
        add_match(data, shortWords[i], match2);
        add_match(data, shortWords[i], match3);
    return data


# In[8]:

def parse_loinc():
    reader = csv.reader(open(loincFilePath, encoding='utf8'))
    index = -1
    loincs = list()
    shortToLong = defaultdict(Counter)
    componentParsed = list()
    systemParsed = list()
    longParsed = list()
    for fields in reader:
        index = index + 1
        if index == 0:
            loincNumInd = fields.index('LOINC_NUM')
            componentInd = fields.index('COMPONENT')
            systemInd = fields.index('SYSTEM')
            shortNameInd = fields.index('SHORTNAME')
            longNameInd = fields.index('LONG_COMMON_NAME')
            classTypeInd = fields.index('CLASSTYPE')
            continue
        loincNum = fields[loincNumInd]
        component = fields[componentInd].upper()
        system = fields[systemInd].upper();
        shortName = fields[shortNameInd].upper();
        longName = fields[longNameInd].upper();
        classType = fields[classTypeInd];
        if classType != "1" and classType != "2": continue  #only keep the lab and clinical class types
        loincs.append(loincNum)
        shortWords = clean_terms(shortName)
        longName = re.sub(r"\[([A-Za-z0-9]*\s*\/*)*\]", "", longName)
        longWords = clean_terms(longName)
        componentWords = clean_terms(component)
        systemWords = clean_terms(system)
        componentParsed.append(" ".join(componentWords))
        systemParsed.append(" ".join(systemWords))
        longParsed.append(" ".join(longWords))
        shortToLong = expand_words(shortToLong, shortWords, longWords)
    short_to_long_df = pd.DataFrame(data=[[outer_key, inner_key, shortToLong[outer_key][inner_key]] 
        for outer_key in shortToLong for inner_key in shortToLong[outer_key]], columns=['Token',
        'TokenMap', 'Count'])
    parsed_loinc_fields_df = pd.DataFrame(data=list(zip(loincs, componentParsed, systemParsed, longParsed)),
        columns=['LOINC', 'Component', 'System', 'LongName'], dtype=object)
    if config.write_file_loinc_parsed:
        short_to_long_df.to_csv(config.out_dir + "LOINC_Name_Map.csv", sep="|", index=False)
        parsed_loinc_fields_df.to_csv(config.out_dir + "LOINC_Parsed_Component_System_Longword.csv", sep="|", index=False)
    return [short_to_long_df, parsed_loinc_fields_df]


# In[ ]:



