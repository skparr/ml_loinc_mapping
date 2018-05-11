
# coding: utf-8

# In[1]:

import re
import config
import pandas as pd
from collections import defaultdict, Counter


# In[2]:

REJECTION_THRESHOLD = config.rejection_threshold


# In[3]:

def import_source_data():
    if config.print_status == 'Y':
        print('Importing source data')
    testNameList = defaultdict(list)
    specimenList = defaultdict(list)
    
    reader = open(config.in_file, 'r')
    index = -1
    for line in reader:
        fields = line.split(config.delim)
        if index == -1:
            siteIdentifierCol = fields.index(config.site)
            testCol = fields.index(config.test_col)
            specimenCol = fields.index(config.spec_col)
        index = index + 1
        if index == 0: continue
        if (fields[siteIdentifierCol].upper().strip() not in testNameList.keys() or
            (fields[siteIdentifierCol].upper().strip() in testNameList.keys() and
            fields[testCol].upper().strip() not in testNameList[fields[siteIdentifierCol]])):
            testNameList[fields[siteIdentifierCol].upper().strip()].append(fields[testCol].upper().strip())
        if (fields[siteIdentifierCol].upper().strip() not in specimenList.keys() or 
            (fields[siteIdentifierCol].upper().strip() in specimenList.keys() and
             fields[specimenCol].upper().strip() not in specimenList[fields[siteIdentifierCol]])):
            specimenList[fields[siteIdentifierCol].upper().strip()].append(fields[specimenCol].upper().strip())
    
    cleanedTests = clean_terms(testNameList, 'testNames')
    cleanedSpecimens = clean_terms(specimenList, 'specimenNames')
    return [cleanedTests, cleanedSpecimens]


# In[4]:

def clean_terms(sourceData, dataType):
    if config.print_status == 'Y':
        print('Cleaning source data')
    insigWords = ["IN", "FROM", "ON", "OR", "OF", "BY", "AND", "&", "TO", "BY", "", " "]
    siteWordCount = defaultdict(Counter)
    siteTotalWordCount = defaultdict(int)
    cleanedList = defaultdict(lambda: defaultdict(list))
    discardedTerms = defaultdict(list)
    for siteKey in sourceData.keys():
        for term in sourceData[siteKey]:
            modTerm = (term.replace("'", "").replace(",", " ").replace(".", " ")                 .replace(":", " ").replace('\t', " ").replace("^", " ").replace("+", " ")                 .replace("*", " ").replace("~", " ").replace("(", " ").replace(")", " ")                 .replace("!",  " ").replace("[", " ").replace("]", " ")                 .replace("_", " ").replace("|", " ").replace('"', " ").split(" "))

            i = 0
            while i < len(modTerm):
                modTerm[i] = re.sub(r"\d{1,2}[\/-]\d{1,4}([\/-]\d{2,4})*|\d{6}", "", modTerm[i])
                if modTerm[i] != None and len(modTerm[i]) > 0:
                    i = i + 1
                else:
                    modTerm.remove(modTerm[i])

            j = 0
            nameSplit = list()
            while j < len(modTerm):
                splits = modTerm[j].replace("/", " ").replace("\\", " ").replace("-", " ").split(" ")
                k = 0
                while ((k < len(splits)) and (len(splits[k]) > 0) and (splits[k] not in insigWords)):
                    newWord = splits[k].strip()
                    nameSplit.append(newWord)
                    siteWordCount[siteKey][newWord] += 1
                    k = k + 1
                j = j + 1

            if siteKey not in cleanedList.keys():
                cleanedList[siteKey][term] = nameSplit
            if term not in cleanedList[siteKey].keys():
                cleanedList[siteKey][term] = nameSplit
    
    for site in siteWordCount.keys():
        siteTotalWordCount[site] = sum(siteWordCount[site].values())
        
    if dataType == "testNames":
        if REJECTION_THRESHOLD is not None:
            filter_out_frequent_tokens(cleanedList, siteWordCount, siteTotalWordCount, discardedTerms)
        cleanedList = convert_to_df(cleanedList, dataType)
        if config.write_file_source_data_cleaning == 'Y':
            cleanedList.to_csv(config.out_dir + "\\Cleaned_Lab_Names.csv", sep='|', index=False)
            write_word_ct_csv(config.out_dir + "\\By_Site_Lab_Word_Count.csv", siteWordCount)
            if len(discardedTerms) > 0:
                write_discarded_terms(config.out_dir + "\\Discarded_Lab_Names.csv", discardedTerms)
    if dataType == 'specimenNames':
        cleanedList = convert_to_df(cleanedList, dataType)
        if config.write_file_source_data_cleaning == 'Y':
            cleanedList.to_csv(config.out_dir + "\\Cleaned_Specimen_Names.csv", sep='|', index=False)
            write_word_ct_csv(config.out_dir + "\\By_Site_Specimen_Word_Count.csv", siteWordCount)
    return cleanedList


# In[5]:

def convert_to_df(cleanedList, dataType):
    if dataType == "testNames":
        cols = ['Site', 'OriginalTestName', 'CleanedTestName']
    else:
        cols = ['Site', 'OriginalSpecimen', 'CleanedSpecimen']
    return pd.DataFrame(([outer_key, inner_key, 
            " ".join(cleanedList[outer_key][inner_key])] for outer_key in cleanedList.keys() 
            for inner_key in cleanedList[outer_key].keys()),
            columns=cols)


# In[6]:

def filter_out_frequent_tokens(cleanedTestNameList, siteWordCount, siteTotalWordCount, discardedTerms):
    for site in siteWordCount.keys():
        for token in siteWordCount[site].keys():
            siteWordCtPct = 100.0 * siteWordCount[site][token] / siteTotalWordCount[site]
            if (siteWordCtPct > REJECTION_THRESHOLD) and (token != "%") and (token != "#"):
                for key in cleanedTestNameList[site].keys():
                    if token in cleanedTestNameList[site][key]:
                        cleanedTestNameList[site][key].remove(token)
                if ((site not in discardedTerms.keys()) or (token not in discardedTerms[site])):
                    discardedTerms[site].append(token)


# In[7]:

def write_cleaned_terms(pathName, data):
    data.to_csv(pathName, sep='|')


# In[8]:

def write_word_ct_csv(pathName, data):
    with open(pathName, 'w') as out_file:
        out_file.write("Site|Term|Count|Percent\n")
        for site in data.keys():
            total_num = sum(data[site].values())
            for word in data[site].keys():
                count = data[site][word]
                percent = 100.0 * data[site][word] / total_num
                out_file.write(site + "|" + word + "|" + str(count) + "|" + str(percent) + "\n")


# In[9]:

def  write_discarded_terms(pathName, discardedTerms):
    with open(pathName, 'w') as out_file:
        out_file.write("Site|DiscardedName\n")
        for site in discardedTerms.keys():
            for term in discardedTerms[site]:
                out_file.write(site + "|" + term + "\n")

