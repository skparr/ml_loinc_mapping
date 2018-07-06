# ml_loinc_mapping

The ml_loinc_mapping module was created to map laboratory tests to their respective LOINC codes by building a model that leverages noisy labels in existing laboratory data.

## Dependencies
### Python
- Python (>= 3.5)
- numpy (>= 1.11.0)
- scikit-learn (>= 0.17.1)
- networkx=1.11
- rpy2
- hyperopt

### R 
- stringdist (used in Python code for string distance matching via Python rpy2 package)

### UMLS
- User will need to have a UMLS username/password login to obtain a UMLS api key (detailed in configuration file)

### LOINC
- User will need to download a copy of the publicly-available loinc.csv file from https://loinc.org/downloads/loinc/ into the working directory.

## Configuration
Please carefully review the config.py file. 
This program assumes that you have collected raw laboratory data and aggregated it, grouping by [Site, Lab Test Name, Specimen Type, Units (missingness allowed), and locally-recorded LOINC code (missingness allowed)].
Within the above grouping, the user is expected to have created summary measures of the numeric test results within the aggregated groups. The summary measures used as features in the algorithm include [Mean, Minimum, Maximum, 5th Percentile, 25th Percentile, Median, 75th Percentile, 95th percentile, and Count]
Because laboratory test names may contain commas, it is recommended that the user ensure that the delimiter in their raw data file is unique (i.e. '|') 
The aggregate source data should be stored in a .txt, .csv, or .xls file
### The following fields in the config.py file require the user to enter information:
- out_dir
- in_file
- loinc_file_path
- lib_loc (for R strindist package)
- test_col (see config.py file)
- spec_col (see config.py file)
- units (see config.py file)
- loinc_col (see config.py file)
- min_col (see config.py file)
- max_col (see config.py file)
- mean_col (see config.py file)
- perc_5 (see config.py file)
- perc_25 (see config.py file)
- median_col (see config.py file)
- perc_75 (see config.py file)
- perc_95 (see config.py file)
- count (see config.py file)
- site (see config.py file)
- missing (see config.py file)
- api_key

### The following fields in the config.py file may optionally be configured:
- delim
- write_file_source_data_cleaning
- write_file_loinc_parsed
- write_file_umls_cuis
- rejection_threshold
- num_cuis
- min_sites_per_loinc_key
- min_tests_per_loinc_group
- min_row_count_per_loinc_group
- run_cv
- n_splits
- tuning_evals
- max_features 
- max_depth
- min_samples_split
- n_estimators 



## Following Configuration:
Once all dependencies are installed and configurations are completed in the config.py file, the user needs only to execute the **Analysis.py** file, which fits the machine learning model, makes predictions on both labeled and unlabeled data, and outputs results to the **labeled_data_predictions.csv** and **unlabeled_data_predictions.csv** files.
