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

### R (R stringdist package is used in Python code for string distance matching via Python rpy2 package)
- stringdist
