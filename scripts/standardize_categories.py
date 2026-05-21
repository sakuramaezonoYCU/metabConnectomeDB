"""
standardize_categories.py
-------------------------
A modular utility script to unify and standardize chemical superclass and subclass categories 
(the values in the 'Super_Class' and 'Sub_Class' columns) across the MetabConnectomeDB pipeline.
Dynamically loads standard category mappings from sister CSV files.
"""

import os
import pandas as pd
import numpy as np

# Dynamically load the mapping dictionary for Super_Class
SUPERCLASS_UNIFICATION_MAP = {}
superclass_csv_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'input', 'superclass_mapping.csv')

try:
    if os.path.exists(superclass_csv_path):
        df_map = pd.read_csv(superclass_csv_path)
        # Store clean lowercase keys for robust case-insensitive lookup
        for _, row in df_map.iterrows():
            if pd.notna(row['Raw_Category']) and pd.notna(row['Standard_Category']):
                raw_val = str(row['Raw_Category']).strip().lower()
                std_val = str(row['Standard_Category']).strip()
                if raw_val and std_val:
                    SUPERCLASS_UNIFICATION_MAP[raw_val] = std_val
except Exception as e:
    print(f"Warning: Failed to load superclass mapping CSV from {superclass_csv_path}: {e}")


# Dynamically load the mapping dictionary for Sub_Class
SUBCLASS_UNIFICATION_MAP = {}
subclass_csv_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'input', 'subclass_mapping.csv')

try:
    if os.path.exists(subclass_csv_path):
        df_map_sub = pd.read_csv(subclass_csv_path)
        # Store clean lowercase keys for robust case-insensitive lookup
        for _, row in df_map_sub.iterrows():
            if pd.notna(row['Raw_Category']) and pd.notna(row['Standard_Category']):
                raw_val = str(row['Raw_Category']).strip().lower()
                std_val = str(row['Standard_Category']).strip()
                if raw_val and std_val:
                    SUBCLASS_UNIFICATION_MAP[raw_val] = std_val
except Exception as e:
    print(f"Warning: Failed to load subclass mapping CSV from {subclass_csv_path}: {e}")


def standardize_superclass(val):
    """
    Standardizes a chemical superclass taxonomy category name string.
    Handles null values, whitespace cleaning, and composite (piped) category strings.
    """
    if pd.isna(val) or not isinstance(val, str):
        return np.nan
        
    val_clean = val.strip()
    if not val_clean:
        return np.nan
        
    # Split composite values (e.g. "Lipids | Organic nitrogen compounds")
    parts = [p.strip().lower() for p in val_clean.split('|') if p.strip()]
    
    if not parts:
        return np.nan
        
    # Attempt to resolve the first recognized category component
    for p in parts:
        if p in SUPERCLASS_UNIFICATION_MAP:
            return SUPERCLASS_UNIFICATION_MAP[p]
            
    # Fallback to the first capitalized clean part if not explicitly mapped
    first_part = parts[0]
    # Sentence-case fallback
    if len(first_part) > 1:
        return first_part[0].upper() + first_part[1:]
    return first_part.upper()


def standardize_subclass(val):
    """
    Standardizes a chemical subclass taxonomy category name string.
    Handles null values, whitespace cleaning, and composite (piped) category strings.
    """
    if pd.isna(val) or not isinstance(val, str):
        return np.nan
        
    val_clean = val.strip()
    if not val_clean:
        return np.nan
        
    # Split composite values (e.g. "Fatty acids | Amines")
    parts = [p.strip().lower() for p in val_clean.split('|') if p.strip()]
    
    if not parts:
        return np.nan
        
    # Attempt to resolve the first recognized category component
    for p in parts:
        if p in SUBCLASS_UNIFICATION_MAP:
            return SUBCLASS_UNIFICATION_MAP[p]
            
    # Fallback to the first capitalized clean part if not explicitly mapped
    first_part = parts[0]
    # Sentence-case fallback
    if len(first_part) > 1:
        return first_part[0].upper() + first_part[1:]
    return first_part.upper()
