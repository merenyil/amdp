# -*- coding: utf-8 -*-
"""
geomag_utils_preferences.py
Created on Mon Jul 27 20:12:01 2026

@author: LaszloM
"""

import yaml
import misc_utils as mu

magnetic_comp_plot_colors = None
FIGURE_SIZE = (10, 3)   # (width, height)
TIMELABEL_PREC = 5
nT_PREC = 2
nT_RESID_PREC = 2 
INCL_DECL_DEG_PREC = 5
# MAG_DATA_DEFAULT_BASE_PATH = None
UI_DATETIME_FORMAT_SEC = None
UI_DATETIME_FORMAT_MIN = None


def load_preferences(_geomag_utils_preferences_yaml_path):
    with open(_geomag_utils_preferences_yaml_path, "r", encoding="utf-8") as _f:
        _preferences = yaml.safe_load(_f)
    
    global magnetic_comp_plot_colors, FIGURE_SIZE, TIMELABEL_PREC, nT_RESID_PREC, nT_PREC, INCL_DECL_DEG_PREC, MAG_DATA_DEFAULT_BASE_PATH, UI_DATETIME_FORMAT_SEC, UI_DATETIME_FORMAT_MIN
    magnetic_comp_plot_colors = _preferences.pop("magnetic_comp_plot_colors")
    
    FIGURE_SIZE = tuple(_preferences.pop("FIGURE_SIZE"))  # (width, height)
    
    TIMELABEL_PREC      = _preferences.pop("TIMELABEL-DECIMAL-DIGITS")      # Number of decimal digits in float type ('Excel encoded') time labels, when printing/saving them
    nT_RESID_PREC       = _preferences.pop("nT-RESID-DECIMAL-DIGITS")       # Number of decimal digits of residuals (in nT), when printing/saving them
    nT_PREC             = _preferences.pop("nT-DECIMAL-DIGITS")             # Number of decimal digits of magnetic values (X/Y/Z) (in nT), when printing/saving them
    INCL_DECL_DEG_PREC  = _preferences.pop("INCL-DECL-DECIMAL-DIGITS")  # Number of decimal digits of Inclination/declination
    # M_SCALE_PREC        = 6
    
    UI_DATETIME_FORMAT_SEC = _preferences.pop("TIMELABEL-FORMAT_SEC","%Y.%m.%d %H:%M:%S")
    UI_DATETIME_FORMAT_MIN = _preferences.pop("TIMELABEL-FORMAT_MIN","%Y.%m.%d %H:%M:%S")
    
    mu.DEFAULT_DATETIMESTR_FORMATS = [UI_DATETIME_FORMAT_SEC, UI_DATETIME_FORMAT_MIN ]
   
    return _preferences