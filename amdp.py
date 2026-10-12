# -*- coding: utf-8 -*-
"""

admp.py -  Absolute Magnetic Data Processing

Python software package to process ADS/DIM measurements data.
@author: Laszlo Merenyi and AdAm Domjan, Mingeo Ltd., Hungary

"""

import os
import argparse

import yaml
import dim_data_process
import geomag_utils
import geomag_obs_dataset
import geomag_utils_preferences
import dim_instr
import dim_azimuth
import dat_file_reader
import misc_utils as mu
import cli_utils as cu

program_version = "0.3"
program_name = "amdp - Absolute Magnetic Data Processing Software v.{}".format(program_version)

print(program_name)       
print("Copyright (c) 2026 Laszlo Merenyi, Adam Domjan, Mingeo Ltd.")
print("Licensed under the MIT License. See LICENSE file for details.")
print("")


# --- Optional command-line overrides for the YAML configuration files ---
# All are optional; if not given on the command line, the previous
# hard-coded default locations are used, so existing usage (running
# amdp.py with no arguments) keeps working unchanged.

arg_parser = argparse.ArgumentParser(description=program_name)
arg_parser.add_argument("main_config_file_path", nargs="?", default = None)
arg_parser.add_argument("--run-parameters", dest="run_config_file", default=None,
                         help="Path to the conversion settings YAML file (default: amdp-run-parameters.yaml)")
arg_parser.add_argument("--piers-data", dest="piers_data_file", default=None,
                         help="Path to the piers-data YAML file "
                              "(default: piers-data.yaml next to amdp.py)")
arg_parser.add_argument("--observatory-datasets", dest="observatory_datasets_file", default=None,
                         help="Path to the observatory-datasets YAML file (default: observatory-datasets.yaml)")
arg_parser.add_argument("--ui-preferences", dest="ui_preferences_file", default=None,
                         help="Path to the ui-preferences YAML file (default: ui-preferences.yaml)")
cli_args = arg_parser.parse_args()

config_files_def_dir = "" # The current directory, or, for the program directory, set it to "os.path.dirname(__file__)"
config_files_def_dir = cli_args.main_config_file_path or ""

run_config_file = cli_args.run_config_file or os.path.join( config_files_def_dir , "amdp-run-parameters.yaml"  )

def get_config_file_with_path(  file_name, cli_arg  ):
    if cli_arg is not None:
        return cli_arg  # Config file explicitly defined in a program argument
    else:
        result = os.path.join(config_files_def_dir, file_name)
        # If user supplied a custom path, but this config file is not there, then we revert to the programs dir:        
        if config_files_def_dir != "" and (not os.path.exists(result)):
            result = os.path.join("", file_name)
        return result
    
def data_files_relative_to_absolute_path(path):
    if (len(path) > 2) and (path[0] == '.') and (path[1] in ['/','\\'] ):
        # dat file's directory spec starting with "." will expand to the directory of the YAML file, and not to the current directory!
        global run_config_file_path
        return os.path.join( run_config_file_path,  path[2:])  
    
    else:
        return path # No conversion
    
    
piers_data_file = get_config_file_with_path( "piers-data.yaml", cli_args.piers_data_file)
observatory_datasets_file = get_config_file_with_path( "observatory-datasets.yaml", cli_args.observatory_datasets_file)
ui_preferences_file = get_config_file_with_path("ui-preferences.yaml", cli_args.ui_preferences_file)

geomag_utils.pillar_data_init( piers_data_file )
geomag_obs_dataset.read_obs_datasets( observatory_datasets_file )
preferences_yaml = geomag_utils_preferences.load_preferences( ui_preferences_file )

# Optional override of the DAT-file cache directory (dat_file_reader.CACHE_DIR),
# read from "data-cache-dir" in ui-preferences.yaml. If the key is not present
# (or is empty), then cache dir is not used.
dat_file_reader.CACHE_DIR =  preferences_yaml.pop("data-cache-dir", None)
    
assert os.path.isfile(run_config_file), "Program-run configuration file \"{}\" does not exsit!".format(run_config_file)

run_config_file_path =  os.path.dirname( os.path.abspath(run_config_file) )

with open(run_config_file, 'r') as f:
    run_config = yaml.safe_load(f)
    
title = run_config.pop("title", None)
if title is not None: 
    cu.print_framed(title,'-',3)
    print()
    

# We accept both "pier"  and "pillar" names
if "pier" in run_config.keys():
    assert not "pillar" in run_config.keys()
    run_config["pillar"] = run_config["pier"]
    
meas_pos = geomag_utils.measurement_pos_class(  run_config.pop("site"), run_config.pop("pillar") ) 
mag_obs_dataset = geomag_obs_dataset.mag_obs_dataset( run_config.pop("mag-obs-dataset") ) 



DIM_instr = dim_instr.DIM_instrument_model_class(run_config.pop("DIM-instr"))

DIM_data_files_format = run_config.pop("DIM-data-files-format", "ADS-raw-dat-files")
DIM_data_files = run_config.pop("DIM-data-files")
DIM_data_files_path, DIM_data_files_filenames = os.path.split(DIM_data_files)

DIM_data_files_path = data_files_relative_to_absolute_path(DIM_data_files_path)



# Base configuration of azimuth mark measurement data source
azimuth_data_files = run_config.pop("DIM-azimuth-readings-data-files", "$ADS-data-files")
if azimuth_data_files ==  "$ADS-data-files":
    def_azimark_meas_data_source_def = dim_azimuth.azimark_measured_data_source_def_class(dim_azimuth.azimark_meas_data_source_type_enum.dim_dat_files,  None, config_dict = run_config)
else:
    azimuth_data_files=data_files_relative_to_absolute_path(azimuth_data_files) 
    
    def_azimark_meas_data_source_def = dim_azimuth.azimark_measured_data_source_def_class(dim_azimuth.azimark_meas_data_source_type_enum.separate_dat_files, azimuth_data_files, config_dict = run_config, azimark_name = run_config.pop("azimuth-mark")   )
   
calc_type =  run_config.pop("calculation-type","spot-baselines")
assert calc_type == "spot-baselines", f'Unsupported calculation method: {calc_type}! (This version supports only "spot-baselines" calculation method!)'

calc_params_list = []
not_used_calculation_params_indeces = []
for list_item_index, calc_params_dict in enumerate(run_config["calculation-parameters"]):
   if calc_params_dict.pop("use", True): 
       calc_params = dim_data_process.spot_baselines_calculation_parameters_class( 
               calc_params_dict,
               meas_pos = meas_pos,            
               DIM_instr = DIM_instr,
               # baseline_calc_approach = baseline_calc_approach,
               refmark_measured_data_sources = [ def_azimark_meas_data_source_def ],               
               cycle_point_filters=[],  # [ ["CYCLE_POINT_FILTER__EXCLUDE_POS_IDS", (26,51)] ] # To process null method withouth the scale check point, and withouth any extra point
             ) 
       calc_params_list +=  [ calc_params ]
   else:
       not_used_calculation_params_indeces += [ list_item_index ] 

for list_item_do_remove in  reversed(not_used_calculation_params_indeces):
    del run_config["calculation-parameters"][list_item_do_remove]
       
assert len(calc_params_list)  > 0, "There is no active calculation parameter setting!"    


class spot_baselines_common_calculation_parameters_class:
     variometer_time_shift_sec = None

     F_time_shift_sec = None  # ToDo: move it to 'spot_baselines_calculation_parameters_class'
     variometer_slopecalc_halfwindow_sec = 0
     dI_pillar_diff = 0
     dD_pillar_diff = 0
     F_pillar_diff_variometer = 0  # "F_absolute_pilar - F_variometer_pilar"  If positive, then F in the absolute pillar is higher than F in absolute pillar
     do_self_checks = True
     
     def __init__(self, params):
         
        self.DIM_time_shift_sec = None
        self.DIM_time_shift_sec = None
        self.data_selection_filters = params.pop("data-selection-filters", None ) # If not defined, then it should be an empty dictionary
        if self.data_selection_filters is None: self.data_selection_filters = {} #  If key is exist in YAML file, but its items are commented out, then it will set to None. In that case, this should be an empty dict instead.
        self.F_pillar_diff_scalar = params.pop("F-pillar-diff[nT]",0.0)
        self.remove_means = params.pop("remove-means", None)     # None or list of keywords: $baselines
                                                    # Used to set baseline values (and offset, scale, etc. curvews) will be around 0.
                                                    # It can be used to easily compare baseline and instrument parameter variatoins at different settings.
        if self.remove_means is not None: self.remove_means = mu.ensure_list(self.remove_means)
        self.post_processing_steps = params.pop("post-processings")
        assert type(self.post_processing_steps) == list        
        
        self.ext_magnetometer_data_files_reader = None

        if "DIM-data-files-overrides" in params.keys():
            assert type(params["DIM-data-files-overrides"]) is dict
        #    params["DIM-data-files-overrides"] = params["DIM-data-files-overrides"][0]
            ext_magnetometer_data_files_def = params["DIM-data-files-overrides"].pop("use-magnetometer-data-from-files", None)
            if ext_magnetometer_data_files_def:
                # DIM M data is coming from a separate, external file
                if (type(ext_magnetometer_data_files_def) is dict):
                    # external M data file is defined by a dict
                    if ext_magnetometer_data_files_def.get("components"):
                        # ToDo: check if only one components is defined in YAML file
                        pass
                    else:
                        # Defining "components" in YAML is not mandatory 
                        ext_magnetometer_data_files_def["components"] = [ "M-ext" ]
                    colindex = ext_magnetometer_data_files_def.pop("M-colindex")
                    ext_magnetometer_data_files_def["column-indeces"] = [colindex]
                self.ext_magnetometer_data_files_reader = dat_file_reader.dat_file_reader_class_from_config_def(ext_magnetometer_data_files_def)
        
        # MNY_temp_reg = np.genfromtxt("p:/dat/MNY/env-mon/2025-temp-B_fixcols.dat", names = True, usecols = [0, 5])
        # MNY_temp_reg = MNY_temp_reg[  ~np.isnan(MNY_temp_reg["T_diddsensor"]) ]

        
        self.remaining_user_params = params
        self.baseline_calc_approach = params.pop("baseline-calc-approach")
        self.resid_std_acceptance_limit = params.pop("residual-stddev-limit[nT]", 10.0)
        

            


         
calc_params_common = spot_baselines_common_calculation_parameters_class(run_config)

processing_window_results, pointwise_results, saved_file_names, warn_logs  = dim_data_process.dim_data_process(DIM_data_files_format, DIM_data_files_path, DIM_data_files_filenames, meas_pos, DIM_instr, mag_obs_dataset, calc_params_common, calc_params_list, preferences_yaml, aux_data = None )

# Listing the not used YAML entries.
# But first, removing empty subitems:

def remove_empty_items(data):
    
    def _is_empty(value):
        return isinstance(value, (dict, list)) and len(value) == 0

    """Recursively remove empty lists and dicts from parsed YAML data."""
    if isinstance(data, dict):
        cleaned = {k: remove_empty_items(v) for k, v in data.items()}
        return {k: v for k, v in cleaned.items() if not _is_empty(v)}
    elif isinstance(data, list):
        cleaned = [remove_empty_items(v) for v in data]
        return [v for v in cleaned if not _is_empty(v)]
    else:
        return data
    
run_config =   remove_empty_items(run_config)  
            
if len(run_config) > 0:
    print("Not used calculation settings:")
    for key, value in run_config.items():
        print(f'- {key}: {value}')
        
if len(preferences_yaml) > 0:
    print("Not used UI preferences:")
    for key, value in preferences_yaml.items():
        print(f'- {key}: {value}')
        
            
warn_logs += mag_obs_dataset.get_warnings()
        
if len(warn_logs)>0:
    print()
    print("Warning{}:".format( mu.plural_suffix(len(warn_logs)) )  )
    for warn_log in warn_logs:
        print("- "+warn_log)
  
if saved_file_names and len(saved_file_names)>0:    
    print()
    print("Results saved in file{}:".format( mu.plural_suffix(len(saved_file_names)) )  )    
    for file_name in saved_file_names:
        print(f'- "{file_name}"')
      