# -*- coding: utf-8 -*-
"""

dim_measurement_loader.py  

@author: Merényi László

Reads and pre-processes (separates into processing windows) DIM data files.
Supported data files:
    - ADS "raw" data files saved by "ADS-control"
    - Mingeo XLSX files
    - MagPY v1.0 files
    

"""

import geomag_utils as gu
import misc_utils as mu
import cli_utils  as cu
import plot_utils
import matplotlib.pyplot as plt
import matplotlib
import re
import os
import numpy as np
import sys
import tabulate
import datetime

import dim_instr

import numpy.lib.recfunctions as rcf
import geomag_utils_preferences as preferences
import geomag_obs_dataset
import dim_azimuth
import absmeas_xlsx_import
import dim_utils as du


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

SUPPORTED_DIM_DATA_FILE_FOMRATS = [ "ADS-raw-dat-files", "MagPy-autodetect" ]

ALL_MEASUREMENTS_TOGETHER = -1  # Constant: this can be used as a parameter in some functions, to request all cycles points

# ToDo: ket csoportra venni a filtereket: point-filter (pl. mag min) / cycle-filter
# FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS = True

# FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN = True
# FILTER_VC_TICK_MIN = 200000
# FILTER_VC_TICK_MAX = 200000
# FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE = 2000  # Close to vertical position, the data processing is subject to higher error
# FILTER_VC_DEG_MIN = 50.0
# FILTER_VC_TICK_MIN = 18000
# FILTER_VC_DEG_MAX = 180.0
# FILTER_HC_DEG_MIN = 20
# FILTER_HC_DEG_MAX = 340
# FILTER_MAG1_MIN=-10000  # -23000 "real nT" to remove wrong MAG2 at direct ObsDaq CH
# FILTER_MAG1_MAX=19000
# FILTER_CYCLE_ID_MIN = 1
# FILTER_CYCLE_ID_MOD100_MIN = 3 # Used for joined cycles, where each new measurements-raw.dat starts with a new '100'
# FILTER_CYCLE_ID_MAX = 5
# FILTER_POS_ID_MAX = 22
# FILTER_CYCLE_ID_EXCLUDE_LIST=  [404]
# FILTER_CYCLE_ID_INCLUDE_LIST=  [403,404]
# FILTER_CYCLE_ID_MOD2 = 1
# FILTER_CYCLE_ID_REMOVE_MODULO = [2,1]
# FILTER_SKIP_FIRST_SAMPLES_ON_EACH_CYCLE = 5
# FILTER_METHOD_STR_INCLUDE = [ b"nullmeth*" ]
# FILTER_METHOD_STR_INCLUDE = [ b"sinus_VC=*" ]
# FILTER_METHOD_STR_INCLUDE = [ b"sinus_VC=*"  ]
# FILTER_METHOD_STR_EXCLUDE = [ b"NM-fast" ] # [ b"sinus_VC=*" , b"grid*", b"*-BWD" ]
# FILTER_METHOD_EXTRA_INFO_STR_EXCLUSIVE= [ b"*-BWD" ]
# FILTER_METHOD_EXTRA_INFO_STR_EXCLUSIVE= [ b"intermed" ]
# FILTER_EXCLUDE_BEFORE=mu.datestr_2_exceltime("2024.07.10 00:00")
# FILTER_EXCLUDE_AFTER=mu.datestr_2_exceltime("2024.06.01 00:00")
# FILTER_EXCLUDE_DATE_RANGES=[    ]

MIN_MAG_MEAS_POINTS_PER_ABMEAS_PROCESSING_WINDOW=8  # Non-repeating residual meausrement without scale check is taken by 8 points

# NULLMETHOD_METHOD_STRING_IDS_BYTEARRAY = [ b"nullmeth", b"NM", b"NM-fast" ]
MEAS_TYPE_STR__NULLMETH_1     = "NM"
MEAS_TYPE_STR__NULLMETH_2     = "ABSM" # From 2026.
MEAS_TYPE_STR__NULLMETH_DEF   = MEAS_TYPE_STR__NULLMETH_2
MEAS_TYPE_STR__MIRA_READING   = "RP"
DAT_FILE_NULLMETHOD_METHOD_STRING_IDS = [ "nullmeth", "NM", "NM-fast", "NM-slow" ]

ADS_DAT_FILE__REFPOS_DETERMINATION_TYPES_STR_LIST = [b"SD-new", b"SU-new", b"SU-user", b"SD-user", b"SD-last",b"SU-last",b"SD-inv",b"SU-inv" ]
ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID =  -91  # pos_id of reference mark reading within absolute meaurement points
COMPNAME__IS_M_RECORD = "isM"

MAX_FILE_PROGRAM_NAME_LEN=40
MAX_FPILLAR_NAME_LEN = 30

MEAS_FLAGS_EMPTY_FIELD      = "_" # This character is used for empty "meas_flag" component
MEAS_FLAGS_EMPTY_FIELD_b    = b"_"
MEAS_FLAGS_SEPARATOR        = "|"
MEAS_FLAGS_KEY_VALUE_SEPARATOR = ":"
# MEAS_FLAGS_KEYS__CYCLE_NAME = "CN"  Deprecated "meas_flag" key ("CycleName"), used for XLSX-to-ADS-dat-file conversion
MEAS_FLAGS_KEYS__PILLAR = "P"
# MEAS_FLAGS_KEYS__SENSOR_SN = "MSN"

COMPNAME__PROGRAMSTEP_TYPE_FULL_STR = "prg_step_type_str"  # A type for a record. ~COMPNAME__MEAS_TYPE_FULL_STR~
PROGRAMSTEP_TYPE_FULL_STR_MAX_LEN       = 50  # ToDo: check max len!

COMPNAME__TIMESTAMP          = "datetime"
COMPNAME__HC_TICK            = "hctick"
COMPNAME__VC_TICK            = "vctick"
COMPNAME__MEASURED_HC_RAD    = "meas-hc"   # Angles from XLSX and other (non ADS) dat files
COMPNAME__MEASURED_VC_RAD    = "meas-vc"   # Angles from XLSX and other (non ADS) dat files

COMPNAME__MEAS_COUNT_SINCESTART = "measurement_count"   # Integer
COMPNAME__POS_ID             = "pos_id"   # Integer
COMPNAME__MAG                = "MnT"
COMPNAME__F                  = "PPM_MnT"
# COMPNAME__PROGRAM_DEF_IDX__NOTUSED    = "PDI"  # Index of 'file_absmeas_window_data' array
COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX      = "cycle_id"  # Index from file
COMPNAME__MEAS_FLAGS         = "meas_flags"
COMPNAME__T_ELECTRONICS      = "T_electronics°C"
COMPNAME__T_ROOM             = "T_room°C"
COMPNAME__MOVETOPOINT_TIME   = "meas_durations"
COMPNAME__MAG_NOISE         = "M_avg_stddevnT"
COMPNAME__MAG_SLOPE          = "M_avg_slopenTsec"
MEAS_FLAGS_FULL_STR_LEN      = 40  # ToDo: check and use max "meas_flags" lengtth!

COMPNAME__FILE_ABMEAS_WIN_INDEX = "FI"
COMPNAME__PROGRAM_MAINTYPE_INDEX = "PMI" # Index of 'program_maintype_str_list' (refers to e.g. 'NM')
DATAWINDOW_COMPNAME__M_RECORD_INDECES = "M-rec-indeces"
DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES =  { "FL":"__rf_FL", "FR":"__rf_FR" }

COMPNAME__TELESCOPE_CARDINAL_DIRECTION = '__card_dir'  # North/South/East/West
COMPNAME__TELESCOPE_FACE_POSITION     = '__tel_pos'   # FL/FR
COMPTYPE__TELESCOPE_CARDINAL_DIRECTION = "U1"
COMPTYPE__TELESCOPE_FACE_POSITION = "U1"
    
    

# --- Created components --
# Component names enbding with "_" have calc-variation specific values

#COMPNAME__PRECORR_HC_RAD_             = "hcrad"   # Pre-corrected angle, for each calcvar
#COMPNAME__PRECORR_VC_RAD_             = "vcrad"   # Pre-corrected angle, for each calcvar


# COMPNAME__PROGRAM_MAINTYPE_STR = "PMTS"  #
# COMPNAME__PROGRAM_MAINTYPE_STR_MAX_LEN   = 12
# COMPTYPE__PROGRAM_MAINTYPE_STR   = f"U{COMPNAME__PROGRAM_MAINTYPE_STR_MAX_LEN}"

# COMPNAME__PROGRAM_SUBTYPE_IDX   = "PSTI" # index of 'program_mainsubtype_str_list' (e.g "NM/8*3")

COMPNAME__ADS_FILE__RECORD_WITHIN_PROGRAM  = "win_rec_nr"

COMPNAME__PROGRAM_MAINTYPE_IDX = "meas_type"   # In ADS raw data files, it is an index of the program name in self.program_maintype_str_list
                                               # E.g.: "ABSM"
COMPNAME__PROGRAM_PARAM1_IDX   = "meas_type_param1"  # In ADS raw data files, it is an index in 'self.meas_type_param1_str_list'. E.g. "E-FL-2" in case of mag records
COMPNAME__PROGRAM_PARAM2_IDX   = "meas_type_param2"
# COMPNAME__CYCLE_SEQUENCE_ID  = "seq_id"
# COMPNAME__EXTRA_FLAG         = "extra_flag"

# COMPNAME__REFPOS_HC_THEO_SU_TICK_ = "SU_hc_theo_tick"
# COMPNAME__REFPOS_HC_THEO_SD_TICK_ = "SD_hc_theo_tick"
# COMPNAME__REFPOS_HC_RAD_THEO_SU = "SU_hc_theo"  # Note: if calc.var. index is added, then it will be in form of COMPNAME__REFPOS_HC_THEO_SU + CALC_VAR_COMPNAME_EXTENSION + calcvarnr.
# COMPNAME__REFPOS_HC_RAD_THEO_SD = "SD_hc_theo"
COMPNAME__REFPOS_HC_WGS84_SU_  = "SU_hc_NED"  # The reference pos angle (in other name: WGS84)
COMPNAME__REFPOS_HC_WGS84_SD_  = "SD_hc_NED"

# COMPNAME__PRECALC_HC_NED_ =  "hc_NED"
# COMPNAME__PRECALC_VC_NED_ =  "vc_NED"

# COMPNAME__PROCESSING_CYCLE_INDEX_ = "proccycleidx"  # For each calc variation. Has '-1' if not filled, and -2, if the record belongs to more processing cycles
COMPNAME__INPUT_FILE_INDEX   = "file-index"
# COMPNAME__RAW_DATA_INDEX   = "rawdataindex"
# CALC_VAR_COMPNAME_EXTENSION = "_calcvar"
COMPNAME__MISC_USE  = "misc"
from dim_azimuth import COMPNAME__AZIMARK_INDEX
INV_REFMARK_INDEX = 0xFFFFFFFF # It is an uin32 type
COMPNAME__PILLAR_INDEX = "pillar_index"
COMPNAME__AZIMARK_WGS84 = "azidir"


FILE_ABSMEAS_WINDOW_DATA_DTYPE =  [ 
    ( COMPNAME__FILE_ABMEAS_WIN_INDEX, int), 
    ( COMPNAME__PROGRAM_MAINTYPE_INDEX, int), # Index of 'program_maintype_str_list'
    # ('M_start_record_index',int) ,   ('M_end_record_index',int),  # Used for separation in plots -> using DATAWINDOW_COMPNAME__M_RECORD_INDECES[0] / DATAWINDOW_COMPNAME__M_RECORD_INDECES[-1] instead
    # ('variation_idx', int),   # ToDo: remove 
    # ('name', "|S{}".format(MAX_FILE_PROGRAM_NAME_LEN) ),   # 
    ("pillar-index",int), 
     # ('index',int),  ('pillar_name', "|S{}".format(MAX_FPILLAR_NAME_LEN)
    (COMPNAME__INPUT_FILE_INDEX,int)  ,
    (DATAWINDOW_COMPNAME__M_RECORD_INDECES, "O"),   # it will be a list
    (DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FL"], "O") ,  # it will be a list
    (DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FR"], "O")   # it will be a list
    # Note: in case of XLSX input files, the components are extended in abmeas_xlsx_import.py
    ]    


class processing_window_data_class:

    def __init__(self):
        self.M_record_indeces = None  # A list of record indeces, including reference pos records, for each calculation variation
        # self.refpos_record_indeces = None
        # self.mag_only_record_indeces = None  # 'M_record_indeces' minus the 'refpos_record_indeces'
        # self.cycle_own_refpos_record_indeces = None
        # self:date_limits = None
        self.complete_azi_meas_data = None
        self.pillar_name = None
        self.name, self.first_M_timelabel, self.last_M_timelabel = None, None, None
        
    def get_problem(self):
        if self.excluded_by_data_selection_filter: return "date-exl."
        if self.complete_azi_meas_data is None: return "azi-dat"
        if self.missing_variometer_data: return "var-dat"
        if self.missing_F_data: return "F-dat"
        if self.excluded_by_sample_number: return "smpl-nr"
        return None

    # Function used in error messages, to identify processing cycle:
    #def get_identification_info_str(self):
    #    cycle_first_timelabel =  raw_data[ self.M_record_indeces[0]][COMPNAME__TIMESTAMP]
    #    return "started at {}".format(cycle_first_timelabel)


# def ensure_no_filter_selected():
#     for var in globals():
#         if var[0:7]=="FILTER_":
#             # if (var == "FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS") and (FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS==False): continue # We accept this filter with this setting, as it does not filter out anything
#             if (var == "FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN") and (FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN==False): continue
#             cu.warning( "There is an active ADS data file filter: {}".format(var))


# class ADS_data_absmeas_window_determination_conditions_class:

#     def __init__(self, at_pos_id_cont_change = False, at_new_measurement_on_the_same_pos = False, same_pos_delta_tick_limit = 100):
#         self.at_file_meas_window_index_change = True  # In  case of ADS data files: New window begins in records where input file has a cycle id change
#             # For other files: window index is from the manual measurement count, and the default separation methods is based in this index
#         self.at_measrement_count_backjump = True # file_format == "ADS-raw-dat-files"  # ADS data files have a measurement_count componet
#         self.at_pos_id_jump_back = True # file_format == "ADS-raw-dat-files"
#         self.at_pos_id_cont_change = at_pos_id_cont_change  # New measurement begins in records where input file has pos_id continuity change
#         self.at_new_measurement_on_the_same_pos = at_new_measurement_on_the_same_pos
#         self.same_pos_delta_tick_limit = same_pos_delta_tick_limit
#         self.at_big_time_jump = True
#         self.time_jump_limit_minutes = 400.0  # Changed from 20.0 To 400.0 in 2025.02.27, for some special MNY measuremnts

       
# Add components (if not yet exist) to a structured array
# def extend_structured_array(a, desired_dtype):
#     """
#     a: existing structured ndarray
#     desired_dtype: list of (name, format[, shape]) tuples describing all desired fields
#     """
#     existing_names = set(a.dtype.names)
#     missing = [field for field in desired_dtype if field[0] not in existing_names]

#     if not missing:
#         return a  # nothing to add

#     new_names = [f[0] for f in missing]
#     new_dtypes = [f[1:] if len(f) > 1 else "f8" for f in missing]  # keep format(+shape)

#     # default fill value 0 for each new field; append_fields needs actual data arrays
#     fill_data = [np.zeros(len(a), dtype=dt) for dt in new_dtypes]

#     b = rcf.append_fields(a, new_names, fill_data, dtypes=new_dtypes, usemask=False)
#     return b


# def warn_user_for_data_correction():
#    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#    print("!!!! USING CORRECTED MEASUREMENT DATA!!!!!!!!!!!")
#    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


# def correct_magdat_with_cosinus( omega, phase, ampl, offset, compindex_rotating_angle_in_deg ):
#     warn_user_for_data_correction()
#     for i in range(0, len(file_records[:,0])):
#         file_records[i,compindex_Mag] =   file_records[i,compindex_Mag]  + ( offset  + ampl * math.cos ( omega * math.radians(file_records[i,compindex_rotating_angle_in_deg]) + phase)  )

# def correct_degdat_with_plate1_corr():
#     warn_user_for_data_correction()
#     hc_corr_diff_arcsec = np.zeros(len(file_records[:,compindex_recalc_hc_deg]))
#     for i in range(0, len(file_records[:,0])):
#         orig_hc_deg = file_records[i,compindex_recalc_hc_deg]
#         file_records[i,compindex_recalc_hc_deg] =  dim_instr.plate1_corr_rawtick_to_deg(file_records[i,compindex_hc_raw_tick], False)
#         hc_corr_diff_arcsec[i] = 3600 * (orig_hc_deg  - file_records[i,compindex_recalc_hc_deg])
#         file_records[i,compindex_recalc_vc_deg] =  dim_instr.plate1_corr_rawtick_to_deg(file_records[i,compindex_vc_raw_tick], True)

#     plt.xlabel('hc (°)')
#     plt.ylabel('angle-corr (arcsec)')
#     plt.plot(  file_records[:,compindex_recalc_hc_deg], hc_corr_diff_arcsec, 'o', color='green');
#     plt.show()


# def correct_rotdat_with_coserr(rotation_axis, ampl_arcmin, phase_deg, omega = 1.0):
#     ampl_deg = ampl_arcmin / 60.0
#     for i in range(0, len(file_records[:,0])):
#         if rotation_axis == ROTATION_PLANE_HORIZONTAL:
#             file_records[i,compindex_recalc_hc_deg] =  file_records[i,compindex_recalc_hc_deg] + ampl_deg * math.cos ( math.radians(file_records[i,compindex_recalc_hc_deg]*omega + phase_deg) )
#         elif rotation_axis == ROTATION_PLANE_VERTICAL:
#             file_records[i,compindex_recalc_vc_deg] =  file_records[i,compindex_recalc_vc_deg] + ampl_deg * math.cos ( math.radians(file_records[i,compindex_recalc_vc_deg]*omega + phase_deg) )



# In ADS data files, a text column indicate the cycle's measurement method (sin/Null/etc.) and optionally one or two point-specific paramters, separated by "/".
# This function parses this field, and converts some old notations to new ones.
def separate_full_meas_type_str(full_meas_type_str):
    # assert type(full_meas_type_str == <class 'numpy.bytes_>'  # From files, string is read as bytes
    # Converting all "-" seaprator of pseudo 0 nT points (used untily 2023-02-23) of the ADS scanning Null method to the new separator of "/".
    if  (full_meas_type_str[0:14] == b"nullmeth-slow-") or (full_meas_type_str[0:14] == b"nullmeth-fast-"): full_meas_type_str=full_meas_type_str[0:13]+b"/"+full_meas_type_str[14::]

    # To make less plot groups in 3D process, we convert the grid information to extra information:
    if (full_meas_type_str[0:10] == b"grid(4x4)_"): full_meas_type_str= full_meas_type_str[:9] +  b"/" + full_meas_type_str[10::]

    fields = full_meas_type_str.split(b'/')


    if len(fields)==1:
        fields = fields + [ b"", b""]
    elif len(fields)==2:
        fields = fields + [ b""]
    else:
        assert len(fields) in [3,4], "Problem with meas-type field of \"{}\"".format(full_meas_type_str)

    # Reference pos definition: 3 items
    #    RP/RP-ID/<determination-flag>
    # Reference mark definition: 4 itemx
    #    RM/<REFPOS-ID>/<REFMARK-ID>/<determination-flag>

    # Quick conversion from "refmark" system to "refpos" system:
    if len(fields)==4:
        assert fields[0]==b"RM"
        fields[0]=MEAS_TYPE_STR__MIRA_READING # 'b"RP"
        fields[1]=fields[2]
        fields[2]=fields[3]

    assert len(fields) in [3,4]

    if fields[0]  in [ b'ED', b'WD', b'EU', b'WU', b'NU', b'SD', b'ND', b'SU', b'CH']:
        # Old (probably deprecated) notation of null method coordinates:
        fields[0]= MEAS_TYPE_STR__NULLMETH_DEF
    elif fields[0] in DAT_FILE_NULLMETHOD_METHOD_STRING_IDS:
        fields[0] = MEAS_TYPE_STR__NULLMETH_DEF
    # if fields[0] == b"NM" and  fields[1] == b"Scale-check": fields[1] = "" # We classify scale check points as normal null method point

    return fields[0].decode("utf-8"), fields[1].decode("utf-8"), fields[2].decode("utf-8")


# Parses the "meas_flag" component of ADS data records
# "meas_flags" can have one or more fields, separated by MEAS_FLAGS_SEPARATOR character
# A field has a format of  "<param_name>:<param_value>" .
# Result is a dictionaly
def separate_meas_flags(meas_flags, row_nr = -1):
    meas_flags = meas_flags.decode("utf-8")
    result_dict = {}
    if meas_flags == MEAS_FLAGS_EMPTY_FIELD: return result_dict
    for flag in meas_flags.split(MEAS_FLAGS_SEPARATOR):
        # print(">{}<".format(flag))
        if (flag[0] in ['+','-']) and (not ":" in flag):
            # MIRA measurement records can have Mira determination flags in the format of "-starting-value"  or "+20240810-1502 "
            result_dict['REFMARKDET_RES'] =flag
        elif flag=="HA":
            # An unknown reason flag in MIRA records NAQ dat files ()
            result_dict['UNKNOWN'] = flag
        else:
            s = flag.split(MEAS_FLAGS_KEY_VALUE_SEPARATOR)
            assert len(s) == 2,"{} {}".format(meas_flags,row_nr)
            [key,value] = s
            result_dict[key] = value
    return result_dict



class DIM_dat_file_reader_class:


    def __init__(self):
        # self.file_cycle_names = []
        # self.cycle_plot_color = []
        self.vc0_from_file = None
        self.pillars = []
        # self.processing_window_data = None  # It is moved under calc_vars
        self.mag_sensor_names           = []
        # self.processing_window_M_date_limits  = []  # Date limites of all calc. variations.
        self.load_main_obs_data_set_data = None
        self.load_main_obs_data_set_data_missing_record_indeces = []
        self.load_main_obs_data_set_data_with_window = []
        self.warn_logs = []
        self.warn_log_count = {}
        

    # Splits a line (passed as string) from the data file (the header line or a data line) to fields (as strings)
    def split_record_str(self, record_as_string):
        if self.record_separator==" ":
            # To treat more spaces as one separatore, this trick is required:
                return re.split(r'\s+', record_as_string)
        return record_as_string.split(self.record_separator)
    

    def apply_record_deselection_filter(self, remove_record_condittion,selection_filter_name,record_index):
        if remove_record_condittion:
            if not record_index in self.filtered_out_record_indeces:
                self.filtered_out_record_indeces = self.filtered_out_record_indeces + [record_index]
                if (selection_filter_name!=""):
                    if not selection_filter_name in self.applied_remove_filters:
                        self.applied_remove_filters[selection_filter_name] = 1
                    else:
                        self.applied_remove_filters[selection_filter_name] += 1
        return remove_record_condittion
    


    # Main function of the class
    
    def load(self, file_format, base_path, filenames, DIM_instr, calc_params_common, calc_vars, meas_pos):  
    
        assert file_format in SUPPORTED_DIM_DATA_FILE_FOMRATS  + list(absmeas_xlsx_import.supported_input_file_types), f'Unsupported DIM data file format: "{file_format}"'        
        self.file_format = file_format
        self.calc_params_common = calc_params_common
        self.calc_vars = calc_vars
        self.instr = DIM_instr
        self.meas_pos = meas_pos
        
        base_path = os.path.join(base_path, '')
        assert os.path.isdir(base_path), "DIM data file directory {:} does not exist".format(base_path)

        # if MIN_MAG_MEAS_POINTS_PER_ABMEAS_PROCESSING_WINDOW < 0: MIN_MAG_MEAS_POINTS_PER_ABMEAS_PROCESSING_WINDOW = DEF_MIN_RECORDS_PER_CYCLE
        # if max_nr_of_mag_sensors is None:  max_nr_of_mag_sensors = min_nr_of_mag_sensors
        # assert max_nr_of_mag_sensors >= min_nr_of_mag_sensors
        # if ADS_data_window_determination_conditions_class == None:
        
        # Depending on the YAML contents, 'calc_params_common.data_selection_filters '  a list of dictionaries, or just a dictionaly
        # Recently only the latter is supported
        data_selection_filters = calc_params_common.data_selection_filters
        assert type(data_selection_filters) == dict
        # Converting keys to uppercase:
        data_selection_filters = {k.upper(): v for k, v in data_selection_filters.items()}
            
        
        date_selection_filters = {}
        
        if 'FROM-DATE' in data_selection_filters.keys():
            from_date = mu.datetimestr_2_datetime(data_selection_filters.pop('FROM-DATE'))
            date_selection_filters["FROM-DATE"] = (from_date,  mu.datetime_2_exceltime(from_date))
        if 'TO-DATE' in data_selection_filters.keys():
            to_date = mu.datetimestr_2_datetime(data_selection_filters.pop('TO-DATE'))
            date_selection_filters["TO-DATE"] = (to_date,   mu.datetime_2_exceltime(to_date))
        date_selection_filters['EXCLUDE-DATE-RANGES'] = []            
        date_selection_filters['INCLUDE-DATE-RANGES'] = []   

        self.exclude_windows_with_dates = []
        if 'EXCLUDE-DATES' in data_selection_filters.keys():
            excude_dates = data_selection_filters.pop('EXCLUDE-DATES')
            excude_dates = mu.ensure_list(excude_dates)
            for exclude_date in excude_dates:
                if type(exclude_date) != float:  # Both float format ('46270.5166') and string format is accepted
                    exclude_date = mu.datetimestr_2_exceltime(exclude_date)
                self.exclude_windows_with_dates += [ exclude_date ]
        
         
        if 'EXCLUDE-DATE-RANGES' in data_selection_filters.keys():
            for date_range in data_selection_filters['EXCLUDE-DATE-RANGES']:
                assert (len(date_range) == 2) and isinstance(date_range[0], float) and isinstance(date_range[1], float)
                assert date_range[1] > date_range[0]
                date_selection_filters['EXCLUDE-DATE-RANGES'] += ( date_range[0], date_range[1])
        if 'INCLUDE-DATE-RANGES' in data_selection_filters.keys():
            for date_range in data_selection_filters['INCLUDE-DATE-RANGES']:
                assert (len(date_range) == 2) and isinstance(date_range[0], float) and isinstance(date_range[1], float)
                assert date_range[1] > date_range[0]
                date_selection_filters['INCLUDE-DATE-RANGES'] += ( date_range[0], date_range[1])
                
        self.nr_of_M_records_outside_of_file_absmeas_windows = 0
        
        filename_masks = mu.ensure_list(filenames)
        self.applied_remove_filters = dict()  # Dictionary with key/value pairs of "filter-name" / "removed-records-due-to-filter"        
    
        if file_format != "ADS-raw-dat-files":
            # In case of 'ADS-raw-dat-files' format, 'date_selection_filters' is parsed and applied inside 'self.read_in_ads_raw_dat_files()' function
            date_from = date_selection_filters.get("FROM-DATE", [datetime.datetime(1970,1,1)]) [0]  # The Excel-time number is on the first element of the list, belonging to 'FROM-DATE" key.
            date_to   = date_selection_filters.get("TO-DATE",   [datetime.datetime(2999,1,1)]) [0]
        
        match file_format:
           case "ADS-raw-dat-files":
                ads_raw_data = self.read_in_ads_raw_dat_files(base_path, filename_masks, date_selection_filters)
                assert len(self.input_file_data)> 0  # 'self.input_file_data' is build by the 'read_in_ads_raw_dat_files()' function
                self.file_absmeas_window_data, self.M_readings, refmark_readings_data_from_DIM_dat_files = self.separate_ads_raw_dat(ads_raw_data)
           case "Mingeo-manual-XLSX":
                # Excel reading function needs from/to date?
                self.input_file_data,  self.file_absmeas_window_data, self.M_readings, refmark_readings_data_from_DIM_dat_files, self.azimuth_mark_list, self.program_maintype_str_list = absmeas_xlsx_import.read_absmeas_XLS_files(meas_pos, date_from, date_to, base_path = base_path, filename_mask = filename_masks, input_file_types = [file_format], DIM_instr = DIM_instr )
                # Addings some ADS specific columns:
                assert len(self.program_maintype_str_list) == 1
           case "MagPy-autodetect":
                import magpy_wrapper
                self.input_file_data,  self.file_absmeas_window_data, self.M_readings, refmark_readings_data_from_DIM_dat_files, self.azimuth_mark_list, self.program_maintype_str_list = magpy_wrapper.load_DI_measurement_files(meas_pos.pillars, date_from, date_to, base_path, filename_masks)
              
           case _ :
                assert False, f'DIM file format "{file_format}" is not supported'
                  
        if not COMPNAME__PROGRAM_MAINTYPE_IDX in self.M_readings.dtype.names:
            # 'COMPNAME__PROGRAM_MAINTYPE_IDX' is read only from ADS data files
            self.M_readings = rcf.append_fields(self.M_readings, COMPNAME__PROGRAM_MAINTYPE_IDX, np.zeros( len(self.M_readings) ) , np.uint32, usemask=False  )
            
        # ----------------------------------------------------------
        # ----  Pre-Processing DIM data  ---------------------------
        # ----------------------------------------------------------

        if calc_params_common.DIM_time_shift_sec is not None:
            self.M_readings[COMPNAME__TIMESTAMP] += float(calc_params_common.DIM_time_shift_sec)/86400.0
            for FLFR in ["FL","FR"]:
                refmark_readings_data_from_DIM_dat_files[FLFR][COMPNAME__TIMESTAMP] += calc_params_common.DIM_time_shift_sec / 86400.0

        if calc_params_common.ext_magnetometer_data_files_reader:
            # Replacing original DIM mag data with mag data read from external files (for two-sensor ADS tests)
            errmsg, aux_mag_timelabels, aux_mag, aux_mag_noncomplete_record_indeces, calc_stat = calc_params_common.ext_magnetometer_data_files_reader.get_data(self.M_readings[COMPNAME__TIMESTAMP]) 
            assert errmsg is None 
            aux_mag = aux_mag.squeeze(axis=1)
            aux_mag_nans = np.where( np.isnan(aux_mag) )[0]
            aux_mag_rangeouts = np.where( np.abs(aux_mag)   > 2000 ) [0]
            
            if len(aux_mag_nans)>0:
                for nan_idx in aux_mag_nans:
                    print(  mu.exceltime_2_datetimestr(self.M_readings[nan_idx][COMPNAME__TIMESTAMP]), self.M_readings[nan_idx][COMPNAME__MAG] )
                
            if len(aux_mag_rangeouts)>0:
                aux_mag[aux_mag_rangeouts] = np.nan
                
            self.M_readings[COMPNAME__MAG][:] = aux_mag[:]
        
        self.nr_of_mag_sensors = 1

        for imagsensor in range(self.nr_of_mag_sensors):
             # self.compindeces_mag  += [len(file_comps_def)]
             self.mag_sensor_names +=  [ "M-{}".format(imagsensor+1) ]
             
             # Scale correction from 2026.05.28
             if DIM_instr.mag_sensors[imagsensor].M_scale_rawdata_corr is not None:
                 self.M_readings[COMPNAME__MAG] *= DIM_instr.mag_sensors[imagsensor].M_scale_rawdata_corr

        # M_scale temperature correction:
        # raw_data[:,self.compindeces_mag[0]] = raw_data[:,self.compindeces_mag[0]] * ( 0* 0.9111762908221671  + raw_data[:,self.compindex_T_electronics] * 0.0004328101072152061)
        # raw_data[COMPNAME__MAG] = raw_data[COMPNAME__MAG] * 20 * 1.13

        # M_offset temperature correction:
        # raw_data[:,self.compindeces_mag[0]] -= ( raw_data[:,self.compindex_T_electronics] - np.mean(raw_data[:,self.compindex_T_electronics])  ) * 1.95

        # Note: conversion with 'mag_sensor[0].M_scale_rawdata_corr' is implemented above, in 2026.05.28


        # -------------------------------
        # --   Configuring 'instr' ------
        # -------------------------------
        if self.instr.vc0 is None:
            assert isinstance( self.vc0_from_file, float), "vc0 is not defined!"
            self.instr.vc0 = self.vc0_from_file
        else:
            if self.vc0_from_file is not None:
                if not np.isclose( self.instr.vc0, self.vc0_from_file, atol = 1e-1):
                    self.add_warning(f"Predefined vc0 ({self.instr.vc0}) differs from DAT file's vc0 ({self.vc0_from_file}). Using the first one")
            
        # ---------------------------------------------------
        # --  Removing records with "all-in-one" steps ------
        # ---------------------------------------------------
        
        # Storing original record indeces before the record removals
        # raw_data[COMPNAME__MISC_USE] = np.arange( len(raw_data) )

        # Note: in case of Excel files, records are already selected by date filters.
        #       In case of ADS files, also 
        

        for date_range in date_selection_filters["EXCLUDE-DATE-RANGES"]:
            assert False
            # self.remove_deselected_records( raw_data, 'EXCLUDE-DATE-RANGES', (date_range[0] <= raw_data[COMPNAME__TIMESTAMP]) & (raw_data[COMPNAME__TIMESTAMP] <= date_range[1]) )
        for date_range in date_selection_filters["INCLUDE-DATE-RANGES"]:            
            assert False
            # self.remove_deselected_records( raw_data, 'INCLUDE-DATE-RANGES', ~(  (date_range[0] <= raw_data[COMPNAME__TIMESTAMP]) & (raw_data[COMPNAME__TIMESTAMP] <= date_range[1]) ) )  # Double parentheses is needed!

        # if 'VC-MIN[TICK]' in data_selection_filters.keys():
        #     remove_deselected_records('FILTER_VC_TICK_MIN',raw_data[COMPNAME__VC_TICK] >= data_selection_filters["VC-MIN[TICK]"])
        # if 'VC-MAX[TICK]' in data_selection_filters.keys():
        #     remove_deselected_records('FILTER_VC_TICK_MAX',raw_data[COMPNAME__VC_TICK] <= data_selection_filters["VC-MAX[TICK]"])
        # if 'FILTER_CYCLE_ID_MIN' in globals():
        #     remove_deselected_records("FILTER_CYCLE_ID_MIN", raw_data[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] < FILTER_CYCLE_ID_MIN )
        # if 'FILTER_CYCLE_ID_MAX' in globals():
        #     remove_deselected_records("FILTER_CYCLE_ID_MAX", raw_data[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] > FILTER_CYCLE_ID_MAX )
        # if "FILTER_POS_ID_MAX"  in globals():
        #     remove_deselected_records("FILTER_POS_ID_MAX", raw_data[COMPNAME__POS_ID] > FILTER_POS_ID_MAX )

        # if DIM_instr.ext_mag_sensor_datfile_reader is not None:
        #    # Reading MAG data from an other file
        #    # Reading the external mag sensor data before MAG sensor related filters:
        #    _timelabel,M_from_other_file, sink = base_instr.ext_mag_sensor_datfile_reader.get_data("EXACT_TIME", raw_data[COMPNAME__TIMESTAMP], in_case_of_missing_data = ["none","drop-error"][1])
        #    raw_data[COMPNAME__MAG][:] = M_from_other_file[...,0]

        # if ('FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN' in globals()) and FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN:
        #     assert self.nr_of_mag_sensors == 1
        #     # for imag in range(min_nr_of_mag_sensors):
        #     # Record is not a reference position measurement (where Mag should be NAN)
        #     remove_deselected_records('FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN', np.logical_and(  np.isnan( raw_data[COMPNAME__MAG] )  , raw_data[COMPNAME__POS_ID] != ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID) )


        # if ('FILTER_MAG1_MIN' in globals()): # NOTE: Mag data is NOT YET corrected by 'AUTEO_DAT_FILE_MAG_SCALE_PRECORR' at this point
        #     remove_deselected_records('FILTER_MAG1_MIN', np.logical_and( ~np.isnan(raw_data[COMPNAME__MAG]), raw_data[COMPNAME__MAG]<FILTER_MAG1_MIN ) )
        # if ('FILTER_MAG1_MAX' in globals()):
        #     remove_deselected_records('FILTER_MAG1_MAX', np.logical_and( ~np.isnan(raw_data[COMPNAME__MAG]), raw_data[COMPNAME__MAG]>FILTER_MAG1_MAX ) )

        assert len(self.M_readings) > 0, "No record left after all-in-one record filtering"

        # assert np.all( raw_data[COMPNAME__PROCESSING_CYCLE_INDEX_ + calc_vars[0].id] == -1 )
        
        # Checking if timelabels are monolithically increasing:
        check_timestamps = [ [self.M_readings ] , [ "magnetometer" ] ]
        if refmark_readings_data_from_DIM_dat_files is not None:
            check_timestamps[0] += [ refmark_readings_data_from_DIM_dat_files["FL"], refmark_readings_data_from_DIM_dat_files["FR"] ]
            check_timestamps[1] += [  "reference-mark-FL", "reference-mark-FR"]
            
        for readings, readings_tpye_name in zip( check_timestamps[0], check_timestamps[1] ):
            if len(readings) == 0: continue # Can happen if there is not refmark_readings_data_from_DIM_dat_files
            d_timeindeces = np.diff( readings[COMPNAME__TIMESTAMP] )
            accepted_time_jumps = (d_timeindeces >= 0)  # Pseudo refpos measurements may have the same time label inaccuracy
            nr_of_time_jumps = np.count_nonzero( ~accepted_time_jumps)  # - nof_read_NM_calculated_zeropos_records
            if (nr_of_time_jumps!=0):
                print("There is/are {} decreasing time labels in {} readings!".format(nr_of_time_jumps, readings_tpye_name) )
                print("Records with time back jumps: ")
                not_accepted_time_jump_indeces = np.flatnonzero(~accepted_time_jumps)  #  [i for i in range(len(accepted_time_jumps)) if not accepted_time_jumps[i]]
                for i in not_accepted_time_jump_indeces:
                    print("timelabels:{:.6f} => {:.6f}   jump:{:d} sec".format(readings[COMPNAME__TIMESTAMP][i], readings[COMPNAME__TIMESTAMP][i+1], round(d_timeindeces[i]*86400) ))
                    # Other possible report:                 assert False, "Time labels are not monolitically increasing! Time back jump:{:.3f} s at record nr. {}\nTiem labels: {} => {} ".format( drectime[min_index]* 86400.0, min_index,mu.exceltime_2_datetime( raw_data[:, COMPNAME__TIMESTAMP][min_index] ),
                # assert False
                cu.warning( "Time jumps in DIM data records, see above", serious = True);
        
        # Calculating some derived data:
        for calc_var in calc_vars:
            calc_var.M_reading_theo_angles = {}
            calc_var.main_refmark_reading_angles = None if refmark_readings_data_from_DIM_dat_files is None else { "FL": {}, "FR":{} }
            for hcvc_enum, hcvc_str in zip(dim_instr.axis_enum, ["hc","vc"]):
                if DIM_instr.axis_angle_data_source == "digital":
                    # ADS eccentricity correction:
                    source_compname = COMPNAME__HC_TICK if hcvc_str=="hc" else COMPNAME__VC_TICK
                    calc_var.M_reading_theo_angles[hcvc_str] = self.instr.axis[hcvc_enum].tick_to_rad( self.M_readings[source_compname], calc_var.angle_precorr_funcs[hcvc_str], 'calcvar-M' ) 
                    if refmark_readings_data_from_DIM_dat_files is not None:
                        calc_var.main_refmark_reading_angles["FL"][hcvc_str] = self.instr.axis[hcvc_enum].tick_to_rad( refmark_readings_data_from_DIM_dat_files["FL"][source_compname], calc_var.angle_precorr_funcs[hcvc_str] )
                        calc_var.main_refmark_reading_angles["FR"][hcvc_str] = self.instr.axis[hcvc_enum].tick_to_rad( refmark_readings_data_from_DIM_dat_files["FR"][source_compname], calc_var.angle_precorr_funcs[hcvc_str] )
                else:
                    source_compname = COMPNAME__MEASURED_HC_RAD if hcvc_str=="hc" else COMPNAME__MEASURED_VC_RAD
                    calc_var.M_reading_theo_angles[hcvc_str] = self.M_readings[source_compname]
                    if refmark_readings_data_from_DIM_dat_files is not None:
                        calc_var.main_refmark_reading_angles["FL"][hcvc_str] = refmark_readings_data_from_DIM_dat_files["FL"][source_compname]
                        calc_var.main_refmark_reading_angles["FR"][hcvc_str] = refmark_readings_data_from_DIM_dat_files["FR"][source_compname]

        # --------------------------------------------------------------------------
        # --   Processing reference mark measurements embedded in DIM data file -----
        # --------------------------------------------------------------------------
        # self.refpos_data_calculated_for_each_record = False
        # self.refpos_data_calculated_for_each_processing_cycle = False


        # del invalidated_refpos_measurement_flag, SU_pos_flag,  refpos_measurements

        # Selectring refpos readings with real or user forced measurements (param2 is "SU-new" or "SD-new" or "SU-user" or "SD-user")
        # real_refpos_measurement_flag = np.empty( (len(self.all_refpos_measurement)), dtype=bool)
        # real_refpos_measurement_flag[:] = False
        # for accepted_refmark_reading_type in ADS_DAT_FILE__REFPOS_DETERMINATION_TYPES_STR_LIST[:4]:
        #     if accepted_refmark_reading_type in self.meas_type_param2_str_list:
        #         index = self.meas_type_param2_str_list.index(accepted_refmark_reading_type)
        #         real_refpos_measurement_flag[ self.all_refpos_measurement[ COMPNAME__PROGRAM_PARAM2_IDX ] == index ] = True


        # all_real_refpos_measurement = all_refpos_measurement[ real_refpos_measurement_flag ]

        # -----------------------------------------------------
        # -- Determining measurement (program) variations   ---
        # -----------------------------------------------------
        # for file_program_def in self.file_absmeas_window_data:
        #     file_cycle_records = slice(file_program_def['start_record_index'], file_program_def['end_record_index']+1)
        #     mag_records = raw_data[file_cycle_records]
        #     mag_records = mag_records[  mag_records[COMPNAME__POS_ID] !=  ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID ]
        #     # assert np.all( mag_records[COMPNAME__FILE_ABMEAS_WIN_INDEX] == mag_records[0][COMPNAME__FILE_ABMEAS_WIN_INDEX] )
            
        #     abmeas_window_data = self.file_absmeas_window_data [ mag_records[0][COMPNAME__FILE_ABMEAS_WIN_INDEX] ]  
        #     maintype_list_index = abmeas_window_data[COMPNAME__PROGRAM_MAINTYPE_INDEX]
        #     variation_str =  self.program_maintype_str_list[  maintype_list_index  ] # .decode()
        #     variation_str += "/{}".format( len(mag_records) )
        #     if not variation_str in self.program_variation_str_list:
        #         self.program_variation_str_list += [variation_str]
        #         self.pos_id_list += [ list()  ]  # We maintain a pos-id list for each program variation type
        #         for pos_id in mag_records[COMPNAME__POS_ID]:
        #             self.pos_id_list[-1].append(pos_id)

        #     variation_idx = self.program_variation_str_list.index(variation_str)
        #     file_program_def['variation_idx'] =  variation_idx

        #     # ToDo: visszatenni 2025.09.2x után:
        #     # for idx,pos_id in enumerate(mag_records[COMPNAME__POS_ID]):
        #     #    if not  pos_id in self.pos_id_list[variation_idx]:
        #     #        assert False, "Unknown pos_id: {} at time {}".format(pos_id, mag_records[idx][COMPNAME__TIMESTAMP])


        # ---------------------------------------------
        # --   End of loading and pre-processing ------
        # ---------------------------------------------

        # self.number_of_absmeas_processing_windows = len(self.file_absmeas_window_data)

        print()
        nof_files = len(self.input_file_data)
        refmark_readings_count = "n.a" if refmark_readings_data_from_DIM_dat_files is None else len(refmark_readings_data_from_DIM_dat_files)
        print(f"{nof_files} DIM data file{'s' if nof_files > 1 else ''} loaded")
        print(f"\tStored records: mag:{len(self.M_readings)} reference-marks: {refmark_readings_count}")
        print(f"\tNot-used records: {self.nr_of_M_records_outside_of_file_absmeas_windows}")
        print(f"\tLast not completed cycle's records: {self.nr_of_M_records_outside_of_file_absmeas_windows}")
        print(f"\tfile-abmeas-windows: {len(self.file_absmeas_window_data)}")
        # print("\tempty/remark lines: {}\n", nan_lines)
        print("\tMeasurement main type{}: {}".format( mu.plural_suffix(self.program_maintype_str_list), mu.list_to_string(self.program_maintype_str_list)))
        #print("Measurement sub type(s) : "+ mu.list_to_string(self.program_variation_str_list))
        #print("Measurement parameters-1: "+ mu.list_to_string(self.meas_type_param1_str_list))
        #print("Measurement parameters-2: "+ mu.list_to_string(self.meas_type_param2_str_list))
        #print("Azimuth mark{}:   {}".format(  "s" if len(self.azimuth_mark_list)>1 else "",  mu.list_to_string(self.azimuth_mark_list)))
        print("")
        

        self.refmark_readings = []
        if refmark_readings_data_from_DIM_dat_files is not None:
            self.refmark_readings += [ refmark_readings_data_from_DIM_dat_files ]

        # ----------------------------------------------------------
        # ---  Reading azimuth reading data from other files   -----
        # ----------------------------------------------------------
        for calc_var in calc_vars:
            for refmark_data_source in calc_var.refmark_measured_data_sources:
                new_azi_marks_from_files = []                
                match refmark_data_source.source_type:
                    case dim_azimuth.azimark_meas_data_source_type_enum.dim_dat_files: 
                        assert refmark_readings_data_from_DIM_dat_files is not None, "DIM data files do not have reference mark readings data"
                        refmark_data_source.refmark_readings_array_index = 0
                    case dim_azimuth.azimark_meas_data_source_type_enum.spreadsheet_table: 
                        spreadsheet_table_refmark_data, new_azi_marks_from_files = dim_azimuth.read_refmark_readings_from_spreadsheet(refmark_data_source) 
                        self.refmark_readings += [spreadsheet_table_refmark_data]
                        refmark_data_source.refmark_readings_array_index = len(self.refmark_readings)-1
                    case dim_azimuth.azimark_meas_data_source_type_enum.separate_dat_files:   
                        assert refmark_data_source.azimark_name is not None  # 'azimark_name' should come from YAML config file
                        if not refmark_data_source.azimark_name in self.azimuth_mark_list:
                            self.azimuth_mark_list += [refmark_data_source.azimark_name]
                        azimark_index = self.azimuth_mark_list.index(refmark_data_source.azimark_name)
                        pillar_index = 0 # Only one pillar is allowed
                        file_azimuth_data = dim_azimuth.read_refmark_readings_from_dat_files(refmark_data_source, pillar_index, azimark_index, self.M_readings[COMPNAME__TIMESTAMP], meas_pos, self.instr, COMPNAME__PILLAR_INDEX)
                        self.refmark_readings += [file_azimuth_data]
                        # In case of ADS, 'file_azimuth_data' contains raw tick values.
                        # It will be eccentricity corrected later
                        refmark_data_source.refmark_readings_array_index = len(self.refmark_readings)-1
                    case _:
                        assert False, "Unsupported reference data source"
                for new_azimark in new_azi_marks_from_files:
                    if not new_azimark in self.azimuth_mark_list: self.azimuth_mark_list += new_azimark
                
        

        # if ( self.number_of_absmeas_processing_windows > 1):
        #     # Checking if points in cycles are changing or not (last cycle is checked later)
        #     if ( self.number_of_absmeas_processing_windows > 2):
        #         for i in range(0,self.number_of_absmeas_processing_windows-2):
        #             if (points_per_cycle[i] != points_per_cycle[i+1]):
        #                 print("Note: number of points in cycles {} and {} are different!".format ( i, i+1 ) )

        #     # Checking last cycle (if last cycle in incomplete, we are cutting it:
        #     if ( points_per_cycle[self.number_of_absmeas_processing_windows-1] < points_per_cycle[self.number_of_absmeas_processing_windows-2] ):
        #         print("Last cycle is incomplete (has {} points) - cutting it !". format( points_per_cycle[self.number_of_absmeas_processing_windows-1]) )
        #         self.number_of_absmeas_processing_windows = self.number_of_absmeas_processing_windows -1
        #     elif ( points_per_cycle[self.number_of_absmeas_processing_windows-1] > points_per_cycle[self.number_of_absmeas_processing_windows-2] ):
        #         print("Note: number of points in the last cycle is bigger!" )



        # Deprecated feature to add mag sensors data from external files:
        # for i in range(1,self.nr_of_mag_sensors):
        #     # Adding new column to file_records[] array
        #     self.compindeces_mag =  self.compindeces_mag + [ len(self.file_records[0])  ]
        #     self.file_records = np.append(self.file_records, np.empty([n,1]), axis = 1)

        #     # Reading 2nd mag data from maglin/ObsDaq data file
        #     TIME_DIFF_SEC = -1
        #     # orig_sensor_M = np.concatenate( (orig_sensor_M, mag) )
        #     for i,itimelabel in enumerate(self.file_records[:,  COMPNAME__TIMESTAMP]):
        #         self.file_records[i,self.compindeces_mag[-1]] = dat_file_reader.second_mag_sensor.get_data(itimelabel + TIME_DIFF_SEC / 86400.0)[0]
        #     # second_sensor_M = np.concatenate( (second_sensor_M, mag) )
        #     # second_sensor_M = np.copy(self.file_records[:,self.compindex_Mag])
        #     print("================ Comparing sensors ==================")
        #     offset,scale = du.LSQ(self.file_records[:,  self.compindeces_mag[0]], self.file_records[:,  self.compindeces_mag[-1]])
        #     print("Without TempCo: mag_2ndsensor = {:.3f}  + {:.5f} * mag".format(offset,scale))
        #     M_diff = (self.file_records[:,  self.compindeces_mag[0]] * scale + offset) - self.file_records[:,  self.compindeces_mag[-1]]
        #     plt.plot(self.file_records[:,  COMPNAME__TIMESTAMP],M_diff)
        #     plt.show()
        #     # plt.plot(hc_tick,M_diff)
        #     # plt.show()
        return True


    # Returns a list of file record indeces between 'file_record_index_start' and 'file_record_index_end', filtered by 'filter_defs'
    # It is used to create "processing cycles" and "file cycles"
    # def get_filtered_file_record_indeces(self, file_record_index_start, file_record_index_end, filter_defs):
    #     assert file_record_index_end>=0
    #     assert file_record_index_end<raw_data.shape[0]
    #     file_record_indeces = list(range(file_record_index_start, file_record_index_end+1))

    #     # Applying cycle point filters, if there is any
    #     for filter_def in filter_defs:  # Filter defs may contain any, one or more filter definitions
    #         assert(len(filter_def)==2)  # Each filter definition has two parameters: first is the filter identification name, second is filter params
    #         filter_name = filter_def[0]
    #         filter_params = filter_def[1]
    #         if filter_name == "CYCLE_POINT_FILTER__EXCLUDE_POS_IDS":
    #             exclude_pos_IDs = filter_params
    #             for i in range( len(file_record_indeces)-1, -1, -1 ):
    #                 if raw_data[COMPNAME__POS_ID][ file_record_indeces[i]  ] in exclude_pos_IDs:
    #                     # print( len(file_record_indeces), i )
    #                     file_record_indeces.pop(i)
    #         elif filter_name == "CYCLE_POINT_FILTER__INCLUDE_POS_IDS":
    #             include_pos_IDs = filter_params
    #             for i in range( len(file_record_indeces)-1, -1, -1 ):
    #                 if not raw_data[COMPNAME__POS_ID][ file_record_indeces[i]  ] in include_pos_IDs:
    #                     # print( len(file_record_indeces), i )
    #                     file_record_indeces.pop(i)

    #         else:
    #             assert False,"Unsupported point filter: {}".format(filter_name)

    #     return file_record_indeces
    
    
    def add_warning(self, warnmsg):
        if warnmsg in self.warn_logs:
            self.warn_log_count[warnmsg] += 1
        else:
            self.warn_logs += [ warnmsg ]
            self.warn_log_count[warnmsg] = 1
            
    def get_warnings(self):
        result = []
        for warnmsg in self.warn_logs:
            if self.warn_log_count[warnmsg] == 1:
                result += [ warnmsg]
            else:
                result += [ "[{}x] {}".format(self.warn_log_count[warnmsg],  warnmsg)]
        return result
                
    
    def read_in_ads_raw_dat_files(self,  base_path, filename_masks, date_selection_filters):
        
        user_set_pillar_setting_used = False

        # This sub-function extract some important ADS parameters from the ADS dat file's remark header lines.
        # ADS raw data files have some information in the header. This sub-function extract the important fields:
        # - vc0 (ADS vc axis offset)
        # - pillar name
        def parse_remark_line(line, file_index, filename):
            words = re.split(r'\s+', line)
            nonlocal new_input_file_data, user_set_pillar_setting_used
            for word in words:
                if word[:4] == "vc0=":
                    vc0_str = word[4:]
                    new_vc0=np.radians(float(vc0_str))
                    if (self.vc0_from_file is not None) and (new_vc0 != self.vc0_from_file):
                        self.add_warning("New vc0 value read from {} file is changed from previous file's {:.3} deg value to {:} deg".format(filename, np.degrees(self.vc0_from_file), vc0_str))
                    self.vc0_from_file = new_vc0
                    new_input_file_data["vc0"] = self.vc0_from_file
                elif word[:7] == "pillar=":
                    pillar_name_from_header = word[7:]
                    
                    if "DIM-data-files-overrides" in self.calc_params_common.remaining_user_params.keys():
                        assert type(self.calc_params_common.remaining_user_params["DIM-data-files-overrides"]) == dict
                        if "set-pillar" in self.calc_params_common.remaining_user_params["DIM-data-files-overrides"].keys():
                            pillar_name_from_header = self.calc_params_common.remaining_user_params["DIM-data-files-overrides"]["set-pillar"]
                            user_set_pillar_setting_used = True
                    
                    assert pillar_name_from_header in self.meas_pos.pillars, f'Pillar "{pillar_name_from_header}" from ADS dat file is not within actual pillars list ("{self.meas_pos.pillars}")'
                    pillar_index = self.meas_pos.pillars.index(pillar_name_from_header)
                    if not "pillar_index" in new_input_file_data.keys():
                        new_input_file_data["pillar_index"] = pillar_index
                    else:
                        assert (new_input_file_data["pillar_index"] == pillar_index), "Pillar name changed in file header data"
                        # A pillar is already read from the header (can happen, ir more files are read)
                        #if self.pillar_name_from_files[file_index].decode() != pillar_name_from_header:
                        #    self.pillar_name_from_files[file_index] = "more_pillars"
                        
        def remove_deselected_records(selection_filter_name, delete_records):
            nonlocal raw_data
            if np.any(delete_records):
                raw_data = np.delete(raw_data, delete_records, axis = 0)
                if selection_filter_name is not None:
                    if selection_filter_name in self.applied_remove_filters:
                       self.applied_remove_filters[selection_filter_name] += np.count_nonzero(delete_records)
                    else:
                       self.applied_remove_filters[selection_filter_name] = np.count_nonzero(delete_records)
                return True
            else:
                return False
                                
                        
        # def how_may_times_this_pos_have_been_measured_before(rec_index, delta_tick):
        #     n = 0
        #     for i in range(0,rec_index):
        #        if not np.isnan(raw_data[i][COMPNAME__HC_TICK]):
        #            if  abs( raw_data[i][COMPNAME__HC_TICK]-raw_data[rec_index][COMPNAME__HC_TICK]   ) <= delta_tick \
        #                and \
        #                abs( raw_data[i][COMPNAME__VC_TICK]-raw_data[rec_index][COMPNAME__VC_TICK]   ) <= delta_tick :
        #                    n += 1
        #     return n
                        
        filenames, used_filenamemasks = mu.filenamemask_to_filenames(base_path, filename_masks, date_selection_filters)
        assert len(filenames) > 0, "No ADS raw data file found with date mask(s) \"{}\" (file mask:\"{}\") in directory \"{}\"".format(filename_masks, used_filenamemasks[0], base_path)
        
        self.record_separator=" "  # Field separator of lines        
        self.input_file_data = []
        
        nr_of_files = len(filenames)
        print("Loading ADS raw data file{} from directory \"{}\"...".format( mu.plural_suffix(nr_of_files), base_path ))

        for file_index, filename in enumerate(filenames):
            print(f"  {file_index+1:3d}/{nr_of_files}  {filename} : ", end="")
            
            
            full_filename = os.path.join(base_path, filename)
            assert os.path.isfile(full_filename), "DAT file {:} does not exist".format(filename)

            # Reading in first lines: 1) header line + 2) comment lines.
            # 1) Reading  header.
            text_file = open(full_filename,"r")
            # self.file_records_str = text_file.readlines()
            header_line = text_file.readline()

            header_line_fields = self.split_record_str(header_line)

            # 1) Reaing command lines
            new_input_file_data = {}
            nr_of_comment_header_lines = 0
            while True:
                new_line = text_file.readline()
                if new_line[0]=='#':
                    nr_of_comment_header_lines += 1                    
                    parse_remark_line( new_line, file_index, filename  )
                    # self.comment_lines_in_header = self.comment_lines_in_header + [ new_line ]
                else:
                    break

            text_file.close()

            # Parsing first header line, containing components
            dtypes=[]
            names=[]
            colindeces=[]
            number_of_string_comps = 0
            for col_index,comp_name in enumerate(header_line_fields):
                if comp_name in ['date','time','doy','']: # ,'PPM_M[nT]'
                    # Skipping this columns
                    continue
                comp_type = float
                if comp_name == 'meas_part_id':
                    comp_type = "|S{}".format(PROGRAMSTEP_TYPE_FULL_STR_MAX_LEN)
                    number_of_string_comps += 1
                    comp_name = COMPNAME__PROGRAMSTEP_TYPE_FULL_STR  # Renaming this component to "_str", as we will have meas_part IDs!
                elif comp_name == 'meas_flags':
                    comp_type = "|S{}".format(MEAS_FLAGS_FULL_STR_LEN)
                    number_of_string_comps += 1
                elif comp_name in [ "measurement_count", "cycle_id", "pos_id"  ]:
                    comp_type = np.int32
                elif "Â°C" in comp_name:   # T[°C]" is read as "TÂ°C"
                    comp_name = comp_name.replace("Â°C","°C")

                dtypes += [ comp_type ]
                names += [ comp_name ]
                colindeces += [col_index]

            if (file_index == 0):
                prev_file_names = names
            else:
                assert names == prev_file_names, "Can not process more files with different contents (column setup)!"


            # -------------------------------------------------------------------------
            # ----   Reading in the data from one file --------------------------------
            # -------------------------------------------------------------------------
            try:
                new_file_ads_raw_data = np.genfromtxt(full_filename, skip_header= 1 + nr_of_comment_header_lines, names=names, usecols=colindeces, dtype=dtypes, comments='#' )
            except  Exception as e : # (IOError, ValueError, OSError)
                assert False, "Error loading file \"{}\": {}".format(filename, repr(e))
                
            orig_nr_of_records = len(new_file_ads_raw_data)
            if "FROM-DATE" in date_selection_filters.keys():
                new_file_ads_raw_data = new_file_ads_raw_data[ new_file_ads_raw_data[COMPNAME__TIMESTAMP] >= date_selection_filters["FROM-DATE"][1] ]
            if "TO-DATE" in date_selection_filters.keys():
                new_file_ads_raw_data = new_file_ads_raw_data[ new_file_ads_raw_data[COMPNAME__TIMESTAMP] <= date_selection_filters["TO-DATE"][1] ]
            nr_of_records_after_date_selection = len(new_file_ads_raw_data)
            
            if nr_of_records_after_date_selection == orig_nr_of_records:
                print("{} records".format(orig_nr_of_records))
            elif nr_of_records_after_date_selection == 0:
                print("file not used due to date selection")
                continue
            else:
                print("{} records, after date selection: {} records".format(orig_nr_of_records, nr_of_records_after_date_selection))
                
                
            new_file_ads_raw_data =  rcf.append_fields( new_file_ads_raw_data, COMPNAME__INPUT_FILE_INDEX, [0], np.uint32, usemask=False )
            new_file_ads_raw_data[COMPNAME__INPUT_FILE_INDEX][:] = file_index
            
            if file_index==0:
                raw_data = new_file_ads_raw_data
            else:
                assert new_file_ads_raw_data.dtype.names == raw_data.dtype.names, "ADS input files fields are different! This situation is not supported" 
                raw_data = np.concatenate( (raw_data, new_file_ads_raw_data) )
                del new_file_ads_raw_data
            
            new_input_file_data["file_type"] = "ADS-raw-dat-files"
            new_input_file_data["full_file_name"] = full_filename
            new_input_file_data["file-name"] = filename
            new_input_file_data["start-date"] = raw_data[0][COMPNAME__TIMESTAMP]
            new_input_file_data["end-date"] = raw_data[-1][COMPNAME__TIMESTAMP]
            # "instr_name":instr_name, "pillar_index": pillars.index(pillar), "date_from_filename":None }  ]                                                        
            self.input_file_data +=  [ new_input_file_data ] 

        # ---- All files are read  ---
        
        if user_set_pillar_setting_used:
            del self.calc_params_common.remaining_user_params["DIM-data-files-overrides"]["set-pillar"]
        
        # -------------------------------------------------------------------------
        # ----   Some conversions for old files -----------------------------------
        # -------------------------------------------------------------------------
        
        assert len(self.input_file_data) > 0, "No DIM data file was read"

        global COMPNAME__MAG
        if not COMPNAME__MAG in raw_data.dtype.names:
            assert "M1nT" in raw_data.dtype.names  # Files in 2023 has "M1[nT]" component name for mag
            raw_data = rcf.rename_fields(raw_data, { "M1nT": COMPNAME__MAG} ) 

        # Some ADS system does not have room temperature sensor, but has '"T_room[]" column in data files, filled with invalid values.
        # Deleting this column:
        if ( COMPNAME__T_ROOM in raw_data.dtype.names ):
            if np.all( raw_data[COMPNAME__T_ROOM]<-99  ) :
                raw_data =  rcf.drop_fields( raw_data, COMPNAME__T_ROOM, usemask=False)
            else:
               # In periods when temp. sensor is not attached, room temperature is -99
               raw_data[COMPNAME__T_ROOM] [  raw_data[COMPNAME__T_ROOM]<-90 ]  = np.nan
            
        # In some old ADS data files, angle data can be empty
        remove_deselected_records(None,  np.isnan(raw_data[COMPNAME__HC_TICK]) | np.isnan(raw_data[COMPNAME__VC_TICK]) )
        
        # Removing the pseudo 0 nT records:
        # if ('FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS' in globals()) and FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS:
        #if remove_deselected_records('FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS', (raw_data[COMPNAME__POS_ID] < 0)  &  (raw_data[COMPNAME__POS_ID] > ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID) ):
        #        nof_removed_NM_calculated_zeropos_records = self.applied_remove_filters['FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS']
        pseudo_0nT_records = (raw_data[COMPNAME__POS_ID] < 0)  &  (raw_data[COMPNAME__POS_ID] > ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID)
        remove_deselected_records(None, pseudo_0nT_records)

        if "FROM-DATE" in date_selection_filters.keys():
            remove_deselected_records(None, raw_data[COMPNAME__TIMESTAMP] < date_selection_filters["FROM-DATE"][1])  

        if "TO-DATE" in date_selection_filters.keys():
            remove_deselected_records(None, raw_data[COMPNAME__TIMESTAMP] > date_selection_filters["TO-DATE"][1])  
            
        # raw_data = extend_structured_array(raw_data, rawdata_base_dtype) 
            
       
        # Extending 'raw_data' structured array with some derived components
        fill_val_uint = np.zeros( raw_data.shape[0], dtype = np.uint32  )
        fill_val_char = np.empty( raw_data.shape[0], dtype = 'U1'  )
        fill_val_char[:] = '?'
        # fill_val_uint[:]=999999
        # fill_val_int_minus1 = np.zeros( raw_data.shape[0], dtype = np.int32  )
        # fill_val_int_minus1[:]=-1
        # fill_val_nan = np.empty(raw_data.shape[0] )
        # fill_val_nan[:] = np.nan
        

        added_comp_names = [ COMPNAME__IS_M_RECORD  ]
        added_comp_dypes = [ "?" ]   # Boolean field
        added_comp_fillvals = [ [0] ]
        
        if np.all(np.isnan(raw_data[COMPNAME__F])):
            raw_data = rcf.drop_fields(raw_data, "COMPNAME__F", usemask=False)

        #if not COMPNAME__PROGRAM_MAINTYPE_INDEX in raw_data.dtype.names:
        #    added_comp_names    += [ COMPNAME__PROGRAM_MAINTYPE_INDEX ] 
        #    added_comp_dypes    += [ int  ]        
        #    added_comp_fillvals += [  fill_val_uint ]  # This will be converted to string "999999"
        
        # added_comp_names    += [ COMPNAME__PROGRAM_DEF_IDX__NOTUSED ]
        # added_comp_dypes    += [ np.int32  ]        
        # added_comp_fillvals += [ fill_val_int_minus1 ]  # This will be converted to string "999999"

        # Adding more "uint32" type fields at once:        
        uint32_fields = [ COMPNAME__ADS_FILE__RECORD_WITHIN_PROGRAM, COMPNAME__FILE_ABMEAS_WIN_INDEX, COMPNAME__PROGRAM_MAINTYPE_IDX, COMPNAME__PROGRAM_PARAM1_IDX, COMPNAME__PROGRAM_PARAM2_IDX, COMPNAME__AZIMARK_INDEX ] # , COMPNAME__RAW_DATA_INDEX, COMPNAME__EXTRA_FLAG
        for field in uint32_fields:
            added_comp_names    += [ field ]
            added_comp_dypes    += [ np.uint32 ]  
            added_comp_fillvals += [ fill_val_uint ] 
            
        char_fields = [  COMPNAME__TELESCOPE_CARDINAL_DIRECTION, COMPNAME__TELESCOPE_FACE_POSITION ] 
        for field in char_fields:
            added_comp_names    += [ field ]
            added_comp_dypes    += [ COMPTYPE__TELESCOPE_CARDINAL_DIRECTION  if field == COMPNAME__TELESCOPE_CARDINAL_DIRECTION else COMPTYPE__TELESCOPE_FACE_POSITION ]  
            added_comp_fillvals += [ fill_val_char ] 
        
        # for i in range(len(added_comp_names)):
        #   print(i," ",added_comp_names[i]," ", added_comp_dypes[i],"  ",added_comp_fillvals[i][0])
        #   raw_data =  rcf.append_fields( raw_data, names=[added_comp_names[i]], data=[added_comp_fillvals[i]], dtypes= [added_comp_dypes[i]], usemask=False )
        #   print("done: ",raw_data[added_comp_names[i]][0])
        raw_data =  rcf.append_fields( raw_data, names = added_comp_names, data = added_comp_fillvals, dtypes = added_comp_dypes, usemask=False )
        
        raw_data[COMPNAME__IS_M_RECORD] = raw_data[COMPNAME__POS_ID] != ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID
        
        range_out_M_records_selection = raw_data[COMPNAME__IS_M_RECORD]  & ( ( raw_data[COMPNAME__MAG] < self.instr.mag_measuring_range_nT[0] ) | ( raw_data[COMPNAME__MAG] > self.instr.mag_measuring_range_nT[1] ) )
        
        n_out_of_range = np.count_nonzero(range_out_M_records_selection)
        if n_out_of_range > 0:
            warnmsg = "{} magnetometer reading record{} out of the magnetometer measuring range of {}..{} nT: {} excluded from calculations".format(n_out_of_range, mu.plural_suffix(n_out_of_range," is", "s are"),  self.instr.mag_measuring_range_nT[0], self.instr.mag_measuring_range_nT[1], mu.plural_suffix(n_out_of_range,"it is", "they are"))
            self.add_warning(warnmsg)
            raw_data =  raw_data[~range_out_M_records_selection]
            

        # ----------------------------------------------------------------------------------
        # --- Applying data include / exlude filters that can be checked within a loop ----
        # --  Converting method name parts into IDs --
        # ----------------------------------------------------------------------------------
        self.pos_id_list = [   ]    # We collect pos id separatlely for each measurement type. So this will be a '2D' list.
        self.azimuth_mark_list = []   # Collecting measured reference marks' names. If azimuth readings are in a separate files, then this list will be filled later, when readings that files
        self.program_variation_str_list  = []
        self.meas_type_param1_str_list = []
        self.meas_type_param2_str_list = []
        
        prev_absmeas_win_index = -9999
        prev_measurement_count = -1
        prev_maintype_str = None
        # ads_raw_data_index_after_filtering = 0

        convert_old_type_names =  raw_data[COMPNAME__TIMESTAMP][0] < mu.datetimestr_2_exceltime("2024.08.01 00:00", "%Y.%m.%d %H:%M")
        record_within_meas_window = 0
        nr_of_M_data_within_cycle = 0
        noncomplete_meas_record_indeces = []
        self.program_maintype_str_list = []

        for file_records_index,raw_file_rec in enumerate(raw_data):

            program_maintype_str, meas_type_param1_str, meas_type_param2_str = separate_full_meas_type_str(raw_file_rec[COMPNAME__PROGRAMSTEP_TYPE_FULL_STR])

            if convert_old_type_names:
                if program_maintype_str in [ b'nullmeth-slow', b'absmeas' ]:
                    program_maintype_str = MEAS_TYPE_STR__NULLMETH_DEF
                # elif program_maintype_str[:8] == b"sinus_VC":
                #    program_maintype_str = b"SV"

            assert program_maintype_str in [MEAS_TYPE_STR__NULLMETH_1, MEAS_TYPE_STR__NULLMETH_2, MEAS_TYPE_STR__MIRA_READING, b"SV", b'grid(4x4)', b'DIAGONAL-FW', b'DIAGONAL-BW', b'SINGLEAXIS', b"TEST", b"NM-D" ], "Unsupported meas. type: \"{}\"".format(program_maintype_str)

            # Doing some conversions (to be compatible with old files):
            if (program_maintype_str in [ MEAS_TYPE_STR__NULLMETH_1, MEAS_TYPE_STR__NULLMETH_2]) and ( len(meas_type_param1_str) in [1,3] )  :
                # In files before 2026-01: pseudo measruement point of the ADS scanning null method: the determined zero nT position.
                # - "ED"/"WU"/etc means east-down, west-up, etc. calculated points
                # - In files before 2023.02.23, due to a bug, only the first letter of the two letter notation was saved
                # - In later "ED0", "WU0" etc notation is used
                # We already removed these records, see 'pseudo_0nT_records'
                assert False 
                # assert(raw_file_rec[COMPNAME__POS_ID] < 0)  # These points are marked by negativ pos-id
                # method_str_extra_info="nullmeth-calculated-zeropos"
                #if 'FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS' in globals():
                #    if apply_record_deselection_filter( FILTER_EXCLUDE_NULLMETH_CALCULATED_ZEROPOS, "", file_records_index): # Do not logging this type of filter!!!
                #        nof_removed_NM_calculated_zeropos_records += 1
                #        continue
                #    # Using calcualted zeropos record!!!!
                if len(meas_type_param1_str)<2:  meas_type_param1_str = meas_type_param1_str + b"0"
                                          

            filtered_out = False

            # Basic separation of record to cycles (this is only a temporarly cycle separation method for record filter)
            new_meas_program = (
                (raw_file_rec[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] != prev_absmeas_win_index) 
                or (raw_file_rec[COMPNAME__MEAS_COUNT_SINCESTART] != prev_measurement_count + 1)    # Jump of measurement count
                or (program_maintype_str != prev_maintype_str)
            )
            
            # Old method of record separation:
                # new_pos_id = raw_data[COMPNAME__POS_ID][i]
                # assert not  np.isnan(new_pos_id)
                # if nof_records_in_current_cycle == 0:
                #     # A new cycle cycle just starts
                #     assert within_cycle == False
                #     cycle_end = False
                #     cycle_start = True
                # elif (i>0) and absmeas_window_determination_conditions.at_measrement_count_backjump and  (raw_data[COMPNAME__MEAS_COUNT_SINCESTART][i]<=raw_data[COMPNAME__MEAS_COUNT_SINCESTART][i-1]):
                #     # For ADS files only: 'measurement_count' jumps back only in case of program restart!
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_file_meas_window_index_change) and  (prev_absmeas_win_index != raw_data[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX][i]):
                #     # Most important separation: separating when a new measurement starts.
                #     # Basic separation method for Excel files.
                #     # For ADS files, this may not be enough, as index restarts from "0" each time the ADS-control program restarts, wo
                #     #  there can be neaighboring measurements windows with "0" id  
                #     # New cycle id in file: new cycle begins
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_pos_id_jump_back) and (nof_records_in_current_cycle>1) and ( round(new_pos_id) >= 0) and ( round(new_pos_id)  < round(prev_mag_pos_id) ):
                #     # Pos ID resetted to zero
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_pos_id_cont_change) and (nof_records_in_current_cycle>0) and (  abs(new_pos_id - prev_pos_id) != 1) :
                #     # Pos ID jump
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_pos_id_cont_change) and  (nof_records_in_current_cycle>1) and ( (new_pos_id - prev_pos_id) != (prev_pos_id - raw_data[i-2][self.compindex_pos_id]) ) :
                #     # Pos direction change (increasing -> decreasing pos_id, or decreasing -> increasing pos_id)
                #     assert prev_pos_id == raw_data[i-1][self.compindex_pos_id]
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_new_measurement_on_the_same_pos) and  (how_may_times_this_pos_have_been_measured_before(i, absmeas_window_determination_conditions.same_pos_delta_tick_limit ) != same_pos_measurement_count ) :
                #     assert ( how_may_times_this_pos_have_been_measured_before(i) == same_pos_measurement_count  + 1 )
                #     same_pos_measurement_count  += 1
                #     cycle_end = True
                #     cycle_start = True
                # elif (absmeas_window_determination_conditions.at_big_time_jump) and ( i>0  and (raw_data[i][COMPNAME__TIMESTAMP] > raw_data[i-1][COMPNAME__TIMESTAMP] + absmeas_window_determination_conditions.time_jump_limit_minutes / 1440.0 )  ):
                #     cycle_end = True
                #     cycle_start = True
                # else:
                #     cycle_end = False
                #     cycle_start = False
            
            
            if new_meas_program:
                prev_absmeas_win_index = raw_file_rec[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX]
                if (nr_of_M_data_within_cycle < MIN_MAG_MEAS_POINTS_PER_ABMEAS_PROCESSING_WINDOW) and (record_within_meas_window > 0):
                    noncomplete_meas_record_indeces += [ list(range( file_records_index - record_within_meas_window, file_records_index)) ]
                record_within_meas_window = 0
                nr_of_M_data_within_cycle = 0
                                     
            raw_file_rec[COMPNAME__ADS_FILE__RECORD_WITHIN_PROGRAM] = record_within_meas_window
            record_within_meas_window += 1
            prev_measurement_count = raw_file_rec[COMPNAME__MEAS_COUNT_SINCESTART]
            prev_maintype_str = program_maintype_str
            # if ('FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN' in globals()) and FILTER_EXCLUDE_ANY_MAG_SENSOR_IS_NAN:
            #     assert False
            #     if not raw_file_rec[self.compindex_pos_id] in REFPOS_MEASUREMENT_RECORDS_POS_IDS: # Record is not a reference position measurement (where Mag should be NAN)
            #         for imag in range(min_nr_of_mag_sensors):
            #             mag_nan_filtered_out = apply_record_deselection_filter(np.isnan( raw_file_rec[self.compindeces_mag[imag]],file_records_index), "M{:}-nan record".format(imag+1))
            #             if mag_nan_filtered_out: continue
            # # We can apply text based filter at this stage of file process
            # if 'FILTER_METHOD_STR_INCLUDE' in globals():
            #     assert False, "Old function needs to be revised"
            #     mask_fits = False
            #     for filter_mask in FILTER_METHOD_STR_INCLUDE:
            #         # Filter can be defined by mask (starting or ending asterix), or by a string
            #         if mu.compare_with_string_mask( filter_mask, program_maintype_str ):
            #             mask_fits = True
            #             break
            #     if apply_record_deselection_filter(not mask_fits, "FILTER_METHOD_STR_INCLUDE", file_records_index): filtered_out = True
            # if 'FILTER_METHOD_STR_EXCLUDE' in globals():
            #     assert False, "Old function needs to be revised"
            #     mask_fits = False
            #     for filter_mask in FILTER_METHOD_STR_EXCLUDE:
            #         # Filter can be defined by mask (starting or ending asterix), or by a string
            #         if mu.compare_with_string_mask( filter_mask, program_maintype_str ):
            #             mask_fits = True
            #             break
            #     if apply_record_deselection_filter(mask_fits, "FILTER_METHOD_STR_EXCLUDE", file_records_index): filtered_out = True
            # if 'FILTER_METHOD_EXTRA_INFO_STR_EXCLUSIVE' in globals():
            #     assert False, "Old function needs to be revised"
            #     mask_fits = False
            #     for filter_mask in FILTER_METHOD_EXTRA_INFO_STR_EXCLUSIVE:
            #         # Filter can be defined by mask (starting or ending asterix), or by a string
            #         if mu.compare_with_string_mask( filter_mask, method_str_extra_info ):
            #             mask_fits = True
            #             break
            #     if apply_record_deselection_filter(mask_fits, "FILTER_METHOD_EXTRA_INFO_STR_EXCLUSIVE", file_records_index): filtered_out = True
            # if 'FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE' in globals():
            #     assert False
            #     _vc_tick = raw_data[i][COMPNAME__VC_TICK]
            #     apply_record_deselection_filter(i, (_vc_tick<FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE) or ( (_vc_tick > dim_instr.ENCODER_DIVISION / 2 - FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE) and (_vc_tick < dim_instr.ENCODER_DIVISION / 2 + FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE)  )   , "FILTER_VC_VERTICAL_POSITIONS_TICK_RANGE")
            # if 'FILTER_VC_DEG_MIN' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[i][self.compindex_file_vc_deg]<FILTER_VC_DEG_MIN, "FILTER_VC_DEG_MIN")
            # if 'FILTER_VC_DEG_MAX' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[i][self.compindex_file_vc_deg]>FILTER_VC_DEG_MAX, "FILTER_VC_DEG_MAX")
            # if 'FILTER_HC_DEG_MIN' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[i][self.compindex_file_hc_deg]<FILTER_HC_DEG_MIN, "FILTER_HC_DEG_MIN")
            # if 'FILTER_HC_DEG_MAX' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[i][self.compindex_file_hc_deg]>FILTER_HC_DEG_MAX, "FILTER_HC_DEG_MAX")
            # if 'FILTER_CYCLE_ID_MOD100_MIN' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX][i] % 100<FILTER_CYCLE_ID_MOD100_MIN, "FILTER_CYCLE_ID_MOD100_MIN")
            # if 'FILTER_CYCLE_ID_EXCLUDE_LIST' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, raw_data[i][COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] in FILTER_CYCLE_ID_EXCLUDE_LIST, "FILTER_CYCLE_ID_EXCLUDE_LIST")
            # if 'FILTER_CYCLE_ID_INCLUDE_LIST' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, not raw_data[i][COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] in FILTER_CYCLE_ID_INCLUDE_LIST, "FILTER_CYCLE_ID_INCLUDE_LIST")
            # if 'FILTER_CYCLE_ID_MOD2' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, (raw_data[i][COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] % 2) == FILTER_CYCLE_ID_MOD2, "FILTER_CYCLE_ID_MOD2")
            # if 'FILTER_CYCLE_ID_REMOVE_MODULO' in globals():
            #     assert False
            #     apply_record_deselection_filter(i, (raw_data[i][COMPNAME__ADS_FILE_ABMEAS_WIN_INDEX] % FILTER_CYCLE_ID_REMOVE_MODULO[0]) == FILTER_CYCLE_ID_REMOVE_MODULO[1], "FILTER_CYCLE_ID_REMOVE_MODULO")
            # if 'FILTER_SKIP_FIRST_SAMPLES_ON_EACH_CYCLE' in globals():
            #     assert False
            #     # egyes hc és vc forgatásokkor a ciklusok első (induló) mintái kiugróak, val.szin megchanikus kotyogás miatt
            #     # Ezeket itt probaljuk kivenni:
            #     apply_record_deselection_filter(i,record_within_meas_window<=FILTER_SKIP_FIRST_SAMPLES_ON_EACH_CYCLE,"FILTER_SKIP_FIRST_SAMPLES_ON_EACH_CYCLE")
            # #if 'FILTER_EXCLUDE_BEFORE' in globals():
            # #    apply_record_deselection_filter(i, raw_data[i][COMPNAME__TIMESTAMP] < FILTER_EXCLUDE_BEFORE,"FILTER_EXCLUDE_BEFORE")
            # #if 'FILTER_EXCLUDE_AFTER' in globals():
            # #    apply_record_deselection_filter(i, raw_data[i][COMPNAME__TIMESTAMP] > FILTER_EXCLUDE_AFTER,"FILTER_EXCLUDE_AFTER")

            if filtered_out:
                continue

            if not program_maintype_str  in self.program_maintype_str_list:
                self.program_maintype_str_list   += [program_maintype_str]
            if not meas_type_param1_str in self.meas_type_param1_str_list: self.meas_type_param1_str_list += [meas_type_param1_str]
            if not meas_type_param2_str in self.meas_type_param2_str_list: self.meas_type_param2_str_list += [meas_type_param2_str]

            # assert len(program_maintype_str) <= COMPNAME__PROGRAM_MAINTYPE_STR_MAX_LEN, "Field-len:{}  max-len:{}  Field:\"{}\"".format(len(program_maintype_str),COMPNAME__PROGRAM_MAINTYPE_STR_MAX_LEN,program_maintype_str )
            # raw_data[file_records_index][COMPNAME__PROGRAM_MAINTYPE_STR]  =  program_maintype_str

            # raw_data[file_records_index][COMPNAME__RAW_DATA_INDEX] = ads_raw_data_index_after_filtering
            # ads_raw_data_index_after_filtering += 1

            # Setting IDs:
            raw_file_rec[COMPNAME__PROGRAM_MAINTYPE_IDX]  =  self.program_maintype_str_list.index(program_maintype_str)
            raw_file_rec[COMPNAME__PROGRAM_PARAM1_IDX]    =  self.meas_type_param1_str_list.index(meas_type_param1_str)
            raw_file_rec[COMPNAME__PROGRAM_PARAM2_IDX]    =  self.meas_type_param2_str_list.index(meas_type_param2_str)
            
            cardinal_direction_from_file = du.UNKNOWN_CARDINAL_DIRECTION
            face_position_from_file      = du.telescope_face_pos_enum.unknown.value
            
            if (program_maintype_str  == MEAS_TYPE_STR__NULLMETH_2) and (len(meas_type_param1_str)>4) and (meas_type_param1_str[1] == "-") and (meas_type_param1_str[0] in du.TELESCOPE_CARDINAL_DIRECTIONS) :
                # Format from 2027:  "ABSM/E-FL-1/3", from which 'meas_type_param1_str' will be E-FL-1
                cardinal_direction_from_file = meas_type_param1_str[0]
                assert meas_type_param1_str[2] == 'F'
                if meas_type_param1_str[3] == 'L':
                    face_position_from_file = du.telescope_face_pos_enum.FL.value
                elif meas_type_param1_str[3] == 'R':
                    face_position_from_file = du.telescope_face_pos_enum.FR.value
                else:
                    assert False
                    
            raw_file_rec[COMPNAME__TELESCOPE_FACE_POSITION] = face_position_from_file
            raw_file_rec[COMPNAME__TELESCOPE_CARDINAL_DIRECTION] = cardinal_direction_from_file
            
            if raw_file_rec[COMPNAME__IS_M_RECORD]:
                raw_file_rec[COMPNAME__AZIMARK_INDEX] = INV_REFMARK_INDEX
                nr_of_M_data_within_cycle += 1
            else:
                # assert raw_file_rec[COMPNAME__PROGRAM_MAINTYPE_STR] == MEAS_TYPE_STR__MIRA_READING
                azimark_name = meas_type_param1_str
                # if azimark_name== None:
                #    azimark_name = refpos_id_from_method_str
                if not azimark_name in self.azimuth_mark_list: self.azimuth_mark_list += [azimark_name]
                raw_file_rec[COMPNAME__AZIMARK_INDEX] = self.azimuth_mark_list.index(azimark_name)
                
                refpos_determination = meas_type_param2_str
                assert refpos_determination in ADS_DAT_FILE__REFPOS_DETERMINATION_TYPES_STR_LIST, refpos_determination

        # Identifying then removing invalidated refpos records, using Boolean type indexing
        invalidated_refpos_measurement_flag = np.zeros( len(raw_data), dtype=bool) # All elements will starte with "False" value
        for invalidated_refmark_reading_type in ADS_DAT_FILE__REFPOS_DETERMINATION_TYPES_STR_LIST[-2:]:
            assert invalidated_refmark_reading_type[-3:] == b"inv"
            if invalidated_refmark_reading_type in self.meas_type_param2_str_list:
                index = self.meas_type_param2_str_list.index(invalidated_refmark_reading_type)
                invalidated_refpos_measurement_flag[ (refmark_records[COMPNAME__POS_ID] == ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID) & (refmark_records[ COMPNAME__PROGRAM_PARAM2_IDX ] == index) ] = True
        raw_data = raw_data[ ~invalidated_refpos_measurement_flag ]
            
        n = len(raw_data)
        assert n > 0, "No record left after record filters"
        
        return raw_data
    
    
    def separate_ads_raw_dat(self, raw_data):
        
        # self.refpos_measurements = {} # Creating a vocabulary for SU and SD readings
        # self.all_real_refpos_measurement = {} # Creating a vocabulary for SU and SD readings
        # SU_pos_flag = np.zeros( (len(refpos_measurements)), dtype=bool)

        #for refmark_reading_type in ADS_DAT_FILE__REFPOS_DETERMINATION_TYPES_STR_LIST[:-2]:
        #    assert refmark_reading_type[:2] in [ b"SU", b"SD" ]
        #    if (refmark_reading_type[:2] == b"SU") and (refmark_reading_type in  self.meas_type_param2_str_list):
        #        index = self.meas_type_param2_str_list.index(refmark_reading_type)
        #        SU_pos_flag[ refpos_measurements[ COMPNAME__PROGRAM_PARAM2_IDX ] == index  ]   = True

        # self.refpos_measurements["SU"] = refpos_measurements[ SU_pos_flag ]
        # self.refpos_measurements["SD"] = refpos_measurements[ ~SU_pos_flag ]
        # assert False, "refpos_measurements is deprecated!"

        b_has_meas_flags = COMPNAME__MEAS_FLAGS in raw_data.dtype.names # ADS raw files before 2023 did not have 'meas_flags' column        
        # absmeas_window_determination_conditions = ADS_data_absmeas_window_determination_conditions_class()        
        # prev_pos_id = -10
        # prev_mag_pos_id  = -10 # Magnetic pos id's without pseudo null points and MIRA reading
        prev_absmeas_win_index = -1
        # nan_lines = 0
        self.nr_of_M_records_outside_of_file_absmeas_windows = 0
        same_pos_measurement_count = 0 # How many times the same point have been measured. It increasess after a full cycle ends
        file_cycle_start_index = None
        cycle_sequence_id = 0
        n = raw_data.shape[0]
        number_of_file_absmeas_processing_windows  = 0
        file_absmeas_window_data = np.empty( 50, dtype = np.dtype(FILE_ABSMEAS_WINDOW_DATA_DTYPE) )        

        within_cycle = False        
        
        for checked_record_index in range(0,n+1):

            if checked_record_index == n: # We are after the last record
                cycle_end = True
                cycle_start = False
            elif raw_data[checked_record_index][COMPNAME__ADS_FILE__RECORD_WITHIN_PROGRAM]==0:
                cycle_end = within_cycle  # Closing previous cycle
                cycle_start = True

            if cycle_end:
                # Separating mag and reference mark measurement records
                file_window_records = raw_data[file_cycle_start_index:checked_record_index]
                
                mag_meas_record_indeces    =  (np.where(file_window_records[COMPNAME__IS_M_RECORD] ) [0] +  file_cycle_start_index ).tolist()
                refpos_meas_record_indeces =  (np.where(~file_window_records[COMPNAME__IS_M_RECORD]) [0] +  file_cycle_start_index ).tolist()
                
                nr_of_refpos_records = len(refpos_meas_record_indeces)
                
                if nr_of_refpos_records>0:
                    # Separating refmark into a "FL"/"FR" dictionary:
                    if self.instr.get_horizontal_vc_angle("FL", "deg") < 170.0: # Just using 180.0 degree can be uncertain, due to rounding and vc0
                        FL_indeces = file_window_records[refpos_meas_record_indeces][COMPNAME__MEASURED_VC_RAD]< np.pi
                    else:
                        assert self.instr.get_horizontal_vc_angle("FR", "deg") < 170.0
                        FL_indeces = file_window_records[refpos_meas_record_indeces][COMPNAME__MEASURED_VC_RAD] > np.pi
                    FR_indeces = ~FL_indeces
                    refpos_meas_record_indeces_FLFR = {  "FL": refpos_meas_record_indeces[FL_indeces].copy(), "FR": refpos_meas_record_indeces[FR_indeces].copy() }
                else:
                    refpos_meas_record_indeces_FLFR = {  "FL": [], "FR": [] }
                del refpos_meas_record_indeces
                
                # cycle_name = str(cycle_sequence_id)
                cycle_custom_pillar_name = None
                for rec in file_window_records:
                    # In some cases, file measurement window name is in the 'COMPNAME__MEAS_FLAGS' component
                    if b_has_meas_flags and (rec[COMPNAME__MEAS_FLAGS] != MEAS_FLAGS_EMPTY_FIELD_b):
                        meas_flags_dict = separate_meas_flags(rec[COMPNAME__MEAS_FLAGS], rec[COMPNAME__TIMESTAMP]) 
                        #if MEAS_FLAGS_KEYS__CYCLE_NAME in meas_flags_dict:
                        #    cycle_name = meas_flags_dict[MEAS_FLAGS_KEYS__CYCLE_NAME]
                        if MEAS_FLAGS_KEYS__PILLAR in meas_flags_dict:
                            assert (cycle_custom_pillar_name is None) or (cycle_custom_pillar_name == meas_flags_dict[MEAS_FLAGS_KEYS__PILLAR]), "There are more pillar names in a cycle: {}/{}".format(cycle_custom_pillar_name, meas_flags_dict[MEAS_FLAGS_KEYS__PILLAR])
                            cycle_custom_pillar_name = meas_flags_dict[MEAS_FLAGS_KEYS__PILLAR]
                
                if len(mag_meas_record_indeces) >=MIN_MAG_MEAS_POINTS_PER_ABMEAS_PROCESSING_WINDOW:
                    number_of_file_absmeas_processing_windows += 1
                    if file_absmeas_window_data.shape[0] < number_of_file_absmeas_processing_windows:
                        file_absmeas_window_data.resize( (number_of_file_absmeas_processing_windows + 50 ), refcheck = False)
                    
                    act_absmeas_window_data = file_absmeas_window_data[number_of_file_absmeas_processing_windows-1] 
                    act_absmeas_window_data[COMPNAME__FILE_ABMEAS_WIN_INDEX] = number_of_file_absmeas_processing_windows-1
                    file_window_records[COMPNAME__FILE_ABMEAS_WIN_INDEX][:]  = number_of_file_absmeas_processing_windows-1
                    
                    # act_absmeas_window_data['M_start_record_index'] = mag_meas_record_indeces[0]
                    # act_absmeas_window_data['M_end_record_index']   = mag_meas_record_indeces[-1]
                    
                    # act_absmeas_window_data['name'] =  cycle_nam  See notes of the deprecated 'MEAS_FLAGS_KEYS__CYCLE_NAME'
                    act_absmeas_window_data[DATAWINDOW_COMPNAME__M_RECORD_INDECES] = mag_meas_record_indeces
                    act_absmeas_window_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FL"]] = refpos_meas_record_indeces_FLFR["FL"]
                    act_absmeas_window_data[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES["FR"]] = refpos_meas_record_indeces_FLFR["FR"]
                    file_index = file_window_records[0][COMPNAME__INPUT_FILE_INDEX]
                    assert file_index == file_window_records[-1][COMPNAME__INPUT_FILE_INDEX], "A detected cycle comes from more files"                                        
                    pillar_index = self.input_file_data[file_index]["pillar_index"] 
                    if cycle_custom_pillar_name is not None:
                        assert cycle_custom_pillar_name == self.meas_point.pillars[pillar_index]
                    
                    
                    if nr_of_refpos_records>4:
                        # Safety check: we expect no refmark reads at all (in old dat files, or ADS_conrol programs withouth refmark read) or exactly 2x2  reference mark reading (SU+SD at start + SU+SD at end).
                        # Or some refpos records may be uncommented,
                        # Later there can be other type of measuremtns also - so maybe this check can be removed
                        cycle_start_time = raw_data[file_cycle_start_index][COMPNAME__TIMESTAMP]
                        cycle_end_time   = raw_data[file_cycle_end_index][COMPNAME__TIMESTAMP]
                        print("Too many reference position measurements: {}".format(refpos_meas_record_indeces) )
                        print("Cycle's period: {} - {}".format(cycle_start_time, cycle_end_time))
                        mu.list_timelabels("Ref. pos. time labels:",  raw_data[COMPNAME__TIMESTAMP][ refpos_meas_record_indeces  ] , 20)
                        assert False
                    cycle_sequence_id += 1
                else:
                    self.nr_of_M_records_outside_of_file_absmeas_windows += len(mag_meas_record_indeces)
                    
                # Removing last cycle start index:
                file_cycle_start_index = None
                within_cycle = False
                cycle_end = False
                
            if cycle_start:
                assert not within_cycle
                file_cycle_start_index = checked_record_index
                # prev_absmeas_win_index = raw_data[COMPNAME__FILE_ABMEAS_WIN_INDEX][i]
                within_cycle = True
                cycle_start = False
                # prev_pos_id = -1000
                # prev_mag_pos_id = -1000
            #if (within_cycle):
                # nof_records_in_current_cycle = nof_records_in_current_cycle  + 1
                # prev_pos_id = new_pos_id
                #if (new_pos_id>=0): prev_mag_pos_id = new_pos_id
            if not within_cycle:
                self.nr_of_M_records_outside_of_file_absmeas_windows += 1

        assert not within_cycle

        file_absmeas_window_data.resize( number_of_file_absmeas_processing_windows, refcheck = False )
        
        # Note: M_readings can contain not used mag records!
        M_readings       = raw_data[ np.where( raw_data[COMPNAME__IS_M_RECORD])  ].copy()
        refmark_readings = raw_data[ np.where(~raw_data[COMPNAME__IS_M_RECORD])  ].copy()
        
        if len(refmark_readings) == 0:
            refmark_readings_separated = None
        else:
            if self.instr.get_horizontal_vc_angle("FL", "deg") < 170.0:
                FL_record_indeces = raw_data[refmark_readings][COMPNAME__MEASURED_VC_RAD] < np.pi
            else:
                assert self.instr.get_horizontal_vc_angle("FR", "deg") < 170.0
                FL_indeces = raw_data[refmark_readings][COMPNAME__MEASURED_VC_RAD] > np.pi
            FR_indeces = ~FL_indeces
            refmark_readings_separated = { "FL": raw_data[refmark_readings][FL_record_indeces].copy(), "FR": raw_data[refmark_readings][FL_record_indeces].copy()} 
                                            
        del raw_data                                      
        
        
        
        return file_absmeas_window_data, M_readings, refmark_readings_separated

    # Separates file records into "processing windows".
    # A "processing windows" is a subset of "file cycles": it may contain the same or less records than a file cycle.
    # A processing cycle typically has one set of reference mark reading (SU+SD before and after the mangetic part).
    # 'calc_var_filters': array (or tuple or list) of filter definitions for calculation variations.
    # For "ALL-IN-ONE" processing, we associate a refmark data for each record
    def prepare_processing_windows(self, nof_concatenated_file_cycles):
        
        print("Grouping magnetometer measurements into processing window data sets...", end='')
        assert len(self.file_absmeas_window_data) > 0
        # assert len(self.file_cycle_end_index) > 0
        # assert len(self.file_cycle_end_index) ==  self.number_of_absmeas_processing_windows
        self.nr_of_calcvars = len(self.calc_vars)
        assert self.nr_of_calcvars >= 1
        # assert self.processing_window_data == None

        if nof_concatenated_file_cycles == ALL_MEASUREMENTS_TOGETHER:
            assert False
            nof_concatenated_file_cycles = self.number_of_absmeas_processing_windows
            assert self.refpos_data_calculated_for_each_record
            assert not self.refpos_data_calculated_for_each_processing_cycle
            self.all_cycles_together_mode = True
        else:
            # calculate_processing_cycles_refpos_data = self.associate_with_MIRA_data  # The calculation is done later, see below
            self.all_cycles_together_mode = False

        # self.processing_cycles_names = []
        # self.processing_cycles_pillar_names = []
        
        for calc_var in self.calc_vars:
            calc_var.processing_window_data = []

        self.processing_cycles__concatenated_cycles = nof_concatenated_file_cycles

        assert nof_concatenated_file_cycles >= 1
        
        number_of_file_absmeas_windows = len(self.file_absmeas_window_data)
        self.number_of_processing_windows = number_of_file_absmeas_windows - nof_concatenated_file_cycles + 1
        self.longest_proessing_window_name = 0
        
           
        for i_file_cycle in range( self.number_of_processing_windows ):

            file_program_def =  self.file_absmeas_window_data[i_file_cycle]
            window_name = f"AM.#{i_file_cycle+1}"  #  file_program_def['name']
            for joined_file_cycle_index in range( i_file_cycle + 1,  min( i_file_cycle +  nof_concatenated_file_cycles, number_of_file_absmeas_windows+1) ):
                window_name = window_name + "+"+self.file_absmeas_window_data[joined_file_cycle_index]['name']
            if len(window_name)>10: window_name = window_name[:8]+"..."

            processing_cycle_pillar_names = [  self.meas_pos.pillars[ file_program_def['pillar-index'] ] ]
            for joined_file_cycle_index in range( i_file_cycle + 1,  min( i_file_cycle +  nof_concatenated_file_cycles, number_of_file_absmeas_windows+1) ):
                pillar_name = self.meas_po.pillarss[ self.file_absmeas_window_data[joined_file_cycle_index]['pillar-index'] ]
                if not pillar_name in processing_cycle_pillar_names:
                    processing_cycle_pillar_names += [ pillar_name ]
            processing_cycle_pillar_name = '/'.join(processing_cycle_pillar_names)

            for i_calc_var, calc_var in enumerate(self.calc_vars):
                
                new_processing_windows_data = processing_window_data_class()
                calc_var.processing_window_data += [ new_processing_windows_data ]
                
                new_processing_windows_data.name = window_name
                if len(new_processing_windows_data.name)> self.longest_proessing_window_name: self.longest_proessing_window_name = len(new_processing_windows_data.name)
                new_processing_windows_data.window_index = len(calc_var.processing_window_data)-1
                new_processing_windows_data.pillar_name  = processing_cycle_pillar_name
                new_processing_windows_data.M_record_indeces = []  # Can be different in each calc vars
                new_processing_windows_data.DIM_data_refmark_record_indeces = { "FL": [], "FR": [] }  # Can be different in each calc vars
                new_processing_windows_data.refmark_records = []  # Will be filled in the next step
                new_processing_windows_data.plot_color = plot_utils.rainbow_color_get( new_processing_windows_data.window_index, number_of_file_absmeas_windows )

                # Determing all records, filtered by 'calc_var.cycle_point_filters'
                # new_processing_windows_data.all_record_indeces = self.get_filtered_file_record_indeces( file_program_def['start_record_index'], self.file_absmeas_window_data[i_file_cycle + nof_concatenated_file_cycles - 1]['end_record_index'], calc_var.cycle_point_filters )
                M_record_indeces = file_program_def[DATAWINDOW_COMPNAME__M_RECORD_INDECES]
                too_few_M_samples = False
                if calc_var.use_only_sample_indeces:
                    if len(M_record_indeces)<=max(calc_var.use_only_sample_indeces):
                        too_few_M_samples = True
                    else:
                        M_record_indeces = np.array(M_record_indeces)[calc_var.use_only_sample_indeces]
                        M_record_indeces = M_record_indeces.tolist()
                    
                new_processing_windows_data.M_record_indeces       +=   M_record_indeces
                
                new_processing_windows_data.first_M_timelabel = self.M_readings[ new_processing_windows_data.M_record_indeces[0] ][COMPNAME__TIMESTAMP]  # Or: self.processing_window_M_date_limits[processing_window_index][0]
                new_processing_windows_data.last_M_timelabel  = self.M_readings[ new_processing_windows_data.M_record_indeces[-1]][COMPNAME__TIMESTAMP] # Or: self.processing_window_M_date_limits[processing_window_index][1]
                # new_processing_windows_data.start_datetime = self.M_readings[ new_processing_windows_data.M_record_indeces[0]][COMPNAME__TIMESTAMP] 
                # new_processing_windows_data.end_datetime   = self.M_readings[ new_processing_windows_data.M_record_indeces[-1]][COMPNAME__TIMESTAMP] 
                
                if i_file_cycle == 0:
                    calc_var.all_processing_windows_first_M_timelabel = new_processing_windows_data.first_M_timelabel
                    calc_var.all_processing_windows_last_M_timelabel  = new_processing_windows_data.last_M_timelabel
                    calc_var.nr_of_excluded_processing_windows = 0
                else:
                    calc_var.all_processing_windows_first_M_timelabel = min( calc_var.all_processing_windows_first_M_timelabel, new_processing_windows_data.first_M_timelabel)
                    calc_var.all_processing_windows_last_M_timelabel  = max( calc_var.all_processing_windows_last_M_timelabel,  new_processing_windows_data.last_M_timelabel)
                    
                for FLFR in ["FL", "FR"]:
                    new_processing_windows_data.DIM_data_refmark_record_indeces[FLFR] +=  file_program_def[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES[FLFR]] 
                new_processing_windows_data.file_id_str = self.get_cycle_file_id_str(new_processing_windows_data.M_record_indeces)
                
                new_processing_windows_data.excluded_by_data_selection_filter = False
                new_processing_windows_data.excluded_by_sample_number = False
                half_second = 0.5 /86400.
                for exclude_datetime in self.exclude_windows_with_dates:
                    
                    if (new_processing_windows_data.first_M_timelabel <= exclude_datetime+half_second) and (new_processing_windows_data.last_M_timelabel >= exclude_datetime-half_second):
                        new_processing_windows_data.excluded_by_data_selection_filter = True
                        calc_var.nr_of_excluded_processing_windows += 1
                        break
                    
                if too_few_M_samples:
                    new_processing_windows_data.excluded_by_sample_number = True
                    calc_var.nr_of_excluded_processing_windows += 1
                    break
                    
                # ToDo: Now 'calc_var.cycle_point_filters' is not applied for 'cycle_own_refpos_record_indeces'. 'calc_var.cycle_point_filters' is applied only for mag records.
                # It is ok?

                # Pairing file cycle's refpos data with processing cycle's refpos data - for concatenated file cycles we will ensure that the joined file cycles have the same refpos data, see below
                #new_processing_windows_data.cycle_own_refpos_record_indeces = file_program_def[DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES]
                #for joined_file_cycle_index in range( i_file_cycle + 1,  min( i_file_cycle +  nof_concatenated_file_cycles, self.number_of_absmeas_processing_windows+1) ):
                #    for i_refpos_record_index in self.file_absmeas_window_data[joined_file_cycle_index][DATAWINDOW_COMPNAME__REFMARK_RECORD_INDECES]:
                #        if not i_refpos_record_index in new_processing_windows_data.cycle_own_refpos_record_indeces:
                #            if not self.all_cycles_together_mode:
                #                assert False, "More refpos records are joined in a processing cycle"
                #            else:
                #                new_processing_windows_data.cycle_own_refpos_record_indeces += [i_refpos_record_index]

                # Setting 'COMPNAME__PROCESSING_CYCLE_INDEX_' fields for records that belong to this processing cycles:
                #for rec_index in new_processing_windows_data.all_record_indeces:
                #    fieldname = COMPNAME__PROCESSING_CYCLE_INDEX_ + calc_var.id
                #    if raw_data[fieldname][rec_index]==-1:
                #        # This field is not yet filled
                #        raw_data[fieldname][rec_index] = act_processing_cycle_index
                #    else:
                #        # This field is filled: it belongs to more processing cycles:
                #        raw_data[fieldname][rec_index] = -2
                
        
        # Determining 'mag_only_record_indeces' of 'processing_window_data'
        #for cycle_index in range( self.number_of_processing_windows() ) :
        #    for i_calc_var in range(len(self.calc_vars)):
        #        # file_record_indeces_with_refpos_records =
        #        # Removing reference positions measurements' indeces:
        #        self.processing_window_data[i_calc_var][cycle_index].mag_only_record_indeces =  [ i for i in self.processing_window_data[i_calc_var][cycle_index].all_record_indeces if i not in self.processing_window_data[i_calc_var][cycle_index].cycle_own_refpos_record_indeces  ]
        #        for icheck in self.processing_window_data[i_calc_var][cycle_index].mag_only_record_indeces:
        #            assert raw_data[COMPNAME__POS_ID][icheck] != ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID
        
        print("  done") # , number of processing windows: {}".format(self.number_of_processing_windows))

        self.all_calc_vars_min_timelabel = 99999999.
        self.all_calc_vars_max_timelabel = 0.
        
        for calc_var in self.calc_vars:
            self.all_calc_vars_min_timelabel = min(self.all_calc_vars_min_timelabel, calc_var.all_processing_windows_first_M_timelabel) 
            self.all_calc_vars_max_timelabel = max(self.all_calc_vars_max_timelabel, calc_var.all_processing_windows_last_M_timelabel)
            if calc_var.nr_of_excluded_processing_windows > 0:
                print("Note: {} processing window{} removed from further analyzis by a date selection filter from calc. var \"{}\"".format( calc_var.nr_of_excluded_processing_windows, mu.plural_suffix(calc_var.nr_of_excluded_processing_windows,"s is", " are"), calc_var.name   ))
        print()

        # Applying reference mark measurement:        
        print("Assinging azimuth mark readings to processing windows... ",end='')
        nr_of_noncomplete_processing_windows = 0
        
            
        for i_calc_var, calc_var in enumerate(self.calc_vars):            
            calc_var.processing_window_indeces_with_noncomplete_refmark_data = []
           
            for processing_window_index, processing_window_data in enumerate(calc_var.processing_window_data) :

                if processing_window_data.excluded_by_data_selection_filter: continue
                # Using first usable reference positions measured data
                complete_azimuth_measurement_data = None

                for refpos_meas_data_source_def in calc_var.refmark_measured_data_sources:

                    if refpos_meas_data_source_def is None: continue
                
                    assert isinstance(refpos_meas_data_source_def, dim_azimuth.azimark_measured_data_source_def_class)
                    
                    match refpos_meas_data_source_def.source_type:
                        
                        case dim_azimuth.azimark_meas_data_source_type_enum.dim_dat_files: 
                            complete_azimuth_measurement_data, pillar_index = self.determine_processing_window_complete_azimuth_data_from_DIM_azi_data(refpos_meas_data_source_def, processing_window_data, calc_var, True, processing_window_data.first_M_timelabel, processing_window_data.last_M_timelabel)
                        
                        case dim_azimuth.azimark_meas_data_source_type_enum.spreadsheet_table:  # if refpos_meas_data_source_def.type=="external-table":
                            # complete_azimuth_measurement_data = calc_var.true_angle_calculator.get_complete_azimuth_measurement_data_from_ext_table( refpos_meas_data_source_def, calc_var.instr, [from_timelabel,to_timelabel] )
                            # Excel table should define a complete (SU/SD) reading:
                            # if complete_azimuth_measurement_data is not None: assert complete_azimuth_measurement_data.is_complete()
                            complete_azimuth_measurement_data, pillar_index = calc_var.true_angle_calculator.get_complete_azimuth_measurement_data_from_load_data( self.refmark_readings[refpos_meas_data_source_def.refmark_readings_array_index], from_timelabel,  to_timelabel)
                            
                        case dim_azimuth.azimark_meas_data_source_type_enum.separate_dat_files:   # elif refpos_meas_data_source_def.type=="refpos-dat-file" :
                            complete_azimuth_measurement_data, pillar_index = calc_var.true_angle_calculator.get_complete_azimuth_measurement_data_from_load_data( refpos_meas_data_source_def, self.refmark_readings[refpos_meas_data_source_def.refmark_readings_array_index], processing_window_data.first_M_timelabel,  processing_window_data.last_M_timelabel, self.instr, calc_var)
                            
                        case _:
                            assert False, "Unsupported reference data source: {}".format(refpos_meas_data_source_def.source_type)
                            
                    if complete_azimuth_measurement_data is not None: break # First valid ref-pos data is found

                if complete_azimuth_measurement_data is None:
                    # Cycle has now own refpos data (any or not enough MIRA readings): we here allow passing through, but will not process this cycle later
                    # assert False, "No reference position data for processing cycle "+self.processing_window_data[i_calc_var][processing_window_index].get_identification_info_str()
                    nr_of_noncomplete_processing_windows += 1
                    calc_var.processing_window_indeces_with_noncomplete_refmark_data += [processing_window_index]
                    pass
                else:
                    # assert complete_azimuth_measurement_data.WGS84 is not None, "No WGS84 coordinate for cycle {}".format(processing_window_index)
                    processing_window_data.pillar_index = pillar_index                                        
                    assert complete_azimuth_measurement_data.is_complete()
                    
                    # Getting out the values of the selected magnetometer sensor (needs for multi-magnetometer sensor configuration tests)
                    # dim_mag = self.M_readings[ M_record_indeces ][COMPNAME__MAG] 
                    # new_processing_windows_data.NM_position_indeces, measurement_balanced =  du.determine_DIM_measurement_positions(dim_mag, hc_NED, vc_theo_rad, DIM_instr, processing_window_raw_data[dim_dat_file_reader.COMPNAME__TELESCOPE_DIRECTION], processing_window_raw_data[dim_dat_file_reader.COMPNAME__TELESCOPE_POSITION])
                    
                processing_window_data.complete_azi_meas_data = complete_azimuth_measurement_data
                # print(complete_azi_meas_data.refpos_measured_data_source)
            
            if len(calc_var.processing_window_indeces_with_noncomplete_refmark_data)>0:
                print("   nr. of processing windows without ref. mark readings in calc.param \"{}\" : {}".format(calc_var.name, len(calc_var.processing_window_indeces_with_noncomplete_refmark_data)))
       
        
       #  all_calc_vars_all_processing_window_data_M_date_limits = [abs_min_timelabel, abs_max_timelabel ]            

        assert nr_of_noncomplete_processing_windows < len(self.calc_vars) * self.number_of_processing_windows, "None of the processing windows has valid azimuth mark data!"

        self.refpos_data_calculated_for_each_processing_cycle = True
        print(" done")
        print()
        
        self.determine_refpos_data_for_each_record()
        
    def print_processing_windows(self):
        # -------------------------------------------------------------
        # -- Printing processing window statistics in a text table  ---
        # -------------------------------------------------------------
        entire_period = self.all_calc_vars_max_timelabel - self.all_calc_vars_min_timelabel # Note: this period can include processing windows with incomplete azimuth readings and processing windows excluded by dates
        print("DIM measurement data are separated into {} processing windows.".format(  self.number_of_processing_windows))
        print("Entire time period of all windwos: {:.1f} day(s)".format( entire_period  )  )
        print("Processing windows properties:")
        header =  ["Win\nidx","File\nID","Name", "Date/time"  ]
        for icalcvar in range( self.nr_of_calcvars ): header += [ "calcvar#{:}\nazi-mark".format(icalcvar+1) , "calcvar#{:}\nrecords".format(icalcvar+1), "calcvar#{:}\ndur.(min)".format(icalcvar+1), "calcvar#{:}\nproblem".format(icalcvar+1) ]
        table = [ header ]

        # We use a two line header
        # print("idx",end="")
        # for icalcvar in range( self.nr_of_calcvars ): print("\t____calc-var-{}____".format(icalcvar+1),end='')
        # print("")
        # for icalcvar in range( self.nr_of_calcvars ): print("\trecs\tdur.(min)\tMIRA-id",end='')
        # print("")
        # print("========="+"===============" * len(calc_var_filters))
        
        has_valid_window = False

        for window_index in range( self.number_of_processing_windows ):

            # NOTE: These parameters are determined from the first calc.var only!
            file_cycle_ID  = self.get_cycle_file_id_str(  self.calc_vars[0].processing_window_data[window_index].M_record_indeces )
            mid_time = np.mean( self.M_readings[ [ self.calc_vars[0].processing_window_data[window_index].M_record_indeces[0], self.calc_vars[0].processing_window_data[window_index].M_record_indeces[-1] ] ][COMPNAME__TIMESTAMP]  )
            mid_time_str = mu.exceltime_2_datetime(mid_time).strftime( preferences.UI_DATETIME_FORMAT_SEC  )

            record = [ window_index, file_cycle_ID, self.calc_vars[0].processing_window_data[window_index].name , mid_time_str ]
            for i_calcvar, calc_var in enumerate(self.calc_vars):

                processing_window_data = calc_var.processing_window_data[window_index]

                nr_of_mag_only_records = len( processing_window_data.M_record_indeces)
                nr_of_mag_and_azimark_records = nr_of_mag_only_records + (0 if processing_window_data.complete_azi_meas_data is None else  processing_window_data.complete_azi_meas_data.nr_of_data())

                # start_time_with_refmark_meas    = min(start_time_without_refmark_meas, 
                # end_time_with_refmark_meas      = raw_data[ processing_window_data.M_record_indeces[-1] ][COMPNAME__TIMESTAMP]
                mag_only_meas_time =  processing_window_data.last_M_timelabel -  processing_window_data.first_M_timelabel
                if processing_window_data.complete_azi_meas_data is None:
                    mag_and_azimark_meas_time  = mag_only_meas_time
                else:
                    mag_and_azimark_meas_time    = max( processing_window_data.complete_azi_meas_data.last_timelabel(), processing_window_data.last_M_timelabel)  - min(processing_window_data.complete_azi_meas_data.first_timelabel(), processing_window_data.first_M_timelabel)                
                assert mag_and_azimark_meas_time >= mag_only_meas_time

                # print("\t{:}\t{:.1f}\t{:}".format(nr_of_records, (end_time - start_time)* 1440.0, self.processed_cycles_MIRA_ids[window_index]) ,end='')
                if processing_window_data.complete_azi_meas_data is not None:
                    refpos_data_source = processing_window_data.complete_azi_meas_data.meas_data_source_str
                    max_len=10;
                    if len(refpos_data_source)>max_len:
                        refpos_data_source = "..."+refpos_data_source[-max_len:] # Maximizing file name len
                else:
                    refpos_data_source = "n.a."
                record += [refpos_data_source ,  "{:}+{:}".format(nr_of_mag_only_records, nr_of_mag_and_azimark_records-nr_of_mag_only_records) , "{:.1f}/{:.1f}".format( (mag_only_meas_time)* 1440.0,    mag_and_azimark_meas_time * 1440.0 ) ]
                
                problem_str = processing_window_data.get_problem()
                if problem_str:
                    record += [ problem_str ]
                else:
                    record += [ "-" ]
                    has_valid_window = True
            # print("")
            table += [ record  ]
        print(tabulate.tabulate(table, headers='firstrow', tablefmt='simple'))  # fancy_grid
        return has_valid_window
        

    


    # Processes the DIM records that stores the reference position readings, to get a complete reference reading for a processing window.
    # It creates and returns an 'complete_azimuth_measurement_data' class instance, that represents the measurements.
    def determine_processing_window_complete_azimuth_data_from_DIM_azi_data(self, refpos_source_def, processing_window_data, calc_params, auto_complete_with_nearby_measurements, cycle_start_timelabel, cycle_end_timelabel ) :

        assert isinstance(refpos_source_def, dim_azimuth.azimark_measured_data_source_def_class)
        assert isinstance(processing_window_data, processing_window_data_class)
        
        #azimark_id = None
        #hc_readings = {"FL":[], "FR":[]}
        #vc_readings = {"FL":[], "FR":[]}
        #timelabels  = {"FL":[], "FR":[]}
        #used_refpos_reading_indeces = { "FL": [], "FR": [] }
        # def use_measurement(refpos_data_record, FLFR, refpos_data_record_index, pillar_index):
        #     new_azimark_id, refpos_determination, hc_rad, vc_rad, hc_tick, vc_tick, valid_for_sensor_pos = self.process_refpos_raw_record_data(refpos_data_record, FLFR, refpos_data_record_index, calc_params)
        #     nonlocal azimark_id, used_refpos_reading_indeces, hc_readings, vc_readings, timelabels
        #     if azimark_id is None:
        #         azimark_id = new_azimark_id
        #     else:
        #         assert new_azimark_id == azimark_id, 'Azimark-id changed from "{}" to "{}"'.format(azimark_id,new_azimark_id)
        #     hc_readings[valid_for_sensor_pos] += [ hc_rad ]
        #     vc_readings[valid_for_sensor_pos] += [ vc_rad ]
        #     timelabels[valid_for_sensor_pos] += [ refpos_data_record[COMPNAME__TIMESTAMP] ]
        #     assert not refpos_data_record_index in used_refpos_reading_indeces[FLFR]
        #     used_refpos_reading_indeces[FLFR] += [ refpos_data_record_index ]

        # First, using cycle own refpos measurements:
        pillar_index = None
        # WGS84_check = None
        complete_azimuth_measrement = dim_azimuth.complete_azimuth_measurement_data_class()
        
        for FLFR in ["FL","FR"]:
            
            # indeces = processing_window_data.DIM_data_refmark_record_indeces[FLFR]
            # data = self.refmark_readings[ refpos_source_def.refmark_readings_array_index ][FLFR]
            # print(FLFR,":",indeces,  "pillar-index:",   data[ indeces[0] ]["pillar_index"])
            
            azi_meas_array = self.refmark_readings[ refpos_source_def.refmark_readings_array_index ][FLFR]
            
            for azimark_array_index in processing_window_data.DIM_data_refmark_record_indeces[FLFR]:
                azimark_read_record = azi_meas_array[azimark_array_index]
                new_pillar_index  = azimark_read_record[COMPNAME__PILLAR_INDEX]
                azimark_index = azimark_read_record[COMPNAME__AZIMARK_INDEX]
                if pillar_index is None:
                    # First record
                    pillar_index = new_pillar_index
                else:
                    assert pillar_index == new_pillar_index,     "Pillar changed in a measurement window"
                #if WGS84_check is None:                    
                #    WGS84_check = azimark_read_record[COMPNAME__AZIMARK_WGS84]
                #    assert np.isclose(WGS84, WGS84_check, atol = angle_utils.angle_comparison_tolerance["rad"])
                #else:
                #    assert np.isclose(azimark_read_record[COMPNAME__AZIMARK_WGS84], WGS84_check, atol = angle_utils.angle_comparison_tolerance["rad"])
                
                # use_measurement(azimark_read_record, FLFR,  azimark_array_index, pillar_index)
                if self.instr.axis_angle_data_source == "digital":
                    hc_tick = azimark_read_record[COMPNAME__HC_TICK]
                    vc_tick = azimark_read_record[COMPNAME__VC_TICK]
                    hc = self.instr.hc_axis( hc_tick, calc_params.hc_precorr_funcs)
                    vc = self.instr.hc_axis( vc_tick, calc_params.vc_precorr_funcs)
                else:
                    hc_tick = None
                    vc_tick = None
                    hc = azimark_read_record[COMPNAME__MEASURED_HC_RAD]
                    vc = azimark_read_record[COMPNAME__MEASURED_VC_RAD]
                    
                complete_azimuth_measrement.use_measurement(FLFR, azimark_index, azimark_read_record[COMPNAME__TIMESTAMP], hc, vc, hc_tick, vc_tick ) 
                

        # Then, using nearby measurements:
        if not complete_azimuth_measrement.is_complete():
            # SU or SD is missing
            if not auto_complete_with_nearby_measurements: return None # We need a complete refmark measurement

            assert False
            # for SU_SD in ["SU","SD"]:
        	   #  if len( hc_readings[SU_SD] )==0:
            # 		nearest_from_start_index = (np.abs(self.refpos_measurements[SU_SD][COMPNAME__TIMESTAMP] - cycle_start_timelabel)).argmin()
            # 		nearest_from_end_index   = (np.abs(self.refpos_measurements[SU_SD][COMPNAME__TIMESTAMP] - cycle_end_timelabel)).argmin()
            # 		nearest_from_start_timediff = np.abs(self.refpos_measurements[SU_SD][nearest_from_start_index][COMPNAME__TIMESTAMP] - cycle_start_timelabel)
            # 		nearest_from_end_timediff   = np.abs(self.refpos_measurements[SU_SD][nearest_from_end_index][COMPNAME__TIMESTAMP] - cycle_end_timelabel)
            # 		time_diff = min((nearest_from_start_timediff,nearest_from_end_timediff))
            # 		if time_diff > refpos_source_def.max_time_diff_to_use_a_refmark_reading__days: return None
            # 		used_index = nearest_from_start_index if nearest_from_start_timediff  < nearest_from_end_timediff else nearest_from_end_index
            # 		use_measurement(self.refpos_measurements[SU_SD][used_index])

        
        #if result.WGS84 is None:
        #    refpos_id, refpos_determination, hc_rad, vc_rad, hc_tick, vc_tick, valid_for_sensor_pos = self.process_refpos_raw_record_data(refpos_meas_indeces[0], calc_var)
        #    assert False, "No WGS84 data for record indeces {} (first refpos-id:{}, first record datetime:{})".format(refpos_meas_indeces,  refpos_id, raw_data[refpos_meas_indeces[0]]["datetime"] )
        
        

        # data_source = "dat-file[{:}]".format(len(used_refpos_reading_indeces)) + extra_source  # This will extend to like this: "dat-file[10, 11, 30, 31]"
        meas_data_source = "DIM-dat-file"
        #valid_from = self.refmark_readings[0] [min(used_refpos_reading_indeces)][COMPNAME__TIMESTAMP]
        #valid_to   = self.refmark_readings[0] [max(used_refpos_reading_indeces)][COMPNAME__TIMESTAMP]
        complete_azimuth_measrement.end_of_data_feed(meas_data_source)
        # complete_azimuth_measrement.set_aux_params(data_source, valid_from, valid_to) 

        return complete_azimuth_measrement, pillar_index


    # Must be called after 'prepare_processing_cycles()' done its task.
    # Required also for ALL-IN-ONE-PROCESSING.
    # Can also be used to check load data.
    def determine_refpos_data_for_each_record(self):
        print("Using azimuth mark readings data to calculate true angles... ", end='');
        assert self.number_of_processing_windows > 0
        assert self.refpos_data_calculated_for_each_processing_cycle  # This is a pre-requisite for this function

        # Adding the new components to 'raw_data'
        # fill_val_nan = np.empty( len(self.M_recordings) )
        # fill_val_nan[:] = np.nan
        # added_comp_names = []
        # added_comp_fillvals = []
        # added_comp_dypes = []
        # for calc_var in self.calc_vars:
        #     added_comp_names    += [ COMPNAME__REFPOS_HC_THEO_SU_TICK_ + calc_var.id, COMPNAME__REFPOS_HC_THEO_SD_TICK_ + calc_var.id, COMPNAME__REFPOS_HC_WGS84_SU_ + calc_var.id, COMPNAME__REFPOS_HC_WGS84_SD_ + calc_var.id, COMPNAME__PRECALC_HC_NED_ + calc_var.id, COMPNAME__PRECALC_VC_NED_ + calc_var.id]
        #     added_comp_fillvals += [ fill_val_nan,fill_val_nan,fill_val_nan,fill_val_nan,fill_val_nan,fill_val_nan]
        #     added_comp_dypes    += [ float,float,float,float,float,float]
        # raw_data =  rcf.append_fields( raw_data, added_comp_names, added_comp_fillvals, added_comp_dypes, usemask=False )


        # self.used_file_reference_position_ids = []
        # self.used_file_reference_positions = []
        # self.valid_file_SU_refmark_reads = np.empty( ( len(self.calc_vars), self.nof_reference_pos_measurement_records), dtype=[ ("datetime",float), ("refpos_index", int ), ("hc", float),("vc",float)   ]  )
        # self.valid_file_SD_refmark_reads = np.empty( ( len(self.calc_vars), self.nof_reference_pos_measurement_records), dtype=[ ("datetime",float), ("refpos_index", int ), ("hc", float),("vc",float)   ]  )
        # n = len(raw_data)

        for calc_var_index, calc_var in enumerate(self.calc_vars):

            # We iterate through processing cycles, as we may want to use both the start and end reference mark readings of an absolute measurement
            for processing_window_index, processing_window_data in enumerate( calc_var.processing_window_data ):

                complete_azi_meas_data = processing_window_data.complete_azi_meas_data
                if complete_azi_meas_data is None: continue
                assert isinstance(complete_azi_meas_data, dim_azimuth.complete_azimuth_measurement_data_class)
                # assert complete_azi_meas_data.WGS84 is not None, "Missing WGS84 data"
                # print()
                # complete_azi_meas_data.print_out()

                # cycle_first_timelabel =  raw_data[ self.processing_window_data[calc_var_index][processing_window_index].M_record_indeces[0]][COMPNAME__TIMESTAMP]
                
                # print( np.degrees(complete_azi_meas_data.WGS84) )
                
                for azi_meas_data_source in calc_var.refmark_measured_data_sources:
                    if azi_meas_data_source.azimark_name is not None:  
                        # In some cases, azimuth mark name is not read from data files, but defined from YAML config file: ' azimark_name = run_config.pop("azimuth-mark")'
                        # In such cases, azimuth mark name is in 'azi_meas_data_source'
                        azimark_index = complete_azi_meas_data.get_azimark_index()
                        assert azi_meas_data_source.azimark_name == self.azimuth_mark_list[ azimark_index ]
                
                calc_var.true_angle_calculator.use_theo_measurement_data(True, complete_azi_meas_data, self.meas_pos.pillars[ processing_window_data.pillar_index], self.azimuth_mark_list ) # , cycle_first_timelabel = cycle_first_timelabel

                # rec_indeces = processing_window_data.M_record_indeces

                # print(calc_var.id, processing_window_index, self.number_of_processing_windows, raw_data[rec_indeces[0]][COMPNAME__TIMESTAMP])
                processing_window_data.M_reading_NED_angles = {  }
                processing_window_data.M_reading_NED_angles["hc"] = calc_var.true_angle_calculator.h_theo_to_ned(  calc_var.M_reading_theo_angles["hc"] [ processing_window_data.M_record_indeces  ]  )
                vc_theo = calc_var.M_reading_theo_angles["vc"][ processing_window_data.M_record_indeces  ] 
                processing_window_data.M_reading_NED_angles["vc"] = dim_azimuth.v_theo_to_ned(  vc_theo )   
                
                # Checking if sensor-position/telescope-orientation read from file correpond to the expected position/orientation
                # processing_window_data.telescope_pos_data = np.empty ( len( processing_window_data.M_record_indeces ), dtype = [ (COMPNAME__TELESCOPE_DIRECTION, TELESCOPE_DIRECTION__TYPE), (COMPNAME__TELESCOPE_POSITION, TELESCOPE_POSITION__TYPE) ]    )
                
                # calculated_telescope_pos_data = du.calc_telescope_pos_data(processing_window_data.M_reading_NED_angles["hc"], vc_theo, self.instr)  
                
                for i in range(len(processing_window_data.M_record_indeces)):
                    M_reading = self.M_readings[ processing_window_data.M_record_indeces[i]   ]  # To overwrite (fix) to original values, we have to use a direct record reference for the original ndarray. 
                    calc_telescope_cardinal_direction = du.calc_NED_cardinal_direction( processing_window_data.M_reading_NED_angles["hc"][i], vc_theo[i], self.instr )
                    calc_telescope_face_position      = du.telescope_face_pos_enum.FL.value if self.instr.is_telescope_in_reference_direction(vc_theo[i]) else du.telescope_face_pos_enum.FR.value
                    assert dim_azimuth.GEODETIC_REFERENCE_FACE_DIRECTION == "FL"                    
                    if M_reading[COMPNAME__TELESCOPE_CARDINAL_DIRECTION] != calc_telescope_cardinal_direction:
                        if M_reading[COMPNAME__TELESCOPE_CARDINAL_DIRECTION] != du.UNKNOWN_CARDINAL_DIRECTION:
                            self.add_warning(f'Mismatch between file ("{M_reading[COMPNAME__TELESCOPE_CARDINAL_DIRECTION]}") and calculated ("{calc_telescope_cardinal_direction}") telescope directions')
                        M_reading[COMPNAME__TELESCOPE_CARDINAL_DIRECTION]  = calc_telescope_cardinal_direction
                    if M_reading[COMPNAME__TELESCOPE_FACE_POSITION] != calc_telescope_face_position:                        
                        if M_reading[COMPNAME__TELESCOPE_FACE_POSITION] != du.telescope_face_pos_enum.unknown.value:
                            self.add_warning(f'Mismatch between file ("{M_reading[COMPNAME__TELESCOPE_FACE_POSITION]}") and calculated ("{calc_telescope_face_position}") telescope position')
                        M_reading[COMPNAME__TELESCOPE_FACE_POSITION]  = calc_telescope_face_position
                        
                # Saving some azimuth measurement infro into 'processing_window_data'
                processing_window_data.azimark_readings_avg_angle      =  calc_var.true_angle_calculator.act_measured_angles_avg()
                
                processing_window_data.azimark_readings_avg_FL_FR_diff = {}
                processing_window_data.azimark_readings_avg_FL_FR_diff["hc"] = calc_var.true_angle_calculator.last_measured_refmark_theo_angle["hc"]["FL"] - calc_var.true_angle_calculator.last_measured_refmark_theo_angle["hc"]["FR"]
                if processing_window_data.azimark_readings_avg_FL_FR_diff["hc"]<0.0: 
                    processing_window_data.azimark_readings_avg_FL_FR_diff["hc"] += np.pi
                else:
                    processing_window_data.azimark_readings_avg_FL_FR_diff["hc"] -= np.pi

                processing_window_data.vc0_from_azimark_readings = np.nan
                processing_window_data.azimark_readings_avg_FL_FR_diff["vc"] = np.nan
                if calc_var.true_angle_calculator.last_measured_refmark_theo_angle["vc"]["FL"] is not None:
                    processing_window_data.azimark_readings_avg_FL_FR_diff["vc"] = calc_var.true_angle_calculator.last_measured_refmark_theo_angle["vc"]["FL"] - du.get_complement_vc_angle(calc_var.true_angle_calculator.last_measured_refmark_theo_angle["vc"]["FR"])
                    if calc_var.true_angle_calculator.last_measured_refmark_true_angle["vc"]["FL"] is not None:
                        processing_window_data.vc0_from_azimark_readings = calc_var.true_angle_calculator.get_last_measured_vc0()
                        
                        
                if calc_var.true_angle_calculator.nr_of_used_WGS84_dirs==1:
                    processing_window_data.azimark_readings_true_calibrating_angle =  calc_var.true_angle_calculator.last_measured_refmark_true_angle["hc"][ calc_var.true_angle_calculator.reference_telescope_position ]
                    # processing_window_data.azimark_readings_theo_angle =  calc_var.true_angle_calculator.last_measured_refmark_theo_angle["hc"][ calc_var.true_angle_calculator.reference_telescope_position ] # ToDo: use an average of FL+ FR readings
                else:
                    processing_window_data.azimark_readings_true_calibrating_angle =  np.nan  # More WGS84 directions have been used
                    # processing_window_data.azimark_readings_theo_angle =  np.nan 
                
                
                # Consider setting it only on new measurements: if true_angle_calculator.new["SU"] and true_angle_calculator.new["SD"]: SU_SD_diff = 
                
                
                # raw_data[COMPNAME__PRECALC_HC_NED_ + calc_var.id] [ rec_indeces ]  =  calc_var.true_angle_calculator.h_theo_to_ned(  raw_data [ rec_indeces ][ COMPNAME__PRECORR_HC_RAD_ + calc_var.id ]  )
                # raw_data[COMPNAME__PRECALC_VC_NED_ + calc_var.id] [ rec_indeces ]  =  dim_azimuth.v_theo_to_ned( raw_data [ rec_indeces ][ COMPNAME__PRECORR_VC_RAD_ + calc_var.id ] )

                #if calc_var.refpos_corr_data is not None: removed in 2025.11.06
                #    # Applying angle correction, read from user processed refemark read data files.
                #    # 'calc_var.refpos_corr_data' must be dictinoary, containing 'SU' and 'SD' fields. Each of them must be a structured array with these fileds: 'datetime',  'hc_change_arcsec'
                #    hc_corr_SU = np.interp(  raw_data[COMPNAME__TIMESTAMP] [ rec_indeces ], calc_var.refpos_corr_data['SU']['datetime'], calc_var.refpos_corr_data['SU']['hc_change_arcsec'] )
                #    hc_corr_SD = np.interp(  raw_data[COMPNAME__TIMESTAMP] [ rec_indeces ], calc_var.refpos_corr_data['SD']['datetime'], calc_var.refpos_corr_data['SD']['hc_change_arcsec'] )
                #    raw_data[COMPNAME__PRECALC_HC_NED_ + calc_var.id] [ rec_indeces ] += -1.0*mu.arcsec_to_rad((hc_corr_SU + hc_corr_SD)/2.0)

                #    # vc_corr_SU = np.interp(  raw_data[COMPNAME__TIMESTAMP] [ rec_indeces ], calc_var.refpos_corr_data['SU']['datetime'], calc_var.refpos_corr_data['SU']['vc_change_arcsec'] )
                #    # vc_corr_SD = np.interp(  raw_data[COMPNAME__TIMESTAMP] [ rec_indeces ], calc_var.refpos_corr_data['SD']['datetime'], calc_var.refpos_corr_data['SD']['vc_change_arcsec'] )
                #    # raw_data[COMPNAME__PRECALC_VC_NED_ + calc_var.id] [ rec_indeces ] += 1.0*mu.arcsec_to_rad((vc_corr_SU - vc_corr_SD)/2.0) # SU and SD should considered in reversed sign

            # for irecord in range(n):

            #     # Processing ADS made reference positions readings

            #     if raw_data[COMPNAME__POS_ID][irecord] == ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID:
            #         # Determining, if this record is a usable SU or SD reference mark reading
            #         new_refpos_id, refpos_determination, hc_rad, vc_rad, hc_tick, vc_tick, valid_for_sensor_pos = self.process_refpos_raw_record_data(irecord, calc_var)
            #         if valid_for_sensor_pos == "SU":
            #             latest_ads_dat_file_SU_reading_index = irecord
            #         elif valid_for_sensor_pos == "SD":
            #             latest_ads_dat_file_SD_reading_index = irecord

            #     # Using first usable reference positions measured data
            #     complete_azi_meas_data = None
            #     refpos_timelabel = raw_data[COMPNAME__TIMESTAMP][irecord]
            #     for refpos_meas_data_source in calc_var.true_angle_calculator.refmark_measured_data_sources:
            #         if refpos_meas_data_source=="external-table":
            #             complete_azi_meas_data = calc_var.true_angle_calculator.get_complete_azimuth_measurement_data_from_ext_table( refpos_timelabel )
            #             latest_ads_dat_file_SU_reading_index = None
            #             latest_ads_dat_file_SD_reading_index = None
            #         elif refpos_meas_data_source=="dat-file":
            #             if (latest_ads_dat_file_SU_reading_index is not None) and (latest_ads_dat_file_SD_reading_index is not None):
            #                 complete_azi_meas_data = self.refpos_meas_indeces__2__complete_azimuth_measurement_data( [latest_ads_dat_file_SU_reading_index, latest_ads_dat_file_SD_reading_index] , calc_var )
            #         else:
            #             assert False, "Unsupported reference data source: {}".format(refpos_meas_data_source)
            #         if complete_azi_meas_data is not None: break # First valid ref-pos data is found

            #     assert complete_azi_meas_data is not None, "No reference position data for processing cycle!"
            #     assert complete_azi_meas_data.WGS84 is not None, "No WGS84 coordinate in source \"{}\"".format(refpos_meas_data_source)

            #     raw_data[COMPNAME__REFPOS_HC_THEO_SU_TICK+calcvar_index_str][irecord] = complete_azi_meas_data.SU_hc_tick
            #     raw_data[COMPNAME__REFPOS_HC_THEO_SD_TICK+calcvar_index_str][irecord] = complete_azi_meas_data.SD_hc_tick


            #     # Making THEO->NED conversion
            #     calc_var.true_angle_calculator.use_measurement_data( complete_azi_meas_data, calc_var.instr, enable_missing_data = False )

            #     raw_data[COMPNAME__REFPOS_HC_WGS84_SU+calcvar_index_str][irecord] = calc_var.true_angle_calculator.last_measured_refpos_hc_WGS84_SU
            #     raw_data[COMPNAME__REFPOS_HC_WGS84_SD+calcvar_index_str][irecord] = calc_var.true_angle_calculator.last_measured_refpos_hc_WGS84_SD

            #     raw_data[COMPNAME__PRECALC_HC_NED + calcvar_index_str][irecord] = calc_var.true_angle_calculator.h_theo_to_ned( calc_var.instr.tick_to_rad(  raw_data[COMPNAME__HC_TICK][irecord], False ) )
            #     raw_data[COMPNAME__PRECALC_VC_NED + calcvar_index_str][irecord] = dim_azimuth.v_theo_to_ned(calc_var.instr.tick_to_rad(  raw_data[COMPNAME__VC_TICK][irecord], True ))

            #     # ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_THEO_SU] = azimuth_calculator.last_measured_refpos_hc_theo_SU
            #     # ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_THEO_SD] = azimuth_calculator.last_measured_refpos_hc_theo_SD
            #     # ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_WGS84_SU]  = azimuth_calculator.last_measured_refpos_hc_WGS84_SU
            #     # ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_WGS84_SD]  = azimuth_calculator.last_measured_refpos_hc_WGS84_SD



        # self.valid_file_SU_refmark_reads = self.valid_file_SU_refmark_reads[:nof_valid_file_SU_refmarkreads]
        # self.valid_file_SD_refmark_reads = self.valid_file_SD_refmark_reads[:nof_valid_file_SD_refmarkreads]
        print("done");
        print()
        self.refpos_data_calculated_for_each_record = True



    # If called after processing cycle separation, then processing cycle borders will be displayed on plots.
    def plot_load_data(self):

        plot_T_electronics = COMPNAME__T_ELECTRONICS in  self.M_readings.dtype.names
        plot_T_room        = COMPNAME__T_ROOM        in  self.M_readings.dtype.names

        fig_nrows = 3 + (1 if (plot_T_electronics or plot_T_room) else 0) # Point positions (in two plot places, for bigger height) / M-vs-time / Temp-vs-time
        fig_ncols = 1 # or: 'self.nr_of_mag_sensors' sensor may have its own position and M data in a separate column
        fig = plt.figure(figsize=(preferences.FIGURE_SIZE[0] * 0.7, preferences.FIGURE_SIZE[1] * fig_nrows * 0.7 )  )
        # plt.rcParams['figure.constrained_layout.use'] = True
        fig.suptitle('-- ADS data file(s) contents --', fontsize=14,  fontweight="bold")
        # axes_mag_positions = []
        # axes_mag_nT        = []
        # subplot_index = 0
        # for imag in range(self.nr_of_mag_sensors):
        #     axes_mag_positions += [ plt.subplot(fig_nrows, fig_ncols, subplot_index ) ]
        #     subplot_index += 1
        # for imag in range(self.nr_of_mag_sensors):
        #     axes_mag_nT += [ plt.suplot(fig_nrows, fig_ncols, subplot_index ) ]
        #     subplot_index += 1
        axes_point_positions = plt.subplot(fig_nrows, fig_ncols, (1,2))  # Double height plot
        # axes_point_positions.title.pad = 10.0
        axes_mag_nT          = plt.subplot(fig_nrows, fig_ncols, 3)
        if plot_T_electronics or plot_T_room:
            axes_temp            = plt.subplot(fig_nrows, fig_ncols, 4 )

        # Plotting mag positions:
        # b_show_legend = False
        plt.sca(axes_point_positions)
        
        if (self.instr.axis_angle_data_source == "digital"):
            plt.plot(self.M_readings[COMPNAME__HC_TICK],   self.M_readings[COMPNAME__VC_TICK], "x",  color='green', label = 'angle positions of measurements' )
            plt.xlabel("hc positions (tick)")
            plt.ylabel("vc positions (tick)")
            axes_point_positions.set_xlim(0, self.instr.hc_axis.encoder_division )
            axes_point_positions.set_ylim(0, self.instr.vc_axis.encoder_division )
        else:
            plt.plot( np.degrees(self.M_readings[COMPNAME__MEASURED_HC_RAD]),   np.degrees(self.M_readings[COMPNAME__MEASURED_VC_RAD]), "x",  color='green', label = 'angle positions of measurements' )
            plt.xlabel("hc positions (deg)")
            plt.ylabel("vc positions (deg)")
            axes_point_positions.set_xlim(0, 360.0 )
            axes_point_positions.set_ylim(0, 360.0 )
            

        # Creating tick marks for xtick/ytick axes
        # angletick_tick_mark_pos = self.instr.hc_axis.full_range_tick_marks_pos
        # angletick_tick_mark_text = []
        # for t in angletick_tick_mark_pos:
        #    angletick_tick_mark_text += [ "{}".format(round(t)) ]
        #    # angletick_tick_mark_text += [ "{} ({:.1f}°)".format(round(t), dim_instr.tick_to_deg( t, False ) ) ]

        # plt.xticks( angletick_tick_mark_pos, angletick_tick_mark_text )
        # plt.yticks( angletick_tick_mark_pos, angletick_tick_mark_text )
        # plt.grid(visible=True, which="major", axis='both', linestyle='dotted' )

        # Filtered out mag positions:
        # ToDo: realize plotting or some of the outfiltered record!

        #if len(self.filtered_out_record_indeces)>0:
        #    # plot_utils.plot_xy(self.filtered_out_records[:,  COMPNAME__HC_TICK], self.filtered_out_records[:, COMPNAME__VC_TICK], "hc raw tick", "vc raw tick", Xaxis_range = [ 0, dim_instr.ENCODER_DIVISION ]  , Yaxis_range = [ 0, dim_instr.ENCODER_DIVISION ],  color='red', label='out-filtered points', close_plot = False  )
        #    b_show_legend = True
        #    plt.plot(self.filtered_out_records[:,  COMPNAME__HC_TICK], self.filtered_out_records[:, COMPNAME__VC_TICK], "x",  color='red', label='removed measurements positions')

        # Plotting reference positions:
        #for refpos_index, refpos_id in enumerate(self.refpos_id_list):
        #    param1_index = self.meas_type_param1_str_list.index(refpos_id)
        #    refpos_points = raw_data[ raw_data[COMPNAME__PROGRAM_PARAM1_IDX] == param1_index ]
        #    refpos_label = "Ref.pos. \"{}\" points".format(refpos_id)
        #    refpos_color = plot_utils.rainbow_color_get(refpos_index, len(self.refpos_id_list)  )
        #    plt.plot( refpos_points[COMPNAME__HC_TICK], refpos_points[COMPNAME__VC_TICK], "o", color = refpos_color, label = refpos_label )
        #    b_show_legend = True

        # Plotting non filtered-out, but NAN mag positions (typically reference positions reading positions):
        # assert self.nr_of_mag_sensors == 1
        #for imag in range(self.nr_of_mag_sensors):
        #    nan_mag_points = raw_data[ np.isnan( raw_data[COMPNAME__MAG] ) ]
        #    nan_mag_points = nan_mag_points[ nan_mag_points[COMPNAME__POS_ID] != ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID ]
        #    if len(nan_mag_points)>0:
        #        b_show_legend = True
        #        plt.plot(  nan_mag_points[COMPNAME__HC_TICK], nan_mag_points[COMPNAME__VC_TICK], ".",  color='black', label="Positions of MAG-{} 'NAN' points".format(imag+1) )
        # if b_show_legend: plt.legend(framealpha=1, frameon=True, bbox_to_anchor=(1.1, 1.05))

        # Plotting mag nT valus:
        plt.sca(axes_mag_nT)
        mags_nT_array = []
        for imag in range(self.nr_of_mag_sensors):
            mags_nT_array += [ self.M_readings[COMPNAME__MAG] ]
        plot_utils.plot_datetime_y( self.M_readings[COMPNAME__TIMESTAMP],  mags_nT_array,  "M [nT]", fmt=['<','>','o'],  label = self.mag_sensor_names, b_close_plot = False )
        # if self.nr_of_mag_sensors>1: plt.legend(framealpha=1, frameon=True, bbox_to_anchor=(1.1, 1.05))

        # Showing processing windows limits (of the first calc var) with y grid lines (if 'prepare_processing_cycles()' have been already called)
        main_calcvar_index = 0
        for calc_var_index, calc_var in enumerate(self.calc_vars):
            calc_var.processing_windows_plot_spearators = {}  # Processing windows plot separators for different X axes
            limiting_file_record_indeces = []
            for processing_window in calc_var.processing_window_data:
                limiting_file_record_indeces += [ processing_window.M_record_indeces[0], processing_window.M_record_indeces[-1] ]   # Collecting the start file-record indeces of the first calc. var
            calc_var.processing_windows_plot_spearators["record-index"] = limiting_file_record_indeces
            calc_var.processing_windows_plot_spearators["date"] = self.M_readings[COMPNAME__TIMESTAMP][limiting_file_record_indeces]
            calc_var.processing_windows_plot_spearators["date"] = mu.exceltime_to_matplotlib_time( calc_var.processing_windows_plot_spearators["date"] )
            # plt.xticks(processing_cycle_date_plot_spearator, minor = True) not working!
            if calc_var_index == main_calcvar_index:
                axes_mag_nT.xaxis.set_minor_locator( matplotlib.ticker.FixedLocator(calc_var.processing_windows_plot_spearators["date"]) )
                plt.grid(visible=True, which="minor", axis='x', linestyle='dotted' )

        if (plot_T_electronics or plot_T_room):
            plt.sca(axes_temp)
            # Breaking temperature lines at the starts of each cycle

            plotted_temperature_series = [ ]
            plotted_temperature_names  = [ ]
            plotted_temperature_colors = [ ]
            
            # Original situation until 2026.09: separating temperature plots by adding NULL values at file absmeas windows
            #M_start_record_indeces = []  # Originally 'self.file_absmeas_window_data['M_start_record_index']' have been used, but 
            #for file_asmeas_window in self.file_absmeas_window_data:
            #    M_start_record_indeces += [ file_asmeas_window[DATAWINDOW_COMPNAME__M_RECORD_INDECES] ]
            
            if plot_T_electronics:
                # T_el__separated_at_cycle_starts = np.copy( self.M_readings[COMPNAME__T_ELECTRONICS] )
                # T_el__separated_at_cycle_starts[  M_start_record_indeces  ] = np.nan
                # plotted_temperature_series += [ T_el__separated_at_cycle_starts  ]
                plotted_temperature_series += [ self.M_readings[COMPNAME__T_ELECTRONICS]  ]
                plotted_temperature_names  += [ "T-electronics °C)" ]
                plotted_temperature_colors += [ 'red' ]
            if plot_T_room:
                plotted_temperature_series += [ self.M_readings[COMPNAME__T_ROOM]  ]
                plotted_temperature_names  += [ "T-room °C)" ]
                plotted_temperature_colors += [ 'blue' ]

            plot_utils.plot_datetime_y(self.M_readings[COMPNAME__TIMESTAMP],  plotted_temperature_series,  plotted_temperature_names, b_close_plot = False, color = plotted_temperature_colors )

            axes_temp.xaxis.set_minor_locator( matplotlib.ticker.FixedLocator(self.calc_vars[main_calcvar_index].processing_windows_plot_spearators["date"]) )
            plt.grid(visible=True, which="minor", axis='x', linestyle='dotted' )

        fig.tight_layout()
        plt.show()

        # Plotting reference position data: it is done in 'check_load_data()', as it is a more difficult task to select valid reference points
        # ads_refpos_data = raw_data[  raw_data[COMPNAME__POS_ID]<=ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID  ]
        #if len(ads_refpos_data)>0:
        #     plot_utils.plot_datetime_y(ads_refpos_data[COMPNAME_DATETIME], [ hc_theo_SU_deg ], "hc-theo-SU [°]" )

        # -----------------------------------------------------------------------------

    def check_load_data(self, magfield):
        assert  isinstance(magfield, gu.magfield_class)
        # print("Checking data..." , end  = "")
        # M_meas_array = np.empty( (len(self.calc_vars),  len(raw_data))  )
        # M_calc_array = np.empty( (len(self.calc_vars),  len(raw_data))  )

        # for icalvcar, calc_var in enumerate(self.calc_vars):


        #     for i,file_record in enumerate(raw_data):

        #         file_record = raw_data[i]
        #         vc_theo_rad___from_raw_tick = dim_instr.tick_to_rad( file_record[COMPNAME__VC_TICK], True, calc_var.instr.vc_axis0)
        #         hc_theo_rad___from_raw_tick = dim_instr.tick_to_rad( file_record[COMPNAME__HC_TICK], False, calc_var.instr.vc_axis0)
        #         # vc_deg_diff = np.degrees(vc_theo_rad___from_raw_tick) - file_record[COMPNAME__VC_DEG]
        #         # if abs(vc_deg_diff) > 0.5:
        #         #     cu.warning("vc angle calculation error for record {}: vc angle from file: {:.1f}°  vc angle calculated: {:.1f}°".format(i,file_record[self.compindex_file_vc_deg], np.degrees(vc_theo_rad___from_raw_tick) ),37)
        #         assert abs(file_record[COMPNAME__HC_DEG] -   np.degrees(hc_theo_rad___from_raw_tick) ) < 0.5, "Record index: {:}  hc_deg_file:{:.3f}  hc_deg_from_tick:{:.3f}".format(i, file_record[COMPNAME__HC_DEG], np.degrees(hc_theo_rad___from_raw_tick))

        #         if not self.refpos_data_calculated_for_each_record: continue

        #         if file_record[COMPNAME__POS_ID] >= ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID:

        #             # Normal mag measurement
        #             if refpos_index<0: continue # Not yet having refmark reading

        #             # print("hc-tick:{}  vc-tick:{}  MIRA-SU:{:.3f}°({})  MIRA-SD:{:.3f}°({})".format(file_record[COMPNAME__HC_TICK],file_record[COMPNAME__VC_TICK], np.degrees( used_reference_positions[refpos_index].SU_hc ), debug_last_SU_hc_tick, np.degrees( used_reference_positions[refpos_index].SD_hc ), debug_last_SD_hc_tick  ) )


        #             M_measured = file_record[ COMPNAME__MAG ]
        #             assert (magfield.F > 10000) and (magfield.F<60000)
        #             M_calc = amdp.dim_instrument_model(magfield.F, vc_NED, hc_NED, magfield.I, magfield.D, instr.mag_sensors[0].epsilon(), instr.mag_sensors[0].delta(), instr.mag_sensors[0].offset, instr.mag_sensors[0].scale)
        #             M_diff = M_measured - M_calc

        #             # print("{:5.1f}°  {:5.1f}°:  {:9.2f}  {:9.2f}".format( np.degrees(hc_NED), np.degrees(vc_theo_rad___from_raw_tick), M_measured, M_calc) )
        #             M_meas_array[icalvcar,nof_mag_measurements_with_refpos_data] = M_measured
        #             M_calc_array[icalvcar,nof_mag_measurements_with_refpos_data] = M_calc
        #             # NED_points[nof_mag_measurements_with_refpos_data]['hc_NED'] = np.degrees(hc_NED)
        #             # NED_points[nof_mag_measurements_with_refpos_data]['vc_NED'] = np.degrees(vc_NED)
        #             # NED_points[nof_mag_measurements_with_refpos_data]['pos_ID'] = file_record[COMPNAME__POS_ID]
        #             # NED_points[nof_mag_measurements_with_refpos_data]['M_meas'] = M_measured
        #             # NED_points[nof_mag_measurements_with_refpos_data]['M_calc'] = M_calc


        #             nof_mag_measurements_with_refpos_data+=1

        #     if (self.associate_with_MIRA_data) and (nof_mag_measurements_with_refpos_data>0) :
        #         plot_utils.plot_xy(M_meas_array[icalvcar,:nof_mag_measurements_with_refpos_data], M_calc_array[icalvcar,:nof_mag_measurements_with_refpos_data] , "M-measured (nT)", "M-precalculated (nT)")

        #         plot_utils.plot_datetime_y(valid_file_SU_refmark_reads["datetime"][:nof_valid_SU_refmarkreads], [  np.degrees( valid_file_SU_refmark_reads["hc"][:nof_valid_SU_refmarkreads] ) ], "SU reference mark hc angle minus WGS84", b_close_plot = False )
        #         plot_utils.plot_datetime_y(valid_file_SD_refmark_reads["datetime"][:nof_valid_SD_refmarkreads], [  (np.degrees( valid_file_SD_refmark_reads["hc"][:nof_valid_SD_refmarkreads] ) + 180.0 ) % 360.0 ], "SD reference mark hc angle minus WGS84 +180°")
        #         print("Please check that M-precalcualted vs. M-measured plot's rough linearity for calc-var nr {}, and the reference mark angles!".format(icalvcar+1))
        #         mu.wait_for_enter()

        #         # mu.savetxt('d:/NED-points-{}.dat'.format(icalvcar), NED_points[:nof_mag_measurements_with_refpos_data])



    # Returns all load records, with associated reference mark data, for one calc_variation.
    # It returns only magnetic data (reference mark records are excluded)
    # It is used only in "all-in-one" measurement type.
    def get_all_processing_cycle_mag_data_together(self, calc_var_index, parallel_mag_data_index):
        assert self.refpos_data_calculated_for_each_record
        assert self.number_of_processing_windows > 0
        assert parallel_mag_data_index == 0

        mag_only_record_indeces = []
        for i_processing_cycle_index in range(self.number_of_processing_windows):
            mag_only_record_indeces += self.processing_window_data[calc_var_index][i_processing_cycle_index].mag_only_record_indeces
        parallel_mag_data = self.load_main_obs_data_set_data[mag_only_record_indeces]

        assert len ( self.load_main_obs_data_set_data_missing_record_indeces[mag_only_record_indeces] ) == 0

        parallel_mag_data_with_window = None if self.load_main_obs_data_set_data_with_window is None else self.load_main_obs_data_set_data_with_window[parallel_mag_data_index][mag_only_record_indeces]
        return raw_data[mag_only_record_indeces],parallel_mag_data, parallel_mag_data_with_window


    #     n = len(raw_data)

    #     use_ext_table_only = True
    #     for refpos_meas_data_source in calc_var.true_angle_calculator.refmark_measured_data_sources:
    #         if refpos_meas_data_source!="external-table":
    #             use_ext_table_only = False

    #     # last_refpos_reading_SU = None
    #     # last_refpos_reading_SD = None

    #     ads_raw_data_magdata_only_records = raw_data[COMPNAME__POS_ID]>ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID

    #     ads_mag_data_with_refpos = np.copy( raw_data[  ads_raw_data_magdata_only_records ]  )

    #     fill_val = np.empty(len(ads_mag_data_with_refpos))
    #     fill_val[:] = np.nan

    #     ads_mag_data_with_refpos = rcf.append_fields( ads_mag_data_with_refpos, [ COMPNAME__REFPOS_HC_THEO_SU, COMPNAME__REFPOS_HC_THEO_SD, COMPNAME__REFPOS_HC_WGS84_SU, COMPNAME__REFPOS_HC_WGS84_SD], [fill_val,fill_val,fill_val,fill_val], usemask=False )

    #     ads_mag_data_with_refpos__index = 0

    #     azimuth_calculator = dim_azimuth.default_azimuth_def  # Use calc_var.azimuth_calcualtor instead

    #     # assert calc_var.true_angle_calculator.refmark_measured_data_sources[0] == "dat-file", "This function works only for dat-file reference data source"

    #     for i in range(0,n):

    #         if raw_data[COMPNAME__POS_ID][i] == ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID:

    #             assert not use_ext_table_only
    #             # It is a reference position measurement
    #             assert np.isnan(raw_data[COMPNAME__MAG][i])
    #             complete_azi_meas_data = self.refpos_meas_indeces__2__complete_azimuth_measurement_data([i], calc_var)
    #             azimuth_calculator.use_measurement_data(complete_azi_meas_data, instr, enable_missing_data = True)

    #             #if np.isfinite(complete_azi_meas_data.SU_hc):
    #             #    assert np.isfinite(complete_azi_meas_data.SU_vc)
    #             #    assert np.isnan(complete_azi_meas_data.SD_hc) # A point can not be SU and SD at the same time
    #             #    last_refpos_reading_SU =  complete_azi_meas_data
    #             #    azimuth_determinator.use_
    #             #elif np.isfinite(complete_azi_meas_data.SD_hc):
    #             #    assert np.isfinite(complete_azi_meas_data.SD_vc)
    #             #    assert np.isnan(complete_azi_meas_data.SU_hc) # A point can not be SU and SD at the same time
    #             #    last_refpos_reading_SD =  complete_azi_meas_data

    #         else:

    #             if use_ext_table_only or np.isnan(azimuth_calculator.last_measured_refpos_hc_theo_SU):
    #                 timelabel = raw_data[COMPNAME__TIMESTAMP][i]
    #                 complete_azi_meas_data = calc_var.true_angle_calculator.get_complete_azimuth_measurement_data_from_ext_table( [ timelabel-0.001, timelabel+0.001 ] )
    #                 assert complete_azi_meas_data is not None, "No reference position data for processing cycle!"
    #                 azimuth_calculator.use_measurement_data(complete_azi_meas_data, instr)

    #             print(azimuth_calculator.last_measured_refpos_hc_theo_SU)
    #             assert not np.isnan(raw_data[COMPNAME__MAG][i])
    #             assert not np.isnan(azimuth_calculator.last_measured_refpos_hc_theo_SU)
    #             assert not np.isnan(azimuth_calculator.last_measured_refpos_hc_theo_SD)
    #             assert not np.isnan(azimuth_calculator.last_measured_refpos_hc_WGS84_SU)
    #             assert not np.isnan(azimuth_calculator.last_measured_refpos_hc_WGS84_SD)

    #             ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_THEO_SU] = azimuth_calculator.last_measured_refpos_hc_theo_SU
    #             ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_THEO_SD] = azimuth_calculator.last_measured_refpos_hc_theo_SD
    #             ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_WGS84_SU]  = azimuth_calculator.last_measured_refpos_hc_WGS84_SU
    #             ads_mag_data_with_refpos[ads_mag_data_with_refpos__index][COMPNAME__REFPOS_HC_WGS84_SD]  = azimuth_calculator.last_measured_refpos_hc_WGS84_SD

    #             ads_mag_data_with_refpos__index += 1


    #     assert ads_mag_data_with_refpos__index == len(ads_mag_data_with_refpos)
    #     assert not np.any( np.isnan(ads_mag_data_with_refpos[COMPNAME__REFPOS_HC_THEO_SU])   )
    #     assert not np.any( np.isnan(ads_mag_data_with_refpos[COMPNAME__REFPOS_HC_THEO_SD])   )


    #     ref_data = self.load_main_obs_data_set_data[parallel_mag_data_index]
    #     # Using only ref_data for the mag recordings
    #     ref_data = ref_data[ ads_raw_data_magdata_only_records  ]
    #     assert len(ref_data) == len(ads_mag_data_with_refpos)

    #     return ads_mag_data_with_refpos, ref_data


    # This function loads magnetic data for the time of the DIM measureemt points.
    # We calcualte spot baseline values for these measurement
    # It returns and index, representing the stored parallel magnetic data
    def load_parallel_mag_recordings(self, mag_obs_dataset, save_filename = None, time_shift_sec = None, use_subsecond_interpol = False ):

        assert isinstance(mag_obs_dataset, geomag_obs_dataset.mag_obs_dataset)
        
        #if mag_obs_dataset in self.load_parallel_mag_obs_datasets:
        #    
        #    # Using the already load magnetic/aux data
        #    already_load_mag_data_index = self.load_parallel_mag_obs_datasets.index(mag_obs_dataset) # This data def is already loaded
        #    # Checking if there is a new component request:
        #    new_component_request = False
        #    for calc_comp in calc_comps:
        #        if not  calc_comp in self.load_main_obs_data_set_data[already_load_mag_data_index].dtype.names:
        #            new_component_request = True
        #            break
        #    if not new_component_request:
        #        return already_load_mag_data_index

        # -------------------------
        # --- Reading in data  ----
        # -------------------------
        print("Loading observatory variometer data, optionally with scalar data...")
        parallel_mag_data, noncomplete_record_indeces =  mag_obs_dataset.load_data(True, False, self.M_readings[COMPNAME__TIMESTAMP], meas_pos = self.meas_pos, in_case_of_missing_error = "none", save_filename = save_filename, time_shift_sec = time_shift_sec, force_using_time_interpolation = use_subsecond_interpol, info_lines_prefix = "  " )
        # Note: 'noncomplete_record_indeces' can be only-variometer, onlyF, or variometer+F 
        print("...done")
        
        
        parallel_mag_data_with_valid_records = np.delete(parallel_mag_data, noncomplete_record_indeces)
        #for comp in mag_obs_dataset.usage_frame +  ("F" if mag_obs_dataset.has_F() else ""):
        #    assert not np.any( np.isnan(parallel_mag_data_with_valid_records[comp]  ) ) 
        
        self.load_main_obs_data_set_data =  parallel_mag_data 
        self.load_main_obs_data_set_data_missing_record_indeces = noncomplete_record_indeces
        # self.load_parallel_mag_obs_datasets += [ mag_obs_dataset ]
        print()
        
        
        F_data_load_printed = False
        from dim_data_process import F_data_source_enum
        for calc_var in self.calc_vars:
            for processing_window in calc_var.processing_window_data:
                missing_indeces_from_window = list(set(noncomplete_record_indeces) & set( processing_window.M_record_indeces ))
                processing_window.missing_variometer_data = len(missing_indeces_from_window) > 0
            if calc_var.F_source == F_data_source_enum.default:
                calc_var.F_source = F_data_source_enum.main_obs_dataset if mag_obs_dataset.has_F_data() else F_data_source_enum.sites_default_F
                
            match calc_var.F_source:
                
                case F_data_source_enum.main_obs_dataset:
                    assert mag_obs_dataset.has_F_data(), "Main observatory data set does not have F component"
                    calc_var.F_abs_pillar = parallel_mag_data['F'] + self.calc_params_common.F_pillar_diff_scalar
                    for processing_window in calc_var.processing_window_data:
                        processing_window.missing_F_data = np.any( np.isnan(calc_var.F_abs_pillar[processing_window.M_record_indeces]) )
                        # Note: we should set mising_F_data for 'processing_window.missing_variometer_data', but it is possible
                        # to have a data set, where different position missing F data and variometer data are. 
                    
                case F_data_source_enum.DIM_data_file:
                    assert "F" in self.M_readings.dtype.names
                    calc_var.F_abs_pillar = self.M_readings["F"] + self.calc_params_common.F_pillar_diff_scalar
                    for processing_window in calc_var.processing_window_data:
                        processing_window.missing_F_data = np.any( np.isnan(calc_var.F_abs_pillar[processing_window.M_record_indeces]) )
                    
                case F_data_source_enum.sites_default_F:
                    calc_var.F_abs_pillar = np.empty( len(self.M_readings) )
                    calc_var.F_abs_pillar[:] = self.meas_pos.F_typical
                    for processing_window in calc_var.processing_window_data:
                        processing_window.missing_F_data = False
                    
                case F_data_source_enum.separate_obs_dataset:
                    if not F_data_load_printed:
                        print("Loading F data for calc.var {}...".format(calc_var.name))
                        F_data_load_printed = True
                        
                    read_data, F_abs_pillar_missing_record_indeces = calc_var.F_source_dataset.load_data(True, False, self.M_readings[COMPNAME__TIMESTAMP], meas_pos = self.meas_pos, report_for_no_data = np.nan, in_case_of_missing_error = "none", time_shift_sec = self.calc_params_common.F_time_shift_sec, force_using_time_interpolation = False)
                    
                    calc_var.F_abs_pillar = read_data["F"]  + self.calc_params_common.F_pillar_diff_scalar

                    for processing_window in calc_var.processing_window_data:
                        processing_window.missing_F_data = len ( list(set(F_abs_pillar_missing_record_indeces) & set( processing_window.M_record_indeces )) ) >0
                        assert processing_window.missing_F_data == np.any( np.isnan(calc_var.F_abs_pillar[processing_window.M_record_indeces]) )
                    #    F_abs_pillar, F_abs_pillar_within_window, F_abs_pillar_missing_record_indeces = DIM_data_loader.processing_cycle_parallel_mag_data_get(processing_cycle_nr, F_data_parallel_mag_data_index, calc_params_index)
                    #    F_abs_pillar = F_abs_pillar['F']
                    #    F_abs_pillar = F_abs_pillar + calc_params_common.F_pillar_diff_scalar
                    
                case F_data_source_enum.calc_from_variometer:
                    calc_var.F_abs_pillar = mag_obs_dataset.calc_F(parallel_mag_data) 
                    for processing_window in calc_var.processing_window_data:
                        processing_window.missing_F_data = processing_window.missing_variometer_data
                case _:
                    assert False, "Unsupported F-source definition: {}".format(calc_var.F_source)
                  
        if F_data_load_printed: 
            print("...done")
            print()

                
        
        if self.calc_params_common.variometer_slopecalc_halfwindow_sec == 0:
            self.load_main_obs_data_set_data_with_window +=[ None ]
        else:
            print("----- Loading variometer data with extended window {} -----".format(self.data_descr))
            assert False, "Check this part"

            mag_data_in_window = np.empty( ( len(raw_data), self.calc_params_common.variometer_slopecalc_halfwindow_sec*2+1 ), dtype = parallel_mag_data.dtype   )
            for i,timelabel in enumerate(raw_data[COMPNAME__TIMESTAMP]):
                timelabels_in_window = np.linspace(  timelabel - self.calc_params_common.variometer_slopecalc_halfwindow_sec/86400.0, timelabel + calc_params_common.variometer_slopecalc_halfwindow_sec/86400.0, num = 2*calc_params_common.variometer_slopecalc_halfwindow_sec + 1  )
                mag_data_in_window_rec, missing_record_indeces = mag_data_set_def.load_mag_data("AUTO_FOR_SECOND" , timelabels_in_window, data_descr = data_descr, in_case_of_missing_error = "drop-error",   )
                mag_data_in_window[i,:] = mag_data_in_window_rec
                # for compname in qqq.dtype.names:
                #    mag_data_in_window[i,:][compname] = qqq[compname]

            # This check is removed in 2025.12, as fails if parallel_mag_data is read by subsecond time interpolation
            # for i in range( len(raw_data) ):
            #     for compname in parallel_mag_data.dtype.names:
            #         if not np.isnan(parallel_mag_data[i][compname]):
            #             assert parallel_mag_data[i][compname] == mag_data_in_window[i,calc_params_common.variometer_slopecalc_halfwindow_sec][compname], "{}: {} <> {}".format(compname, parallel_mag_data[i][compname], mag_data_in_window[i,calc_params_common.variometer_slopecalc_halfwindow_sec][compname])

            self.load_main_obs_data_set_data_with_window +=[ mag_data_in_window ]


        # return len(self.load_main_obs_data_set_data)-1




    def plot_mag_meas_slopes(self, baselines, variometer_data_def):

        # 2024.10: added plots for checking of magnetometer delays.
        # Plotting average slope for different reference positions (for each mesurement type):

        return # QQQ

        if variometer_data_def is not None:

            assert False
            # Calculating slopes from baseline corrected variometer data for DIM instruments measurement position
            # (for each calc. var)
            assert baselines is not None
            assert len(baselines) == self.nr_of_calcvars

            field_slope = np.empty(  (self.nr_of_calcvars, len(raw_data) )  )
            for icalcvar in range( len(self.calc_vars) ):
                field_slope[icalcvar,:] = np.nan
                for record_index,raw_data_record in enumerate(raw_data):
                    if raw_data_record[COMPNAME__POS_ID]<0: continue
                    assert not np.isnan(raw_data_record[COMPNAME__MAG])
                    variometer_meas_data, missing_record_indeces =  variometer_data_def.load_mag_data("MEAN_WITHIN_PERIOD", raw_data_record[COMPNAME__TIMESTAMP], data_descr = "var. data for DIM slope", in_case_of_missing_error = "none")
                    # _F, _I, _D = gu.variometer_to_abs_data(variometer_meas_data, baselines[icalvar], variometer_data_def)
                    # field_slope[icalcvar,record_index] = gu

        assert  len(self.pos_id_list)  == len( self.program_variation_str_list )
        for program_type_str,  pos_id_list  in zip(self.program_variation_str_list, self.pos_id_list):

            if pos_id_list[0] == ADS_DAT_FILE__REFPOS_MEASUREMENT_RECORDS_POS_ID:
                # Reference position readings have one pos_id
                assert program_type_str in [MEAS_TYPE_STR__MIRA_READING]
                assert len(pos_id_list) == 1
                continue
            assert not program_type_str in [MEAS_TYPE_STR__MIRA_READING]

            # Calculating average slope separateky for each position
            nr_of_pos_ids = len(pos_id_list)
            slope_avg = np.empty( nr_of_pos_ids )
            slope_stddev = np.empty( nr_of_pos_ids )
            slope_avg[:] = np.nan
            slope_stddev[:] = np.nan
            slope_raw_data = []
            NM_meas_part_separators = []  # Vertical dotted lines, separating different parts of a null method measurement

            for pos_id_list_index, pos_id in enumerate(pos_id_list):
                # 'pos_id_list_index'  is also variation index
                assert False
                # variation_file_cycles = self.file_absmeas_window_data['variation_idx'] == pos_id_list_index
                boolean_mask = (raw_data[COMPNAME__POS_ID]  ==  pos_id) & (raw_data[COMPNAME__PROGRAM_MAINTYPE_STR] == program_type_str)
                given_pos_id_within_meas_type__raw_data = raw_data[boolean_mask]
                assert len(given_pos_id_within_meas_type__raw_data)>0
                assert np.all( given_pos_id_within_meas_type__raw_data[COMPNAME__POS_ID] == pos_id)
                slope_raw_data += [  given_pos_id_within_meas_type__raw_data[COMPNAME__MAG_SLOPE] ]
                assert len(slope_raw_data[-1])>0
                slope_avg[pos_id_index] = np.mean(  slope_raw_data[-1] )
                slope_stddev[pos_id_index] = np.std(slope_raw_data[-1])

            if (program_type_str in [MEAS_TYPE_STR__NULLMETH_1,MEAS_TYPE_STR__NULLMETH_2]) and  ((nr_of_pos_ids % 8) in [0,1]):  # must have 8 positions, withouth or with scale check
                nr_of_scan_points = nr_of_pos_ids // 8
                for i in range(1,8):
                    end_of_scan_pos_id = pos_id_list[i*nr_of_scan_points-1]
                    start_of_scan_pos_id = pos_id_list[i*nr_of_scan_points]
                    NM_meas_part_separators += [ (end_of_scan_pos_id + start_of_scan_pos_id)/2.0 ]

            plot_utils.plot_xy( pos_id_list, slope_avg, "{} pos-ID".format(program_type_str), "M_slope avg (nT/sec)" , close_plot = False, plot_width = preferences.FIGURE_SIZE[0], vertical_separators = NM_meas_part_separators, horizontal_separators = [ -0.1, 0.0, +0.1 ] )
            # Plotting also error bars:
            # plt.errorbar(pos_id_list, slope_avg, fmt='.', yerr = slope_stddev )

            # Writing pos id above points
            for i, pos_id in enumerate(pos_id_list):
                # Plotting all rad slope data:
                # plt.plot( [pos_id] * len(slope_raw_data[i]), slope_raw_data[i],  '.')
                plt.gca().annotate( "{:}".format(pos_id)  , (pos_id, slope_avg[i]))
            plt.show()


    def get_aux_comp_names(self, only_temperatures):
        result = []
        if COMPNAME__T_ELECTRONICS in self.M_readings.dtype.names:
            result += ['T_electronics'] 
        if COMPNAME__T_ROOM in self.M_readings.dtype.names:
            result += ['T_room'] 
        return result
        

    # ~processing_cycle_data_get~
    # Selects a processing cycle data for processing, for a calc_var.
    def processing_cycle_data_select(self, instr, processing_window_index, calc_var_index, mag_precorr_functions = [], mag_indeces = [0] ):
        assert not self.all_cycles_together_mode
        assert self.number_of_processing_windows>0
        assert processing_window_index < self.number_of_processing_windows
        assert len(mag_indeces)==1 and (mag_indeces[0] == 0)

        processing_window_data = self.calc_vars[calc_var_index].processing_window_data[processing_window_index]
        M_record_indeces = processing_window_data.M_record_indeces

        assert(len(M_record_indeces)>2)

        mag_data_in_window = self.M_readings[ M_record_indeces ].copy()

        method_info =  self.get_program_type_info_str(mag_data_in_window)

        datetime    = mag_data_in_window[COMPNAME__TIMESTAMP].copy()

        temperatures = dict()
        if COMPNAME__T_ELECTRONICS in self.M_readings.dtype.names:
            temp_electronics =  mag_data_in_window[COMPNAME__T_ELECTRONICS]
            temperatures['T_electronics'] = temp_electronics
        else:
            temp_electronics = None

        if COMPNAME__T_ROOM in self.M_readings.dtype.names:
            temp_room =  mag_data_in_window[COMPNAME__T_ROOM]
            temperatures['T_room'] = temp_room

        if self.instr.axis_angle_data_source =="digital":
            hc_raw_tick, vc_raw_tick = mag_data_in_window[COMPNAME__HC_TICK], mag_data_in_window[COMPNAME__VC_TICK]
            meas_time   = mag_data_in_window[COMPNAME__MOVETOPOINT_TIME]
        else:
            hc_raw_tick, vc_raw_tick =  None, None
            meas_time = None

        # hc_theo_rad = d_for_one_cycle[COMPNAME__PRECORR_HC_RAD_ + self.calc_vars[calc_var_index].id]  #  self.calc_vars[calc_var_index].correct_raw_angles( "hc",  dim_instr.tick_to_rad(hc_raw_tick, False) )
        # vc_theo_rad = d_for_one_cycle[COMPNAME__PRECORR_VC_RAD_ + self.calc_vars[calc_var_index].id]   # self.calc_vars[calc_var_index].correct_raw_angles( "vc",  dim_instr.tick_to_rad(vc_raw_tick, True, instr.vc0) )
        hc_theo_rad = self.calc_vars[calc_var_index].M_reading_theo_angles["hc"][M_record_indeces].copy()
        vc_theo_rad = self.calc_vars[calc_var_index].M_reading_theo_angles["vc"][M_record_indeces].copy()
        
        assert len(hc_theo_rad) == len(mag_data_in_window) == len(vc_theo_rad)

        hc_theo_deg = np.degrees(hc_theo_rad)
        vc_theo_deg = np.degrees(vc_theo_rad)

        if len(mag_precorr_functions)>0:
            assert False
            mag = np.copy(d_for_one_cycle[COMPNAME__MAG])
            for mag_correction_function in mag_precorr_functions:
                if isinstance(mag_correction_function, dict):
                    # nT correction table by pos-ID:
                    for i in range(len(d_for_one_cycle)):
                        if d_for_one_cycle[COMPNAME__POS_ID][i] in mag_correction_function:
                            print(mag_correction_function [ d_for_one_cycle[COMPNAME__POS_ID][i]  ])
                            mag[i] += mag_correction_function [ d_for_one_cycle[COMPNAME__POS_ID][i]  ]
                else:
                    assert  not mag_correction_function.is_angle_correction()
                    # for i in range(0,len(mag)): mag[i] = mag_correction_function.calc(mag[i])
                    mag = mag_correction_function.calc(mag, temp_electronics)
        else:
            mag = mag_data_in_window[COMPNAME__MAG]

        F = self.calc_vars[calc_var_index].F_abs_pillar[M_record_indeces]
        

        if processing_window_data.complete_azi_meas_data is not None:
            hc_NED =  processing_window_data.M_reading_NED_angles["hc"].copy()
            vc_NED =  processing_window_data.M_reading_NED_angles["vc"].copy()
            assert len(hc_NED) == len(mag_data_in_window) 
        else:
            # This cycle has no refpos data (or refpos data is not complete)
            hc_NED = None
            vc_NED = None

        return processing_window_data, method_info, datetime, [mag], F, hc_raw_tick, vc_raw_tick, hc_theo_deg, vc_theo_deg, hc_theo_rad, vc_theo_rad, hc_NED, vc_NED, meas_time, temperatures, mag_data_in_window


    # ~print_new_processing_cycle_begin~
    # This is called after 'processing_cycle_data_select()' was called and cycle data was preprocessed,
    # as the 'first_line_extra_str' and 'second_line' generally can be formatted after cycle data preprocessing!
    def print_processing_windows_data_info(self, processing_window_data, print_in_frame, first_line_extra_str, calc_params_info = "", fig = None, nof_plots_added_later  = -1):

        # file_record_indeces = self.processing_cycles_file_record_indeces_without_refpos_records[processing_window_index][calc_var_index]
        #cycle_start_time_str = starting_record_str[0]
        #cycle_id_str = starting_record_str[2]

        togehter_with_cycle_indeces = "" if self.processing_cycles__concatenated_cycles == 1 else "-{}".format(processing_window_data.cycle_index + self.processing_cycles__concatenated_cycles-1)

        datetime_str = "{:.6f} {:}".format(processing_window_data.first_M_timelabel, mu.exceltime_2_datetimestr(processing_window_data.first_M_timelabel) )

        if first_line_extra_str != "": first_line_extra_str = first_line_extra_str + " "
        main_line = calc_params_info
        if len(main_line)>0: main_line += " "
        main_line += "WIN:{:d}{} FILE-ID:\"{}\" START:{} REC-CNT:{} {}".format(processing_window_data.window_index + 1, togehter_with_cycle_indeces,   processing_window_data.file_id_str, datetime_str, len(processing_window_data.M_record_indeces), first_line_extra_str)
        # Other infos: "NAME:{} ".format(processing_window_data.name.ljust(self.longest_proessing_window_name))

        if print_in_frame:
            cu.print_text_in_frame(main_line,calc_params_info)
        else:
            print(main_line)


        # 'fig' has to be created by the caller: it can contain one or more plots - here we show cycle information
        if fig !=  None:

            assert nof_plots_added_later >= 0

            plot_ncols =  2
            if len(processing_window_data.temperatures)>0: plot_ncols +=1

            # Caller should prepared the figure like 'fig = plt.figure(figsize=(preferences.FIGURE_SIZE[0],3 + nof_plots_added_later * 3))'

            # cycle_file_id_info, method_info, datetime, mags, hc_raw_tick, vc_raw_tick, hc_deg, vc_deg, hc, vc, hc_NED, vc_NED, meas_time, temp = self.processing_cycle_data_get(instr,processing_window_index,calc_var,calc_var_index)

            fig.suptitle('Cycle-file-id:{}  method:{:} start-date:{:}'.format(processing_window_data.file_id_str, processing_window_data.method_info[1], mu.exceltime_2_datetime(processing_window_data.start_datetime).strftime("%Y.%m.%d.%H:%M:%S") ), fontsize=14)
            plt.subplot(1 + nof_plots_added_later, plot_ncols, 1)
            plot_utils.plot_xy(processing_window_data.hc_deg, processing_window_data.vc_deg, 'hc [°]', 'vc [°]', 'o', color='green', close_plot = False)
            plt.subplot(1 + nof_plots_added_later, plot_ncols, 2)
            plot_utils.plot_xy(range(processing_window_data.nr_of_records), processing_window_data.meas_time, 'meas. point', 'meas_time [sec]', 'o', color='green', close_plot = False)
            if len(processing_window_data.temperatures)>0:
                assert isinstance(processing_window_data.temperatures, (dict) )
                plt.subplot(1 + nof_plots_added_later, plot_ncols, 3)
                plotted_temperature_series = []
                plotted_temperature_series_names = []
                for temp_ch_name in processing_window_data.temperatures:
                    plotted_temperature_series += [processing_window_data.temperatures[temp_ch_name] ]
                    plotted_temperature_series_names += [temp_ch_name]

                plot_utils.plot_xy(range(processing_window_data.nr_of_records), plotted_temperature_series , 'meas. point', 'temp[°C]', 'o', color='red')

            # plt.subplot(1, nof_plots, 3)
            # vc_tick_mod_40 = vc_tick % 40
            # # vc_tick_mod_40 = np.empty( len(vc_tick) )
            # # for i in range(0, len(vc_tick)): vc_tick_mod_40[i] = vc_tick[i] % 40
            # plot_utils.plot_xy(range(0, len(vc_tick_mod_40)), vc_tick_mod_40, 'meas. point', 'vc tick mod 40', 'o', color='green')

            # plt.show()
        #else:
        #    plot_ncols = 0

        #return plot_ncols


    # Deprected (unless the special case when all records are processed).
    # Use 'prepare_processing_cycles()' then 'processing_cycle_data_get()' instead.
    # Can only be called after entire data file is load by 'load()' function
    # If 'nof_concatenated_file_cycles == ALL_MEASUREMENTS_TOGETHER',  then all records are returned
    # ~cycle_data_get~
    def file_cycle_data_get(self, file_cycle_index, nof_concatenated_file_cycles, calc_var, instr, mag_precorr_functions = [], conv_to_NED = False, mag_indeces = [0], cycle_point_filters = [], return_data="all" ):

        assert isinstance(return_data, str)
        # Checking if data have already been loaded:
        assert(self.number_of_absmeas_processing_windows > 0)
        # Checking cycle_id
        assert( file_cycle_index < len(self.number_of_absmeas_processing_windows) )

        if nof_concatenated_file_cycles == ALL_MEASUREMENTS_TOGETHER:  nof_concatenated_file_cycles = self.number_of_absmeas_processing_windows
        assert(nof_concatenated_file_cycles>=1)

        file_record_indeces = self.get_filtered_file_record_indeces(self.file_absmeas_window_data[file_cycle_index]['start_record_index'], self.file_absmeas_window_data[file_cycle_index + nof_concatenated_file_cycles - 1]['end_record_index'], cycle_point_filters )

        cycle_file_id_str = self.get_cycle_file_id_str(file_record_indeces)

        d_for_one_cycle = sads_raw_data[ file_record_indeces  ]

        if (return_data=="cycle_file_id_str"):
            return cycle_file_id_str


        datetime    = d_for_one_cycle[:,COMPNAME__TIMESTAMP]

        MIRA_id_start = dim_azimuth.use_MIRA_record_from_table( datetime[0], instr  )
        MIRA_id_end   = dim_azimuth.use_MIRA_record_from_table( datetime[-1], instr, print_report = False )
        if MIRA_id_start != MIRA_id_end:
            resp=mu.ask_user("MIRA changed within the cycle!", [ ["a","abort"], ["s", "use starting MIRA"], ["e", "use ending-MIRA"] ])
            assert False, "Check code below"
            if resp == "a": sys.exit(1)
            elif resp == "s": dim_azimuth.use_MIRA_record_from_table( d_for_one_cycle[0,  COMPNAME__TIMESTAMP], instr )
            elif resp == "e": dim_azimuth.use_MIRA_record_from_table( d_for_one_cycle[-1,  COMPNAME__TIMESTAMP], instr  )
            else:
                assert False

        mag = d_for_one_cycle[:,self.compindeces_mag]

        meas_time   = d_for_one_cycle[COMPNAME__MOVETOPOINT_TIME]
        if self.compindex_T_electronics != None:
            temp =   d_for_one_cycle[:,self.compindex_T_electronics]
        else:
            temp = None

        hc_raw_tick = d_for_one_cycle[:,COMPNAME__HC_TICK]
        vc_raw_tick = d_for_one_cycle[:,COMPNAME__VC_TICK]

        hc_rad = instr.hc_axis.tick_to_rad(hc_raw_tick, calc_var.hc_precorr_funcs) 
        vc_rad = instr.vc_axis.tick_to_rad(vc_raw_tick, calc_var.vc_precorr_funcs) 

        for mag_correction_function in mag_precorr_functions:
            assert  not mag_correction_function.is_angle_correction()
            # for i in range(0,len(mag)): mag[i] = mag_correction_function.calc(mag[i])
            mag[:,0] = mag_correction_function.calc(mag[:,0])

        if conv_to_NED:
            assert False
            hc_rad = dim_azimuth.h_theo_to_ned(hc_rad)
            vc_rad = dim_azimuth.v_theo_to_ned(vc_rad)

        hc_deg = np.degrees(hc_rad)
        vc_deg = np.degrees(vc_rad)

        if self.compindex_F >=0:
            F = d_for_one_cycle[:,self.compindex_F]
        else:
            F = None

        if self.has_camera_data:
            cam_X   = d_for_one_cycle[:,self.compindex_cam_X]
            cam_Y   = d_for_one_cycle[:,self.compindex_cam_Y]
            cam_ref_ID = d_for_one_cycle[:,self.compindex_cam_ref_ID]
            if self.compindex_cam_avg_pixel_value>0:
                cam_avg_pixel_value = d_for_one_cycle[:,self.compindex_cam_avg_pixel_value]
            else:
                cam_avg_pixel_value = np.zeros( len(cam_X) )
            cam_data = [cam_X, cam_Y, cam_ref_ID, cam_avg_pixel_value ]
        else:
            cam_X = np.nan
            cam_Y = np.nan
            cam_ref_ID = -1
            cam_data = []

        # mag_corrected = np.copy(mag)
        # if mag_tick4_corr_phase != None:
        #     assert len(angle_corr.hc_mag_modulo4_corr) == MAG_TICK_MODULO
        #     for i in range( len(mag_corrected) ) :
        #         assert( hc_raw_tick[i] == round(hc_raw_tick[i]) ) # check if raw tick is integer
        #         mag_corrected[i] =   mag_corrected[i]  + angle_corr.hc_mag_modulo4_corr[ round(hc_raw_tick[i] + mag_tick4_corr_phase)  %  MAG_TICK_MODULO  ]


        method_info =  self.get_program_type_info_str(d_for_one_cycle)


        self.records_in_last_read_file_cycle = len(datetime)

        # plot_utils.plot_xy(hc_raw_tick, vc_raw_tick, "hc raw tick", "vc raw tick", Xaxis_range = [ 0, dim_instr.ENCODER_DIVISION ]  , Yaxis_range = [ 0, dim_instr.ENCODER_DIVISION ] )

        # self.read_cycle_vc0 = instr.vc0

        self.last_returned_cycle_data = d_for_one_cycle

        # nr_of_removed = du.remove_nans(mag, [ datetime, hc_rad, vc_rad, hc_deg, vc_deg, hc_raw_tick, vc_raw_tick, hc_tick, vc_tick, meas_duration, program_maintype_id, F, temp ]):
        #    print("Removed nan mag points: {}".format(nr_of_removed))

        return ( cycle_file_id_str, d_for_one_cycle[:,COMPNAME__FILE_ABMEAS_WIN_INDEX] ), datetime, hc_deg, vc_deg, hc_rad, vc_rad, mag, hc_raw_tick, vc_raw_tick, hc_corrected_tick, vc_corrected_tick, meas_time, method_info, F, temp, cam_data


    # Can be used after 'cycle_data_get()' to get some extra data
    def last_cycle_data_extra_info_get(self, requested_component_names ):
        result = []
        for component_name in requested_component_names:
            if component_name == "full_method_str":  # Method string for each of the records:
                full_program_maintype_str_list = []
                for i in range( len(self.last_returned_cycle_data) ):
                    program_maintype_id = int(round(self.last_returned_cycle_data[i][self.compindex_methodstr_id]))
                    program_maintype_str     = self.program_maintype_str_list[program_maintype_id].decode('UTF-8')
                    if self.compindex_method_extrainfo_str_id>0:
                        method_extra_info_id = int(round(self.last_returned_cycle_data[i][self.compindex_method_extrainfo_str_id]))
                        method_extra_str= self.meas_type_extra_info_str_list[ method_extra_info_id ]
                        if method_extra_str != b"":
                            program_maintype_str = program_maintype_str + "/" + method_extra_str.decode('UTF-8')
                    full_program_maintype_str_list = full_program_maintype_str_list + [ program_maintype_str ]
                result = result + [ full_program_maintype_str_list ]
            else:

                compindex = self.get_compindex_by_name(component_name)
                result = result + [ self.last_returned_cycle_data[:,compindex] ]
        return result


    def get_compindex_by_name(self, component_name):
        if component_name=="pos_id":            return self.compindex_pos_id
        if component_name=="movetopoint_time":  return self.compindex_movetopoint_time
        if component_name=="tickcount":         return self.compindex_tickcount
        if component_name=="file_vc_deg":       return self.compindex_file_vc_deg
        if component_name=="temp":              return self.compindex_T_electronics
        if component_name=="M_avg_stddev[:]":   return self.compindeces_M_avg_stddev
        if component_name=="M_avg_slope[:]":    return self.compindeces_M_avg_slope
        assert False


    # ~print_new_cycle_begin~
    # def print_new_file_cycle_begin(self, cycle_count, nof_concatenated_file_cycles, extra_str, fig = None, nof_plots_added_later  = -1):
    #     # cycle_start_time = self.file_records[ self.file_cycle_start_index[cycle_count] , COMPNAME__TIMESTAMP]
    #     if self.file_records_str is not None:
    #         starting_record_str = self.split_record_str(self.file_records_str[ self.file_absmeas_window_data[cycle_count]['start_record_index'] ])
    #     else:
    #         starting_record_str = "QQQ"
    #     cycle_start_time_str = starting_record_str[0]
    #     cycle_id_str = starting_record_str[2]
    #     main_line = "---- CYCLE:COUNT:{}{} ID:{} START:{} RECORDS:{}  {} ---------".format(cycle_count, "" if nof_concatenated_file_cycles == 1 else "-{}".format(cycle_count+nof_concatenated_file_cycles-1), cycle_id_str, cycle_start_time_str,  self.records_in_last_read_file_cycle,extra_str)
    #     print("-" * len(main_line))
    #     print(main_line)
    #     print("-" * len(main_line))


    #     # 'fig' has to be created by the caller: it can contain one or more plots - here we show cycle information
    #     if fig !=  None:

    #         assert nof_plots_added_later >= 0

    #         plot_ncols =  3 if self.compindex_T_electronics >=0 else 2
    #         # Caller should prepared the figure like 'fig = plt.figure(figsize=(preferences.FIGURE_SIZE[0],3 + nof_plots_added_later * 3))'

    #         cycle_file_id_info, datetime, hc_deg, vc_deg, hc, vc, mag, hc_raw_tick, vc_raw_tick, hc_tick, vc_tick, meas_time, method_info, F, temp, cam_data = self.file_cycle_data_get(cycle_count,nof_concatenated_file_cycles)

    #         fig.suptitle('Cycle-file-id:{}  method:{:} start-date:{:}'.format(cycle_file_id_info[0], method_info[1], mu.exceltime_2_datetime(datetime[0]).strftime("%Y.%m.%d.%H:%M:%S") ), fontsize=14)
    #         plt.subplot(1 + nof_plots_added_later, plot_ncols, 1)
    #         plot_utils.plot_xy(hc_deg, vc_deg, 'hc [°]', 'vc [°]', 'o', color='green', close_plot = False)
    #         plt.subplot(1 + nof_plots_added_later, plot_ncols, 2)
    #         plot_utils.plot_xy(range(0, len(hc_deg)), meas_time, 'meas. point', 'meas_time [sec]', 'o', color='green', close_plot = False)
    #         if plot_ncols>2 :
    #             plt.subplot(1 + nof_plots_added_later, plot_ncols, 3)
    #             plot_utils.plot_xy(range(0, len(hc_deg)), temp, 'meas. point', 'temp[°C]', 'o', color='red', close_plot = False)

    #         # plt.subplot(1, nof_plots, 3)
    #         # vc_tick_mod_40 = vc_tick % 40
    #         # # vc_tick_mod_40 = np.empty( len(vc_tick) )
    #         # # for i in range(0, len(vc_tick)): vc_tick_mod_40[i] = vc_tick[i] % 40
    #         # plot_utils.plot_xy(range(0, len(vc_tick_mod_40)), vc_tick_mod_40, 'meas. point', 'vc tick mod 40', 'o', color='green')

    #         # plt.show()
    #     #else:
    #     #    plot_ncols = 0

    #     #return plot_ncols


    # Checks if a processing's cycle points' belongs to only one or to more 'cycle_file_id'.
    # Then returns a string accordingly.
    def get_cycle_file_id_str(self, file_record_indeces ):
        # Original solution (until 2024.04):
        # if np.all( d_for_one_cycle[:,COMPNAME__FILE_ABMEAS_WIN_INDEX] == d_for_one_cycle[0,COMPNAME__FILE_ABMEAS_WIN_INDEX] ) :
        #    cycle_file_id_str = "{}".format(round(d_for_one_cycle[0,COMPNAME__FILE_ABMEAS_WIN_INDEX]))
        #else:
        #    cycle_file_id_str = "more-cycles"

        ids = []
        for file_record_index in file_record_indeces:
            rec_id_str = str(self.M_readings[file_record_index][COMPNAME__FILE_ABMEAS_WIN_INDEX])
            if not rec_id_str in ids:
                ids += [ rec_id_str ]
        if len(ids)==1: return ids[0]
        result="["
        for id in ids: result = result + "," + id
        result = result[:len(result)-1] + "]"
        return result


    # Formats 'program_maintype_id' (integer) and 'full_method_name' (string) for a file record array, stored in 'd' (it is a view of raw_data).
    # Checks if method str changes within cycle / returns one string that is representative to all of the records is 'd'.
    # Checking if there are one or more methods within this cycle.
    def get_program_type_info_str(self,d):  # ~get_method_info_str'
        if np.all( d[COMPNAME__PROGRAM_MAINTYPE_IDX] == d[COMPNAME__PROGRAM_MAINTYPE_IDX][0]  ):
            # There is only one main measurement type
            program_maintype_id =  d[COMPNAME__PROGRAM_MAINTYPE_IDX][0]
            # program_name_str= d[COMPNAME__PROGRAM_MAINTYPE_STR][0].decode('UTF-8')
            program_name_str= self.program_maintype_str_list[program_maintype_id]
        else:
            program_maintype_id = -1
            program_name_str = "__mixed__"

        if self.file_format == "ADS-raw-dat-files":
            if np.all( d[COMPNAME__PROGRAM_PARAM1_IDX] == d[COMPNAME__PROGRAM_PARAM1_IDX][0]  ):
                method_str_extrainfo =  self.meas_type_param1_str_list[ d[COMPNAME__PROGRAM_PARAM1_IDX][0] ]
            else:
                if program_name_str in [MEAS_TYPE_STR__NULLMETH_1, MEAS_TYPE_STR__NULLMETH_2]:
                    method_str_extrainfo =""  # In null method, the 'param1' is the position id
                else:
                    method_str_extrainfo = "__mixed__"
        else:
            method_str_extrainfo = ""
            
        if  method_str_extrainfo=='':
            full_program_name = program_name_str
        else:
            full_program_name =  program_name_str + "/" + method_str_extrainfo

        return  [ program_maintype_id,  full_program_name ]

    def nr_mag_sensors_get(self):
        assert ( self.number_of_absmeas_processing_windows > 0)
        return len(self.compindeces_mag)

    def nr_of_load_records(self):
        return len(self.M_readings)

    def get_timelabel_minmax(self):
        assert ( self.number_of_absmeas_processing_windows > 0)
        return raw_data[COMPNAME__TIMESTAMP][0], raw_data[COMPNAME__TIMESTAMP][-1]

    #def ads_dat_file__get_meas_type_param(self, param_number, ads_raw_data_record):
    #    if param_number == 1:
    #        return self.meas_type_param1_str_list[ ads_raw_data_record[COMPNAME__PROGRAM_PARAM1_IDX] ]
    #    elif param_number == 2:
    #        return self.meas_type_param2_str_list[ ads_raw_data_record[COMPNAME__PROGRAM_PARAM2_IDX] ]
    #    else:
    #        assert False

    def processing_cycle_parallel_mag_data_get(self, processing_window_data):
        file_record_indeces = processing_window_data.M_record_indeces
        missing_record_indeces = np.intersect1d( file_record_indeces, self.load_main_obs_data_set_data_missing_record_indeces )
        returned_mag_data_in_window = None # self.load_main_obs_data_set_data_with_window
        return self.load_main_obs_data_set_data[file_record_indeces], returned_mag_data_in_window, missing_record_indeces




